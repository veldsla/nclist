//! A nested containment list (NClist) is a datastructure that can be queried for elements
//! overlapping intervals. It was invented and published by Alexander V and Alekseyenko Christopher
//! J. Lee in Bioinformatics in
//! 2007 (doi: [10.1093/bioinformatics/btl647](https://doi.org/10.1093/bioinformatics/btl647)).
//!
//! # How it works
//! The `NClist` internals rely on the observation that when a set of intervals, where all are
//! non-contained (based on their interval bounds) in any of the other intervals in the set, are
//! sorted on their start coordinate are also sorted on their end coordinate. If this requirement
//! is fulfilled the items overlapping an interval can be found by a binary search on the query
//! start and returning items until the query end coordinate has been passed, giving a complexity
//! of `O(log N + M)` where N is the size of the set and M is the number of overlaps.
//!
//! The only remaining problem is the intervals that **are** contained in another interval. This was
//! solved by taking out these intervals and storing them in a separate set, linking this set to the
//! original interval. Now when you search for overlaps you check for contained intervals and also
//! search these nested sets. This can be implemented recursively (as shown in the paper) or
//! using a queue (which was used for this crate).
//!
//! The linked article also provides details about an on-disk version that can also be efficiently
//! searched, but this crate implementation is in-memory and stores the items in a (single) `Vec`.
//!
//! # How to use
//! You can create a searchable `NClist<T>` from a `Vec<T>` if you implement the `Interval` trait
//! for `T` The `Interval` trait also requires that `T` is `Ord`. Creating the NClist validates
//! that the end coordinate is greater than start. This means negative and zero-width intervals
//! cannot be used in an `NClist<T>`. 
//!
//! # Example
//! ```
//! use nclist::NClist;
//! // Given a set of `T` where `T` implements the `Interval` trait
//! let v = vec![(10..15), (10..20), (1..8)];
//! // Create the NClist, this consumes v
//! let nc = NClist::from(v);
//! // count overlaps, the query is provided as a reference to a `std::ops::Range`
//! assert_eq!(nc.count_overlaps(&(10..12)), 2);
//! // remember, intervals are half open
//! assert_eq!(nc.count_overlaps(&(20..30)), 0);
//!
//! //or query them using an iterator
//! let mut q = nc.overlaps(&(7..10));
//! assert_eq!(q.next(), Some(&(1..8)));
//! assert_eq!(q.next(), None);
//!
//! ```
//! More examples can be found in the `examples` directory on Github
//!
//! # Recommendations for use
//! The `NClist<T>` is not mutable. Any mutable access to the items in the could invalidate the interval
//! bounds (interior mutability using for example a `RefCell` could solve this). Also insertion and
//! deletion are not supported. I can speculate that an interval tree would also be a better for
//! this type of access. For usage in Bioinformatics where interval data is often provided as
//! (sorted) lists (gff, gtf, bed) the `NClist<T>` is a perfect fit and has very nice ergonomics.
//!
//! Obviously the implemtation works better when nesting depth is limited.
use std::collections::VecDeque;
use std::ops::Range;

use itertools::Itertools;

/// The interval trait needs to be implemented for `T` before you can create an `NClist<T>`.
/// An interval is half-open, inclusive start and exclusive end (like `std::ops::Range<T>`), but 
/// `end > start` must always be true.
pub trait Interval {
    /// The coordinate type of the interval. This type must implement `std::ops::Ord`
    type Coord: Ord;

    /// Return a reference to the start coordinate of the interval (inclusive)
    fn start(&self) -> &Self::Coord;

    /// Return a reference to the end coordinate of the interval (non-inclusive)
    fn end(&self) -> &Self::Coord;
}

/// Interval is already implemented for `std::ops::Range`.
impl<N> Interval for Range<N> where N: Ord{
    type Coord = N;

    #[inline(always)]
    fn start(&self) -> &Self::Coord {
        &self.start
    }

    #[inline(always)]
    fn end(&self) -> &Self::Coord {
        &self.end
    }
}

#[derive(Debug)]
pub struct NClist<T> where T: Interval {
    intervals: Vec<T>,
    contained: Vec<Option<(usize, usize)>>
}

struct SlicedNClist<'a, T> where T: 'a + Interval {
    intervals: &'a [T],
    contained: &'a [Option<(usize, usize)>],
    stop_at: &'a T::Coord
}

pub struct Overlaps<'a, T> where T: 'a + Interval {
    nclist: &'a NClist<T>,
    range:  &'a Range<T::Coord>,
    current_pos: usize,
    current_end: usize,
    sublists: VecDeque<(usize, usize)>,
}

pub struct OrderedOverlaps<'a, T> where T: 'a + Interval {
    nclist: &'a NClist<T>,
    range:  &'a Range<T::Coord>,
    current: SlicedNClist<'a, T>,
    queue: Vec<SlicedNClist<'a, T>>
}

impl<T> NClist<T> where T: Interval {
    fn new() -> NClist<T> {
        NClist { intervals: Vec::new(), contained: vec![Some((0,0))] }
    }

    /// Count the number of elements overlapping the `Range` r. Counting overlaps is slightly
    /// faster than iterating over the overlaps. This method is preferred when only the number of
    /// overlapping elements is required.
    pub fn count_overlaps(&self, r: &Range<T::Coord>) -> usize {
        if r.end <= r.start {
            return 0;
        }
        let mut count = 0;
        let mut queue = VecDeque::new();
        queue.push_back(self.contained[0].unwrap());
        while let Some((start, end)) = queue.pop_front() {
            self.slice(start, end, &r.start, &r.end)
                .for_each(|(_, contained)| {
                    count += 1;
                    if let Some(subrange) = *contained {
                        queue.push_back(subrange);
                    }
                });
        }
        count
    }

    /// Returns an iterator that returns overlapping elements to query `r`. During iteration
    /// contained intervals are pushed to a queue an processed in order after yielding the
    /// non-overlapping regions.
    pub fn overlaps<'a>(&'a self, r: &'a Range<T::Coord>) -> Overlaps<'a , T> {
        let current_slice = self.contained[0].as_ref().unwrap();
        //empty or negative width intervals do not overlap anything
        let start = if r.end > r.start {
            self.bin_search_end(current_slice.0, current_slice.1, &r.start)
        } else {
            current_slice.1
        };

        Overlaps { nclist: self, range: r, current_pos: start, current_end: current_slice.1, sublists: VecDeque::new() }
    }

    /// Returns an iterator that returns overlapping elements to query `r` ordered by start
    /// coordinate. This is less efficient that returning without ordering, but doesn't require
    /// allocating storage for all overlapping elements.
    pub fn overlaps_ordered<'a>(&'a self, r: &'a Range<T::Coord>) -> OrderedOverlaps<'a , T> {
        let &(mut start, end) = self.contained[0].as_ref().unwrap();
        if r.end <= r.start {
            start = end;
        }
        OrderedOverlaps { nclist: self, range: r, current: self.slice(start, end, &r.start, &r.end), queue: Vec::new() }
    }

    /// Return the intervals `Vec`. This will run without allocation and return the intervals in a
    /// different order then provided.
    pub fn into_vec(self) -> Vec<T> {
        self.into()
    }

    #[inline]
    fn slice<'a>(&'a self, mut start: usize, end: usize, q:  &T::Coord, q_end: &'a T::Coord) -> SlicedNClist<'a, T> {
        start += match self.intervals[start..end].binary_search_by(|e| e.end().cmp(q))
        {
            Ok(n) => n + 1,
            Err(n) => n
        };
        SlicedNClist { intervals: &self.intervals[start..end], contained: &self.contained[start+1..end+1], stop_at: q_end }
    }

    #[inline]
    fn bin_search_end(&self, start: usize, end: usize, q: &T::Coord) -> usize {
        match self.intervals[start..end].binary_search_by(|e| e.end().cmp(q)) {
            Ok(n) => n + 1,
            Err(n) => n
        }
    }
}

impl<'a, T> Iterator for SlicedNClist<'a, T> where T: Interval {
    type Item = (&'a T, &'a Option<(usize, usize)>);
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some((i, ref mut intervals)) = self.intervals.split_first() {
            if i.start() >= self.stop_at {
                None
            } else {
                let (c, ref mut contained) = self.contained.split_first().unwrap();
                std::mem::replace(&mut self.intervals, intervals);
                std::mem::replace(&mut self.contained, contained);
                Some((i, c))
            }
        } else {
            None
        }
    }
}

impl<'a, T> Iterator for Overlaps<'a, T> where T: Interval {
    type Item = &'a T;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let remaining = self.current_end - self.current_pos;

        if remaining == 0 || *self.nclist.intervals[self.current_pos].start() >= self.range.end {
            if let Some((mut new_start, new_end)) = self.sublists.pop_front() {
                new_start += self.nclist.bin_search_end(new_start, new_end, &self.range.start);
                self.current_pos = new_start;
                self.current_end = new_end;
                self.next()
            } else {
                None
            }
        } else {
            let pos = self.current_pos;
            self.current_pos += 1;
            if let Some(next_sublist) = self.nclist.contained[self.current_pos] {
                self.sublists.push_back(next_sublist);
            }
            Some(&self.nclist.intervals[pos])
        }
    }
}

impl<'a, T> Iterator for OrderedOverlaps<'a, T> where T: Interval {
    type Item = &'a T;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some((interval, contained)) = self.current.next() {
            if let Some((start, end)) = *contained {
                let mut ns = self.nclist.slice(start, end, &self.range.start, &self.range.end);
                std::mem::swap(&mut self.current, &mut ns);
                self.queue.push(ns);
            }
            Some(interval)
        } else if let Some(sublist) = self.queue.pop() {
            self.current = sublist;
            self.next()
        } else {
            None
        }
    }
}

/// Internal intermediate sublist used for creating `NClist<T>`
struct NClistBuilder<T> {
    intervals: Vec<T>,
    contained_pos: usize,
}

fn build_nclist<T: Interval>(sublists: &mut VecDeque<NClistBuilder<T>>, result: &mut NClist<T>) {
    if let Some(sublist) = sublists.pop_front() {
        //iterate over all ranges and take out contained intervals
        let mut it = sublist.intervals.into_iter().peekable();

        let sublist_start = result.intervals.len();
        while let Some(e) = it.next() {
            let interval_pos = result.intervals.len();
            let contained: Vec<_> = it.peeking_take_while(|n| n.end() < e.end()).collect();
            if !contained.is_empty() {
                sublists.push_back(NClistBuilder {intervals: contained, contained_pos: interval_pos + 1});
            }
            result.intervals.push(e);
            result.contained.push(None);
        }

        //store the position and the length of the sublist
        result.contained[sublist.contained_pos] = Some((sublist_start, result.intervals.len()));
    }
}

/// This is currently the only way to create an `NClist<T>`.
impl<T> From<Vec<T>> for NClist<T> where T: Interval {
    fn from(mut v: Vec<T>) -> Self {
        if v.iter().any(|e| e.end() <= e.start()) {
            panic!("Cannot use intervals with zero or negative width");
        }
        v.sort_by(|a, b| a.start().cmp(b.start())
                  .then(a.end().cmp(b.end()).reverse()));

        let mut list = NClist::new();
        let mut sublists = VecDeque::from(vec![NClistBuilder { intervals: v, contained_pos: 0}]);

        while !sublists.is_empty() {
            build_nclist(&mut sublists, &mut list);
        }
        list
    }
}

/// Return the intervals `Vec`. This will run without allocation and return the intervals in a
/// different order then provided.
impl<T> Into<Vec<T>> for NClist<T> where T: Interval {
    fn into(self) -> Vec<T> {
        self.intervals
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // This test is copied from Rust's stdlib. This software relies on the fact that a binary
    // search returns the last matching element. If the stdlib implementation changes this should
    // be caught.
    //
    // See:
    // https://github.com/rust-lang/rust/blob/975e83a32ad8c2c894391711d227786614d61a50/src/libcore/tests/slice.rs#L68
    #[test]
    fn test_binary_search_implementation_details() {
        let b = [1, 1, 2, 2, 3, 3, 3];
        assert_eq!(b.binary_search(&1), Ok(1));
        assert_eq!(b.binary_search(&2), Ok(3));
        assert_eq!(b.binary_search(&3), Ok(6));
        let b = [1, 1, 1, 1, 1, 3, 3, 3, 3];
        assert_eq!(b.binary_search(&1), Ok(4));
        assert_eq!(b.binary_search(&3), Ok(8));
        let b = [1, 1, 1, 1, 3, 3, 3, 3, 3];
        assert_eq!(b.binary_search(&1), Ok(3));
        assert_eq!(b.binary_search(&3), Ok(8));
    }

    #[test]
    fn from_vec() {
        let list: Vec<Range<u64>> = vec![(10..15), (10..20), (1..8)].into_iter().collect();
        let nclist = NClist::from(list);

        assert_eq!(nclist.intervals.len(), 3);
        assert!(nclist.contained[0].is_some());
        assert!(nclist.contained[1].is_none());
        assert!(nclist.contained[2].is_some());
        assert!(nclist.contained[3].is_none());

        let list: Vec<Range<u64>> = Vec::new();
        let nclist = NClist::from(list);
        assert_eq!(nclist.intervals.len(), 0);
    }

    #[test]
    #[should_panic]
    fn interval_width() {
        let list: Vec<Range<u64>> = vec![(5..20), (19..20), (7..7)].into_iter().collect();
        let _nclist = NClist::from(list);
    }

    #[test]
    fn illegal_width_queries() {
        let list: Vec<Range<u64>> = vec![(5..20), (19..20), (7..8)].into_iter().collect();
        let nclist = NClist::from(list);
        assert_eq!(nclist.count_overlaps(&(7..7)), 0);
        assert_eq!(nclist.count_overlaps(&(8..7)), 0);
        assert_eq!(nclist.count_overlaps(&(19..19)), 0);
        assert_eq!(nclist.overlaps(&(7..7)).count(), 0);
        assert_eq!(nclist.overlaps(&(8..7)).count(), 0);
        assert_eq!(nclist.overlaps(&(19..19)).count(), 0);
        assert_eq!(nclist.overlaps_ordered(&(7..7)).count(), 0);
        assert_eq!(nclist.overlaps_ordered(&(8..7)).count(), 0);
        assert_eq!(nclist.overlaps_ordered(&(19..19)).count(), 0);
    }

    #[test]
    fn count() {
        let list: Vec<Range<u64>> = vec![(10..15), (10..20), (1..8)].into_iter().collect();
        let nclist = NClist::from(list);

        assert_eq!(nclist.intervals.len(), 3);
        assert_eq!(nclist.count_overlaps(&(5..20)), 3);
        assert_eq!(nclist.count_overlaps(&(14..18)), 2);
        assert_eq!(nclist.count_overlaps(&(150..180)), 0);
        assert_eq!(nclist.count_overlaps(&(10..10)), 0);
        assert_eq!(nclist.count_overlaps(&(10..11)), 2);
        assert_eq!(nclist.count_overlaps(&(9..10)), 0);
        assert_eq!(nclist.count_overlaps(&(8..9)), 0);
        assert_eq!(nclist.count_overlaps(&(8..10)), 0);
        assert_eq!(nclist.count_overlaps(&(20..100)), 0);

        let list: Vec<Range<u64>> = Vec::new();
        let nclist = NClist::from(list);
        assert_eq!(nclist.count_overlaps(&(100..200)), 0);

    }
    #[test]
    fn overlaps() {
        let list: Vec<Range<u64>> = vec![(10..15), (10..20), (1..8)].into_iter().collect();
        let nclist = NClist::from(list);

        assert_eq!(nclist.intervals.len(), 3);
        assert_eq!(nclist.overlaps(&(5..20)).count(), 3);

        let mut q = nclist.overlaps_ordered(&(5..20));
        assert_eq!(q.next(), Some(&(1..8)));
        assert_eq!(q.next(), Some(&(10..20)));
        assert_eq!(q.next(), Some(&(10..15)));
        assert_eq!(q.next(), None);

        assert_eq!(nclist.overlaps(&(20..100)).count(), 0);
        assert_eq!(nclist.overlaps(&(8..10)).count(), 0);
        assert_eq!(nclist.overlaps(&(8..9)).count(), 0);
    }

    #[test]
    fn duplicate_intervals() {
        let list: Vec<Range<u64>> = vec![(10..15), (11..13), (10..20), (1..8), (11..13), (16..18)].into_iter().collect();
        let nclist = NClist::from(list);
        println!("{:?}", nclist);

        assert_eq!(nclist.overlaps(&(5..20)).count(), 6);
        assert_eq!(nclist.overlaps(&(11..13)).count(), 4);

        let mut q = nclist.overlaps_ordered(&(11..17));
        assert_eq!(q.next(), Some(&(10..20)));
        assert_eq!(q.next(), Some(&(10..15)));
        assert_eq!(q.next(), Some(&(11..13)));
        assert_eq!(q.next(), Some(&(11..13)));
        assert_eq!(q.next(), Some(&(16..18)));
        assert_eq!(q.next(), None);

        assert_eq!(nclist.overlaps(&(20..100)).count(), 0);
        assert_eq!(nclist.overlaps(&(8..10)).count(), 0);
        assert_eq!(nclist.overlaps(&(8..9)).count(), 0);
    }
}
