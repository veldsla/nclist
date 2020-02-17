# NClist
A nested containment list (NClist) is a datastructure that can be queried for elements
overlapping intervals. It was invented and published by Alexander V and Alekseyenko Christopher
J. Lee in Bioinformatics in
2007 (doi: [10.1093/bioinformatics/btl647](https://doi.org/10.1093/bioinformatics/btl647)).

## How it works
The `NClist` internals rely on the observation that when a set of intervals, where all are
non-contained (based on their interval bounds) in any of the other intervals in the set, are
sorted on their start coordinate are also sorted on their end coordinate. If this requirement
is fulfilled the items overlapping an interval can be found by a binary search on the query
start and returning items until the query end coordinate has been passed, giving a complexity
of `O(log N + M)` where N is the size of the set and M is the number of overlaps.

The only remaining problem is the intervals that **are** contained in another interval. This was
solved by taking out these intervals and storing them in a separate set, linking this set to the
original interval. Now when you search for overlaps you check for contained intervals and also
search these nested sets. This can be implemented recursively (as shown in the paper) or
using a queue (which was used for this crate).

The linked article also provides details about an on-disk version that can also be efficiently
searched, but this crate implementation is in-memory and stores the items in a (single) `Vec`.

## How to use
You can create a searchable `NClist<T>` from a `Vec<T>` if you implement the `Interval` trait
for `T` The `Interval` trait also requires that `T` is `Ord`. Creating the NClist validates
that the end coordinate is greater than start. This means negative and zero-width intervals
cannot be used in an `NClist<T>`. 

## Example
```rust
use nclist::NClist;
// Given a set of `T` where `T` implements the `Interval` trait
let v = vec![(10..15), (10..20), (1..8)];
// Create the NClist, this consumes v
let nc = NClist::from(v);
// count overlaps, the query is provided as a reference to a `std::ops::Range`
assert_eq!(nc.count_overlaps(&(10..12)), 2);
// remember, intervals are half open
assert_eq!(nc.count_overlaps(&(20..30)), 0);

//or query them using an iterator
let mut q = nc.overlaps(&(7..10));
assert_eq!(q.next(), Some(&(1..8)));
assert_eq!(q.next(), None);

```
More examples can be found in the `examples` directory on Github

## Recommendations for use
The `NClist<T>` is not mutable. Any mutable access to the items in the could invalidate the interval
bounds (interior mutability using for example a `RefCell` could solve this). Also insertion and
deletion are not supported. I can speculate that an interval tree would also be a better for
this type of access. For usage in Bioinformatics where interval data is often provided as
(sorted) lists (gff, gtf, bed) the `NClist<T>` is a perfect fit and has very nice ergonomics.

Obviously the implemtation works better when nesting depth is limited.
