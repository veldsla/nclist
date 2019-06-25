#[macro_use]
extern crate criterion;
extern crate itertools;
extern crate nclist;

use std::ops::Range;

use itertools::Itertools;
use criterion::Criterion;

use nclist::NClist;

fn make_ranges(n: usize, width: usize, ) -> NClist<Range<usize>> {
    assert!(n

    let mut v = Vec::new();
    let mut n = 0;
    while n < l {
        width = end - start;
        let this_layer = if (l - n) * range_width < width {
            l - n
        } else {
            width / range_width
        };

        v.extend((0..width).step(range_width).take(this_layer).map(|p| (p..p+range_width)));
        n += this_layer;
        eprintln!("l = {}", v.len());
        start +=1;
        end -= 1;
    }
    NClist::from(v) 
}

fn query_nclist(nclist: &NClist<Range<usize>>, q: &[Range<usize>]) {
    for o in q {
        let res = nclist.overlaps(o).count();
    }
}

fn test_nc(c:  &mut Criterion) {

    let l = make_contained_ranges(10000, 0.8);
    let q: Vec<Range<usize>> = (0..1_000_000).step(1000).map(|s| (s..s+1000)).collect();

    c.bench_function("NClist", move |b| b.iter(|| query_nclist(&l, &q)));
}    


criterion_group!(benches, test_nc);
criterion_main!(benches);


