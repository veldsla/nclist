use std::ops::Range;

use criterion::{Criterion, criterion_main, criterion_group};
use nclist::NClist;
use rand::Rng;

fn make_ranges(n: u32, min_width: u32, max_width: u32, l: u32) -> Vec<Range<u32>> {
    let mut rng = rand::thread_rng();
    (0..n).map(|_| {
        let width = rng.gen_range(min_width, max_width);
        let start = rng.gen_range(0, l-width);
        start..start + width
    }).collect()
}

fn query_nclist<T: Ord>(nclist: &NClist<Range<T>>, q: &[Range<T>]) {
    for o in q {
        let _res = nclist.overlaps(o).count();
    }
}

fn test_nc(c:  &mut Criterion) {
    // 100 ranges of size 100-200 (total size 10-20K in 10K space)
    let ranges = make_ranges(100, 100, 200, 10_000);
    let nclist = NClist::from_vec(ranges).unwrap();

    // Query 100 short ranges
    let q = make_ranges(100,10,50,10_000);

    c.bench_function("NClist", move |b| b.iter(|| query_nclist(&nclist, &q)));
}    

fn test_nc_overlaps(c:  &mut Criterion) {
    let ranges: Vec<_> = (100..2000).step_by(200).map(|p| make_ranges(1000, 50, p, 10_000)).collect();
    let nclists: Vec<_> = ranges.into_iter().map(|r| NClist::from_vec(r).unwrap()).collect();

    let q = make_ranges(100,10,50,10_000);

    c.bench_function_over_inputs("NClist vs fraction overlap", move |b, &&i| b.iter(|| query_nclist(&nclists[i], &q)),
        &[0usize,1,2,3,4,5,6,7,8,9]);
}   
criterion_group!(benches, test_nc, test_nc_overlaps);
criterion_main!(benches);

