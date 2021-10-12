use rosalind::{complement_base, parse_dna, parse_fasta, parse_rna};
use std::collections::HashMap;

////////////////////////////////////

pub fn complement() {
    let filename = "data/revc.txt";
    let string = std::fs::read(filename).unwrap();

    let out = string
        .into_iter()
        .filter_map(|c| complement_base(c))
        .rev()
        .collect::<Vec<_>>();
    println!("{}", String::from_utf8(out).unwrap());
}

pub fn gc() {
    let filename = "data/gc.txt";
    let string = std::fs::read(filename).unwrap();
    let strands = parse_fasta(&string).unwrap();

    let max_gc = strands
        .iter()
        .max_by(|a, b| a.gc().partial_cmp(&b.gc()).unwrap())
        .unwrap();
    println!("{}", max_gc.name);
    println!("{:?}", max_gc.gc() * 100.0);
}

pub fn dna_to_rna() {
    let filename = "data/rna.txt";
    let string = std::fs::read(filename).unwrap();
    let out = string
        .into_iter()
        .map(|c| if c == b'T' { b'U' } else { c })
        .collect::<Vec<_>>();
    println!("{}", String::from_utf8(out).unwrap());
}

pub fn dna() {
    let filename = "data/dna.txt";
    let string = std::fs::read(filename).unwrap();
    let mut counter = HashMap::new();

    for c in string {
        *counter.entry(c).or_insert(0) += 1;
    }

    println!(
        "{} {} {} {}",
        counter.get(&b'A').unwrap_or(&0),
        counter.get(&b'C').unwrap_or(&0),
        counter.get(&b'G').unwrap_or(&0),
        counter.get(&b'T').unwrap_or(&0)
    );
}

pub fn subs() {
    let filename = "data/subs.txt";
    let string = std::fs::read(filename).unwrap();
    let tokens = string.split(|c| c == &b'\n').collect::<Vec<_>>();
    let data = parse_dna(tokens[0]).unwrap();
    let pattern = parse_dna(tokens[1]).unwrap();
    for i in pattern.find(&data) {
        print!("{:?} ", i + 1);
    }
}

pub fn prot() {
    let filename = "data/prot.txt";
    let string = std::fs::read(filename).unwrap();
    let rna = parse_rna(&string[..string.len() - 1]).unwrap();
    let proteins = rna.proteins().unwrap();
    for p in &proteins[..proteins.len() - 1] {
        print!("{:?}", p);
    }
}
fn main() {
    // dna();
    // dna_to_rna();
    // complement();
    // gc();
    // subs();
    prot();
}
