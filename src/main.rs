use rosalind::{complement_base, parse_dna, parse_fasta, parse_prot, parse_rna, DNABase};
use std::collections::HashMap;
use std::time::Instant;

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

pub fn grph() {
    let filename = "data/grph.txt";
    let string = std::fs::read(filename).unwrap();
    let strands = parse_fasta(&string).unwrap();
    for (i, left) in strands.iter().enumerate() {
        for (j, right) in strands.iter().enumerate() {
            if i == j {
                continue;
            }
            if left.strand.data[left.strand.data.len() - 3..] == right.strand.data[..3] {
                println!("{} {}", left.name, right.name);
            }
        }
    }
}

pub fn long() {
    let filename = "data/long.txt";
    let string = std::fs::read(filename).unwrap();
    let strands = parse_fasta(&string).unwrap();

    let superstring = strands.superstring().unwrap();
    for c in superstring.data {
        print!("{:?}", c);
    }
}

pub fn fib() {
    let filename = "data/fib.txt";
    let string = std::fs::read(filename).unwrap();
    let mut string = String::from_utf8(string).unwrap();
    string = string.trim().to_string();
    let tokens = string.split(' ').collect::<Vec<_>>();

    let n: usize = tokens[0].parse().unwrap();
    let k: usize = tokens[1].parse().unwrap();

    let mut a = 1;
    let mut b = 1;

    for _ in 0..n - 2 {
        let c = k * a + b;
        a = b;
        b = c;
    }
    println!("{:?}", b);
}

pub fn fibd() {
    let filename = "data/fibd.txt";
    let string = std::fs::read(filename).unwrap();
    let mut string = String::from_utf8(string).unwrap();
    string = string.trim().to_string();
    let tokens = string.split(' ').collect::<Vec<_>>();

    let n: u64 = tokens[0].parse().unwrap();
    let m: usize = tokens[1].parse().unwrap();

    let mut babies = vec![1u64, 0];
    let mut adults = vec![0u64, 1];
    for _ in 0..n - 2 {
        let dying = *babies.get(babies.len() - m).unwrap_or(&0);
        let last_adults = *adults.last().unwrap();
        let newbabies = last_adults;
        let grownups = *babies.last().unwrap();
        let newadults = last_adults - dying + grownups;

        babies.push(newbabies);
        adults.push(newadults);
    }

    let c = adults[adults.len() - 1] + babies[babies.len() - 1];
    println!("{:?}", c);
}
pub fn mrna() {
    let filename = "data/mrna.txt";
    let string = std::fs::read(filename).unwrap();
    let prot = parse_prot(&string[..string.len() - 1]).unwrap();

    let mut p = 1;
    for aa in prot {
        let n = aa.encodings().len();
        p *= n;
        p = p % 1_000_000;
    }
    println!("{:?}", p);
}
pub fn cons() {
    let filename = "data/cons.txt";
    let string = std::fs::read(filename).unwrap();
    let strands = parse_fasta(&string).unwrap();

    let profile = strands.profile_matrix().unwrap();

    for possibilities in &profile {
        let index = possibilities
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.cmp(b))
            .map(|(index, _)| index)
            .unwrap();
        print!("{:?}", DNABase::bases()[index]);
    }
    println!("");
    for i in 0..4 {
        print!("{:?}: ", DNABase::bases()[i]);
        for p in profile.iter() {
            print!("{:?} ", p[i]);
        }
        println!("");
    }
}

fn main() {
    // dna();
    // dna_to_rna();
    // complement();
    // gc();
    // subs();
    // prot();
    // grph();
    // long();
    // fib();
    // fibd();
    // mrna();
    // cons();
    let start = Instant::now();

    println!("Elapsed {:?}", start.elapsed());
}
