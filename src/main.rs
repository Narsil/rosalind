use itertools::Itertools;
use num_traits::pow::Pow;
use rosalind::individuals::{coupling, Dd, Individual, DD};
use rosalind::{parse, parse_fasta, DNABase, RNABase, AA};
use std::collections::HashMap;
use std::io::Read;
use std::time::Instant;

////////////////////////////////////

pub fn complement() {
    let filename = "data/revc.txt";
    let string = std::fs::read(filename).unwrap();

    let dna = parse::<DNABase>(&string).unwrap();

    let cdna = dna.complement();
    println!("{}", cdna);
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
    let data = parse::<DNABase>(tokens[0]).unwrap();
    let pattern = parse::<DNABase>(tokens[1]).unwrap();
    for i in pattern.find(&data.strand) {
        print!("{:?} ", i + 1);
    }
}

pub fn prot() {
    let filename = "data/prot.txt";
    let string = std::fs::read(filename).unwrap();
    let rna = parse::<RNABase>(&string[..string.len() - 1]).unwrap();
    let protein = rna.protein().unwrap();
    for p in &protein[..protein.len() - 1] {
        print!("{:?}", p);
    }
}

pub fn grph() {
    let filename = "data/grph.txt";
    let string = std::fs::read(filename).unwrap();
    let strands = parse_fasta::<DNABase>(&string).unwrap();
    for (i, left) in strands.iter().enumerate() {
        for (j, right) in strands.iter().enumerate() {
            if i == j {
                continue;
            }
            if left.strand[left.strand.len() - 3..] == right.strand[..3] {
                println!("{} {}", left.name, right.name);
            }
        }
    }
}

pub fn long() {
    let filename = "data/long.txt";
    let string = std::fs::read(filename).unwrap();
    let strands = parse_fasta::<DNABase>(&string).unwrap();

    let superstring = strands.superstring().unwrap();
    for c in superstring {
        print!("{:?}", c);
    }
}

pub fn fib() {
    let filename = "data/fib.txt";
    let string = std::fs::read(filename).unwrap();
    let mut string = String::from_utf8(string).unwrap();
    string = string.trim().to_string();
    let tokens = string.split(' ').collect::<Vec<_>>();

    let generations: usize = tokens[0].parse().unwrap();
    let child_per_generation: usize = tokens[1].parse().unwrap();

    let mut a = 1;
    let mut b = 1;

    for _ in 0..generations - 2 {
        let tmp = child_per_generation * a + b;
        a = b;
        b = tmp;
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
    let prot = parse::<AA>(&string[..string.len() - 1]).unwrap();

    let mut p = 1;
    for aa in prot.strand {
        let n = aa.encodings().len();
        p *= n;
        p %= 1_000_000;
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
    println!();
    for i in 0..4 {
        print!("{:?}: ", DNABase::bases()[i]);
        for p in profile.iter() {
            print!("{:?} ", p[i]);
        }
        println!();
    }
}

fn p_dominate(k: usize, m: usize, n: usize) -> f32 {
    let total = m + n + k;
    let total_ = total - 1;
    let p_n = n as f32 / total as f32;
    let p_n_m = m as f32 / total_ as f32;
    let p_n_n = (n - 1) as f32 / total_ as f32;
    let p_m = m as f32 / total as f32;
    let p_m_m = (m - 1) as f32 / total_ as f32;
    let p_m_n = n as f32 / total_ as f32;

    let p = p_n * p_n_n + p_n * p_n_m * 0.5 + 0.5 * p_m * p_m_n + 0.25 * p_m * p_m_m;
    1.0 - p
}

pub fn iprb() {
    let filename = "data/iprb.txt";
    let string = std::fs::read(filename).unwrap();
    let mut string = String::from_utf8(string).unwrap();
    string = string.trim().to_string();
    let tokens = string.split(' ').collect::<Vec<_>>();

    let k: usize = tokens[0].parse().unwrap();
    let m: usize = tokens[1].parse().unwrap();
    let n: usize = tokens[2].parse().unwrap();

    println!("{:?}", p_dominate(k, m, n));
}
pub fn iev() {
    let filename = "data/iev.txt";
    let string = std::fs::read(filename).unwrap();
    let mut string = String::from_utf8(string).unwrap();
    string = string.trim().to_string();
    let numbers: Vec<f32> = string
        .split(' ')
        .map(|c| c.parse().unwrap())
        .collect::<Vec<_>>();

    let out = 2.0 * numbers[0]
        + 2.0 * numbers[1]
        + 2.0 * numbers[2]
        + 2.0 * 0.75 * numbers[3]
        + numbers[4];
    println!("{:?}", out);
}

pub fn lcsm() {
    let filename = "data/lcsm.txt";
    let string = std::fs::read(filename).unwrap();
    let strands = parse_fasta::<DNABase>(&string).unwrap();

    let longest = strands.longest_common_substrand();
    for base in longest {
        print!("{:?}", base);
    }
}

fn c(n: usize, k: usize) -> f64 {
    if k > n {
        panic!("Invalid cnk parameters");
    }
    if n == 0 {
        return 1.0;
    }
    (n - k + 1..=n).map(|i| i as f64).product::<f64>() / factorial(k)
}

pub fn factorial(num: usize) -> f64 {
    match num {
        0 => 1.0,
        num => (1..=num).map(|i| i as f64).product(),
    }
}

pub fn lia() {
    let filename = "data/lia.txt";
    let string = std::fs::read(filename).unwrap();
    let mut string = String::from_utf8(string).unwrap();
    string = string.trim().to_string();
    let tokens = string.split(' ').collect::<Vec<_>>();

    let k: usize = tokens[0].parse().unwrap();
    let n: usize = tokens[1].parse().unwrap();

    let mut tom = Individual::<2>::from([DD, DD]);
    let anonymous = Individual::<2>::from([Dd, Dd]);

    for _ in 0..k {
        tom = coupling(&tom, &anonymous);
    }

    let p = tom.probability([Dd, Dd]) as f64;
    let n_individuals = 2usize.pow(k as u32);

    let mut final_prob = 0.0;
    for i in n..=n_individuals {
        final_prob +=
            c(n_individuals, i) * p.pow(i as f64) * (1.0 - p).pow((n_individuals - i) as f64);
    }
    println!("{:?}", final_prob);
}

pub fn mprt() {
    let filename = "data/mprt.txt";
    let string = std::fs::read(filename).unwrap();
    let mut string = String::from_utf8(string).unwrap();
    string = string.trim().to_string();
    let ids = string.split('\n').collect::<Vec<_>>();

    for id in ids {
        let url = format!("http://www.uniprot.org/uniprot/{}.fasta", id);

        let resp = ureq::get(&url).call().unwrap();
        assert!(resp.has("Content-Length"));
        let len = resp
            .header("Content-Length")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap();

        let mut bytes: Vec<u8> = Vec::with_capacity(len);
        resp.into_reader()
            .take(100_000_000)
            .read_to_end(&mut bytes)
            .unwrap();

        assert_eq!(bytes.len(), len);

        let strands = parse_fasta::<AA>(&bytes).unwrap();

        for strand in strands.iter() {
            let mut exists = false;
            let mut partial_matches = vec![];
            for (i, aa) in strand.strand.iter().enumerate() {
                partial_matches = partial_matches
                    .iter()
                    .filter_map(|ii| {
                        let d = i - ii;
                        match (d, aa) {
                            (1, AA::P) => None,
                            (1, _) => Some(*ii),
                            (2, AA::S | AA::T) => Some(*ii),
                            (3, AA::P) => None,
                            (3, _) => {
                                if !exists {
                                    exists = true;
                                    println!("{}", id);
                                }
                                print!("{} ", ii + 1);
                                None
                            }
                            (_, _) => None,
                        }
                    })
                    .collect();

                if *aa == AA::N {
                    partial_matches.push(i);
                }
            }
            if exists {
                println!();
            }
        }
    }
}
pub fn prtm() {
    let filename = "data/prtm.txt";
    let string = std::fs::read(filename).unwrap();
    let tokens = string.split(|c| c == &b'\n').collect::<Vec<_>>();
    let prot = parse::<AA>(tokens[0]).unwrap();

    // let _water = 18.01056;
    let weight: f64 = prot.strand.iter().map(|aa| aa.weight()).sum::<f64>();
    println!("{:.10}", weight);
}

pub fn orf() {
    let filename = "data/orf.txt";
    let string = std::fs::read(filename).unwrap();
    let dna = parse_fasta::<DNABase>(&string).unwrap();
    for protein in dna.strands[0].candidate_proteins().strands {
        println!("{}", protein);
    }
}

pub fn perm() {
    let filename = "data/perm.txt";
    let string = std::fs::read(filename).unwrap();
    let string = String::from_utf8(string).unwrap();
    let string = string.trim();
    let n: usize = string.parse().unwrap();

    println!("{:?}", factorial(n) as usize);
    let perms = (0..n).permutations(n);
    for perm in perms {
        print!("{:?}", perm[0] + 1);
        for a in &perm[1..] {
            print!(" {:?}", a + 1);
        }
        println!();
    }
}

pub fn revp() {
    let filename = "data/revp.txt";
    let string = std::fs::read(filename).unwrap();
    let strands = parse_fasta::<DNABase>(&string).unwrap();

    for (start, length) in strands.strands[0].find_palindromes(4, 12) {
        println!("{:?} {:?}", start + 1, length);
    }
}
pub fn splc() {
    let filename = "data/splc.txt";
    let string = std::fs::read(filename).unwrap();
    let strands = parse_fasta::<DNABase>(&string).unwrap();

    let proteins = strands.splice();
    for aa in &proteins[0][..proteins[0].len() - 1] {
        print!("{:?}", aa);
    }
    println!();
}
pub fn lexf() {
    let filename = "data/lexf.txt";
    let string = std::fs::read(filename).unwrap();
    let tokens = string.split(|c| c == &b'\n').collect::<Vec<_>>();
    let letters = tokens[0]
        .split(|c| c == &b' ')
        .map(|le| le[0])
        .collect::<Vec<_>>();
    let n: usize = String::from_utf8(tokens[1].to_vec())
        .unwrap()
        .parse()
        .unwrap();
    for permutation in (0..n)
        .map(|_| letters.iter().cloned())
        .multi_cartesian_product()
    {
        println!("{}", String::from_utf8(permutation).unwrap());
    }
}

fn increasing_subsequence<const INCREASING: bool>(sequence: &[usize]) -> Vec<usize> {
    let n = sequence.len();
    let mut predecessors = vec![0; n];
    let mut partial = vec![0; n + 1];
    let mut l = 0;

    for i in 0..n {
        // Binary search for the largest positive j ≤ L
        // such that X[M[j]] < X[i]
        let mut lo = 1;
        let mut hi = l + 1;
        let mut mid;
        while lo < hi {
            mid = lo + (hi - lo) / 2;

            match (INCREASING, sequence[partial[mid]].cmp(&sequence[i])) {
                (true, std::cmp::Ordering::Less) => {
                    lo = mid + 1;
                }
                (false, std::cmp::Ordering::Greater) => {
                    lo = mid + 1;
                }
                _ => {
                    hi = mid;
                }
            }
        }

        // After searching, lo is 1 greater than the
        // length of the longest prefix of X[i]
        let new_l = lo;

        // The predecessor of X[i] is the last index of
        // the subsequence of length newL-1
        predecessors[i] = partial[new_l - 1];
        partial[new_l] = i;

        if new_l > l {
            // If we found a subsequence longer than any we've
            // found yet, update L
            l = new_l;
        }
    }

    // Reconstruct the longest increasing subsequence
    let mut result = vec![0; l];
    let mut k = partial[l];
    for i in 0..l {
        let index = l - i - 1;
        result[index] = sequence[k];
        k = predecessors[k];
    }

    result
}

pub fn lgis() {
    let filename = "data/lgis.txt";
    let string = std::fs::read(filename).unwrap();
    let lines = string.split(|c| c == &b'\n').collect::<Vec<_>>();
    let _n: usize = String::from_utf8(lines[0].to_vec())
        .unwrap()
        .parse()
        .unwrap();
    let sequence: Vec<usize> = lines[1]
        .split(|c| c == &b' ')
        .map(|c| String::from_utf8(c.to_vec()).unwrap().parse().unwrap())
        .collect::<Vec<_>>();

    let increasing = increasing_subsequence::<true>(&sequence);
    let decreasing = increasing_subsequence::<false>(&sequence);

    for i in increasing {
        print!("{:?} ", i);
    }
    println!();
    for i in decreasing {
        print!("{:?} ", i);
    }
    println!();
}

fn main() {
    //dna();
    //dna_to_rna();
    //complement();
    //gc();
    //subs();
    //prot();
    //grph();
    //long();
    //fib();
    //fibd();
    //mrna();
    //cons();
    //iprb();
    //iev();
    //lia();
    //lcsm();
    //mprt();
    //orf();
    //perm();
    //prtm();
    //revp();
    //splc();
    //lexf();
    let start = Instant::now();
    lgis();
    println!("Elapsed {:?}", start.elapsed());
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0), 1.0);
        assert_eq!(factorial(1), 1.0);
        assert_eq!(factorial(2), 2.0);
        assert_eq!(factorial(3), 6.0);
        assert_eq!(factorial(4), 24.0);
    }

    #[test]
    fn test_cnk() {
        assert_eq!(c(0, 0), 1.0);
        assert_eq!(c(1, 0), 1.0);
        assert_eq!(c(1, 1), 1.0);
        assert_eq!(c(2, 0), 1.0);
        assert_eq!(c(2, 1), 2.0);
        assert_eq!(c(2, 2), 1.0);
        assert_eq!(c(3, 0), 1.0);
        assert_eq!(c(3, 1), 3.0);
        assert_eq!(c(3, 2), 3.0);
        assert_eq!(c(3, 3), 1.0);
    }
}
