use itertools::Itertools;
use std::collections::VecDeque;
use std::collections::{HashMap, HashSet};
use thiserror::Error;

pub mod individuals;

pub trait Base:
    std::hash::Hash + Eq + Clone + Copy + std::convert::TryFrom<u8> + std::fmt::Debug
{
    fn terminal() -> Option<Self> {
        None
    }
}

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("IoError")]
    IoError(#[from] std::io::Error),
    #[error("utf-8 error")]
    Utf8Error(#[from] std::string::FromUtf8Error),
    #[error("UnexpectedData in FAST")]
    UnexpectedData,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DNABase {
    A,
    C,
    G,
    T,
}
impl std::convert::TryFrom<u8> for DNABase {
    type Error = ParseError;

    fn try_from(c: u8) -> Result<Self, ParseError> {
        match c {
            b'A' => Ok(DNABase::A),
            b'C' => Ok(DNABase::C),
            b'G' => Ok(DNABase::G),
            b'T' => Ok(DNABase::T),
            _ => Err(ParseError::UnexpectedData),
        }
    }
}
impl Base for DNABase {}

impl DNABase {
    pub fn index(&self) -> usize {
        match self {
            DNABase::A => 0,
            DNABase::C => 1,
            DNABase::G => 2,
            DNABase::T => 3,
        }
    }

    pub fn bases() -> [DNABase; 4] {
        [DNABase::A, DNABase::C, DNABase::G, DNABase::T]
    }

    pub fn complement(&self) -> DNABase {
        match self {
            DNABase::A => DNABase::T,
            DNABase::T => DNABase::A,
            DNABase::C => DNABase::G,
            DNABase::G => DNABase::C,
        }
    }

    pub fn rna(&self) -> RNABase {
        match self {
            DNABase::A => RNABase::A,
            DNABase::T => RNABase::U,
            DNABase::C => RNABase::C,
            DNABase::G => RNABase::G,
        }
    }

    /// 0 -> purine
    /// 1 -> pyrimidine
    pub fn type_(&self) -> usize {
        match self {
            DNABase::A | DNABase::G => 0,
            DNABase::C | DNABase::T => 1,
        }
    }
}

pub struct StrandIterator<'a, 'b, T: Base> {
    pattern_strand: &'a [T],
    data_strand: &'b [T],
    matches: Vec<(usize, usize)>,
    index: usize,
}

fn hamming_distance<T: Base>(left: &[T], right: &[T]) -> usize {
    left.iter()
        .zip(right.iter())
        .map(|(l, r)| if l == r { 0 } else { 1 })
        .sum()
}

impl<'a, 'b, T: Base> Iterator for StrandIterator<'a, 'b, T> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        while self.index < self.data_strand.len() {
            let c = self.data_strand[self.index];

            let mut matched = None;
            self.matches = self
                .matches
                .iter()
                .filter_map(|(data_index, pattern_index)| {
                    if self.pattern_strand.get(pattern_index + 1) == Some(&c) {
                        // Matching
                        if self.pattern_strand.len() == pattern_index + 2 {
                            // Final match
                            matched = Some(*data_index);
                        }
                        Some((*data_index, pattern_index + 1))
                    } else {
                        None
                    }
                })
                .collect();
            if c == self.pattern_strand[0] {
                if self.pattern_strand.len() == 1 {
                    // Final match
                    assert!(self.matches.is_empty());
                    matched = Some(self.index);
                } else {
                    self.matches.push((self.index, 0));
                }
            }
            self.index += 1;
            if matched.is_some() {
                return matched;
            }
        }
        None
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RNABase {
    A,
    C,
    G,
    U,
}
impl std::convert::TryFrom<u8> for RNABase {
    type Error = ParseError;

    fn try_from(c: u8) -> Result<Self, ParseError> {
        match c {
            b'A' => Ok(RNABase::A),
            b'C' => Ok(RNABase::C),
            b'G' => Ok(RNABase::G),
            b'U' => Ok(RNABase::U),
            _ => Err(ParseError::UnexpectedData),
        }
    }
}
impl Base for RNABase {}

impl RNABase {
    pub fn bases() -> [RNABase; 4] {
        [RNABase::A, RNABase::C, RNABase::G, RNABase::U]
    }
}

type Codon = (RNABase, RNABase, RNABase);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AA {
    F,
    L,
    S,
    Y,
    Stop,
    C,
    W,
    P,
    H,
    Q,
    R,
    I,
    M,
    T,
    N,
    K,
    V,
    A,
    D,
    E,
    G,
}
impl std::convert::TryFrom<u8> for AA {
    type Error = ParseError;

    fn try_from(c: u8) -> Result<Self, ParseError> {
        match c {
            b'F' => Ok(AA::F),
            b'L' => Ok(AA::L),
            b'S' => Ok(AA::S),
            b'Y' => Ok(AA::Y),
            b'C' => Ok(AA::C),
            b'W' => Ok(AA::W),
            b'P' => Ok(AA::P),
            b'H' => Ok(AA::H),
            b'Q' => Ok(AA::Q),
            b'R' => Ok(AA::R),
            b'I' => Ok(AA::I),
            b'M' => Ok(AA::M),
            b'T' => Ok(AA::T),
            b'N' => Ok(AA::N),
            b'K' => Ok(AA::K),
            b'V' => Ok(AA::V),
            b'A' => Ok(AA::A),
            b'D' => Ok(AA::D),
            b'E' => Ok(AA::E),
            b'G' => Ok(AA::G),
            _ => Err(ParseError::UnexpectedData),
        }
    }
}
impl Base for AA {
    fn terminal() -> Option<Self> {
        Some(AA::Stop)
    }
}

impl AA {
    pub fn weight(&self) -> f64 {
        match self {
            AA::A => 71.03711,
            AA::C => 103.00919,
            AA::D => 115.02694,
            AA::E => 129.04259,
            AA::F => 147.06841,
            AA::G => 57.02146,
            AA::H => 137.05891,
            AA::I => 113.08406,
            AA::K => 128.09496,
            AA::L => 113.08406,
            AA::M => 131.04049,
            AA::N => 114.04293,
            AA::P => 97.05276,
            AA::Q => 128.05858,
            AA::R => 156.10111,
            AA::S => 87.03203,
            AA::T => 101.04768,
            AA::V => 99.06841,
            AA::W => 186.07931,
            AA::Y => 163.06333,
            AA::Stop => 0.0,
        }
    }
    pub fn encodings(&self) -> Vec<Codon> {
        match self {
            AA::F => vec![
                (RNABase::U, RNABase::U, RNABase::U),
                (RNABase::U, RNABase::U, RNABase::C),
            ],
            AA::L => vec![
                (RNABase::U, RNABase::U, RNABase::A),
                (RNABase::U, RNABase::U, RNABase::G),
                (RNABase::C, RNABase::U, RNABase::U),
                (RNABase::C, RNABase::U, RNABase::C),
                (RNABase::C, RNABase::U, RNABase::A),
                (RNABase::C, RNABase::U, RNABase::G),
            ],
            AA::S => vec![
                (RNABase::U, RNABase::C, RNABase::U),
                (RNABase::U, RNABase::C, RNABase::C),
                (RNABase::U, RNABase::C, RNABase::A),
                (RNABase::U, RNABase::C, RNABase::G),
                (RNABase::A, RNABase::G, RNABase::U),
                (RNABase::A, RNABase::G, RNABase::C),
            ],
            AA::Y => vec![
                (RNABase::U, RNABase::A, RNABase::U),
                (RNABase::U, RNABase::A, RNABase::C),
            ],
            AA::Stop => vec![
                (RNABase::U, RNABase::A, RNABase::A),
                (RNABase::U, RNABase::A, RNABase::G),
                (RNABase::U, RNABase::G, RNABase::A),
            ],
            AA::C => vec![
                (RNABase::U, RNABase::G, RNABase::U),
                (RNABase::U, RNABase::G, RNABase::C),
            ],
            AA::W => vec![(RNABase::U, RNABase::G, RNABase::G)],
            AA::P => vec![
                (RNABase::C, RNABase::C, RNABase::U),
                (RNABase::C, RNABase::C, RNABase::C),
                (RNABase::C, RNABase::C, RNABase::A),
                (RNABase::C, RNABase::C, RNABase::G),
            ],
            AA::H => vec![
                (RNABase::C, RNABase::A, RNABase::U),
                (RNABase::C, RNABase::A, RNABase::C),
            ],
            AA::Q => vec![
                (RNABase::C, RNABase::A, RNABase::A),
                (RNABase::C, RNABase::A, RNABase::G),
            ],
            AA::R => vec![
                (RNABase::C, RNABase::G, RNABase::U),
                (RNABase::C, RNABase::G, RNABase::C),
                (RNABase::C, RNABase::G, RNABase::A),
                (RNABase::C, RNABase::G, RNABase::G),
                (RNABase::A, RNABase::G, RNABase::A),
                (RNABase::A, RNABase::G, RNABase::G),
            ],
            AA::I => vec![
                (RNABase::A, RNABase::U, RNABase::U),
                (RNABase::A, RNABase::U, RNABase::C),
                (RNABase::A, RNABase::U, RNABase::A),
            ],
            AA::M => vec![(RNABase::A, RNABase::U, RNABase::G)],
            AA::T => vec![
                (RNABase::A, RNABase::C, RNABase::U),
                (RNABase::A, RNABase::C, RNABase::C),
                (RNABase::A, RNABase::C, RNABase::A),
                (RNABase::A, RNABase::C, RNABase::G),
            ],
            AA::N => vec![
                (RNABase::A, RNABase::A, RNABase::U),
                (RNABase::A, RNABase::A, RNABase::C),
            ],
            AA::K => vec![
                (RNABase::A, RNABase::A, RNABase::A),
                (RNABase::A, RNABase::A, RNABase::G),
            ],
            AA::V => vec![
                (RNABase::G, RNABase::U, RNABase::U),
                (RNABase::G, RNABase::U, RNABase::C),
                (RNABase::G, RNABase::U, RNABase::A),
                (RNABase::G, RNABase::U, RNABase::G),
            ],
            AA::A => vec![
                (RNABase::G, RNABase::C, RNABase::U),
                (RNABase::G, RNABase::C, RNABase::C),
                (RNABase::G, RNABase::C, RNABase::A),
                (RNABase::G, RNABase::C, RNABase::G),
            ],
            AA::D => vec![
                (RNABase::G, RNABase::A, RNABase::U),
                (RNABase::G, RNABase::A, RNABase::C),
            ],
            AA::E => vec![
                (RNABase::G, RNABase::A, RNABase::A),
                (RNABase::G, RNABase::A, RNABase::G),
            ],
            AA::G => vec![
                (RNABase::G, RNABase::G, RNABase::U),
                (RNABase::G, RNABase::G, RNABase::C),
                (RNABase::G, RNABase::G, RNABase::A),
                (RNABase::G, RNABase::G, RNABase::G),
            ],
        }
    }
}

#[derive(Debug)]
pub enum AAError {
    UnexpectedRNALength,
}

impl Strand<RNABase> {
    pub fn protein(&self) -> Result<Vec<AA>, AAError> {
        let protein: Result<Vec<_>, _> = self
            .strand
            .chunks(3)
            .map(|chunk| match chunk {
                [RNABase::U, RNABase::U, RNABase::U] => Ok(AA::F),
                [RNABase::U, RNABase::U, RNABase::C] => Ok(AA::F),
                [RNABase::U, RNABase::U, RNABase::A] => Ok(AA::L),
                [RNABase::U, RNABase::U, RNABase::G] => Ok(AA::L),
                [RNABase::U, RNABase::C, RNABase::U] => Ok(AA::S),
                [RNABase::U, RNABase::C, RNABase::C] => Ok(AA::S),
                [RNABase::U, RNABase::C, RNABase::A] => Ok(AA::S),
                [RNABase::U, RNABase::C, RNABase::G] => Ok(AA::S),
                [RNABase::U, RNABase::A, RNABase::U] => Ok(AA::Y),
                [RNABase::U, RNABase::A, RNABase::C] => Ok(AA::Y),
                [RNABase::U, RNABase::A, RNABase::A] => Ok(AA::Stop),
                [RNABase::U, RNABase::A, RNABase::G] => Ok(AA::Stop),
                [RNABase::U, RNABase::G, RNABase::U] => Ok(AA::C),
                [RNABase::U, RNABase::G, RNABase::C] => Ok(AA::C),
                [RNABase::U, RNABase::G, RNABase::A] => Ok(AA::Stop),
                [RNABase::U, RNABase::G, RNABase::G] => Ok(AA::W),
                [RNABase::C, RNABase::U, RNABase::U] => Ok(AA::L),
                [RNABase::C, RNABase::U, RNABase::C] => Ok(AA::L),
                [RNABase::C, RNABase::U, RNABase::A] => Ok(AA::L),
                [RNABase::C, RNABase::U, RNABase::G] => Ok(AA::L),
                [RNABase::C, RNABase::C, RNABase::U] => Ok(AA::P),
                [RNABase::C, RNABase::C, RNABase::C] => Ok(AA::P),
                [RNABase::C, RNABase::C, RNABase::A] => Ok(AA::P),
                [RNABase::C, RNABase::C, RNABase::G] => Ok(AA::P),
                [RNABase::C, RNABase::A, RNABase::U] => Ok(AA::H),
                [RNABase::C, RNABase::A, RNABase::C] => Ok(AA::H),
                [RNABase::C, RNABase::A, RNABase::A] => Ok(AA::Q),
                [RNABase::C, RNABase::A, RNABase::G] => Ok(AA::Q),
                [RNABase::C, RNABase::G, RNABase::U] => Ok(AA::R),
                [RNABase::C, RNABase::G, RNABase::C] => Ok(AA::R),
                [RNABase::C, RNABase::G, RNABase::A] => Ok(AA::R),
                [RNABase::C, RNABase::G, RNABase::G] => Ok(AA::R),
                [RNABase::A, RNABase::U, RNABase::U] => Ok(AA::I),
                [RNABase::A, RNABase::U, RNABase::C] => Ok(AA::I),
                [RNABase::A, RNABase::U, RNABase::A] => Ok(AA::I),
                [RNABase::A, RNABase::U, RNABase::G] => Ok(AA::M),
                [RNABase::A, RNABase::C, RNABase::U] => Ok(AA::T),
                [RNABase::A, RNABase::C, RNABase::C] => Ok(AA::T),
                [RNABase::A, RNABase::C, RNABase::A] => Ok(AA::T),
                [RNABase::A, RNABase::C, RNABase::G] => Ok(AA::T),
                [RNABase::A, RNABase::A, RNABase::U] => Ok(AA::N),
                [RNABase::A, RNABase::A, RNABase::C] => Ok(AA::N),
                [RNABase::A, RNABase::A, RNABase::A] => Ok(AA::K),
                [RNABase::A, RNABase::A, RNABase::G] => Ok(AA::K),
                [RNABase::A, RNABase::G, RNABase::U] => Ok(AA::S),
                [RNABase::A, RNABase::G, RNABase::C] => Ok(AA::S),
                [RNABase::A, RNABase::G, RNABase::A] => Ok(AA::R),
                [RNABase::A, RNABase::G, RNABase::G] => Ok(AA::R),
                [RNABase::G, RNABase::U, RNABase::U] => Ok(AA::V),
                [RNABase::G, RNABase::U, RNABase::C] => Ok(AA::V),
                [RNABase::G, RNABase::U, RNABase::A] => Ok(AA::V),
                [RNABase::G, RNABase::U, RNABase::G] => Ok(AA::V),
                [RNABase::G, RNABase::C, RNABase::U] => Ok(AA::A),
                [RNABase::G, RNABase::C, RNABase::C] => Ok(AA::A),
                [RNABase::G, RNABase::C, RNABase::A] => Ok(AA::A),
                [RNABase::G, RNABase::C, RNABase::G] => Ok(AA::A),
                [RNABase::G, RNABase::A, RNABase::U] => Ok(AA::D),
                [RNABase::G, RNABase::A, RNABase::C] => Ok(AA::D),
                [RNABase::G, RNABase::A, RNABase::A] => Ok(AA::E),
                [RNABase::G, RNABase::A, RNABase::G] => Ok(AA::E),
                [RNABase::G, RNABase::G, RNABase::U] => Ok(AA::G),
                [RNABase::G, RNABase::G, RNABase::C] => Ok(AA::G),
                [RNABase::G, RNABase::G, RNABase::A] => Ok(AA::G),
                [RNABase::G, RNABase::G, RNABase::G] => Ok(AA::G),
                _ => Err(AAError::UnexpectedRNALength),
            })
            .collect();
        protein
    }
}

#[derive(Debug, Hash, Eq, PartialEq)]
pub struct Strand<T> {
    pub name: String,
    pub strand: Vec<T>,
}

impl<T: Base> Strand<T> {
    pub fn find<'a, 'b>(&'a self, data_strand: &'b [T]) -> StrandIterator<'a, 'b, T> {
        StrandIterator {
            pattern_strand: &self.strand,
            data_strand,
            matches: vec![],
            index: 0,
        }
    }
}

impl<T: Base> std::fmt::Display for Strand<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let n = self.strand.len();
        let last = if T::terminal().is_some() { n - 1 } else { n };
        for base in &self.strand[..last] {
            write!(f, "{:?}", base)?;
        }
        Ok(())
    }
}

fn candidate_proteins(rna: &[RNABase]) -> Vec<Vec<AA>> {
    let mut patterns = HashMap::new();
    for encoding in AA::M.encodings() {
        patterns.insert(encoding, AA::M);
    }
    for encoding in AA::Stop.encodings() {
        patterns.insert(encoding, AA::Stop);
    }
    let mut outputs = vec![];
    let mut starts = [vec![], vec![], vec![]];

    for (i, window) in rna.windows(3).enumerate() {
        let bucket = i % 3;
        if let [a, b, c] = window {
            if let Some(aa) = patterns.get(&(*a, *b, *c)) {
                match aa {
                    AA::Stop => {
                        let current_starts = &starts[bucket];
                        for st in current_starts {
                            let protein = Strand {
                                name: "".to_string(),
                                strand: rna[*st..=i + 2].to_vec(),
                            }
                            .protein()
                            .unwrap();
                            outputs.push(protein);
                        }
                        starts[bucket].clear();
                    }
                    AA::M => {
                        starts[bucket].push(i);
                    }

                    _ => unreachable!(),
                }
            }
        }
    }
    outputs
}

impl Strand<DNABase> {
    pub fn gc(&self) -> f32 {
        let sum = self
            .strand
            .iter()
            .filter(|c| matches!(c, DNABase::C | DNABase::G))
            .count();
        sum as f32 / self.strand.len() as f32
    }

    pub fn complement(&self) -> Self {
        let out = &self
            .strand
            .iter()
            .map(|c| c.complement())
            .rev()
            .collect::<Vec<_>>();
        Strand {
            name: "".to_string(),
            strand: out.to_vec(),
        }
    }

    pub fn rna(&self) -> Strand<RNABase> {
        let mut strand = Vec::with_capacity(self.strand.len());
        for dbase in &self.strand {
            strand.push(dbase.rna());
        }
        Strand {
            name: self.name.clone(),
            strand,
        }
    }
    pub fn candidate_proteins(&self) -> Strands<AA> {
        let mut left = candidate_proteins(&self.rna().strand);
        left.extend(candidate_proteins(&self.complement().rna().strand));
        let deduplicate: HashSet<_> = left.into_iter().collect();
        let strands = deduplicate.into_iter();
        Strands {
            strands: strands
                .map(|prot| Strand {
                    name: "".to_string(),
                    strand: prot,
                })
                .collect(),
        }
    }

    pub fn find_palindromes(&self, n: usize, m: usize) -> Vec<(usize, usize)> {
        let mut positions = vec![];
        let mut trie = STrie::new();
        let bases = DNABase::bases();
        for i in (n..=m).step_by(2) {
            let half = i / 2;
            for first in (0..half).map(|_| bases.iter()).multi_cartesian_product() {
                let complement = first
                    .iter()
                    .map(|base| base.complement())
                    .rev()
                    .collect::<Vec<_>>();

                let mut full = Vec::with_capacity(m);
                full.extend(first);
                full.extend(complement);
                trie.add(&full, i);
            }
        }
        for (pos, i) in trie.matches(&self.strand) {
            positions.push((pos, i));
        }
        // positions.sort();
        positions
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct PNode<T: Base> {
    subnodes: HashMap<T, PNode<T>>,
    ids: HashSet<usize>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct SNode<T: Base> {
    subnodes: HashMap<T, SNode<T>>,
    leafs: HashSet<usize>,
}

impl<T: Base> PNode<T> {
    fn new() -> Self {
        Self {
            subnodes: HashMap::new(),
            ids: HashSet::new(),
        }
    }
}
impl<T: Base> SNode<T> {
    fn new() -> Self {
        Self {
            subnodes: HashMap::new(),
            leafs: HashSet::new(),
        }
    }
}

#[derive(Debug, Eq, PartialEq)]
struct PTrie<T: Base> {
    root: PNode<T>,
}

impl<T: Base> PTrie<T> {
    fn new() -> Self {
        Self { root: PNode::new() }
    }

    fn add(&mut self, data: &[T], id: usize) {
        let mut node = &mut self.root;
        node.ids.insert(id);

        for base in data {
            let subnode = node.subnodes.entry(*base).or_insert_with(PNode::new);
            node = subnode;
            node.ids.insert(id);
        }
        // node.leafs.insert(id);
    }

    fn from(strands: &[Strand<T>]) -> Self {
        let mut trie = PTrie::new();
        for (i, strand) in strands.iter().enumerate() {
            trie.add(&strand.strand, i);
        }
        trie
    }

    fn find_longest(&self, pattern: &[T]) -> Vec<T> {
        let mut node = &self.root;
        for (i, base) in pattern.iter().enumerate() {
            if let Some(subnode) = node.subnodes.get(base) {
                node = subnode;
            } else {
                return pattern[..i].to_vec();
            }
        }
        pattern.to_vec()
    }
}

#[derive(Debug, Eq, PartialEq)]
struct STrie<T: Base> {
    root: SNode<T>,
}

impl<T: Base> STrie<T> {
    fn new() -> Self {
        Self { root: SNode::new() }
    }

    fn add(&mut self, data: &[T], id: usize) {
        let mut node = &mut self.root;

        for base in data {
            let subnode = node.subnodes.entry(*base).or_insert_with(SNode::new);
            node = subnode;
        }
        node.leafs.insert(id);
    }

    fn matches(&self, long_string: &[T]) -> Vec<(usize, usize)> {
        let mut matches = vec![];
        let mut nodes: Vec<(usize, &SNode<T>)> = vec![];
        for (position, c) in long_string.iter().enumerate() {
            nodes.push((position, &self.root));
            nodes = nodes
                .iter()
                .filter_map(|(p, n)| n.subnodes.get(c).map(|subnode| (*p, subnode)))
                .collect();

            for (p, n) in &nodes {
                for leaf in &n.leafs {
                    matches.push((*p, *leaf));
                }
            }
        }
        matches
    }
}

#[derive(Debug)]
pub struct Strands<T> {
    pub strands: Vec<Strand<T>>,
}

#[derive(Debug)]
pub enum StrandError {
    NoStrands,
    LengthsAreDifferent,
}

impl<T: Base> Default for Strands<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Base> Strands<T> {
    pub fn new() -> Self {
        Self { strands: vec![] }
    }

    pub fn from(strands: Vec<Strand<T>>) -> Self {
        Self { strands }
    }

    pub fn iter(&self) -> std::slice::Iter<Strand<T>> {
        self.strands.iter()
    }

    fn common_substrand(&self, a: &[T], b: &[T]) -> Vec<T> {
        // Small optimization, where we build the trie
        // only for the smaller common substring.
        let (left, right) = if a.len() > b.len() { (b, a) } else { (a, b) };

        // Build a prefix Trie
        let mut trie = PTrie::new();
        for i in 0..left.len() {
            trie.add(&left[i..], 0);
        }
        // Find the longest match within right
        let mut longest = vec![];
        for i in 0..right.len() {
            let current = trie.find_longest(&right[i..]);
            if current.len() > longest.len() {
                longest = current;
            }
        }
        longest
    }

    pub fn longest_common_substrand(&self) -> Vec<T> {
        let mut common = self.strands[0].strand.clone();
        for strand in &self.strands {
            common = self.common_substrand(&common, &strand.strand);
        }
        common
    }

    pub fn superstring(&self) -> Option<Vec<T>> {
        let trie = PTrie::from(&self.strands);

        let mut all_states = vec![];
        for strand in &self.strands {
            let mut states = vec![];
            for data in &strand.strand {
                states.push((&trie.root, 0));
                states = states
                    .into_iter()
                    .filter_map(|(state, i)| {
                        if let Some(subnode) = state.subnodes.get(data) {
                            let newstate = subnode;
                            Some((newstate, i + 1))
                        } else {
                            None
                        }
                    })
                    .collect();
            }
            all_states.push(states);
        }

        let mut match_: Option<Vec<T>> = None;
        for (i, strand) in self.strands.iter().enumerate() {
            let mut current_data = strand.strand.clone();
            let mut index = i;
            let mut visited = HashSet::new();
            visited.insert(i);
            loop {
                let states = &all_states[index];
                if let Some((state, match_size)) = states.get(1) {
                    if let Some(leaf_index) = state
                        .ids
                        .iter()
                        .find(|leaf_index| !visited.contains(leaf_index))
                    {
                        let new_strand = &self.strands[*leaf_index].strand;
                        current_data.extend(new_strand.iter().skip(*match_size));
                        visited.insert(*leaf_index);
                        index = *leaf_index;
                    } else {
                        break;
                    }
                } else {
                    // finalize match
                    break;
                }
            }
            if visited.len() == self.strands.len() {
                if let Some(current_match) = &match_ {
                    if current_match.len() > current_data.len() {
                        match_ = Some(current_data);
                    }
                } else {
                    match_ = Some(current_data);
                }
            }
        }
        match_
    }
}

impl Strands<DNABase> {
    pub fn error_correction(&self) -> Vec<(&Strand<DNABase>, &Strand<DNABase>, bool)> {
        let mut counter: HashMap<&Vec<DNABase>, HashSet<usize>> = HashMap::new();
        let complements: Vec<_> = self.strands.iter().map(|s| s.complement()).collect();
        let mut gcs_contents = vec![];
        for (i, (strand, complement)) in self.strands.iter().zip(&complements).enumerate() {
            counter
                .entry(&complement.strand)
                .or_insert(HashSet::new())
                .insert(i);
            counter
                .entry(&strand.strand)
                .or_insert(HashSet::new())
                .insert(i);
            gcs_contents.push((strand.gc(), i));
        }

        gcs_contents.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut originals = vec![0; gcs_contents.len()];
        for (i, (_, id)) in gcs_contents.iter().enumerate() {
            originals[*id] = i;
        }

        let mut visited = HashSet::new();
        let mut corrections = vec![];
        for (_, set) in counter.iter() {
            if set.len() == 1 {
                let id = *set.iter().next().unwrap();
                if visited.contains(&id) {
                    continue;
                }
                let strand = &self.strands[id].strand;
                let mut i = 1;
                let gc_index = originals[id];
                loop {
                    if gc_index < i && gc_index + i >= gcs_contents.len() {
                        panic!("No error correction found");
                    }
                    if gc_index >= i {
                        let index = gc_index - i;
                        let next_id = gcs_contents[index].1;
                        let lstrand = &self.strands[next_id];
                        if counter.get(&lstrand.strand).unwrap().len() > 1 {
                            let distance = hamming_distance(strand, &lstrand.strand);
                            if distance == 1 {
                                corrections.push((&self.strands[id], lstrand, false));
                                break;
                            }
                            let cstrand = self.strands[next_id].complement();
                            let distance = hamming_distance(strand, &cstrand.strand);
                            if distance == 1 {
                                corrections.push((&self.strands[id], lstrand, true));
                                break;
                            }
                        }
                    }
                    if gc_index + i < gcs_contents.len() {
                        let index = gc_index + i;
                        let next_id = gcs_contents[index].1;

                        let lstrand = &self.strands[next_id];
                        if counter.get(&lstrand.strand).unwrap().len() > 1 {
                            let distance = hamming_distance(strand, &lstrand.strand);
                            if distance == 1 {
                                corrections.push((&self.strands[id], lstrand, false));
                                break;
                            }
                            let cstrand = self.strands[next_id].complement();
                            let distance = hamming_distance(strand, &cstrand.strand);
                            if distance == 1 {
                                corrections.push((&self.strands[id], lstrand, true));
                                break;
                            }
                        }
                    }
                    i += 1;
                }
                visited.insert(id);
            }
        }
        corrections
    }
    pub fn profile_matrix(&self) -> Result<Vec<[usize; 4]>, StrandError> {
        let n = self
            .strands
            .first()
            .ok_or(StrandError::NoStrands)?
            .strand
            .len();
        let mut profile_matrix = vec![[0; 4]; n];

        for strand in &self.strands {
            if strand.strand.len() != n {
                return Err(StrandError::LengthsAreDifferent);
            }

            for (index, base) in strand.strand.iter().enumerate() {
                profile_matrix[index][base.index()] += 1;
            }
        }

        Ok(profile_matrix)
    }

    pub fn splice(&self) -> Vec<Vec<AA>> {
        let reference = &self.strands[0];
        let mut trie = STrie::new();
        for strand in &self.strands[1..] {
            trie.add(&strand.strand, strand.strand.len())
        }

        let mut skip = 0;
        let mut matches: VecDeque<_> = trie.matches(&reference.strand).into_iter().collect();

        let mut all = vec![];
        while let Some((mut pos, mut len)) = matches.pop_front() {
            let mut transcribed = vec![];
            let mut leftovers = vec![];
            for (i, c) in reference.strand.iter().enumerate() {
                if skip > 0 {
                    skip -= 1;
                    if skip > 0 {
                        continue;
                    }
                }
                match i.cmp(&pos) {
                    std::cmp::Ordering::Less => {
                        transcribed.push(*c);
                    }
                    std::cmp::Ordering::Equal => {
                        skip = len;
                        loop {
                            if let Some((p, l)) = matches.pop_front() {
                                if p < i + skip {
                                    leftovers.push((p, l));
                                } else {
                                    pos = p;
                                    len = l;
                                    break;
                                }
                            } else {
                                pos = usize::MAX;
                                len = 0;
                                break;
                            }
                        }
                    }
                    _ => {
                        unreachable!();
                    }
                }
            }
            for (p, l) in leftovers {
                matches.push_back((p, l));
            }

            let strand = Strand {
                name: "".to_string(),
                strand: transcribed,
            };
            if let Ok(protein) = strand.rna().protein() {
                all.push(protein);
            }
        }
        all
    }
}

pub enum ParseState {
    NAME,
    DATA,
}

pub fn parse_fasta<T: Base>(string: &[u8]) -> Result<Strands<T>, ParseError> {
    let mut name = vec![];
    let mut data = vec![];
    let mut state = ParseState::NAME;

    let mut strands = vec![];
    for c in string {
        let newstate = match (c, &state) {
            (b'>', _) => {
                if !name.is_empty() {
                    let sname = String::from_utf8(name)?;
                    if let Some(terminal) = T::terminal() {
                        data.push(terminal);
                    }
                    strands.push(Strand {
                        name: sname,
                        strand: data,
                    });
                    name = vec![];
                    data = vec![];
                }
                Some(ParseState::NAME)
            }
            (b'\n', ParseState::NAME) => Some(ParseState::DATA),
            (c, ParseState::NAME) => {
                name.push(*c);
                None
            }
            (b'\n', ParseState::DATA) => None,
            (c, ParseState::DATA) => {
                let base: T = (*c).try_into().map_err(|_| ParseError::UnexpectedData)?;
                data.push(base);
                None
            }
        };
        if let Some(newstate) = newstate {
            state = newstate;
        }
    }
    if let Some(terminal) = T::terminal() {
        data.push(terminal);
    }
    if !name.is_empty() {
        let sname = String::from_utf8(name)?;
        strands.push(Strand {
            name: sname,
            strand: data,
        });
    }
    Ok(Strands::from(strands))
}

pub fn parse<T: Base>(string: &[u8]) -> Result<Strand<T>, ParseError> {
    let mut parsed: Result<Vec<_>, _> = string
        .iter()
        .map(|c| {
            let base: Result<T, _> = (*c).try_into().map_err(|_| ParseError::UnexpectedData);
            base
        })
        .collect();
    if let Some(terminal) = T::terminal() {
        if let Ok(ref mut parsed) = parsed {
            parsed.push(terminal);
        }
    }
    match parsed {
        Ok(parsed) => Ok(Strand {
            name: "".to_string(),
            strand: parsed,
        }),
        Err(err) => Err(err),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use maplit::{hashmap, hashset};

    #[test]
    fn test_find_loci() {
        let data = parse::<DNABase>(b"ATGCTTCAGAAAGGTCTTACG").unwrap();
        let pattern = parse::<DNABase>(b"T").unwrap();

        assert_eq!(
            pattern.find(&data.strand).collect::<Vec<_>>(),
            vec![1, 4, 5, 14, 16, 17]
        );
    }

    #[test]
    fn test_prot() {
        let data =
            parse::<RNABase>(b"AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA").unwrap();

        let expected = parse::<AA>(b"MAMAPRTEINSTRING").unwrap().strand;
        assert_eq!(data.protein().unwrap(), expected);
    }

    #[test]
    fn test_trie_add() {
        let mut trie = PTrie::new();
        trie.add(&vec![DNABase::A], 0);
        assert_eq!(
            trie,
            PTrie {
                root: PNode {
                    subnodes: hashmap! {
                        DNABase::A => PNode{
                            subnodes: hashmap!{},
                            ids: hashset!{0},
                        }
                    },
                    ids: hashset! {0},
                }
            }
        );
    }

    #[test]
    fn test_strie_add() {
        let mut trie = STrie::new();
        trie.add(&vec![DNABase::A], 0);
        assert_eq!(
            trie,
            STrie {
                root: SNode {
                    subnodes: hashmap! {
                        DNABase::A => SNode{
                            subnodes: hashmap!{},
                            leafs: hashset!{0},
                        }
                    },
                    leafs: hashset! {},
                }
            }
        );
    }

    #[test]
    fn test_matches() {
        let data = parse::<DNABase>(b"ATGCTTCAGAAAGGTCTTACG").unwrap();
        let pattern = parse::<DNABase>(b"TT").unwrap();
        let pattern2 = parse::<DNABase>(b"ATG").unwrap();
        let pattern3 = parse::<DNABase>(b"TGC").unwrap();
        let mut trie = STrie::new();
        trie.add(&pattern.strand, 10);
        trie.add(&pattern2.strand, 11);
        trie.add(&pattern3.strand, 12);

        assert_eq!(
            trie.matches(&data.strand),
            vec![(0, 11), (1, 12), (4, 10), (16, 10)]
        );
    }
}
