use thiserror::Error;

pub fn complement_base(c: u8) -> Option<u8> {
    match c {
        b'A' => Some(b'T'),
        b'T' => Some(b'A'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        _ => None,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DNABase {
    A,
    C,
    G,
    T,
}

#[derive(Debug)]
pub struct DNA {
    pub data: Vec<DNABase>,
}

pub struct DNAIterator<'a, 'b> {
    pattern_strand: &'a DNA,
    data_strand: &'b DNA,
    matches: Vec<(usize, usize)>,
    index: usize,
}

impl<'a, 'b> Iterator for DNAIterator<'a, 'b> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        while self.index < self.data_strand.data.len() {
            let c = self.data_strand.data[self.index];

            let mut matched = None;
            self.matches = self
                .matches
                .iter()
                .filter_map(|(data_index, pattern_index)| {
                    if self.pattern_strand.data.get(pattern_index + 1) == Some(&c) {
                        // Matching
                        if self.pattern_strand.data.len() == pattern_index + 2 {
                            // Final match
                            matched = Some(*data_index);
                        }
                        Some((*data_index, pattern_index + 1))
                    } else {
                        None
                    }
                })
                .collect();
            if c == self.pattern_strand.data[0] {
                if self.pattern_strand.data.len() == 1 {
                    // Final match
                    assert!(self.matches.len() == 0);
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

impl DNA {
    pub fn new() -> Self {
        Self { data: vec![] }
    }
    pub fn find<'a, 'b>(&'a self, data_strand: &'b DNA) -> DNAIterator<'a, 'b> {
        DNAIterator {
            pattern_strand: self,
            data_strand,
            matches: vec![],
            index: 0,
        }
    }
}

impl FromIterator<DNABase> for DNA {
    fn from_iter<I: IntoIterator<Item = DNABase>>(iter: I) -> Self {
        let mut c = DNA::new();

        for i in iter {
            c.data.push(i);
        }

        c
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RNABase {
    A,
    C,
    G,
    U,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Protein {
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
    U,
}

#[derive(Debug)]
pub enum ProteinError {
    UnexpectedRNALength,
}

pub struct RNA {
    pub data: Vec<RNABase>,
}

impl RNA {
    pub fn new() -> Self {
        Self { data: vec![] }
    }

    pub fn proteins(&self) -> Result<Vec<Protein>, ProteinError> {
        let protein: Result<Vec<_>, _> = self
            .data
            .chunks(3)
            .map(|chunk| match &chunk[..] {
                [RNABase::U, RNABase::U, RNABase::U] => Ok(Protein::F),
                [RNABase::C, RNABase::U, RNABase::U] => Ok(Protein::L),
                [RNABase::A, RNABase::U, RNABase::U] => Ok(Protein::I),
                [RNABase::G, RNABase::U, RNABase::U] => Ok(Protein::V),
                [RNABase::U, RNABase::U, RNABase::C] => Ok(Protein::F),
                [RNABase::C, RNABase::U, RNABase::C] => Ok(Protein::L),
                [RNABase::A, RNABase::U, RNABase::C] => Ok(Protein::I),
                [RNABase::G, RNABase::U, RNABase::C] => Ok(Protein::V),
                [RNABase::U, RNABase::U, RNABase::A] => Ok(Protein::L),
                [RNABase::C, RNABase::U, RNABase::A] => Ok(Protein::L),
                [RNABase::A, RNABase::U, RNABase::A] => Ok(Protein::I),
                [RNABase::G, RNABase::U, RNABase::A] => Ok(Protein::V),
                [RNABase::U, RNABase::U, RNABase::G] => Ok(Protein::L),
                [RNABase::C, RNABase::U, RNABase::G] => Ok(Protein::L),
                [RNABase::A, RNABase::U, RNABase::G] => Ok(Protein::M),
                [RNABase::G, RNABase::U, RNABase::G] => Ok(Protein::V),
                [RNABase::U, RNABase::C, RNABase::U] => Ok(Protein::S),
                [RNABase::C, RNABase::C, RNABase::U] => Ok(Protein::P),
                [RNABase::A, RNABase::C, RNABase::U] => Ok(Protein::T),
                [RNABase::G, RNABase::C, RNABase::U] => Ok(Protein::A),
                [RNABase::U, RNABase::C, RNABase::C] => Ok(Protein::S),
                [RNABase::C, RNABase::C, RNABase::C] => Ok(Protein::P),
                [RNABase::A, RNABase::C, RNABase::C] => Ok(Protein::T),
                [RNABase::G, RNABase::C, RNABase::C] => Ok(Protein::A),
                [RNABase::U, RNABase::C, RNABase::A] => Ok(Protein::S),
                [RNABase::C, RNABase::C, RNABase::A] => Ok(Protein::P),
                [RNABase::A, RNABase::C, RNABase::A] => Ok(Protein::T),
                [RNABase::G, RNABase::C, RNABase::A] => Ok(Protein::A),
                [RNABase::U, RNABase::C, RNABase::G] => Ok(Protein::S),
                [RNABase::C, RNABase::C, RNABase::G] => Ok(Protein::P),
                [RNABase::A, RNABase::C, RNABase::G] => Ok(Protein::T),
                [RNABase::G, RNABase::C, RNABase::G] => Ok(Protein::A),
                [RNABase::U, RNABase::A, RNABase::U] => Ok(Protein::Y),
                [RNABase::C, RNABase::A, RNABase::U] => Ok(Protein::H),
                [RNABase::A, RNABase::A, RNABase::U] => Ok(Protein::N),
                [RNABase::G, RNABase::A, RNABase::U] => Ok(Protein::D),
                [RNABase::U, RNABase::A, RNABase::C] => Ok(Protein::U),
                [RNABase::C, RNABase::A, RNABase::C] => Ok(Protein::H),
                [RNABase::A, RNABase::A, RNABase::C] => Ok(Protein::N),
                [RNABase::G, RNABase::A, RNABase::C] => Ok(Protein::D),
                [RNABase::U, RNABase::A, RNABase::A] => Ok(Protein::Stop),
                [RNABase::C, RNABase::A, RNABase::A] => Ok(Protein::Q),
                [RNABase::A, RNABase::A, RNABase::A] => Ok(Protein::K),
                [RNABase::G, RNABase::A, RNABase::A] => Ok(Protein::E),
                [RNABase::U, RNABase::A, RNABase::G] => Ok(Protein::Stop),
                [RNABase::C, RNABase::A, RNABase::G] => Ok(Protein::Q),
                [RNABase::A, RNABase::A, RNABase::G] => Ok(Protein::K),
                [RNABase::G, RNABase::A, RNABase::G] => Ok(Protein::E),
                [RNABase::U, RNABase::G, RNABase::U] => Ok(Protein::C),
                [RNABase::C, RNABase::G, RNABase::U] => Ok(Protein::R),
                [RNABase::A, RNABase::G, RNABase::U] => Ok(Protein::S),
                [RNABase::G, RNABase::G, RNABase::U] => Ok(Protein::G),
                [RNABase::U, RNABase::G, RNABase::C] => Ok(Protein::C),
                [RNABase::C, RNABase::G, RNABase::C] => Ok(Protein::R),
                [RNABase::A, RNABase::G, RNABase::C] => Ok(Protein::S),
                [RNABase::G, RNABase::G, RNABase::C] => Ok(Protein::G),
                [RNABase::U, RNABase::G, RNABase::A] => Ok(Protein::Stop),
                [RNABase::C, RNABase::G, RNABase::A] => Ok(Protein::R),
                [RNABase::A, RNABase::G, RNABase::A] => Ok(Protein::R),
                [RNABase::G, RNABase::G, RNABase::A] => Ok(Protein::G),
                [RNABase::U, RNABase::G, RNABase::G] => Ok(Protein::W),
                [RNABase::C, RNABase::G, RNABase::G] => Ok(Protein::R),
                [RNABase::A, RNABase::G, RNABase::G] => Ok(Protein::R),
                [RNABase::G, RNABase::G, RNABase::G] => Ok(Protein::G),
                _ => Err(ProteinError::UnexpectedRNALength),
            })
            .collect();
        protein
    }
}

impl FromIterator<RNABase> for RNA {
    fn from_iter<I: IntoIterator<Item = RNABase>>(iter: I) -> Self {
        let mut c = RNA::new();

        for i in iter {
            c.data.push(i);
        }

        c
    }
}

#[derive(Debug)]
pub struct Strand {
    pub name: String,
    pub strand: DNA,
}

impl Strand {
    pub fn gc(&self) -> f32 {
        let sum = self
            .strand
            .data
            .iter()
            .filter(|c| match c {
                DNABase::C | DNABase::G => true,
                _ => false,
            })
            .count();
        sum as f32 / self.strand.data.len() as f32
    }
}

pub type Strands = Vec<Strand>;

pub enum ParseState {
    NAME,
    DATA,
}

pub fn parse_fasta(string: &[u8]) -> Result<Strands, ParseError> {
    let mut name = vec![];
    let mut data = vec![];
    let mut state = ParseState::NAME;

    let mut strands = vec![];
    for c in string {
        let newstate = match (c, &state) {
            (b'>', _) => {
                if name.len() != 0 {
                    let sname = String::from_utf8(name)?;
                    strands.push(Strand {
                        name: sname,
                        strand: DNA { data },
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
            (b'A', ParseState::DATA) => {
                data.push(DNABase::A);
                None
            }
            (b'C', ParseState::DATA) => {
                data.push(DNABase::C);
                None
            }
            (b'G', ParseState::DATA) => {
                data.push(DNABase::G);
                None
            }
            (b'T', ParseState::DATA) => {
                data.push(DNABase::T);
                None
            }
            (b'\n', ParseState::DATA) => None,
            (_, ParseState::DATA) => {
                return Err(ParseError::UnexpectedData);
            }
        };
        if let Some(newstate) = newstate {
            state = newstate;
        }
    }
    if name.len() != 0 {
        let sname = String::from_utf8(name)?;
        strands.push(Strand {
            name: sname,
            strand: DNA { data },
        });
    }
    Ok(strands)
}

pub fn parse_dna(string: &[u8]) -> Result<DNA, ParseError> {
    string
        .into_iter()
        .map(|c| match c {
            b'A' => Ok(DNABase::A),
            b'C' => Ok(DNABase::C),
            b'G' => Ok(DNABase::G),
            b'T' => Ok(DNABase::T),
            _ => Err(ParseError::UnexpectedData),
        })
        .collect()
}

pub fn parse_rna(string: &[u8]) -> Result<RNA, ParseError> {
    string
        .into_iter()
        .map(|c| match c {
            b'A' => Ok(RNABase::A),
            b'C' => Ok(RNABase::C),
            b'G' => Ok(RNABase::G),
            b'U' => Ok(RNABase::U),
            _ => Err(ParseError::UnexpectedData),
        })
        .collect()
}

pub fn parse_prot(string: &[u8]) -> Result<Vec<Protein>, ParseError> {
    let mut result = string
        .into_iter()
        .map(|c| match c {
            b'F' => Ok(Protein::F),
            b'L' => Ok(Protein::L),
            b'S' => Ok(Protein::S),
            b'Y' => Ok(Protein::Y),
            b'C' => Ok(Protein::C),
            b'W' => Ok(Protein::W),
            b'P' => Ok(Protein::P),
            b'H' => Ok(Protein::H),
            b'Q' => Ok(Protein::Q),
            b'R' => Ok(Protein::R),
            b'I' => Ok(Protein::I),
            b'M' => Ok(Protein::M),
            b'T' => Ok(Protein::T),
            b'N' => Ok(Protein::N),
            b'K' => Ok(Protein::K),
            b'V' => Ok(Protein::V),
            b'A' => Ok(Protein::A),
            b'D' => Ok(Protein::D),
            b'E' => Ok(Protein::E),
            b'G' => Ok(Protein::G),
            b'U' => Ok(Protein::U),
            _ => Err(ParseError::UnexpectedData),
        })
        .collect::<Result<Vec<_>, _>>()?;
    result.push(Protein::Stop);
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_loci() {
        let data = parse_dna(b"ATGCTTCAGAAAGGTCTTACG").unwrap();
        let pattern = parse_dna(b"T").unwrap();

        assert_eq!(
            pattern.find(&data).collect::<Vec<_>>(),
            vec![1, 4, 5, 14, 16, 17]
        );
    }

    #[test]
    fn test_prot() {
        let data = parse_rna(b"AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA").unwrap();
        let expected = parse_prot(b"MAMAPRTEINSTRING").unwrap();
        assert_eq!(data.proteins().unwrap(), expected);
    }
}
