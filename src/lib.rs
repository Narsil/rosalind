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
            println!("{:?}", c);

            let mut matched = None;
            self.matches = self
                .matches
                .iter()
                .filter_map(|(data_index, pattern_index)| {
                    println!("{:?}"
                    if self.pattern_strand.data.get(pattern_index + 1) == Some(&c) {
                        // Matching
                        if self.pattern_strand.data.len() == pattern_index + 1 {
                            // Final match
                            matched = Some(*data_index);
                        }
                        println!("second match");
                        Some((*data_index, pattern_index + 1))
                    } else {
                        None
                    }
                })
                .collect();
            if c == self.pattern_strand.data[0] {
                println!("First match");
                if self.pattern_strand.data.len() == 1 {
                    // Final match
                    assert!(self.matches.len() == 0);
                    matched = Some(self.index);
                } else {
                    self.matches.push((self.index, 1));
                }
            }
            self.index += 1;
            if matched.is_some() {
                println!("MATCH");
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
    pub fn find_locus<'a, 'b>(&'a self, data_strand: &'b DNA) -> DNAIterator<'a, 'b> {
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

#[derive(Debug)]
pub struct Strand {
    pub name: String,
    pub data: Vec<DNABase>,
}

impl Strand {
    pub fn gc(&self) -> f32 {
        let sum = self
            .data
            .iter()
            .filter(|c| match c {
                DNABase::C | DNABase::G => true,
                _ => false,
            })
            .count();
        sum as f32 / self.data.len() as f32
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
                    strands.push(Strand { name: sname, data });
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
        strands.push(Strand { name: sname, data });
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
