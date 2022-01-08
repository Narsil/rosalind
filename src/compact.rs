use crate::{Base, ParseError, ParseState, Strand, Strands};
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DNACompact {
    AAA,
    AAT,
    AAC,
    AAG,
    AAN,

    ACA,
    ACT,
    ACC,
    ACG,
    ACN,

    AGA,
    AGT,
    AGC,
    AGG,
    AGN,

    ATA,
    ATT,
    ATC,
    ATG,
    ATN,

    ANA,
    ANT,
    ANC,
    ANG,
    ANN,

    CAA,
    CAT,
    CAC,
    CAG,
    CAN,

    CCA,
    CCT,
    CCC,
    CCG,
    CCN,

    CGA,
    CGT,
    CGC,
    CGG,
    CGN,

    CTA,
    CTT,
    CTC,
    CTG,
    CTN,

    CNA,
    CNT,
    CNC,
    CNG,
    CNN,

    GAA,
    GAT,
    GAC,
    GAG,
    GAN,

    GCA,
    GCT,
    GCC,
    GCG,
    GCN,

    GGA,
    GGT,
    GGC,
    GGG,
    GGN,

    GTA,
    GTT,
    GTC,
    GTG,
    GTN,

    GNA,
    GNT,
    GNC,
    GNG,
    GNN,

    NAA,
    NAT,
    NAC,
    NAG,
    NAN,

    NCA,
    NCT,
    NCC,
    NCG,
    NCN,

    NGA,
    NGT,
    NGC,
    NGG,
    NGN,

    NTA,
    NTT,
    NTC,
    NTG,
    NTN,

    NNA,
    NNT,
    NNC,
    NNG,
    NNN,
}
impl std::convert::TryFrom<[u8; 3]> for DNACompact {
    type Error = ParseError;

    fn try_from(c: [u8; 3]) -> Result<Self, ParseError> {
        match c {
            [b'A', b'A', b'A'] => Ok(DNACompact::AAA),
            [b'A', b'A', b'T'] => Ok(DNACompact::AAT),
            [b'A', b'A', b'C'] => Ok(DNACompact::AAC),
            [b'A', b'A', b'G'] => Ok(DNACompact::AAG),
            [b'A', b'A', _] => Ok(DNACompact::AAN),

            [b'A', b'C', b'A'] => Ok(DNACompact::ACA),
            [b'A', b'C', b'T'] => Ok(DNACompact::ACT),
            [b'A', b'C', b'C'] => Ok(DNACompact::ACC),
            [b'A', b'C', b'G'] => Ok(DNACompact::ACG),
            [b'A', b'C', _] => Ok(DNACompact::ACN),

            [b'A', b'G', b'A'] => Ok(DNACompact::AGA),
            [b'A', b'G', b'T'] => Ok(DNACompact::AGT),
            [b'A', b'G', b'C'] => Ok(DNACompact::AGC),
            [b'A', b'G', b'G'] => Ok(DNACompact::AGG),
            [b'A', b'G', _] => Ok(DNACompact::AGN),

            [b'A', b'T', b'A'] => Ok(DNACompact::ATA),
            [b'A', b'T', b'T'] => Ok(DNACompact::ATT),
            [b'A', b'T', b'C'] => Ok(DNACompact::ATC),
            [b'A', b'T', b'G'] => Ok(DNACompact::ATG),
            [b'A', b'T', _] => Ok(DNACompact::ATN),

            [b'A', _, b'A'] => Ok(DNACompact::ANA),
            [b'A', _, b'T'] => Ok(DNACompact::ANT),
            [b'A', _, b'C'] => Ok(DNACompact::ANC),
            [b'A', _, b'G'] => Ok(DNACompact::ANG),
            [b'A', _, _] => Ok(DNACompact::ANN),

            [b'C', b'A', b'A'] => Ok(DNACompact::CAA),
            [b'C', b'A', b'T'] => Ok(DNACompact::CAT),
            [b'C', b'A', b'C'] => Ok(DNACompact::CAC),
            [b'C', b'A', b'G'] => Ok(DNACompact::CAG),
            [b'C', b'A', _] => Ok(DNACompact::CAN),

            [b'C', b'C', b'A'] => Ok(DNACompact::CCA),
            [b'C', b'C', b'T'] => Ok(DNACompact::CCT),
            [b'C', b'C', b'C'] => Ok(DNACompact::CCC),
            [b'C', b'C', b'G'] => Ok(DNACompact::CCG),
            [b'C', b'C', _] => Ok(DNACompact::CCN),

            [b'C', b'G', b'A'] => Ok(DNACompact::CGA),
            [b'C', b'G', b'T'] => Ok(DNACompact::CGT),
            [b'C', b'G', b'C'] => Ok(DNACompact::CGC),
            [b'C', b'G', b'G'] => Ok(DNACompact::CGG),
            [b'C', b'G', _] => Ok(DNACompact::CGN),

            [b'C', b'T', b'A'] => Ok(DNACompact::CTA),
            [b'C', b'T', b'T'] => Ok(DNACompact::CTT),
            [b'C', b'T', b'C'] => Ok(DNACompact::CTC),
            [b'C', b'T', b'G'] => Ok(DNACompact::CTG),
            [b'C', b'T', _] => Ok(DNACompact::CTN),

            [b'C', _, b'A'] => Ok(DNACompact::CNA),
            [b'C', _, b'T'] => Ok(DNACompact::CNT),
            [b'C', _, b'C'] => Ok(DNACompact::CNC),
            [b'C', _, b'G'] => Ok(DNACompact::CNG),
            [b'C', _, _] => Ok(DNACompact::CNN),

            [b'G', b'A', b'A'] => Ok(DNACompact::GAA),
            [b'G', b'A', b'T'] => Ok(DNACompact::GAT),
            [b'G', b'A', b'C'] => Ok(DNACompact::GAC),
            [b'G', b'A', b'G'] => Ok(DNACompact::GAG),
            [b'G', b'A', _] => Ok(DNACompact::GAN),

            [b'G', b'C', b'A'] => Ok(DNACompact::GCA),
            [b'G', b'C', b'T'] => Ok(DNACompact::GCT),
            [b'G', b'C', b'C'] => Ok(DNACompact::GCC),
            [b'G', b'C', b'G'] => Ok(DNACompact::GCG),
            [b'G', b'C', _] => Ok(DNACompact::GCN),

            [b'G', b'G', b'A'] => Ok(DNACompact::GGA),
            [b'G', b'G', b'T'] => Ok(DNACompact::GGT),
            [b'G', b'G', b'C'] => Ok(DNACompact::GGC),
            [b'G', b'G', b'G'] => Ok(DNACompact::GGG),
            [b'G', b'G', _] => Ok(DNACompact::GGN),

            [b'G', b'T', b'A'] => Ok(DNACompact::GTA),
            [b'G', b'T', b'T'] => Ok(DNACompact::GTT),
            [b'G', b'T', b'C'] => Ok(DNACompact::GTC),
            [b'G', b'T', b'G'] => Ok(DNACompact::GTG),
            [b'G', b'T', _] => Ok(DNACompact::GTN),

            [b'G', _, b'A'] => Ok(DNACompact::GNA),
            [b'G', _, b'T'] => Ok(DNACompact::GNT),
            [b'G', _, b'C'] => Ok(DNACompact::GNC),
            [b'G', _, b'G'] => Ok(DNACompact::GNG),
            [b'G', _, _] => Ok(DNACompact::GNN),

            [_, b'A', b'A'] => Ok(DNACompact::NAA),
            [_, b'A', b'T'] => Ok(DNACompact::NAT),
            [_, b'A', b'C'] => Ok(DNACompact::NAC),
            [_, b'A', b'G'] => Ok(DNACompact::NAG),
            [_, b'A', _] => Ok(DNACompact::NAN),

            [_, b'C', b'A'] => Ok(DNACompact::NCA),
            [_, b'C', b'T'] => Ok(DNACompact::NCT),
            [_, b'C', b'C'] => Ok(DNACompact::NCC),
            [_, b'C', b'G'] => Ok(DNACompact::NCG),
            [_, b'C', _] => Ok(DNACompact::NCN),

            [_, b'G', b'A'] => Ok(DNACompact::NGA),
            [_, b'G', b'T'] => Ok(DNACompact::NGT),
            [_, b'G', b'C'] => Ok(DNACompact::NGC),
            [_, b'G', b'G'] => Ok(DNACompact::NGG),
            [_, b'G', _] => Ok(DNACompact::NGN),

            [_, b'T', b'A'] => Ok(DNACompact::NTA),
            [_, b'T', b'T'] => Ok(DNACompact::NTT),
            [_, b'T', b'C'] => Ok(DNACompact::NTC),
            [_, b'T', b'G'] => Ok(DNACompact::NTG),
            [_, b'T', _] => Ok(DNACompact::NTN),

            [_, _, b'A'] => Ok(DNACompact::NNA),
            [_, _, b'T'] => Ok(DNACompact::NNT),
            [_, _, b'C'] => Ok(DNACompact::NNC),
            [_, _, b'G'] => Ok(DNACompact::NNG),
            [_, _, _] => Ok(DNACompact::NNN),
        }
    }
}

impl Base for DNACompact {}

pub fn parse_fasta_compact<T: Base + std::convert::TryFrom<[u8; 3]>>(
    string: &[u8],
) -> Result<Strands<T>, ParseError> {
    let mut name = vec![];
    let mut data = vec![];
    let mut state = ParseState::NAME;

    let mut strands = vec![];
    let mut pack = [0u8; 3];
    let mut index = 0;
    for c in string {
        let newstate = match (c, &state) {
            (b'>', _) => {
                if !name.is_empty() {
                    let sname = String::from_utf8(name)?;
                    if let Some(terminal) = T::terminal() {
                        data.push(terminal);
                    }
                    data.shrink_to_fit();
                    // println!("Strand {:?} {:?}", data.len(), data.capacity());
                    let strand = Strand {
                        name: sname,
                        strand: data,
                    };
                    strands.push(strand);
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
                if index < 3 {
                    pack[index] = *c;
                    index += 1;
                } else {
                    let base: T = pack
                        .try_into()
                        .map_err(|_| ParseError::UnexpectedData(*c))?;
                    data.push(base);
                    index = 0;
                }
                None
            }
        };
        if let Some(newstate) = newstate {
            if newstate != ParseState::DATA && state == ParseState::DATA && index != 0 {
                for tmp_index in index..3 {
                    pack[tmp_index] = b'N';
                }
                let base: T = pack
                    .try_into()
                    .map_err(|_| ParseError::UnexpectedData(*c))?;
                data.push(base);
                index = 0;
            }
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
