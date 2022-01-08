use rosalind::{parse_fasta, Base, DNABase, Strand, Strands};
use std::env;

#[derive(Debug)]
struct Mutation<T: Base> {
    locus: usize,
    before: Option<T>,
    after: Option<T>,
}

#[derive(Debug)]
struct MutationError;

impl<T: Base> Mutation<T> {
    fn from(string: &[u8]) -> Result<Mutation<T>, MutationError> {
        let n = string.len();
        let before: Option<T> = match string[0] {
            b'-' => None,
            c => Some(c.try_into().map_err(|_| MutationError)?),
        };
        let after: Option<T> = match string[n - 1] {
            b'-' => None,
            c => Some(c.try_into().map_err(|_| MutationError)?),
        };
        let locus = std::str::from_utf8(&string[1..n - 1])
            .map_err(|_| MutationError)?
            .parse::<usize>()
            .map_err(|_| MutationError)?;
        Ok(Self {
            locus,
            before,
            after,
        })
    }
}

fn parse_mut<T: Base>(mutation: &[u8]) -> Result<Vec<Mutation<T>>, MutationError> {
    let mut mutation_ = vec![];
    let mut mutations = vec![];
    for c in mutation {
        match c {
            b',' => {
                mutations.push(Mutation::from(&mutation_)?);
                mutation_ = vec![];
            }
            b' ' => {}
            c => mutation_.push(*c),
        }
    }
    Ok(mutations)
}

fn apply<T: Base>(
    strand: &Strand<T>,
    mutations: &[Mutation<T>],
) -> Result<Strands<T>, MutationError> {
    let mut mutations_iter = mutations.iter();
    let mut mutation = mutations_iter.next();
    let mut result = vec![];
    for (i, base) in strand.strand.iter().enumerate() {
        if let Some(mut _mutation) = mutation {
            let locus = i + 21563;
            // println!("locus {:?}", locus);
            while locus > _mutation.locus {
                mutation = mutations_iter.next();
                _mutation = mutation.unwrap();
            }
            if locus == _mutation.locus {
                // TODO Support insertions
                assert!(_mutation.before == Some(*base));
                if let Some(after) = _mutation.after {
                    // Replacement
                    result.push(after);
                } else {
                    // Deletion
                }
                mutation = mutations_iter.next();
            } else {
                result.push(*base);
            }
        }
    }
    Ok(Strands::from(vec![Strand {
        name: strand.name.clone(),
        strand: result,
    }]))
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let fasta_filename = &args[1];
    let string = std::fs::read(fasta_filename).unwrap();
    let strands = parse_fasta::<DNABase>(&string).unwrap();

    let mut_filename = &args[2];
    let string = std::fs::read(mut_filename).unwrap();
    let mutations = parse_mut::<DNABase>(&string).unwrap();

    let out = apply(&strands.strands[0], &mutations).unwrap();

    println!("{}", out);
}
