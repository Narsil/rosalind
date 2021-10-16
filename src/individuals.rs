use std::ops::Index;

#[derive(Debug, Clone)]
pub struct Individual<const N: usize> {
    probabilities: [[f32; 3]; N],
}

impl<const N: usize> Individual<N> {
    pub fn new(probabilities: [[f32; 3]; N]) -> Self {
        Self { probabilities }
    }
    pub fn from(allels: [Allel; N]) -> Self {
        let mut probabilities = [[0.0; 3]; N];
        for i in 0..N {
            probabilities[i] = match allels[i] {
                Allel::DD => [1.0, 0.0, 0.0],
                Allel::Dd => [0.0, 1.0, 0.0],
                Allel::dd => [0.0, 0.0, 1.0],
            };
        }
        Self { probabilities }
    }

    pub fn probability(&self, allels: [Allel; N]) -> f32 {
        let mut p = 1.0;
        for (probas, allel) in self.probabilities.iter().zip(allels) {
            p *= match allel {
                Allel::DD => probas[0],
                Allel::Dd => probas[1],
                Allel::dd => probas[2],
            };
        }
        p
    }
}

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy)]
pub enum Allel {
    DD,
    Dd,
    dd,
}

impl Index<Allel> for [f32; 3] {
    type Output = f32;

    fn index(&self, allel: Allel) -> &Self::Output {
        let i = allel as usize;
        &self[i]
    }
}

impl std::ops::IndexMut<Allel> for [f32; 3] {
    fn index_mut(&mut self, allel: Allel) -> &mut Self::Output {
        let i = allel as usize;
        &mut self[i]
    }
}
pub static DD: Allel = Allel::DD;
#[allow(non_upper_case_globals)]
pub static Dd: Allel = Allel::Dd;
#[allow(non_upper_case_globals)]
pub static dd: Allel = Allel::dd;

pub fn coupling<const N: usize>(left: &Individual<N>, right: &Individual<N>) -> Individual<N> {
    let mut probabilities = [[0.0; 3]; N];

    for ((ps, lp), rp) in probabilities
        .iter_mut()
        .zip(left.probabilities)
        .zip(right.probabilities)
    {
        ps[DD] = lp[DD] * rp[DD]
            + 0.5 * lp[Dd] * rp[DD]
            + lp[DD] * 0.5 * rp[Dd]
            + 0.5 * lp[Dd] * 0.5 * rp[Dd];
        ps[Dd] = lp[DD] * 0.5 * rp[Dd]
            + lp[DD] * rp[dd]
            + 2.0 * 0.5 * lp[Dd] * 0.5 * rp[Dd]
            + 0.5 * lp[Dd] * rp[dd]
            + lp[dd] * 0.5 * rp[DD]
            + lp[dd] * 0.5 * rp[Dd];
        ps[dd] = lp[dd] * rp[dd]
            + 0.5 * lp[Dd] * rp[DD]
            + lp[dd] * 0.5 * rp[Dd]
            + 0.5 * lp[Dd] * 0.5 * rp[Dd];
    }
    Individual::<N>::new(probabilities)
}
