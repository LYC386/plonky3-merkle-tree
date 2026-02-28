use p3_bn254::Bn254;
use p3_bn254::Poseidon2Bn254;
use p3_field::PrimeCharacteristicRing;
use p3_symmetric::{PseudoCompressionFunction, TruncatedPermutation};
use rand::rngs::ThreadRng;

// WIDTH=3: state = [left, right, 0(capacity)] → permute → output state[0]
// CHUNK=1: each input/output is a single Bn254 element
// N=2: two inputs compressed to one output
type Poseidon2Compress = TruncatedPermutation<Poseidon2Bn254<3>, 2, 1, 3>;

#[derive(Debug, Clone)]
pub struct MerkleTree {
    pub root: Bn254,
    pub depth: usize,
    pub nodes: Vec<Vec<Bn254>>,
}

impl MerkleTree {
    pub fn new(depth: usize, leaves: Option<Vec<Bn254>>) -> Self {
        let mut rng = ThreadRng::default();
        let perm = Poseidon2Bn254::<3>::new_from_rng(8, 22, &mut rng);
        let compress = Poseidon2Compress::new(perm);

        let zero = Bn254::ZERO;
        // Build each level bottom-up; start with the zero digest
        let mut current_level = match leaves {
            Some(leaves) => {
                let total_leaves_len = 1 << depth;
                if leaves.len() >= total_leaves_len {
                    panic!("Too many leaves");
                }
                let mut current_level = leaves.clone();
                current_level.append(&mut vec![zero; total_leaves_len - leaves.len()]);
                current_level
            }
            None => vec![zero; 1 << depth],
        };
        let mut nodes = vec![current_level.clone()];

        for _ in 0..depth - 1 {
            current_level = current_level
                .chunks(2)
                .map(|pair| compress.compress([[pair[0]], [pair[1]]])[0])
                .collect();
            nodes.push(current_level.clone());
        }

        let root = nodes.last().unwrap()[0];
        Self { root, depth, nodes }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let leaves: Vec<Bn254> = vec![Bn254::from_u32(1), Bn254::ZERO, Bn254::from_u32(2)];
        let mt = MerkleTree::new(3, Some(leaves));
        println!("{:?}", mt);
    }
}
