use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_field::PrimeCharacteristicRing;
use p3_symmetric::{PseudoCompressionFunction, TruncatedPermutation};

// WIDTH=3: state = [left, right, 0(capacity)] → permute → output state[0]
// CHUNK=1: each input/output is a single Bn254 element
// N=2: two inputs compressed to one output
type Poseidon2Compress = TruncatedPermutation<Poseidon2BabyBear<16>, 2, 1, 16>;
#[derive(Debug, Clone)]
pub struct MerkleTree {
    pub root: BabyBear,
    pub depth: usize,
    pub nodes: Vec<Vec<BabyBear>>,
}

impl MerkleTree {
    pub fn new(depth: usize, leaves: Option<Vec<BabyBear>>) -> Self {
        let zero = BabyBear::ZERO;
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

        for _ in 0..depth {
            current_level = current_level
                .chunks(2)
                .map(|pair| poseidon2_two_to_one(pair[0], pair[1]))
                .collect();
            nodes.push(current_level.clone());
        }

        let root = nodes.last().unwrap()[0];
        Self { root, depth, nodes }
    }

    pub fn gen_merkle_proof(&self, index: usize) -> Vec<BabyBear> {
        let mut current_index = index;
        let mut proof_nodes = vec![];
        for i in 0..self.depth {
            if current_index & 1 == 0 {
                proof_nodes.push(self.nodes[i][current_index + 1]);
            } else {
                proof_nodes.push(self.nodes[i][current_index - 1])
            }
            current_index = current_index >> 1;
        }
        proof_nodes
    }

    pub fn verify_merkle_proof(&self, proof: Vec<BabyBear>, index: usize) -> bool {
        let mut current_index = index;
        let mut result = self.nodes[0][index];
        for node in proof {
            if current_index & 1 == 0 {
                result = poseidon2_two_to_one(result, node);
            } else {
                result = poseidon2_two_to_one(node, result);
            }
            current_index = current_index >> 1;
        }
        result == self.root
    }
}

pub fn poseidon2_two_to_one(a: BabyBear, b: BabyBear) -> BabyBear {
    // let mut rng = ThreadRng::default();
    // let perm = Poseidon2BabyBear::<16>::new_from_rng(8, 22, &mut rng);
    // let compress = Poseidon2Compress::new(perm);
    let perm = p3_baby_bear::default_babybear_poseidon2_16();
    let compress = Poseidon2Compress::new(perm);
    compress.compress([[a], [b]])[0]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_proof() {
        let leaves: Vec<BabyBear> = vec![
            BabyBear::from_u32(1),
            BabyBear::ZERO,
            BabyBear::from_u32(2),
            BabyBear::from_u32(46),
            BabyBear::from_u32(78),
        ];
        let mt = MerkleTree::new(3, Some(leaves));
        let pf = mt.gen_merkle_proof(4);
        println!("tree: {:?}", mt);
        println!("proof: {:?}", pf);
        println!("result: {}", mt.verify_merkle_proof(pf, 4))
    }
}
