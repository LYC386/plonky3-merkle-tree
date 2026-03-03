use p3_air::{Air, AirBuilder, AirBuilderWithContext, BaseAir};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;

use p3_baby_bear::{BabyBear, GenericPoseidon2LinearLayersBabyBear};
use p3_challenger::{HashChallenger, SerializingChallenger32};
use p3_circle::CirclePcs;
use p3_commit::ExtensionMmcs;
use p3_field::extension::BinomialExtensionField;
use p3_keccak::Keccak256Hash;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_mersenne_31::Mersenne31;
use p3_uni_stark::{StarkConfig, prove, verify};

use p3_poseidon2_air::{Poseidon2Air, Poseidon2Cols, RoundConstants};
use std::borrow::Borrow;

use p3_uni_stark::SubAirBuilder;

// Consider making it generic in MerkleTreeAir struct
const WIDTH: usize = 16;
const SBOX_DEGREE: u64 = 7;
const SBOX_REGISTERS: usize = 1;
const HALF_FULL_ROUNDS: usize = 4;
const PARTIAL_ROUNDS: usize = 13;

type Poseidon2AirType = Poseidon2Air<
    BabyBear,
    GenericPoseidon2LinearLayersBabyBear,
    WIDTH,
    SBOX_DEGREE,
    SBOX_REGISTERS,
    HALF_FULL_ROUNDS,
    PARTIAL_ROUNDS,
>;

//TODO - put root, depth, leaf_value in public_values of AirBuilder instead of hard codeing like now.
pub struct MerkleTreeAir {
    pub root: u32,
    pub depth: usize,
    pub leaf_value: u32,
    pub poseidon2_air: Poseidon2Air<
        BabyBear,
        GenericPoseidon2LinearLayersBabyBear,
        WIDTH,
        SBOX_DEGREE,
        SBOX_REGISTERS,
        HALF_FULL_ROUNDS,
        PARTIAL_ROUNDS,
    >, // pub node: u32,
}

impl<F: Field> BaseAir<F> for MerkleTreeAir {
    fn width(&self) -> usize {
        // is_odd_index + Poseidon2Col2(includes current_node and sibling of the merkle proof)
        self.poseidon2_air.width() + 1
    }
}

impl<AB: AirBuilder<F = BabyBear>> Air<AB> for MerkleTreeAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main.row_slice(0).unwrap();
        let next = main.row_slice(1).unwrap();

        // create Poseidon2 sub air to evaluate Poseidon2 hash
        // given input as column range [1..]
        let p2_col_count = self.poseidon2_air.width();
        let mut sub: SubAirBuilder<AB, Poseidon2AirType, AB::Var> =
            SubAirBuilder::new(builder, 1..1 + p2_col_count);
        self.poseidon2_air.eval(&mut sub);

        // Use helper function from Poseidon2Cols to get the hash result
        // Cast the slice into Poseidon2Cols (starting at col 1)
        let p2_cols: &Poseidon2Cols<
            _,
            WIDTH,
            SBOX_DEGREE,
            SBOX_REGISTERS,
            HALF_FULL_ROUNDS,
            PARTIAL_ROUNDS,
        > = local[1..].borrow();
        let hash_output = p2_cols.ending_full_rounds[HALF_FULL_ROUNDS - 1].post[0];

        // Constrain: first left (col[1]) or right (col[2]) input should be the leaf_value based on
        // is_odd_index (col[0])
        builder.when_first_row().assert_zero(
            (local[0] * (AB::Expr::from_u32(self.leaf_value) - local[2]))
                + ((AB::Expr::ONE - local[0]) * (AB::Expr::from_u32(self.leaf_value) - local[1])),
        );

        // Constrain: hash output == next row's left (col[1]) or right (col[2]) input based on
        // is_odd_index (col[0])
        builder.when_transition().assert_zero(
            (next[0] * (hash_output - next[2]))
                + ((AB::Expr::ONE - next[0]) * (hash_output - next[1])),
        );

        // Constrain: final output == root
        let root = AB::Expr::from_u32(self.root);
        builder.when_last_row().assert_eq(root, hash_output);

        // The above constrain should cover both below
        // Constrain: hash output == next row's current node
        // builder.when_transition().assert_eq(hash_output, next[0]);

        // Constrain: hash input == (sibling, current_node) or (current_node, sibling) based on is_odd_index
        // builder.when_transition().assert_zeros([]);

        // Constrain: is_odd should be bool
        builder.when_first_row().assert_bool(local[0].clone());
        builder.when_transition().assert_bool(local[0].clone());
        builder.when_last_row().assert_bool(local[0].clone());
    }
}
