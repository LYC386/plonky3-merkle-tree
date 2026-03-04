use p3_air::{Air, AirBuilder, BaseAir};
use p3_baby_bear::{BabyBear, GenericPoseidon2LinearLayersBabyBear};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;
use p3_poseidon2_air::{Poseidon2Air, Poseidon2Cols, RoundConstants};
use std::borrow::Borrow;

use p3_uni_stark::SubAirBuilder;

// Consider making it generic in MerkleTreeAir struct
pub const WIDTH: usize = 16;
pub const SBOX_DEGREE: u64 = 7;
pub const SBOX_REGISTERS: usize = 1;
pub const HALF_FULL_ROUNDS: usize = 4;
pub const PARTIAL_ROUNDS: usize = 13;
pub const P2_ROUND_CONSTANTS: RoundConstants<BabyBear, WIDTH, HALF_FULL_ROUNDS, PARTIAL_ROUNDS> =
    RoundConstants::<BabyBear, WIDTH, HALF_FULL_ROUNDS, PARTIAL_ROUNDS>::new(
        p3_baby_bear::BABYBEAR_RC16_EXTERNAL_INITIAL,
        p3_baby_bear::BABYBEAR_RC16_INTERNAL,
        p3_baby_bear::BABYBEAR_RC16_EXTERNAL_FINAL,
    );

type Poseidon2AirType = Poseidon2Air<
    BabyBear,
    GenericPoseidon2LinearLayersBabyBear,
    WIDTH,
    SBOX_DEGREE,
    SBOX_REGISTERS,
    HALF_FULL_ROUNDS,
    PARTIAL_ROUNDS,
>;

//TODO - put root, leaf_value in public_values of AirBuilder instead of hard codeing like now.
pub struct MerkleTreeAir {
    pub root: u32,
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
        // Transmute the slice into Poseidon2Cols (starting at col 1)
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

pub fn generate_merkle_proof_trace(
    leaf_index: u32,
    leaf_value: BabyBear,
    merkle_proof: Vec<BabyBear>,
) -> RowMajorMatrix<BabyBear> {
    // let round_constants = RoundConstants::<BabyBear, WIDTH, HALF_FULL_ROUNDS, PARTIAL_ROUNDS>::new(
    //     p3_baby_bear::BABYBEAR_RC16_EXTERNAL_INITIAL,
    //     p3_baby_bear::BABYBEAR_RC16_INTERNAL,
    //     p3_baby_bear::BABYBEAR_RC16_EXTERNAL_FINAL,
    // );
    let mut current_node = leaf_value;
    let mut current_index = leaf_index;
    let mut mt_trace_vec = vec![];

    // The number of trace rows needs to be a power of two for the STARK prover.
    // Might require some padding depends on merkle_proof.len()
    for sibling in merkle_proof {
        let is_odd = current_index & 1 == 1;
        let input = if is_odd {
            let mut input = vec![[BabyBear::ZERO; WIDTH]];
            input[0][0] = sibling;
            input[0][1] = current_node;
            input
        } else {
            let mut input = vec![[BabyBear::ZERO; WIDTH]];
            input[0][0] = current_node;
            input[0][1] = sibling;
            input
        };
        // generate trace for poseidon2 air
        let mut p2_trace = p3_poseidon2_air::generate_trace_rows::<
            BabyBear,
            GenericPoseidon2LinearLayersBabyBear,
            WIDTH,
            SBOX_DEGREE,
            SBOX_REGISTERS,
            HALF_FULL_ROUNDS,
            PARTIAL_ROUNDS,
        >(input, &P2_ROUND_CONSTANTS, 0);
        // Use helper function from Poseidon2Cols to get the hash result
        // Transmute the slice into Poseidon2Cols (starting at col 1)
        let p2_cols: &Poseidon2Cols<
            _,
            WIDTH,
            SBOX_DEGREE,
            SBOX_REGISTERS,
            HALF_FULL_ROUNDS,
            PARTIAL_ROUNDS,
        > = p2_trace.values[..].borrow();
        let output = p2_cols.ending_full_rounds[HALF_FULL_ROUNDS - 1].post[0];
        // push is_odd + p2_trace_vec to mt_trace_vec
        mt_trace_vec.push(BabyBear::from_bool(is_odd));
        mt_trace_vec.append(&mut p2_trace.values);
        //update current_node, current_index for next round
        current_node = output;
        current_index >>= 1;
    }
    // mt_trace_width = p2_trace_width + 1;
    let mt_trace_width = p3_poseidon2_air::num_cols::<
        WIDTH,
        SBOX_DEGREE,
        SBOX_REGISTERS,
        HALF_FULL_ROUNDS,
        PARTIAL_ROUNDS,
    >() + 1;
    let mt_trace = RowMajorMatrix::new(mt_trace_vec, mt_trace_width);
    mt_trace
}
