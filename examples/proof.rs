use core::fmt::Debug;

use p3_baby_bear::BabyBear;
use p3_challenger::{HashChallenger, SerializingChallenger32};
use p3_commit::ExtensionMmcs;
use p3_field::extension::BinomialExtensionField;
use p3_field::{PrimeCharacteristicRing, PrimeField32};
use p3_fri::{HidingFriPcs, create_benchmark_fri_params_zk};
use p3_keccak::{Keccak256Hash, KeccakF};
use p3_merkle_tree::MerkleTreeHidingMmcs;
use p3_poseidon2_air::Poseidon2Air;
use p3_symmetric::{CompressionFunctionFromHasher, PaddingFreeSponge, SerializingHasher};
use p3_uni_stark::{StarkConfig, prove, verify};
use plonky3_merkle_tree::merkle_tree::{
    air::{MerkleTreeAir, P2_ROUND_CONSTANTS, generate_merkle_proof_trace},
    merkle_tree::MerkleTree,
};
use rand::SeedableRng;
use rand::rngs::SmallRng;
#[cfg(target_family = "unix")]
use tracing_forest::ForestLayer;
use tracing_forest::util::LevelFilter;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};

type Dft = p3_dft::Radix2DitParallel<BabyBear>;
// Settings taken from Plonky3/poseidon2-air/examples/prove_poseidon2_baby_bear_keccak_zk.rs
fn main() -> Result<(), impl Debug> {
    let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();

    Registry::default()
        .with(env_filter)
        .with(ForestLayer::default())
        .init();

    type Val = BabyBear;
    type Challenge = BinomialExtensionField<Val, 4>;

    type ByteHash = Keccak256Hash;
    let byte_hash = ByteHash {};

    type U64Hash = PaddingFreeSponge<KeccakF, 25, 17, 4>;
    let u64_hash = U64Hash::new(KeccakF {});

    type FieldHash = SerializingHasher<U64Hash>;
    let field_hash = FieldHash::new(u64_hash);

    type MyCompress = CompressionFunctionFromHasher<U64Hash, 2, 4>;
    let compress = MyCompress::new(u64_hash);

    type ValMmcs = MerkleTreeHidingMmcs<
        [Val; p3_keccak::VECTOR_LEN],
        [u64; p3_keccak::VECTOR_LEN],
        FieldHash,
        MyCompress,
        SmallRng,
        4,
        4,
    >;
    let rng = SmallRng::seed_from_u64(1);
    let val_mmcs = ValMmcs::new(field_hash, compress, 0, rng);

    type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());

    type Challenger = SerializingChallenger32<Val, HashChallenger<u8, ByteHash, 32>>;
    let challenger = Challenger::from_hasher(vec![], byte_hash);

    // create a test merkle tree
    let leaves: Vec<BabyBear> = vec![
        BabyBear::from_u32(1),
        BabyBear::ZERO,
        BabyBear::from_u32(2),
        BabyBear::from_u32(46),
        BabyBear::from_u32(78),
    ];
    let mt = MerkleTree::new(4, Some(leaves));
    let pf = mt.gen_merkle_proof(3);

    let p2_air = Poseidon2Air::new(P2_ROUND_CONSTANTS);
    let mt_air = MerkleTreeAir {
        root: mt.root.as_canonical_u32(),
        leaf_value: 46,
        poseidon2_air: p2_air,
    };

    let fri_params = create_benchmark_fri_params_zk(challenge_mmcs);

    let mt_trace = generate_merkle_proof_trace(3, mt.nodes[0][3], pf);

    let dft = Dft::default();

    type Pcs = HidingFriPcs<Val, Dft, ValMmcs, ChallengeMmcs, SmallRng>;
    let pcs = Pcs::new(dft, val_mmcs, fri_params, 4, SmallRng::seed_from_u64(1));

    type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;
    let config = MyConfig::new(pcs, challenger);

    let proof = prove(&config, &mt_air, mt_trace, &[]);

    verify(&config, &mt_air, &proof, &[])
}
