# Reference genome library for Rust
This is a library that will load a reference genome into memory for quick lookups.
This is largely an in-memory wrapper for a FASTA file.

Basic example:
```
let reference_fn = "./test_data/test_reference.fa";
let simple_reference_fn: PathBuf = PathBuf::from(reference_fn);
let reference_genome = ReferenceGenome::from_fasta(&simple_reference_fn).unwrap();

assert_eq!(reference_genome.contig_keys(), &[
    "chr1".to_string(),
    "chr2".to_string()
]);

//chr1 = ACGTACGT
let chr1_string: Vec<u8> = "ACGTACGT".as_bytes().to_vec();
assert_eq!(reference_genome.get_slice(&"chr1", 0, 8), &chr1_string);
```
