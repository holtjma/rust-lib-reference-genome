
use bio::io::fasta;
use flate2::bufread::MultiGzDecoder;
use log::{debug, warn};
use rustc_hash::FxHashMap as HashMap;
use simple_error::{bail, SimpleError};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Wrapper structure for a reference genome
pub struct ReferenceGenome {
    /// The filename we loaded 
    filename: PathBuf,
    /// Contains the keys in order of the reference load
    contig_keys: Vec<String>,
    /// Map where keys are contig names and value is ASCII formatted sequence
    contig_map: HashMap<String, Vec<u8>>
}

impl ReferenceGenome {
    /// Creates an empty reference genome, which we can be populated through `add_contig(...)`
    pub fn empty_reference() -> Self {
        Self {
            filename: PathBuf::from(""),
            contig_keys: vec![],
            contig_map: Default::default()
        }
    }

    /// Loads a reference genome from a given FASTA file
    /// # Arguments
    /// * `fasta_fn` - the FASTA filename, gzip is allowed
    /// # Errors
    /// This will pass through any error detected from loading the provided FASTA file.
    /// This includes file reading and/or record reading errors.
    pub fn from_fasta(fasta_fn: &Path) -> Result<ReferenceGenome, Box<dyn std::error::Error>> {
        debug!("Loading {:?}...", fasta_fn);
        let mut contig_keys: Vec<String> = Default::default();
        let mut contig_map: HashMap<String, Vec<u8>> = Default::default();
        
        // needletail can technically read FASTA and FASTQ, not sure we can check for that easy though
        let fasta_file: std::fs::File = std::fs::File::open(fasta_fn)?;
        let file_reader = BufReader::new(fasta_file);
        let fasta_reader: fasta::Reader<Box<dyn BufRead>> = if fasta_fn.extension().unwrap_or_default() == "gz" {
            debug!("Detected gzip extension, loading reference with MultiGzDecoder...");
            let gz_decoder = MultiGzDecoder::new(file_reader);
            let bufreader = BufReader::new(gz_decoder);
            fasta::Reader::from_bufread(Box::new(bufreader))
        } else {
            debug!("Loading reference as plain-text file...");
            fasta::Reader::from_bufread(Box::new(file_reader))
        };

        for entry in fasta_reader.records() {
            let record: fasta::Record = entry?;
            let seq_id: String = record.id().to_string();
            let sequence: Vec<u8> = record.seq().to_ascii_uppercase();

            contig_keys.push(seq_id.clone());
            contig_map.insert(seq_id, sequence);
        }
        debug!("Finished loading {} contigs.", contig_map.len());

        Ok(ReferenceGenome {
            filename: fasta_fn.to_path_buf(),
            contig_keys,
            contig_map
        })
    }

    /// Adds a new contig to the reference genome
    /// # Arguments
    /// * `contig_key` - the name of the contig
    /// * `contig_sequence` - the sequence to add; all sequence is automatically upper-cased
    pub fn add_contig(&mut self, contig_key: String, contig_sequence: &str) -> Result<(), SimpleError> {
        if self.contig_map.contains_key(&contig_key) {
            bail!("Contig key \"{contig_key}\" is already in the reference genome");
        }

        // create the uppercase byte form
        let byte_form = contig_sequence.to_ascii_uppercase().into_bytes();
        
        // save everything
        self.contig_keys.push(contig_key.clone());
        self.contig_map.insert(contig_key, byte_form);
        Ok(())
    }

    pub fn filename(&self) -> &Path {
        &self.filename
    }

    pub fn contig_keys(&self) -> &[String] {
        &self.contig_keys
    }

    /// Retrieves a reference slice from a given 0-based coordinates.
    /// If `start` or `end` goes past the full contig length, it will be truncated to the full contig length.
    /// # Arguments
    /// * `chromosome` - the chromosome to slice from
    /// * `start` - the 0-based start index (included)
    /// * `end` - the 0-based end index (excluded)
    /// # Panics
    /// * if `chromosome` was not in the FASTA file
    /// * if `start` > `end`
    pub fn get_slice(&self, chromosome: &str, start: usize, end: usize) -> &[u8] {
        let full_contig = self.contig_map.get(chromosome).expect("a chromosome from the reference file");
        assert!(start <= end, "start > end: {start} > {end}");
        let truncated_start = if start <= full_contig.len() { start } else {
            warn!("Received get_slice({:?}, {}, {}), truncated start to {}", chromosome, start, end, full_contig.len());
            full_contig.len()
        };
        let truncated_end = if end <= full_contig.len() { end } else {
            warn!("Received get_slice({:?}, {}, {}), truncated end to {}", chromosome, start, end, full_contig.len());
            full_contig.len()
        };
        &full_contig[truncated_start..truncated_end]
    }

    /// Retrieves a full chromosome by name
    /// # Arguments
    /// * `chromosome` - the chromosome to slice from
    /// # Panics
    /// * if `chromosome` was not in the FASTA file
    pub fn get_full_chromosome(&self, chromosome: &str) -> &[u8] {
        let full_contig = self.contig_map.get(chromosome).expect("a chromosome from the reference file");
        full_contig
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    #[test]
    fn test_simple_reference() {
        let references = vec![
            "./test_data/test_reference.fa",
            "./test_data/test_reference.fa.gz"
        ];
        for &reference_fn in references.iter() {
            let simple_reference_fn: PathBuf = PathBuf::from(reference_fn);
            let reference_genome = ReferenceGenome::from_fasta(&simple_reference_fn).unwrap();

            assert_eq!(reference_genome.contig_keys(), &[
                "chr1".to_string(),
                "chr2".to_string()
            ]);

            //chr1 = ACGTACGT
            let chr1_string: Vec<u8> = "ACGTACGT".as_bytes().to_vec();
            for i in 0..8 {
                assert_eq!(reference_genome.get_slice(&"chr1", i, 8), &chr1_string[i..]);
            }

            //chr2 = ACCATGTA
            let chr1_string: Vec<u8> = "ACCATGTA".as_bytes().to_vec();
            assert_eq!(reference_genome.get_slice(&"chr2", 0, 8), chr1_string);
        }
    }

    #[test]
    fn test_add_contig() {
        let mut reference_genome = ReferenceGenome::empty_reference();
        reference_genome.add_contig("test".to_string(), "Acgt").unwrap();
        reference_genome.add_contig("test2".to_string(), "TGNA").unwrap();

        assert_eq!(reference_genome.contig_keys(), &["test".to_string(), "test2".to_string()]);
        assert_eq!(reference_genome.get_full_chromosome("test"), b"ACGT");
        assert_eq!(reference_genome.get_full_chromosome("test2"), b"TGNA");
    }
}