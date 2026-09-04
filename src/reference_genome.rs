
use anyhow::ensure;
use bio::io::fasta;
use flate2::bufread::MultiGzDecoder;
use indexmap::IndexMap;
use log::{debug, warn};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Wrapper structure for a reference genome
pub struct ReferenceGenome {
    /// The filename we loaded 
    filename: PathBuf,
    /// Map where keys are contig names and value is ASCII formatted sequence
    contig_map: IndexMap<String, Vec<u8>>
}

impl ReferenceGenome {
    /// Creates an empty reference genome, which we can be populated through `add_contig(...)`
    pub fn empty_reference() -> Self {
        Self {
            filename: PathBuf::from(""),
            contig_map: Default::default()
        }
    }

    /// Loads a reference genome from a given FASTA file
    /// # Arguments
    /// * `fasta_fn` - the FASTA filename, gzip is allowed
    /// # Errors
    /// * if the FASTA file cannot be opened or read
    /// * if a FASTA record cannot be parsed
    /// * if a contig key is duplicated in the FASTA file
    pub fn from_fasta(fasta_fn: &Path) -> anyhow::Result<ReferenceGenome> {
        debug!("Loading {:?}...", fasta_fn);
        
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

        // create a new reference genome where we will add contigs as we read them
        let mut reference_genome = Self {
            filename: fasta_fn.to_path_buf(),
            contig_map: Default::default()
        };

        // read the FASTA file and add contigs to the reference genome
        for entry in fasta_reader.records() {
            let record: fasta::Record = entry?;
            let seq_id = record.id().to_string();
            let sequence = record.seq().to_vec();
            reference_genome.add_contig(seq_id, sequence)?;
        }
        debug!("Finished loading {} contigs.", reference_genome.contig_map.len());

        Ok(reference_genome)
    }

    /// Adds a new contig to the reference genome
    /// # Arguments
    /// * `contig_key` - the name of the contig
    /// * `contig_sequence` - the sequence to add; all sequence is automatically upper-cased in place
    /// # Errors
    /// * if `contig_key` is already in the reference genome
    pub fn add_contig(&mut self, contig_key: String, mut contig_sequence: Vec<u8>) -> anyhow::Result<()> {
        // check for a duplicate contig key
        ensure!(
            !self.contig_map.contains_key(&contig_key),
            "Duplicate contig key detected: {contig_key}"
        );

        // make the sequence uppercase in place and then add it to the lookup map
        contig_sequence.make_ascii_uppercase();
        self.contig_map.insert(contig_key, contig_sequence);
        Ok(())
    }

    /// Adds a new contig from a UTF-8 string by converting it to bytes and calling [`Self::add_contig`].
    /// # Arguments
    /// * `contig_key` - the name of the contig
    /// * `contig_sequence` - the sequence to add; all sequence is automatically upper-cased
    /// # Errors
    /// * if `contig_key` is already in the reference genome
    pub fn add_contig_string(&mut self, contig_key: String, contig_sequence: String) -> anyhow::Result<()> {
        self.add_contig(contig_key, contig_sequence.into_bytes())
    }

    /// Returns the filename of the reference genome.
    pub fn filename(&self) -> &Path {
        &self.filename
    }

    /// Returns an iterator over the contig keys in the reference genome.
    /// The order of the contig keys is guaranteed to be the same as the order as they were added to the reference genome.
    pub fn contig_keys(&self) -> impl Iterator<Item = &String> {
        self.contig_map.keys()
    }

    /// Retrieves a reference slice from a given 0-based coordinates.
    /// If `start` or `end` goes past the full contig length, it will be truncated to the full contig length.
    /// # Arguments
    /// * `chromosome` - the chromosome to slice from
    /// * `start` - the 0-based start index (included)
    /// * `end` - the 0-based end index (excluded)
    /// # Errors
    /// * if `chromosome` was not in the FASTA file
    /// * if `start` > `end`
    pub fn get_slice(&self, chromosome: &str, start: usize, end: usize) -> anyhow::Result<&[u8]> {
        // check the inputs
        let full_contig = self.contig_map.get(chromosome)
            .ok_or_else(|| anyhow::anyhow!("chromosome {chromosome:?} was not in the reference file"))?;
        ensure!(start <= end, "start > end: {start} > {end}");

        // truncate the start and end if they go past the full contig length
        let truncated_start = if start <= full_contig.len() { start } else {
            warn!("Received get_slice({:?}, {}, {}), truncated start to {}", chromosome, start, end, full_contig.len());
            full_contig.len()
        };
        let truncated_end = if end <= full_contig.len() { end } else {
            warn!("Received get_slice({:?}, {}, {}), truncated end to {}", chromosome, start, end, full_contig.len());
            full_contig.len()
        };
        Ok(&full_contig[truncated_start..truncated_end])
    }

    /// Retrieves a full chromosome by name
    /// # Arguments
    /// * `chromosome` - the chromosome to slice from
    /// # Errors
    /// * if `chromosome` was not in the FASTA file
    pub fn get_full_chromosome(&self, chromosome: &str) -> anyhow::Result<&[u8]> {
        self.contig_map.get(chromosome)
            .map(|s| s.as_slice())
            .ok_or_else(|| anyhow::anyhow!("chromosome {chromosome:?} was not in the reference file"))
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

            assert_eq!(
                reference_genome.contig_keys().map(String::as_str).collect::<Vec<_>>(),
                ["chr1", "chr2"]
            );

            //chr1 = ACGTACGT
            let chr1_string: Vec<u8> = "ACGTACGT".as_bytes().to_vec();
            for i in 0..8 {
                assert_eq!(reference_genome.get_slice(&"chr1", i, 8).unwrap(), &chr1_string[i..]);
            }

            //chr2 = ACCATGTA
            let chr2_string: Vec<u8> = "ACCATGTA".as_bytes().to_vec();
            assert_eq!(reference_genome.get_slice(&"chr2", 0, 8).unwrap(), chr2_string);
        }
    }

    #[test]
    fn test_add_contig() {
        let mut reference_genome = ReferenceGenome::empty_reference();
        reference_genome.add_contig_string("test".to_string(), "Acgt".to_string()).unwrap();
        reference_genome.add_contig_string("test2".to_string(), "TGNA".to_string()).unwrap();

        assert_eq!(
            reference_genome.contig_keys().map(String::as_str).collect::<Vec<_>>(),
            ["test", "test2"]
        );
        assert_eq!(reference_genome.get_full_chromosome("test").unwrap(), b"ACGT");
        assert_eq!(reference_genome.get_full_chromosome("test2").unwrap(), b"TGNA");
    }

    #[test]
    fn test_lookup_errors() {
        let mut reference_genome = ReferenceGenome::empty_reference();
        reference_genome.add_contig("chr1".to_string(), b"ACGT".to_vec()).unwrap();

        assert!(reference_genome.get_full_chromosome("chrX").is_err());
        assert!(reference_genome.get_slice("chrX", 0, 1).is_err());
        assert!(reference_genome.get_slice("chr1", 3, 1).is_err());
        assert!(reference_genome.add_contig("chr1".to_string(), b"AAAA".to_vec()).is_err());
    }
}
