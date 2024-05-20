#![allow(incomplete_features)]
#![allow(clippy::type_complexity)]

use bincode::{DefaultOptions, Options};
use cbl::kmer::Kmer;
use cbl::CBL;
use needletail::parse_fastx_file;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};

type T = u64;
const K: usize = 21;

// deserialize a given CBL
pub fn deserialize_cbl(input_filename: &str) -> CBL<K, T> {
    //let input_filename = format!("{}/{}.cbl", output_dir, input_index);
    let index =
        File::open(input_filename).unwrap_or_else(|_| panic!("Failed to open {}", input_filename));
    let reader = BufReader::new(index);
    DefaultOptions::new()
        .with_varint_encoding()
        .reject_trailing_bytes()
        .deserialize_from(reader)
        .unwrap()
}

pub fn serialize_cbl(cbl: &CBL<K, T>, output_filename: &str) {
    let _ = fs::remove_file(&output_filename);
    let output = File::create(output_filename).unwrap();
    let mut writer = BufWriter::new(output);
    DefaultOptions::new()
        .with_varint_encoding()
        .reject_trailing_bytes()
        .serialize_into(&mut writer, cbl)
        .unwrap();
}

pub fn create_cbl_from_fasta(input_filename: &str) -> CBL<K, T> {
    let mut reader = parse_fastx_file(input_filename).unwrap();
    let mut cbl = CBL::<K, T>::new();
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid record");
        cbl.insert_seq(&seqrec.seq());
    }
    cbl
}

pub fn cbl_printer(cbl: &CBL<K, T>, output_path: &str) -> std::io::Result<()> {
    if cbl.is_empty() {
        println!("Empty solution.");
        return Ok(());
    }
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    for (index, kmer) in cbl.iter().enumerate() {
        writeln!(writer, ">kmer{}", index)?;
        writer.write_all(&kmer.to_nucs())?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
