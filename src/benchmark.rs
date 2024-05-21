#![feature(generic_const_exprs)]

use anti_reindeer::utils::{
    create_cbl_from_fasta, deserialize_cbl, read_fof_file_csv, serialize_cbl,
};
use std::env;
use std::fs;
use std::path::Path;
use std::time::Instant;

use cbl::CBL;

type T = u64;
const K: usize = 21;

fn merge_cbls_in_batches(file_paths: &[String], batch_size: usize) -> CBL<K, T> {
    let mut input_iter = file_paths.chunks(batch_size);
    let mut global_cbl = if let Some(input_filename_chunk) = input_iter.next() {
        let mut cbls_chunk: Vec<_> = input_filename_chunk
            .iter()
            .map(|input_filename| deserialize_cbl(input_filename))
            .collect();
        CBL::<K, T>::merge(cbls_chunk.iter_mut().collect())
    } else {
        panic!("No CBL files to merge");
    };

    for input_filename_chunk in input_iter {
        let mut cbls_chunk: Vec<_> = input_filename_chunk
            .iter()
            .map(|input_filename| deserialize_cbl(input_filename))
            .collect();
        global_cbl |= &mut CBL::<K, T>::merge(cbls_chunk.iter_mut().collect());
    }

    global_cbl
}

fn intersect_cbls_in_batches(file_paths: &[String], batch_size: usize) -> CBL<K, T> {
    let mut input_iter = file_paths.chunks(batch_size);
    let mut global_cbl = if let Some(input_filename_chunk) = input_iter.next() {
        let mut cbls_chunk: Vec<_> = input_filename_chunk
            .iter()
            .map(|input_filename| deserialize_cbl(input_filename))
            .collect();
        CBL::<K, T>::intersect(cbls_chunk.iter_mut().collect())
    } else {
        panic!("No CBL files to intersect");
    };

    for input_filename_chunk in input_iter {
        let mut cbls_chunk: Vec<_> = input_filename_chunk
            .iter()
            .map(|input_filename| deserialize_cbl(input_filename))
            .collect();
        global_cbl &= &mut CBL::<K, T>::intersect(cbls_chunk.iter_mut().collect());
    }

    global_cbl
}



fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "usage: {} <input_file_list> <mode> <max_batch_size>",
            args[0]
        );
        std::process::exit(1);
    }
    let input_file_list = &args[1];
    let mode = &args[2];
    let batch_size: usize = args[3].parse().expect("max batch size must be a number");
    if mode != "union" && mode != "intersection" {
        eprintln!("Invalid mode. Use 'union' or 'intersection'.");
        std::process::exit(1);
    }

    let output_dir = "results_benchmark";
    let output_path = Path::new(output_dir);
    if !output_path.exists() {
        fs::create_dir_all(output_dir).unwrap();
    }

    println!("Loading files and writing CBLs...");

    let (to_load, col_nb) = read_fof_file_csv(&input_file_list).unwrap();
    let indices: Vec<usize> = (0..col_nb).collect();
    let mut file_paths: Vec<String> = vec![];
    for (i, input_filename) in to_load.iter().enumerate() {
        let output_filename = format!("{}/{}u.cbl", output_dir, indices[i]);
        let output_file_path = Path::new(&output_filename);
        if !output_file_path.exists() {
            let cbl = create_cbl_from_fasta(input_filename);
            let start_serialize = Instant::now();
            serialize_cbl(&cbl, &output_filename);
            let duration_serialize = start_serialize.elapsed();
            println!("Created CBL for file: {} in {:?}", input_filename, duration_serialize);
        } else {
            println!("CBL already exists for file: {}", input_filename);
        }

        file_paths.push(output_filename);
    }
    let cbl = merge_cbls_in_batches(&file_paths, batch_size);
    let kmer_total = cbl.count();
    println!("Total number of k-mers: {}", kmer_total);
    println!("Processing to benchmark");
    for bs in 1..=batch_size {
        println!("Batch size {}", bs);
        if mode == "union" {
            let start = Instant::now();
            let merged_cbl = merge_cbls_in_batches(&file_paths, bs);
            let duration = start.elapsed();
            println!("Time taken for union of batch size {}: {:?}", bs, duration);
            let duration_secs = duration.as_secs_f64();
            let kmer_per_sec = kmer_total as f64 / duration_secs;
            println! {"K-mers per second: {}", kmer_per_sec}
            println!("Total number of k-mers in result: {}", merged_cbl.count());
        } else if mode == "intersection" {
            let start = Instant::now();
            let inter_cbl = intersect_cbls_in_batches(&file_paths, bs);
            let duration = start.elapsed();
            println!(
                "Time taken for intersection of batch size {}: {:?}",
                bs, duration
            );
            let duration_secs = duration.as_secs_f64();
            let kmer_per_sec = kmer_total as f64 / duration_secs;
            println! {"K-mers per second: {} ", kmer_per_sec}
            println!("Total number of k-mers in result: {}", inter_cbl.count());
        }
    }
}
