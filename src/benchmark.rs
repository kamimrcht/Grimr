#![feature(generic_const_exprs)]

use anti_reindeer::utils::{
    create_cbl_from_fasta, deserialize_cbl, read_fof_file_csv, serialize_cbl,
};
use std::env;
use std::fs;
use std::path::Path;
use std::time::Instant;
use std::io;

use cbl::CBL;

use serde::Serialize;
#[derive(Serialize)]
struct BenchmarkData {
    batch_size: usize,
    time: f64,
    nb_input_kmer: usize,
    nb_output_kmer: usize,
    input_kmer_per_sec: f64,
    output_kmer_per_sec: f64,
}
#[derive(Serialize)]
struct SerializationData {
    filename: String,
    time: f64,
}

#[derive(Serialize)]
struct DeserializationData {
    time: f64,
    nb_kmers: usize,
}

type T = u64;
const K: usize = 21;

#[derive(Serialize)]
struct GlobalIOs {
    ser_time: f64,
    deser_time: f64,
    nb_kmers: usize,
}

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

fn difference_cbls_in_batches(global_cbl: &mut CBL<K, T>, input_filenames: &[String], batch_size: usize) {
        let mut input_iter = input_filenames.chunks(batch_size);
        let mut local_cbl = if let Some(input_filename_chunk) = input_iter.next() {
            let mut cbls_chunk: Vec<_> = input_filename_chunk
                .iter()
                .map(|input_filename| deserialize_cbl(input_filename))
                .collect();
            CBL::<K, T>::intersect(cbls_chunk.iter_mut().collect())
        } else {
            unreachable!()
        };

        for input_filename_chunk in input_iter {
            let mut cbls_chunk: Vec<_> = input_filename_chunk
                .iter()
                .map(|input_filename| deserialize_cbl(input_filename))
                .collect();
            local_cbl &= &mut CBL::<K, T>::intersect(cbls_chunk.iter_mut().collect());
        }

        *global_cbl -= &mut local_cbl;
}

fn write_csv<T: Serialize>(data: &[T], filename: &str) -> io::Result<()> {
    let mut wtr = csv::Writer::from_path(filename)?;
    for record in data {
        wtr.serialize(record)?;
    }
    wtr.flush()?;
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: {} <input_file_list> <mode> <max_batch_size> [deser]", args[0]);
        std::process::exit(1);
    }
    let input_file_list = &args[1];
    let mode = &args[2];
    let max_batch_size: usize = args[3].parse().expect("Max batch size must be a number");
    let do_deserialize = args.get(4).map_or(false, |arg| arg == "deser");

    if mode != "union" && mode != "intersection" && mode != "difference" {
        eprintln!("Invalid mode. Use 'union', 'intersection', or 'difference'.");
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

	let mut global_cbl = CBL::<K, T>::new();


    let mut serialization_data: Vec<SerializationData> = Vec::new();
    for (i, input_filename) in to_load.iter().enumerate() {
        let output_filename = format!("{}/{}b.cbl", output_dir, indices[i]);
        let output_file_path = Path::new(&output_filename);

        if !output_file_path.exists() {
            let cbl = create_cbl_from_fasta(input_filename);
            let start_serialize = Instant::now();
            serialize_cbl(&cbl, &output_filename);
            let duration_serialize = start_serialize.elapsed().as_secs_f64();
            serialization_data.push(SerializationData {
                filename: input_filename.clone(),
                time: duration_serialize,
            });
            println!("Created CBL for file: {} in {:?}", input_filename, duration_serialize);
        } else {
            println!("CBL already exists for file: {}", input_filename);
        }

        file_paths.push(output_filename);
    }

    println!("Calculating total number of kmers...");
    let total_kmer_cbl = merge_cbls_in_batches(&file_paths, file_paths.len());
    let total_kmers = total_kmer_cbl.count();
    println!("Total number of kmers: {}", total_kmers);

    let mut benchmark_data: Vec<BenchmarkData> = Vec::new();

    for batch_size in 1..=max_batch_size {
        println!("Testing batch size {}", batch_size);
        let start = Instant::now();
        let result_cbl = match mode.as_str() {
            "union" => merge_cbls_in_batches(&file_paths, batch_size),
            "intersection" => intersect_cbls_in_batches(&file_paths, batch_size),
            "difference" => {
                let mut global_cbl = total_kmer_cbl.clone();
                difference_cbls_in_batches(&mut global_cbl, &file_paths, batch_size);
                global_cbl
            }
            _ => unreachable!(),
        };
        let duration = start.elapsed();
        let duration_secs = duration.as_secs_f64();
        let input_kmer_per_sec = total_kmers as f64 / duration_secs;
        let output_kmer_total = result_cbl.count();
        let output_kmer_per_sec = output_kmer_total as f64 / duration_secs;

        benchmark_data.push(BenchmarkData {
            batch_size,
            time: duration_secs,
            nb_input_kmer: total_kmers,
            nb_output_kmer: output_kmer_total,
            input_kmer_per_sec,
            output_kmer_per_sec,
        });

        println!(
            "Time taken for {} of batch size {}: {:?}, input_kmers: {}, output_kmers: {}, input_kmers/sec: {}, output_kmers/sec: {}",
            mode, batch_size, duration, total_kmers, output_kmer_total, input_kmer_per_sec, output_kmer_per_sec
        );
        global_cbl = result_cbl.clone();
    }

    let csv_filename = format!("{}_benchmark.csv", mode);
    write_csv(&benchmark_data, &csv_filename).expect("Failed to write CSV file");

	if !serialization_data.is_empty() {
        write_csv(&serialization_data, "serialization_times.csv").expect("Failed to write serialization CSV file");
    }
   if do_deserialize {
        println!("Benchmarking deserialization of the global dataset...");
        let global_cbl_filename = format!("{}/global.cbl", output_dir);
        
        let start_serialize_global = Instant::now();
        serialize_cbl(&global_cbl, &global_cbl_filename);
        let duration_serialize_global = start_serialize_global.elapsed().as_secs_f64();

        let start_deserialize = Instant::now();
        let deserialized_cbl = deserialize_cbl(&global_cbl_filename);
        let duration_deserialize = start_deserialize.elapsed().as_secs_f64();
        
        let global_data = vec![GlobalIOs {
            ser_time: duration_serialize_global,
            deser_time: duration_deserialize,
            nb_kmers: deserialized_cbl.count(),
        }];
        
        println!("Serialized the global dataset in {:?}", duration_serialize_global);
        println!("Deserialized the global dataset in {:?}", duration_deserialize);

        write_csv(&global_data, "global_serialization_deserialization_times.csv").expect("Failed to write global serialization/deserialization CSV file");
    }
}
