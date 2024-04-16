#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![allow(clippy::type_complexity)]

use bincode::{deserialize, serialize};
use cbl::CBL;
use needletail::parse_fastx_file;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write, Result};
use std::path::Path;
use cbl::kmer::Kmer;

type T = u64;
const K: usize = 21;

// parse input query file
fn parse_line(line: &str) -> Vec<i32> {
    line.split('\t').last().unwrap_or("").trim_matches(|p| p == '[' || p == ']').split(',').filter_map(|s| s.trim().parse::<i32>().ok()).collect() 
}

fn parse_list_of_lists(line: &str) -> Vec<Vec<i32>> {
    line.split('\t').last().unwrap_or("").trim_matches(|p| p == '[' || p == ']').split("],[").map(parse_line).collect() 
}

fn parse_file<P: AsRef<Path>>(
    path: P,
) -> io::Result<(Vec<i32>, Vec<Vec<i32>>, Vec<Vec<i32>>, Vec<i32>)> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let list1 = lines
        .next()
        .map_or(Ok(vec![]), |l| l.map(|line| parse_line(&line)))?;
    let list_of_lists1 = lines
        .next()
        .map_or(Ok(vec![]), |l| l.map(|line| parse_list_of_lists(&line)))?;
    let list_of_lists2 = lines
        .next()
        .map_or(Ok(vec![]), |l| l.map(|line| parse_list_of_lists(&line)))?;
    let list2 = lines
        .next()
        .map_or(Ok(vec![]), |l| l.map(|line| parse_line(&line)))?;
    Ok((list1, list_of_lists1, list_of_lists2, list2))
}

// reads fastas from a file of file
fn read_fof_file_csv(file_path: &str) -> io::Result<(Vec<String>, usize)> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut file_paths = Vec::new(); 
    let mut color_number = 0;

    for line in reader.lines() {
        let line = line?;
        if let Some(first_column) = line.split_whitespace().next() {
            file_paths.push(first_column.to_string());
            color_number += 1;
        }
    }

    Ok((file_paths, color_number))
}
// create and serialize cbls
fn create_and_serialize_cbls(input_files: Vec<String>, output_dir: &str) {
    // dir where serialized cbls are stored
    fs::create_dir_all(output_dir).unwrap();

    // for each fasta i of the fof, create a cbl
    for (i, input_filename) in input_files.iter().enumerate() {
        let mut reader = parse_fastx_file(input_filename).unwrap();
        let mut cbl = CBL::<K, T>::new();
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid record");
            cbl.insert_seq(&seqrec.seq());
        }

        // serialize the cbl and save it to a file
        let serialized_cbl = serialize(&cbl).unwrap();
        let output_filename = format!("{}/{}.cbl", output_dir, i);
        fs::write(Path::new(&output_filename), &serialized_cbl).unwrap();
    }
}

// deserialize a given CBL
fn deserialize_cbl(input_index: i32, output_dir: &str) -> CBL<K, T> {
    let input_path = format!("{}/{}.cbl", output_dir, input_index);
    let data = fs::read(Path::new(&input_path)).unwrap();
    let deserialized_cbl: CBL<K, T> = deserialize(&data).unwrap();
    deserialized_cbl
}

// smallest vec for Bstar
fn find_smallest_vec_and_index(list_of_vecs: &[Vec<i32>]) -> (usize, Vec<i32>) {
    let mut smallest_index = 0;
    let mut smallest_vec = Vec::new();

    if !list_of_vecs.is_empty() {
        smallest_vec.clone_from(&list_of_vecs[0]);
        let mut smallest_len = smallest_vec.len();

        for (index, vec) in list_of_vecs.iter().enumerate() {
            if vec.len() < smallest_len {
                smallest_index = index;
                smallest_vec.clone_from(vec);
                smallest_len = vec.len();
            }
        }
    }

    (smallest_index, smallest_vec)
}

fn query_cbls(Acup: Vec<i32>, Bstar: Vec<Vec<i32>>, Cstar: Vec<Vec<i32>>, Dcup: Vec<i32>, output_dir: &str,) -> CBL<K, T> {
    // load cbls and build union for A cup and D cup

    let mut global_cbl = CBL::<K, T>::new();
    let mut returned_cbl = CBL::<K, T>::new();

    let mut Bstar_work = Bstar.clone();
    if Acup.is_empty() {
        if Bstar.is_empty() {
            println!("Warning: a putatively large union operation is needed for this query.");
            // union all
            for index in Acup {
                let mut cbl_a = deserialize_cbl(index, output_dir);
                global_cbl |= &mut cbl_a;
            }
            for C in &Cstar {
                for index in C {
                    let mut cbl_c = deserialize_cbl(*index, output_dir);
                    global_cbl |= &mut cbl_c;
                }
            }
            for index in &Dcup {
                let mut cbl_d = deserialize_cbl(*index, output_dir);
                global_cbl |= &mut cbl_d;
            }
        } else {
            let (ind, smallest_vecB) = find_smallest_vec_and_index(&Bstar_work);
            Bstar_work.remove(ind); // remove smallest B from B*
            let mut local_cbl = CBL::<K, T>::new();
            for index in smallest_vecB {
                let mut cbl_b = deserialize_cbl(index, output_dir);
                local_cbl |= &mut cbl_b;
            }
            global_cbl |= &mut local_cbl; // clone
        }
    } else {
        let mut local_cbl = CBL::<K, T>::new();
        for (i, index) in Acup.iter().enumerate() {
            if i == 0 {
                local_cbl = deserialize_cbl(*index, output_dir);
            } else {
                let mut cbl_a = deserialize_cbl(*index, output_dir);
                local_cbl &= &mut cbl_a;
            }
        }
        global_cbl |= &mut local_cbl; // clone
    }
    for C in &Cstar {
        let mut local_cbl = CBL::<K, T>::new();
        local_cbl |= &mut global_cbl;
        for index in C {
            let mut cbl_c = deserialize_cbl(*index, output_dir);
            local_cbl &= &mut cbl_c;
        }
        global_cbl -= &mut local_cbl;
    }
    for index in &Dcup {
        let mut cbl_d = deserialize_cbl(*index, output_dir);
        global_cbl -= &mut cbl_d;
    }
    for B in Bstar_work {
        let mut local_cbl = CBL::<K, T>::new();
        for index in B {
            let mut cbl_b = deserialize_cbl(index, output_dir);
            let mut cbl_tmp = CBL::<K, T>::new();
            cbl_tmp |= &mut global_cbl; // "clone"
            cbl_tmp &= &mut cbl_b;
            local_cbl |= &mut cbl_tmp;
        }
        returned_cbl |= &mut local_cbl;
    }
    returned_cbl
}


fn cbl_printer(cbl: &CBL<K, T>) -> std::io::Result<()> {
    if cbl.is_empty() {
		println!("Empty solution.");
        return Ok(());  
    }

    let file = File::create("output_anti_reindeer.fasta")?;
    let mut writer = BufWriter::new(file);
	for (index, kmer) in cbl.iter().enumerate() {
        let nuc_array = kmer.to_nucs();
        let slice_nuc = nuc_array.as_slice(); 
        let nucs = String::from_utf8_lossy(slice_nuc); // @Igor, meilleure manière de faire?
        writeln!(writer, ">kmer{}", index)?;
        writeln!(writer, "{}", nucs)?;
    }
    Ok(())
}

fn main() {
    // parse args
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: {} <mode> <input_file_list> <labels_list>", args[0]);
        std::process::exit(1);
    }
    let mode = args[1].clone();
    let input_file_list = args[2].clone();
    let label_file_list = args[3].clone();
    let output_dir = "serialized_cbls";

    if mode == "index" {
        // read the fof
        let (input_files, _col_nb) = read_fof_file_csv(&input_file_list).unwrap(); // use of col_nb?
                                                                               // create and serialize CBLs
        create_and_serialize_cbls(input_files, output_dir);
    } else if mode == "query" {
        let labels = parse_file(label_file_list).unwrap();
        let (Acup, Bstar, Cstar, Dcup) = labels;
        let cbl = query_cbls(Acup, Bstar, Cstar, Dcup, output_dir);
        cbl_printer(&cbl);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_parse_line() {
        assert_eq!(parse_line("A	ALL	[1,2,3]"), vec![1, 2, 3]);
        assert_eq!(parse_line("A	ALL	[ 1 , 2 , 3 ]"), vec![1, 2, 3]);
        assert_eq!(parse_line("A	ALL	[]"), vec![]);
        assert_eq!(parse_line("B	ANY	[[],[]]"), vec![]);
        assert_eq!(parse_line("A	ALL	[100]"), vec![100]);
    }

    #[test]
    fn test_read_fof_file_csv() {
        let result = read_fof_file_csv("test_files/metadata.csv").unwrap();
        assert_eq!(
            result.0,
            vec![
                "test_files/test1.fa".to_string(),
                "test_files/test2.fa".to_string(),
                "test_files/test3.fa".to_string(),
                "test_files/test4.fa".to_string(),
                "test_files/test5.fa".to_string(),
                "test_files/test6.fa".to_string()
            ]
        );
        assert_eq!(result.1, 6);
    }

    #[test]
    fn test_smallest1() {
        let list_of_vecs = &[vec![1, 2, 3]];
        let (index, smallest_vec) = find_smallest_vec_and_index(list_of_vecs);
        assert_eq!(index, 0);
        assert_eq!(smallest_vec, vec![1, 2, 3]);
    }

    #[test]
    fn test_smallest2() {
        let list_of_vecs = &[vec![1], vec![2], vec![3]];
        let (index, smallest_vec) = find_smallest_vec_and_index(list_of_vecs);
        assert_eq!(index, 0);
        assert_eq!(smallest_vec, vec![1]);
    }

    #[test]
    fn test_smallest3() {
        let list_of_vecs = &[vec![1, 2, 3, 4], vec![5, 6], vec![7, 8, 9]];
        let (index, smallest_vec) = find_smallest_vec_and_index(list_of_vecs);
        assert_eq!(index, 1);
        assert_eq!(smallest_vec, vec![5, 6]);
    }

    #[test]
    fn test_smallest4() {
        let list_of_vecs = &[vec![1, 2], vec![3, 4], vec![5, 6]];
        let (index, smallest_vec) = find_smallest_vec_and_index(list_of_vecs);
        assert_eq!(index, 0);
        assert_eq!(smallest_vec, vec![1, 2]);
    }
}
