#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![allow(clippy::type_complexity)]

use bincode::{DefaultOptions, Options};
use cbl::kmer::Kmer;
use cbl::CBL;
use needletail::parse_fastx_file;
use serde_json::from_str;
use std::collections::HashSet;
use std::convert::TryInto;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

type T = u64;
const K: usize = 21;

// parse labels and obtain files for all, any, not all, not any
fn parse_label_file<P: AsRef<Path>>(
    path: P,
) -> io::Result<(Vec<i32>, Vec<Vec<i32>>, Vec<Vec<i32>>, Vec<i32>)> {
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);

    let mut vec_all = Vec::new();
    let mut vec_any = Vec::new();
    let mut vec_not_all = Vec::new();
    let mut vec_not_any = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() != 3 {
            continue; // Skip lines that do not conform to the expected format.
        }

        let typ = parts[1];
        let data_str = parts[2];

        match typ {
            "ALL" | "NOT-ANY" => {
                let vec: Vec<i32> = from_str(data_str)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
                if typ == "ALL" {
                    vec_all.extend(vec);
                } else {
                    vec_not_any.extend(vec);
                }
            }
            "ANY" | "NOT-ALL" => {
                let vec_of_vec: Vec<Vec<i32>> = from_str(data_str)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
                if typ == "ANY" {
                    vec_any.extend(vec_of_vec);
                } else {
                    vec_not_all.extend(vec_of_vec);
                }
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Unknown type encountered",
                ))
            }
        }
    }

    Ok((vec_all, vec_any, vec_not_all, vec_not_any))
}

// reads fastas from a file of file (csv), only needed files are processed
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

// select files necessary to load in cbls and serialize
fn select_files_to_load(
    input_files: &[String],
    a_cup: &[i32],
    b_star: &[Vec<i32>],
    c_star: &[Vec<i32>],
    d_cup: &[i32],
) -> Vec<String> {
    let mut to_load = Vec::new();

    if a_cup.is_empty() && b_star.is_empty() {
        // Corrected logical AND
        // a_cup and b_star are empty, load everything
        to_load.extend(input_files.iter().cloned());
    } else {
        // load only necessary datasets
        let mut to_load_id = create_unique_vec(
            a_cup.to_vec(),
            b_star.to_vec(),
            c_star.to_vec(),
            d_cup.to_vec(),
        );
        for id in to_load_id {
            if id < input_files.len() as i32 {
                // Make sure the ID is a valid index
                let input_filename = &input_files[id as usize];
                to_load.push(input_filename.to_string());
            }
        }
    }
    to_load
}

// create and serialize cbls
fn create_and_serialize_cbls(
    input_files: Vec<String>,
    output_dir: &str,
    a_cup: Vec<i32>,
    b_star: Vec<Vec<i32>>,
    c_star: Vec<Vec<i32>>,
    d_cup: Vec<i32>,
) {
    // dir where serialized cbls are stored
    fs::create_dir_all(output_dir).unwrap();

    //  create cbls only if needed (all if a, b empty, else, only indexes that appear)
    let to_load = select_files_to_load(&input_files, &a_cup, &b_star, &c_star, &d_cup);
    if to_load.is_empty() {}
    for (i, input_filename) in to_load.iter().enumerate() {
        let mut reader = parse_fastx_file(input_filename).unwrap();
        let mut cbl = CBL::<K, T>::new();
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid record");
            cbl.insert_seq(&seqrec.seq());
        }

        // serialize the cbl and save it to a file
        let output_filename = format!("{}/{}.cbl", output_dir, i);
        let output = File::create(output_filename).unwrap();
        let mut writer = BufWriter::new(output);
        DefaultOptions::new()
            .with_varint_encoding()
            .reject_trailing_bytes()
            .serialize_into(&mut writer, &cbl)
            .unwrap();
    }
}

// deserialize a given CBL
fn deserialize_cbl(input_index: i32, output_dir: &str) -> CBL<K, T> {
    let input_filename = format!("{}/{}.cbl", output_dir, input_index);
    let index = File::open(input_filename.as_str())
        .unwrap_or_else(|_| panic!("Failed to open {}", input_filename.as_str()));
    let reader = BufReader::new(index);
    DefaultOptions::new()
        .with_varint_encoding()
        .reject_trailing_bytes()
        .deserialize_from(reader)
        .unwrap()
}

// smallest vec for b_star
fn find_smallest_vec_and_index(list_of_vecs: &[Vec<i32>]) -> (usize, Vec<i32>) {
    let mut smallest_index = 0;
    let mut smallest_vec = Vec::new();
    let mut smallest_len = usize::MAX;

    for (index, vec) in list_of_vecs.iter().enumerate() {
        if !vec.is_empty() && vec.len() < smallest_len {
            smallest_index = index;
            smallest_vec = vec.clone();
            smallest_len = vec.len();
        }
    }

    // if smallest_vec remains empty, there were no non-empty vectors in the input
    if smallest_vec.is_empty() {
        return (0, Vec::new());
    }

    (smallest_index, smallest_vec)
}

fn create_unique_vec(
    a_cup: Vec<i32>,
    b_star: Vec<Vec<i32>>,
    c_star: Vec<Vec<i32>>,
    d_cup: Vec<i32>,
) -> Vec<i32> {
    let mut set = HashSet::new();
    for num in a_cup {
        set.insert(num);
    }
    for vec in b_star {
        for num in vec {
            set.insert(num);
        }
    }
    for vec in c_star {
        for num in vec {
            set.insert(num);
        }
    }
    for num in d_cup {
        set.insert(num);
    }
    let mut result: Vec<i32> = set.into_iter().collect();
    result
}

fn query_cbls(
    a_cup: Vec<i32>,
    b_star: Vec<Vec<i32>>,
    c_star: Vec<Vec<i32>>,
    d_cup: Vec<i32>,
    output_dir: &str,
) -> CBL<K, T> {
    // load cbls and build union for A cup and D cup

    let mut global_cbl = CBL::<K, T>::new();
    let mut returned_cbl = CBL::<K, T>::new();

    let mut b_star_work = b_star.clone();
    if a_cup.is_empty() {
        if b_star.is_empty() {
            // union all
            let all_datasets =
                create_unique_vec(a_cup.clone(), b_star.clone(), c_star.clone(), d_cup.clone());
            for index in all_datasets {
                let mut cbl_i = deserialize_cbl(index, output_dir);
                global_cbl |= &mut cbl_i;
            }
        } else {
            let (ind, smallest_vec_b) = find_smallest_vec_and_index(&b_star_work);
            b_star_work.remove(ind); // remove smallest b from b*
            let mut local_cbl = CBL::<K, T>::new();
            for (i, index) in smallest_vec_b.iter().enumerate() {
                if i == 0 {
                    local_cbl = deserialize_cbl(ind.try_into().unwrap(), output_dir);
                    // todo vérifier avec Florian ligne 11 algo
                } else {
                    let mut cbl_b = deserialize_cbl(*index, output_dir);
                    local_cbl |= &mut cbl_b;
                }
            }
            global_cbl |= &mut local_cbl; // clone
            let c = global_cbl.count();
        }
    } else {
        let mut local_cbl = CBL::<K, T>::new();
        for (i, index) in a_cup.iter().enumerate() {
            if i == 0 {
                local_cbl = deserialize_cbl(*index, output_dir);
            } else {
                let mut cbl_a = deserialize_cbl(*index, output_dir);
                local_cbl &= &mut cbl_a;
            }
        }
        global_cbl |= &mut local_cbl; // clone
    }
    for c in &c_star {
        if !c.is_empty() {
            // todo voir avec florian
            let mut local_cbl = CBL::<K, T>::new();
            local_cbl |= &mut global_cbl;
            for index in c {
                let mut cbl_c = deserialize_cbl(*index, output_dir);
                local_cbl &= &mut cbl_c;
            }
            global_cbl -= &mut local_cbl;
            let count = global_cbl.count();
        }
    }
    for index in &d_cup {
        let mut cbl_d = deserialize_cbl(*index, output_dir);
        global_cbl -= &mut cbl_d;
        let count = global_cbl.count();
    }
    if b_star_work.iter().all(|vec| vec.is_empty()) {
        //check if the rest of B star is only empty vecs
        returned_cbl |= &mut global_cbl;
    } else {
        for b in b_star_work {
            // todo voir avec Florian ligne 32, si vide
            if !b.is_empty() {
                let mut local_cbl = CBL::<K, T>::new();
                for index in b {
                    let mut cbl_b = deserialize_cbl(index, output_dir);
                    let mut cbl_tmp = CBL::<K, T>::new();
                    cbl_tmp |= &mut global_cbl; // "clone"
                    cbl_tmp &= &mut cbl_b;
                    local_cbl |= &mut cbl_tmp;
                }
                returned_cbl |= &mut local_cbl;
            }
        }
    }
    returned_cbl
}

fn cbl_printer(cbl: &CBL<K, T>, output_path: &str) -> std::io::Result<()> {
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

fn main() {
    // parse args
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: {} <mode> <input_metadata>", args[0]);
        std::process::exit(1);
    }
    let mode = args[1].clone();
    let input_file_list = args[2].clone();
    let label_file_list = args[3].clone();
    let output_dir = "serialized_cbls";
    let labels = parse_label_file(label_file_list).unwrap(); //todo test
    let (a_cup, b_star, c_star, d_cup) = labels;
    if mode == "index" {
        // read the fof
        let (input_files, _col_nb) = read_fof_file_csv(&input_file_list).unwrap(); // use of col_nb?
                                                                                   // create and serialize CBLs
        create_and_serialize_cbls(input_files, output_dir, a_cup, b_star, c_star, d_cup);
    } else if mode == "query" {
        let cbl = query_cbls(a_cup, b_star, c_star, d_cup, output_dir);
        let output_path = "output_anti_reindeer_query.txt";
        let _ = fs::remove_file(output_path);
        cbl_printer(&cbl, output_path).expect("Failed to print CBL");
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
        assert_eq!(parse_line("b	ANY	[[],[]]"), vec![]);
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

    #[test]
    fn test_smallest_empty() {
        let list_of_vecs = &[vec![], vec![3, 4], vec![5, 6], vec![7]];
        let (index, smallest_vec) = find_smallest_vec_and_index(list_of_vecs);
        assert_eq!(index, 3);
        assert_eq!(smallest_vec, vec![7]);
    }

    #[test]
    fn test_create_unique_vec1() {
        let a_cup = vec![1, 2, 3];
        let b_star = vec![vec![4, 5], vec![6]];
        let c_star = vec![vec![7, 8], vec![9]];
        let d_cup = vec![10, 11];
        let mut result = create_unique_vec(a_cup, b_star, c_star, d_cup);
        result.sort();
        assert_eq!(result, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
    }
    #[test]
    fn test_create_unique_vec2() {
        let a_cup = vec![1, 2, 3];
        let b_star = vec![];
        let c_star = vec![vec![4, 5], vec![6]];
        let d_cup = vec![7, 8];
        let mut result = create_unique_vec(a_cup, b_star, c_star, d_cup);
        result.sort();
        assert_eq!(result, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }
    #[test]
    fn test_create_unique_vec3() {
        let a_cup = vec![1, 2, 2];
        let b_star = vec![vec![2], vec![3, 4]];
        let c_star = vec![vec![3, 4, 4], vec![5]];
        let d_cup = vec![1, 6];
        let mut result = create_unique_vec(a_cup, b_star, c_star, d_cup);
        result.sort();
        assert_eq!(result, vec![1, 2, 3, 4, 5, 6]);
    }
    #[test]
    fn test_create_unique_vec4() {
        let a_cup = vec![];
        let b_star = vec![vec![0], vec![]];
        let c_star = vec![vec![], vec![]];
        let d_cup = vec![];
        let mut result = create_unique_vec(a_cup, b_star, c_star, d_cup);
        result.sort();
        assert_eq!(result, vec![0]);
    }
    #[test]
    fn test_printer() {
        let cbl_a = deserialize_cbl(0, "test_files");
        let output_path = "test_files/test_printer.fa";
        let _ = fs::remove_file(output_path);
        cbl_printer(&cbl_a, output_path).unwrap();
        assert!(Path::new(output_path).exists(), "The file was not created.");
        let file = File::open(output_path).expect("Failed to open file");
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let expected_lines = [
            ">kmer0",
            "CTAAAAAACCGTCAATGTGAA",
            ">kmer1",
            "TAAAAAACCGTCAATGTGAAA",
            ">kmer2",
            "AAACTGGAACGGTTAGAGAAA",
            ">kmer3",
            "AAAAAGACGGACAAGAAGCGA",
            ">kmer4",
            "TGCAGTTAAAAAGCTCGTAGT",
            ">kmer5",
            "GCAGTTAAAAAGCTCGTAGTT",
        ];
        for expected_line in expected_lines.iter() {
            if let Some(Ok(actual_line)) = lines.next() {
                assert_eq!(
                    actual_line, *expected_line,
                    "Mismatch in file content at expected line: {}",
                    expected_line
                );
            } else {
                panic!("File has fewer lines than expected");
            }
        }
    }
    #[test]
    fn test_load_all_files_1empty() {
        let input_files = vec![
            "file1.txt".to_string(),
            "file2.txt".to_string(),
            "file3.txt".to_string(),
        ];
        let a_cup: Vec<i32> = vec![];
        let b_star: Vec<Vec<i32>> = vec![];
        let c_star: Vec<Vec<i32>> = vec![vec![1]];
        let d_cup: Vec<i32> = vec![1];

        let loaded_files = select_files_to_load(&input_files, &a_cup, &b_star, &c_star, &d_cup);
        assert_eq!(loaded_files, input_files);
    }

    #[test]
    fn test_load_all_files2() {
        let input_files = vec![
            "file1.txt".to_string(),
            "file2.txt".to_string(),
            "file3.txt".to_string(),
        ];
        let a_cup: Vec<i32> = vec![1];
        let b_star: Vec<Vec<i32>> = vec![vec![2]];
        let c_star: Vec<Vec<i32>> = vec![vec![1]];
        let d_cup: Vec<i32> = vec![1];

        let mut loaded_files = select_files_to_load(&input_files, &a_cup, &b_star, &c_star, &d_cup);
        loaded_files.sort();
        assert_eq!(
            loaded_files,
            vec!["file2.txt".to_string(), "file3.txt".to_string()]
        );
    }
}
