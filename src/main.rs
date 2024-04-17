#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![allow(clippy::type_complexity)]

use bincode::{DefaultOptions, Options};
use cbl::kmer::Kmer;
use cbl::CBL;
use needletail::parse_fastx_file;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

type T = u64;
const K: usize = 21;

// parse input query file
fn parse_line(line: &str) -> Vec<i32> {
    line.split('\t')
        .last()
        .unwrap_or("")
        .trim_matches(|p| p == '[' || p == ']')
        .split(',')
        .filter_map(|s| s.trim().parse::<i32>().ok())
        .collect()
}

fn parse_list_of_lists(line: &str) -> Vec<Vec<i32>> {
    line.split('\t')
        .last()
        .unwrap_or("")
        .trim_matches(|p| p == '[' || p == ']')
        .split("],[")
        .map(parse_line)
        .collect()
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
            println!("Warning: a putatively large union operation is needed for this query.");
            // union all
            // todo correct
            for index in a_cup {
                let mut cbl_a = deserialize_cbl(index, output_dir);
                global_cbl |= &mut cbl_a;
            }
            for c in &c_star {
                for index in c {
                    let mut cbl_c = deserialize_cbl(*index, output_dir);
                    global_cbl |= &mut cbl_c;
                }
            }
            for index in &d_cup {
                let mut cbl_d = deserialize_cbl(*index, output_dir);
                global_cbl |= &mut cbl_d;
            }
        } else {
			println!("here");
            let (ind, smallest_vec_b) = find_smallest_vec_and_index(&b_star_work);
            println!("{}", ind);
            b_star_work.remove(ind); // remove smallest b from b*
            let mut local_cbl = CBL::<K, T>::new();
            for index in smallest_vec_b {
				println!("{}", index);
                let mut cbl_b = deserialize_cbl(index, output_dir);
                local_cbl |= &mut cbl_b;
            }
            global_cbl |= &mut local_cbl; // clone
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
        let mut local_cbl = CBL::<K, T>::new();
        local_cbl |= &mut global_cbl;
        for index in c {
            let mut cbl_c = deserialize_cbl(*index, output_dir);
            local_cbl &= &mut cbl_c;
        }
        global_cbl -= &mut local_cbl;
    }
    for index in &d_cup {
        let mut cbl_d = deserialize_cbl(*index, output_dir);
        global_cbl -= &mut cbl_d;
    }
    for b in b_star_work {
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
        let (a_cup, b_star, c_star, d_cup) = labels;
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
}
