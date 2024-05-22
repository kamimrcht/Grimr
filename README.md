# anti reindeer

```sh
git clone https://github.com/kamimrcht/anti_reindeer.git
cd anti_reindeer
```

Install Rust nightly: [rustup.rs](https://rustup.rs/) then `rustup install nightly`.

```sh
cargo +nightly build --release
cargo +nightly test
```
## Creating queries with Python

Make sure you have a recent version of Python on your machine.

In `label_editor.py`, at line 83, make sure to indicate the proper metadata file. If you wish to do so, the file also provide tools to create new tags from the ones present in the metadata. Then, run the script, or use command line:

```sh
python label_editor.py
```

This will generate a pickle file containing a dictionary associating each tag to its list of indexes of references, as well as a `map.txt` file containing the association between the indexes and the references themselves (in the same order as the metadate file).

Then, in `query_builder.py`, after line 62, write your query using pre-existing tags. Then run the script or use command line:

```sh
python query_builder.py
```
This will generate a file `query.txt` that you can use with the subsequent Rust part of this tool.

## Index mode

```sh
cd test_files
```

```sh
cargo +nightly run --bin anti_reindeer --release -- index test_files/metadata.csv test_files/query2.txt
```

## Query mode

```sh
cd test_files
```

```sh
cargo +nightly run --bin anti_reindeer --release -- query test_files/metadata.csv test_files/query2.txt
```

## Useful commands

Update Rust:
```sh
rustup update
```

Update CBL:
```sh
cargo update cbl
```

Format code:
```sh
cargo fmt
```

Show clippy warnings:
```sh
cargo clippy
```

Fix warnings automatically:
```sh
cargo clippy --fix
```
