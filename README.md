# anti reindeer

```sh
git clone https://github.com/kamimrcht/anti_reindeer.git
cd anti_reindeer
```

```sh
cargo +nightly build --release
cargo +nightly test
```

To install Rust nightly: [rustup.rs](https://rustup.rs/) then `rustup install nightly`.

## Index mode

```sh
cargo +nightly run --release -- index test_files/fof.txt test_files/input_florian.txt
```

## Query mode

```sh
cargo +nightly run --release -- query test_files/fof.txt test_files/input_florian.txt
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
