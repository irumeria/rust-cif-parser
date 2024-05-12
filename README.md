## Cif (Crystallographic Information File) Parser of Rust

This library is for parsing the cif file into a **Structure** struct. No grantee for now. Welcome for new contribution or issue.


### Usage

```rust
use cif_parser::cif::parse_cif;

let cif_path = Path::new("assets/ABW-0.cif"); // read the cif 
let structure = parse_cif(cif_path);
```

The structure struct containing the members:
```rust
pub struct Structure {
    pub lattice: Array2<f64>,
    pub lattice_parameters: [f64; 6],
    pub atom_numbers: Vec<u32>,
    pub positions: Array2<f64>,
    pub frac_positions: Array2<f64>,
    pub atom_symbols: Vec<String>,
}
```
Noted that the `Array`s is from the [rust-ndarray](https://github.com/rust-ndarray/ndarray).

### Interact with spglib

We can analysis the spacegroup and symmetric operations by using the rust wapper of spglib [spglib-rs](https://github.com/spglib/spglib-rs).

```rust
extern crate spglib;

fn array2_to_vec(array: Array2<f32>) -> Vec<[f32; 3]> {
    array.map_axis(ndarray::Axis(1), |row| {
        let arr: [f32; 3] = [row[0], row[1], row[2]];
        arr
    }).to_vec()
}

fn main(){

    ...

    let atom_numbers = structure.atom_numbers;
    let atom_numbers_i32: &[i32] =
        &atom_numbers.iter().map(|x| *x as i32).collect::<Vec<i32>>();
    let lattice = structure.lattice;
    let lattice = [
        [lattice[[0, 0]], lattice[[0, 1]], lattice[[0, 2]]],
        [lattice[[1, 0]], lattice[[1, 1]], lattice[[1, 2]]],
        [lattice[[2, 0]], lattice[[2, 1]], lattice[[2, 2]]],
    ];
    let positions = structure.frac_positions;
    let positions = array2_to_vec(positions);
    let dataset = spglib::dataset::Dataset::new(
        &mut spglib::cell::Cell::new(&lattice, &positions, atom_numbers_i32),
        1e-5,
    );
    println!("{:?}", dataset.equivalent_atoms);
    println!("{:?}", dataset.spacegroup_number);
}
```


### Known issue

+ The lines inside `_symmetry_equiv_pos_as_xyz` block **must** be surrounded by `'`, like : `'1/2+x,1/2+y,1/2+z'`

### Acknowledgement

+ Parts of the code is inspired or based on the excellent project [Pymatgen](https://github.com/materialsproject/pymatgen). 