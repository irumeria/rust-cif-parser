use crate::structure::Structure;
use crate::sym::Symmetry;
use ndarray::{array, Array2, Axis};
use ndarray::{stack, Array1, ArrayView1, Zip};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

fn in_coord_list_pbc(fcoord_list: &Vec<[f64; 3]>, fcoord: Array1<f64>, atol: f64) -> bool {
    if fcoord_list.is_empty() {
        return false;
    }

    let fcoords_views: Vec<ArrayView1<_>> = vec![fcoord.view(); fcoord_list.len()];
    let fcoords = stack(Axis(0), &fcoords_views).unwrap();
    let mut fdist = Array2::from_shape_vec(
        (fcoord_list.len(), 3),
        fcoord_list.clone().into_iter().flatten().collect(),
    )
    .unwrap();

    Zip::from(&mut fdist)
        .and(&fcoords)
        .for_each(|a, &b| *a -= b);

    fdist.mapv_inplace(|x| x - x.round());
    // println!("fdisk {:?}", fdist);

    let mask = fdist.mapv(|x| x.abs() < atol);
    let ret = mask.axis_iter(Axis(0))
        .any(|row| row.iter().all(|&x| x));
    // println!("mask {:?}", mask);
    // println!("ret {:?}", ret);
    ret
}

fn restore_equiv_atoms(
    symops: &Vec<Symmetry>,
    positions: &Vec<[f64; 3]>,
    atol: f64,
) -> (Vec<[f64; 3]>, Vec<usize>) {
    let mut coords_out: Vec<[f64; 3]> = Vec::new();
    let mut mapping_to_primitive: Vec<usize> = Vec::new();
    for (i, tmp_coord) in positions.iter().enumerate() {
        for op in symops {
            let mut coord = op.operate(array![tmp_coord[0], tmp_coord[1], tmp_coord[2]]);
            coord.mapv_inplace(|x| x - x.floor());
            if !in_coord_list_pbc(&coords_out, coord.clone(), atol) {
                coords_out.push(coord.to_vec().as_slice().try_into().unwrap());
                mapping_to_primitive.push(i);
            }
        }
    }


    (coords_out, mapping_to_primitive)
}

fn parse_lattice_line(line: &str) -> f64 {
    let iter = line.split_whitespace();
    let last_word = iter.last().unwrap();
    let mut split = last_word.split('(');
    let value: f64 = split.next().unwrap().parse().unwrap();
    let _unit: i32 = if last_word.contains('(') {
        split.next().unwrap().trim_end_matches(')').parse().unwrap()
    } else {
        0
    };
    // value * 10f64.powi(unit)
    value
}


pub fn parse_cif(cif_path: &Path) -> Structure {
    let cif_file = File::open(cif_path).unwrap();
    let cif_reader = io::BufReader::new(cif_file);
    let mut abc = [0.0, 0.0, 0.0];
    let mut angles = [0.0, 0.0, 0.0];
    let mut positions = vec![];
    let mut primitive_atom_symbols = vec![];
    let mut atom_site_flag = false;
    let mut atom_site_p:[usize; 4] = [0, 0, 0, 0];
    let mut atom_site_p_counter = 0;
    let mut xyz_strs = vec![];
    let mut xyz_flag = false;
    let mut xyz_p = [0];
    let mut xyz_p_counter = 0;
    for line in cif_reader.lines() {
        let line = line.unwrap();
        let line = line.trim();
        if line.starts_with("_cell_length_a") {
            abc[0] = parse_lattice_line(&line);
        } else if line.starts_with("_cell_length_b") {
            abc[1] = parse_lattice_line(&line);
        } else if line.starts_with("_cell_length_c") {
            abc[2] = parse_lattice_line(&line);
        } else if line.starts_with("_cell_angle_alpha") {
            angles[0] = parse_lattice_line(&line);
        } else if line.starts_with("_cell_angle_beta") {
            angles[1] = parse_lattice_line(&line);
        } else if line.starts_with("_cell_angle_gamma") {
            angles[2] = parse_lattice_line(&line);
        } else if line.trim().starts_with("_atom_site"){
            atom_site_flag = true;
        } else if line.trim().starts_with("_symmetry_equiv") {
            xyz_flag = true;
        }
        if xyz_flag {
            if line.is_empty() || line.starts_with("loop_") {
                xyz_flag = false;
                continue;
            }
            let tiled_trimed = line.trim().to_owned();
            if tiled_trimed.starts_with("_symmetry_equiv"){
                if tiled_trimed== "_symmetry_equiv_pos_as_xyz" {
                    xyz_p[0] = xyz_p_counter;
                }
                xyz_p_counter += 1;
                continue;
            }
            let line_parts = line.split("'").collect::<Vec<&str>>();
            
            let xyz_str = line_parts[1].to_owned(); // todo: make sure of it
            xyz_strs.push(xyz_str);
        }
        if atom_site_flag {
            // println!("{:?}", line);
            if line.is_empty() {
                break;
            }
            let tiled_trimed = line.trim().to_owned();
            if tiled_trimed.starts_with("_atom_site"){
                if tiled_trimed== "_atom_site_type_symbol" {
                    atom_site_p[0] = atom_site_p_counter;
                } else if tiled_trimed == "_atom_site_fract_x" {
                    atom_site_p[1] = atom_site_p_counter;
                } else if tiled_trimed == "_atom_site_fract_y" {
                    atom_site_p[2] = atom_site_p_counter;
                } else if tiled_trimed == "_atom_site_fract_z" {
                    atom_site_p[3] = atom_site_p_counter;
                }
                atom_site_p_counter += 1;
                continue;
            }
            let line_parts = line.split_whitespace().collect::<Vec<&str>>();
            let atom_symbol = line_parts[atom_site_p[0]].to_owned();
            let pos_vec: Vec<f64> = vec![
                line_parts[atom_site_p[1]].parse().unwrap(),
                line_parts[atom_site_p[2]].parse().unwrap(),
                line_parts[atom_site_p[3]].parse().unwrap(),
            ];
            positions.push([pos_vec[0], pos_vec[1], pos_vec[2]]);
            primitive_atom_symbols.push(atom_symbol);
        }
    }

    let symops = xyz_strs
        .iter()
        .map(|x| Symmetry::from_xyz_str(x))
        .collect::<Vec<Symmetry>>();
    let lattice_ndarray = Structure::cellpar_to_cell(
        &vec![abc[0], abc[1], abc[2], angles[0], angles[1], angles[2]],
        array![0.0, 0.0, 1.0],
        None,
    );
    let mut mapping_atomtype_to_index:HashMap<&String, Vec<usize>> = HashMap::new();
    for (i, atom_symbol) in primitive_atom_symbols.iter().enumerate() {
        let indices = mapping_atomtype_to_index.entry(atom_symbol).or_insert(vec![]);
        indices.push(i);
    }
    let mut atom_coords = vec![];
    let mut atom_symbols = vec![];
    for (_atom_symbol, indices) in mapping_atomtype_to_index.iter() {
        let positions_of_atomtype: Vec<_> = indices.iter().map(|&i| positions[i]).collect();
        let (_atom_coords, mapping_to_primitive) = restore_equiv_atoms(&symops, &positions_of_atomtype, 1e-3);
        for (i, atom_coord) in _atom_coords.iter().enumerate() {
            let index = indices[mapping_to_primitive[i]];
            atom_coords.push(atom_coord.clone());
            atom_symbols.push(primitive_atom_symbols[index].clone());
        }
    }
    let atom_numbers: Vec<u32> = atom_symbols.iter().map(|x| periodic_table(x)).collect();
    let atom_coords_ndarray = Array2::from_shape_vec(
        (atom_coords.len(), 3),
        atom_coords.iter().flatten().cloned().collect(),
    ).unwrap();
    let cart_coords = atom_coords_ndarray.dot(&lattice_ndarray); // todo: make sure of it
    
    let structure = Structure {
        lattice: lattice_ndarray,
        lattice_parameters: [abc[0], abc[1], abc[2], angles[0], angles[1], angles[2]],
        atom_numbers: atom_numbers,
        positions: cart_coords ,
        frac_positions: atom_coords_ndarray,
        atom_symbols: atom_symbols,
    };
    // println!("{:?}", structure);
    structure
}

pub fn periodic_table(atom_symbol: &str) -> u32 {
    let pt = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
        "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
        "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
        "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
        "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
        "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    ];
    let mut index = 0;
    for (i, symbol) in pt.iter().enumerate() {
        if symbol == &atom_symbol {
            index = i;
            break;
        }
    }
    (index + 1) as u32 // +1 because the atomic number starts from 1
}
