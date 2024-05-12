extern crate regex;
extern crate ndarray;
use regex::Regex;
use ndarray::{s, Array1, Array2};

#[derive(Debug, Clone)]
pub struct Symmetry {
    pub rot_matrix: Array2<f64>,
    pub trans_vec: Array1<f64>,
    pub affine_matrix: Array2<f64>,
}

impl Symmetry {
    pub fn from_xyz_str(xyz_str: &str) -> Symmetry {
        let (rot_matrix, trans) = operators_from_xyz_str(xyz_str);
        Symmetry::new(&rot_matrix, &trans)
    }

    pub fn new(rot_matrix: &Array2<f64>, trans_vec: &Array1<f64>) -> Symmetry {
        let mut affine_matrix = Array2::<f64>::eye(4);
        for i in 0..3 {
            for j in 0..3 {
                affine_matrix[[i, j]] = rot_matrix[[i, j]];
            }
            affine_matrix[[i, 3]] = trans_vec[i];
        }
        Symmetry {
            rot_matrix: rot_matrix.clone(),
            trans_vec: trans_vec.clone(),
            affine_matrix: affine_matrix,
        }
    }
    
    pub fn operate(&self, point: Array1<f64>) -> Array1<f64> {
        let mut affine_point = Array1::<f64>::zeros(4);
        for i in 0..3 {
            affine_point[i] = point[i];
        }
        affine_point[3] = 1.0;
        let result = self.affine_matrix.dot(&affine_point);
        result.slice(s![0..3]).to_owned()
    }
}

pub fn operators_from_xyz_str(xyz_str: &str) -> (Array2<f64>, Array1<f64>) {
    let mut rot_matrix = Array2::<f64>::zeros((3, 3));
    let mut trans = Array1::<f64>::zeros(3);
    let binding = xyz_str.trim().replace(" ", "").to_lowercase();
    let tokens: Vec<&str> = binding.split(",").collect();
    let re_rot = Regex::new(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])").unwrap();
    let re_trans = Regex::new(r"([+-]?)([\d\.]*)/?([\d\.]*)([+-]?[x-z])").unwrap();

    for (i, tok) in tokens.iter().enumerate() {
        for cap in re_rot.captures_iter(tok) {
            let mut factor = if &cap[1] == "-" { -1.0 } else { 1.0 };
            if &cap[2] != "" {
                factor *= cap[2].parse::<f64>().unwrap_or(1.0) / cap[3].parse::<f64>().unwrap_or(1.0);
            }
            let j = cap[4].chars().next().unwrap() as usize - 120;
            rot_matrix[[i, j]] = factor;
        }
        for cap in re_trans.captures_iter(tok) {
            if cap.get(4).is_some() {
                let factor = if &cap[1] == "-" { -1.0 } else { 1.0 };
                let num = cap[2].parse::<f64>().unwrap_or(0.0) / cap[3].parse::<f64>().unwrap_or(1.0);
                trans[i] += num * factor;
            }
        }
    }
    (rot_matrix, trans)
}
