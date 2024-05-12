use crate::cif::parse_cif;
use ndarray::{array, stack, Array1, Array2, Axis};
use ndarray_linalg::Norm;
use std::path::Path;
use crate::sym::Symmetry;

fn unit_vector(vector: &Array1<f64>) -> Array1<f64> {
    vector / vector.norm()
}

pub fn cross(a: &Array1<f64>, b: &Array1<f64>) -> Array1<f64> {
    array![
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ]
}

#[derive(Debug, Clone)]
pub struct Structure {
    pub lattice: Array2<f64>,
    pub lattice_parameters: [f64; 6],
    pub atom_numbers: Vec<u32>,
    pub positions: Array2<f64>,
    pub frac_positions: Array2<f64>,
    pub atom_symbols: Vec<String>,
    pub symops: Vec<Symmetry>,
}

impl Structure {
    pub fn from_cif(cif_path: &Path) -> Structure {
        parse_cif(cif_path)
    }
    
    #[allow(unused_variables, unused_mut, unused_assignments)]
    pub fn to_cif(&self, cif_path: &Path) -> ! {
        todo!("Implement this function")
    }
    
    pub fn len(&self) -> usize {
        self.frac_positions.len_of(Axis(0))
    }

    pub fn cellpar_to_cell(
        cellpar: &Vec<f64>,
        ab_normal: Array1<f64>,
        a_direction: Option<Array1<f64>>,
    ) -> Array2<f64> {
        let a_direction = match a_direction {
            Some(dir) => dir,
            None => {
                if cross(&ab_normal, &array![1.0, 0.0, 0.0]).norm() < 1e-5 {
                    array![0.0, 0.0, 1.0]
                } else {
                    array![1.0, 0.0, 0.0]
                }
            }
        };
        let z = unit_vector(&ab_normal);
        let x = unit_vector(&(&a_direction - a_direction.dot(&z) * &z));
        let y = cross(&z, &x);

        let (a, b, c, alpha, beta, gamma) = match cellpar.len() {
            1 => (cellpar[0], cellpar[0], cellpar[0], 90.0, 90.0, 90.0),
            3 => (cellpar[0], cellpar[1], cellpar[2], 90.0, 90.0, 90.0),
            6 => (
                cellpar[0], cellpar[1], cellpar[2], cellpar[3], cellpar[4], cellpar[5],
            ),
            _ => panic!("Invalid cellpar"),
        };

        let eps = 2.0 * f64::EPSILON * 90.0;
        let cos_alpha = if (alpha.abs() - 90.0).abs() < eps {
            0.0
        } else {
            alpha.to_radians().cos()
        };
        let cos_beta = if (beta.abs() - 90.0).abs() < eps {
            0.0
        } else {
            beta.to_radians().cos()
        };
        let (cos_gamma, sin_gamma) = if (gamma - 90.0).abs() < eps {
            (0.0, 1.0)
        } else if (gamma + 90.0).abs() < eps {
            (0.0, -1.0)
        } else {
            (gamma.to_radians().cos(), gamma.to_radians().sin())
        };

        let va = a * array![1.0, 0.0, 0.0];
        let vb = b * array![cos_gamma, sin_gamma, 0.0];
        let cx = cos_beta;
        let cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
        let cz_sqr = 1.0 - cx * cx - cy * cy;
        assert!(cz_sqr >= 0.0);
        let cz = cz_sqr.sqrt();
        let vc = c * array![cx, cy, cz];

        let abc: Array2<f64> = stack(Axis(0), &[va.view(), vb.view(), vc.view()]).unwrap();
        let t: Array2<f64> = stack(Axis(0), &[x.view(), y.view(), z.view()]).unwrap();
        abc.dot(&t)
    }
}
