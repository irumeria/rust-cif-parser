pub mod cif;
pub mod structure;
pub mod sym;


#[cfg(test)]
mod tests {
    use crate::cif::parse_cif;
    use std::path::Path;
    
    #[test]
    fn test_parse_cif() {
        let cif_path = Path::new("assets/ABW-0.cif");
        let structure = parse_cif(cif_path);
        assert_eq!(structure.positions.len(), 72);
    }
}