mod molecule;

pub use molecule::{Molecule, Atom};

#[derive(Debug, PartialEq, Clone)]
pub struct Structure {
    molecules:     Vec<Molecule>,
    lat_info:      String,
    num_molecules: usize,
    energy:        f64,  // not sure that we need it
}

impl Structure {
    pub fn new() -> Self {
        Structure {
            molecules:    Vec::new(),
            lat_info: String::from(""),
            num_molecules: 0,
            energy: 0f64,
        }
    }

    pub fn set_lat_info(&mut self, lat_info: String) {
        self.lat_info = lat_info
    }

    pub fn has_lat_info(&self) -> bool {
        !self.lat_info.is_empty()
    }

    pub fn add_molecule(&mut self, molecule: Molecule) {
        self.molecules.push(molecule)
    }

    pub fn update_num_molecules(&mut self) {
        self.num_molecules = self.molecules.len()
    }

    pub fn len(&self) -> usize {
        self.molecules.len()
    }

    pub fn reset(&mut self) {
        self.molecules.clear();
        self.lat_info.clear();
        self.num_molecules = 0;
    }

    pub fn is_empty(&self) -> bool {
        self.molecules.len() == 0
    }

    pub fn get_oxygen_atoms(&self) -> Vec<Atom> {
        let mut oxygen_atoms: Vec<Atom> = Vec::new();
        for molecule in &self.molecules {
            oxygen_atoms.push(molecule.get_oxygen_atom().unwrap());
        }
        oxygen_atoms
    }

    ///
    /// ## Get distances for list of Atoms
    ///
    ///
    pub fn get_distances(atoms: Vec<Atom>) -> Vec<f64> {
        let mut distances: Vec<f64> = Vec::new();
        for (idx_i, atom_i) in atoms.iter().enumerate() {
            for (idx_j, atom_j) in atoms.iter().enumerate() {
                if idx_i < idx_j {
                    distances.push(atom_i.distance_to(atom_j));
                }
            }
        }

        distances
    }

    /// ## Implementation part of energy: [LJ Energy for oxygen atoms only] 
    /// 
    /// ```
    ///
    ///                 /   A    \ 6     /   B    \ 12
    /// energy_LJ =  - |   ---    |   + |   ---    |     [kJ/mol]
    ///                 \ r_{OO} /       \ r_{OO} /
    ///
    /// ```
    /// 
    /// `A` = 0.37122 (kJ/mol)^{1/6}  * nm => 3.7122 (kJ/mol)^{1/6}  * angstrom
    /// 
    /// `B` = 0.34280 (kJ/mol)^{1/12} * nm => 3.4280 (kJ/mol)^{1/12} * angstrom
    /// 
    /// Assuming the r in angstrom -- [1e-10 m]
    fn get_lj_oo_energy(&self) -> f64 {
        let mut lj_oo_energy = 0f64;
        let oxygen_atoms = self.get_oxygen_atoms();
        let distances    = Structure::get_distances(oxygen_atoms);
        let a = 3.7122; 
        let b = 3.4280;
        for distance in distances {
            lj_oo_energy += -(a/distance).powi(6) + (b/distance).powi(12)
        }
        lj_oo_energy
    }


    /// ## Implementation polarization correction to the electrostatic contribution:
    /// 
    /// ```
    ///                               
    ///               1     (µ - µ_{0})^{2}
    /// correction = --- ∑ ----------------  = 5.22 [kJ / mol]
    ///               2  i       α_{i} 
    ///
    /// ```
    /// 
    /// `µ`     -- electric dipole moment of the effectively polarized water molecule (2.35 D for the SPC/E model)
    /// 
    /// `μ_{0}` -- is the dipole moment of an isolated water molecule (1.85 D from experiment)
    /// 
    /// `α` is an isotropic polarizability constant, with a value of 1.608×10−40 F·m
    ///
    /// Since the charges in the model are constant, this correction just results in adding 1.25 kcal/mol (5.22 kJ/mol) to the total energy.
    /// 
    fn get_polarization_correction(&self) -> f64 {
        let correction: f64 = 5.22; // kJ/mol
        correction
    }

    /// ## Implementation electrostatic contribution to the energy: 
    /// 
    /// ```
    ///                    [A atom] [B atom]
    ///                                    Kc * q_{i} * q_{j}
    /// energy_electro =      ∑        ∑   -------------------    [kJ / mol]
    ///                       i        j           r_{ij}
    ///
    /// ```
    /// `Kc`           = 332.1               [angstrom * kcal / ( mol * e^{2} )]
    /// 
    /// `Kc*`          = 332.1*4.2 = 1394.82 [angstrom * kJ   / ( mol * e^{2} )]
    /// 
    /// `q_{i/j}` in [e]
    /// 
    /// Assuming the r_{ij} in angstrom -- [1e-10 m]
    /// 
    fn get_inter_electrostatic_contrib(&self) -> f64 {
        let mut electro_energy = 0f64;
        let Kc = 1394.82f64; // angstrom * kJ / (mol * e^{2})

        for (idx_mol_i, mol_i) in self.molecules.iter().enumerate() {
            for (idx_mol_j, mol_j) in self.molecules.iter().enumerate() {
                if idx_mol_i < idx_mol_j {

                    for atom_k in mol_i.get_atoms().iter() {
                        for atom_l in mol_j.get_atoms().iter() {
                            let r_ij: f64 = atom_k.distance_to(atom_l);
                            electro_energy += (Kc * atom_k.get_charge() * atom_l.get_charge()) / r_ij;
                        }
                    }
                }
            }
        }
        electro_energy
    }

    /// ## Calculating SPC energy
    ///
    /// `SPC_energy` = `energy_LJ` [kJ/mol] + `energy_electro` [kJ/mol]
    ///
    pub fn calc_energy(&self) {
        let energy: f64;
        let mut intra_electro: f64 = 0.0f64;
        
        // self energy via model
        for mol in self.molecules.iter() { intra_electro += mol.get_self_energy(); }
        // LJ O-O
        let lj_oo         = self.get_lj_oo_energy();
        // Inter molecular interactions
        let inter_electro = self.get_inter_electrostatic_contrib();
        // Polarization corrections
        let polarization_correction: f64 = self.get_polarization_correction();

        energy = lj_oo + inter_electro + intra_electro + polarization_correction;

        println!("Structure: #{} ", self.len());
        println!("            LJ[kJ/mol]: {:>15.6} | hartree: {:>10.6}",   lj_oo,                   lj_oo                   / 2600f64);
        println!("[Inter]Electro[kJ/mol]: {:>15.6} | hartree: {:>10.6}",   inter_electro,           inter_electro           / 2600f64);
        println!("[Intra]Electro[kJ/mol]: {:>15.6} | hartree: {:>10.6}",   intra_electro,           intra_electro           / 2600f64);
        println!("Polarization  [kJ/mol]: {:>15.6} | hartree: {:>10.6}",   polarization_correction, polarization_correction / 2600f64);
        println!("  Total Energy[kJ/mol]: {:>15.6} | hartree: {:>10.6}\n", energy,                  energy                  / 2600f64);

    }
}


#[cfg(test)]
mod tests {
    use super::*;

    
    #[test]
    fn structure_test() {

        let mut str1 = Structure::new();

        assert_eq!(false, str1.has_lat_info());
        str1.set_lat_info(String::from("Lattice..."));
        assert_eq!(true, str1.has_lat_info());

        let mut mol1 = Molecule::new();                  /*   A    A    A  */
        let atom1 = Atom::new("O".to_string(), -0.8476, vec![1.0, 1.0, 1.0], vec![1.0, 1.0, 1.0]);
        let atom2 = Atom::new("H".to_string(), 0.4238,  vec![1.0, 1.0, 2.0], vec![1.0, 1.0, 1.0]);
        let atom3 = Atom::new("H".to_string(), 0.4238,  vec![1.0, 1.0, 0.0], vec![1.0, 1.0, 1.0]);

        mol1.add_atom(atom1.clone());
        mol1.add_atom(atom2);
        mol1.add_atom(atom3);
        mol1.set_name("H2O_1");

        let mut mol2 = Molecule::new();                  /*   A    A    A  */
        let atom4 = Atom::new("O".to_string(), -0.8476, vec![2.0, 2.0, 2.0], vec![1.0, 1.0, 1.0]);
        let atom5 = Atom::new("H".to_string(), 0.4238,  vec![2.0, 2.0, 3.0], vec![1.0, 1.0, 1.0]);
        let atom6 = Atom::new("H".to_string(), 0.4238,  vec![2.0, 2.0, 1.0], vec![1.0, 1.0, 1.0]);

        mol2.add_atom(atom4.clone());
        mol2.add_atom(atom5);
        mol2.add_atom(atom6);
        mol2.set_name("H2O_2");

        str1.add_molecule(mol1);
        str1.add_molecule(mol2);

        assert_eq!(2, str1.get_oxygen_atoms().len());
        assert_eq!(Some(&atom1), str1.get_oxygen_atoms().get(0));
        assert_eq!(Some(&atom4), str1.get_oxygen_atoms().get(1));

        // Distances
        assert_eq!(1, Structure::get_distances(str1.get_oxygen_atoms()).len());
        assert!( (3f64.sqrt() - Structure::get_distances(str1.get_oxygen_atoms()).get(0).unwrap()) < 1e-9  );

        // E: LJ O-O
        assert!( (3515.1980452787766 - str1.get_lj_oo_energy() ).abs() < 1e-9 );

        // E: electro
        assert_eq!(-29.675702726834345f64 + 5.22f64,  str1.get_electrostatic_contrib());

    }
}
