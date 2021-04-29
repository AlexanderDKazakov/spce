mod molecule;

const W1_ENERGY_HARTREE: f64 = -17.164634771108034;

pub use molecule::{Molecule, Atom};

#[derive(Debug, PartialEq, Clone)]
pub struct Structure {
    molecules:     Vec<Molecule>,
    lat_info:      String,
    num_molecules: usize,
    energy:        f64,  
    energy_diff:   f64, 
    dipole_moment: Vec<f64>,
}

impl Structure {
    pub fn new() -> Self {
        Structure {
            molecules:    Vec::new(),
            lat_info: String::from(""),
            num_molecules: 0,
            energy: 0f64,
            energy_diff: 0f64,
            dipole_moment: vec![0.0f64, 0.0f64, 0.0f64],
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

    pub fn _update_num_molecules(&mut self) {
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
    ///                                    kc * q_{i} * q_{j}
    /// energy_electro =      ∑        ∑   -------------------    [kJ / mol]
    ///                       i        j           r_{ij}
    ///
    /// ```
    /// `kc`           = 332.1               [angstrom * kcal / ( mol * e^{2} )]
    /// 
    /// `kc*`          = 332.1*4.2 = 1394.82 [angstrom * kJ   / ( mol * e^{2} )]
    /// 
    /// `q_{i/j}` in [e]
    /// 
    /// Assuming the r_{ij} in angstrom -- [1e-10 m]
    /// 
    fn get_inter_electrostatic_contrib(&self) -> f64 {
        let mut electro_energy = 0f64;
        let kc = 1394.82f64; // angstrom * kJ / (mol * e^{2})

        for (idx_mol_i, mol_i) in self.molecules.iter().enumerate() {
            for (idx_mol_j, mol_j) in self.molecules.iter().enumerate() {
                if idx_mol_i < idx_mol_j {

                    for atom_k in mol_i.get_atoms().iter() {
                        for atom_l in mol_j.get_atoms().iter() {
                            let r_ij: f64 = atom_k.distance_to(atom_l);
                            electro_energy += (kc * atom_k.get_charge() * atom_l.get_charge()) / r_ij;
                        }
                    }
                }
            }
        }
        electro_energy
    }

    /// ## Parse some quantities from DFT calculations
    ///
    /// Examples: 
    ///           energy
    ///           dipole_moment
    ///           ...?
    pub fn parse(&self, what:&str) -> Vec<f64> {
        let mut output: Vec<f64> = Vec::new();

        /*_Lattice="10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" 
         * Properties=species:S:1:pos:R:3:forces:R:3 
         * energy=VAL pbc="T T T" dipole="X Y Z"*/
        let iter: Vec<&str> = self.lat_info.split("=").collect();

        if what == "energy" {
            if let Some(value) = iter.get(3) {
                // VAL pbc
                let energy_p_something: Vec<&str> = value.split_whitespace().collect();
                if let Some (value) = energy_p_something.get(0) {
                    if let Ok(energy_dft) = value.parse::<f64>() {
                        output.push(energy_dft);
                    }
                }
            }
        }
        if what == "dipole_moment" {
            if let Some(value) = iter.get(5) {
                // "VALX VALY VALZ"
                let cutted_string: &str = &value[1..value.len()-1];
                let dipole_moment_xyz: Vec<&str> = cutted_string.split_whitespace().collect();
                for dipole_moment_value in dipole_moment_xyz {
                    if let Ok(dipole_moment_per_axis) = dipole_moment_value.parse::<f64>() {
                        output.push(dipole_moment_per_axis);
                    }
                }
            }
        }

        output
    }

    /// ## Calculating SPC energy
    ///
    /// `SPC_energy` = `energy_LJ` [kJ/mol] + `energy_electro` [kJ/mol]
    ///
    pub fn calc_energy(&mut self) {

        let nmulw1: f64 = self.len() as f64 * W1_ENERGY_HARTREE;

        let energy_dft_hartree: f64;
        match self.parse("energy").get(0) {
            Some(energy) => energy_dft_hartree = energy.clone(),
            None         => energy_dft_hartree = 0.0f64,
        }

        match self.parse("dipole_moment").len() {
            3 => self.dipole_moment = self.parse("dipole_moment"), 
            _ => self.dipole_moment = vec![0.0f64, 0.0f64, 0.0f64],
        }

        // LJ O-O
        let lj_oo:                   f64 = self.get_lj_oo_energy();
        // Inter molecular interactions
        let inter_electro:           f64 = self.get_inter_electrostatic_contrib();
        // Polarization corrections
        let polarization_correction: f64 = self.get_polarization_correction();

        self.energy      = lj_oo + inter_electro + polarization_correction;
        self.energy_diff = (energy_dft_hartree - nmulw1) - self.energy / 2600f64;

        eprintln!("Structure: #{} ", self.len());
        eprintln!("            LJ[kJ/mol]: {:>15.6} | [hartree] : {:>10.6}", lj_oo,                   lj_oo                   / 2600f64);
        eprintln!("[Inter]Electro[kJ/mol]: {:>15.6} | [hartree] : {:>10.6}", inter_electro,           inter_electro           / 2600f64);
        eprintln!("  Polarization[kJ/mol]: {:>15.6} | [hartree] : {:>10.6}", polarization_correction, polarization_correction / 2600f64);
        eprintln!("  Total Energy[kJ/mol]: {:>15.6} | [hartree] : {:>10.6}", self.energy,             self.energy             / 2600f64);
        eprintln!("[hartree]   DFT energy: {:>15.6} | N*W1      : {:>10.6} | [W1] {:10.6}", energy_dft_hartree, nmulw1, W1_ENERGY_HARTREE);
        eprintln!("[hartree] (DFT-N*W1)-SPC/E: {:>11.6} | [DFT-N*W1]: {:>10.6}\n",  self.energy_diff, energy_dft_hartree - nmulw1);
    }

    /// ## Prepare output for NNs
    ///
    /// `exyz` fileformat
    ///
    /// <number_of_atoms>
    /// Lattice="10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" Properties=species:S:1:pos:R:3:forces:R:3 energy=VAL pbc="T T T" dipole="X Y Z"
    /// <ATOM> POS_X POS_Y POS_Z FPOS_X FPOS_Y FPOS_Z
    /// ...
    ///
    pub fn make_exyz(&self) {
        println!("{}", self.len() as u32 * 3u32);
        println!("Lattice=\"10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0\" Properties=species:S:1:pos:R:3:forces:R:3 energy={} pbc=\"T T T\" dipole=\"{}\"",
            self.energy_diff,
            self.dipole_moment.iter().fold(String::new(), |acc, &val| acc + &val.to_string() + " "),
            );
        for molecule in &self.molecules {
            for atom in molecule.get_atoms() {
                println!("{} {} {}", 
                    atom.get_name(), 
                    atom.get_pos().iter().fold(String::new(), |acc, &pos| acc + &pos.to_string() + " "), 
                    atom.get_fpos().iter().fold(String::new(), |acc, &pos| acc + &pos.to_string() + " ")
                    );
            }
        }

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
        assert!( (-29.675702726834345f64 - str1.get_inter_electrostatic_contrib() ).abs() < 1e-9 );

    }
}
