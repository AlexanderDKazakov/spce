pub mod atom;

pub use atom::Atom;

#[derive(Debug, PartialEq, Clone)]
pub struct Molecule {
    pub name:   String,
    pub atoms:  Vec<Atom>,
    pub energy: f64,
}

impl Molecule {

    pub fn new() -> Self {
        Molecule {
            name:   String::from("Molecule"),
            atoms:  Vec::new(),
            energy: 0.0f64,
        }
    }

    pub fn set_name(&mut self, name: &str) {
        self.name = String::from(name)
    }

    pub fn get_name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom)
    }

    pub fn reset(&mut self) {
        self.name.clear();
        self.atoms.clear();
        self.energy = 0.0f64;
    }

    pub fn get_atoms(&self) -> Vec<Atom> {
        self.atoms.clone()
    }

    pub fn get_oxygen_atom(&self) -> Option<Atom> {
        for atom in &self.atoms {
            if atom.get_name() == "O" {
                return Some(atom.clone())
            }
        }
        None        
    }

    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// # Molecule energy corresponds to SPC/E model
    ///
    /// Implementation electrostatic contribution to the energy: 
    /// ```
    ///                    [A atom] [B atom]
    ///                                    Kc * q_{i} * q_{j}
    /// energy_electro =      ∑        ∑   -------------------    [kJ / mol]
    ///                       i        j           r_{ij}
    /// ```
    /// `Kc`           = 332.1               [angstrom * kcal / ( mol * e^{2} )]
    /// 
    /// `Kc*`          = 332.1*4.2 = 1394.82 [angstrom * kJ   / ( mol * e^{2} )]
    /// 
    /// q_{i/j} in [e]
    ///
    /// Get model distance 
    /// ```
    ///  model             TIPS     SPC    TIP3P    SPC/E
    /// NOTE: r(OH), Å    0.9572    1.0    0.9572    1.0 
    ///       α(HOH) deg  104.52   109.47  104.52   109.47 
    ///       r(HH), Å    1.5816   1.6330  1.5816   1.6330
    /// ```
    /// 
    /// Assuming the r_{ij} in angstrom -- [1e-10 m]
    pub fn update_self_energy_with_model(&mut self, model: &str) {
        let mut energy: f64 = 0.0f64;
        let Kc = 1394.82f64; // angstrom * kJ / (mol * e^{2})

        // go through all 3 atoms
        for (idx_i, atom_i) in self.atoms.iter().enumerate() {
            for (idx_j, atom_j) in self.atoms.iter().enumerate() {
                if idx_i < idx_j {
                    if let Some(r_ij) = atom_i.model_distance_to(model, atom_j) {
                        energy += (Kc * atom_i.get_charge() * atom_j.get_charge()) / r_ij;
                    } else {
                        panic!("No model distance for {} {}", atom_i.get_name(), atom_j.get_name());
                    }
                }
            }
        }

        self.energy = energy;
    }

    pub fn get_self_energy(&self) -> f64 {
        self.energy
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    
    #[test]
    fn molecule_distances() {
        let mut mol = Molecule::new();

        let atom1 = Atom::new("O".to_string(), -0.8476, vec![1.0, 1.0, 1.0], vec![1.0, 1.0, 1.0]);
        let atom2 = Atom::new("H".to_string(),  0.4238, vec![1.0, 1.0, 1.0], vec![1.0, 1.0, 1.0]);
        let atom3 = Atom::new("H".to_string(),  0.4238, vec![1.0, 1.0, 1.0], vec![1.0, 1.0, 1.0]);

        mol.add_atom(atom1);
        mol.add_atom(atom2);
        mol.add_atom(atom3);

        mol.set_name("H2O");

        assert!(mol.get_oxygen_atom().is_some());

        mol.update_self_energy_with_model("SPC");
        assert!( (-848.6645422369294 - mol.energy).abs() < 1e-9  );

        mol.update_self_energy_with_model("SPC/E");
        assert!( (-848.6645422369294 - mol.energy).abs() < 1e-9  );
    }
}
