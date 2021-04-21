use std::process;
use std::env;
use std::path::{Path, PathBuf};
use std::fs;
use std::io::{self, BufRead};
use std::str::FromStr;

//extern crate ndarray;
//use ndarray::prelude::*;

#[derive(Debug, PartialEq, Clone)]
struct Structure {
    atoms:     Vec<Atom>,
    lat_info:  String,
    num_atoms: usize,
    energy:    f64,
}

impl Structure {
    fn new() -> Self {
        Structure {
            atoms:    Vec::new(),
            lat_info: String::from(""),
            num_atoms: 0,
            energy: 0f64,
        }
    }
    fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom)
    }

    fn set_lat_info(&mut self, lat_info: String) {
        self.lat_info = lat_info
    }

    fn update_num_atoms(&mut self) {
        self.num_atoms = self.atoms.len()
    }

    fn reset(&mut self) {
        self.atoms.clear();
        self.lat_info = String::from("");
        self.num_atoms = 0;
    }

    fn is_empty(&self) -> bool {
        self.atoms.len() == 0
    }

    fn get_oxygen_atoms(&self) -> Vec<Atom> {
        let mut oxygen_atoms: Vec<Atom> = Vec::new();
        for atom in &self.atoms {
            if atom.get_name() == "O" {
                oxygen_atoms.push(atom.clone());
            }
        }
        oxygen_atoms
    }

    fn get_distances(atoms: Vec<Atom>) -> Vec<f64> {
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

    /// Implementation part of energy: [LJ Energy for oxygen atoms only] 
    ///
    ///                 /   A    \ 6     /   B    \ 12
    /// energy_LJ =  - |   ---    |   + |   ---    |     [kJ/mol]
    ///                 \ r_{OO} /       \ r_{OO} /
    ///
    /// A = 0.37122 (kJ/mol)^{1/6}  * nm
    /// B = 0.34280 (kJ/mol)^{1/12} * nm
    /// Assuming the r in nm -- [1-e9 m]
    fn get_lj_oo_energy(&self) -> f64 {
        let mut lj_oo_energy = 0f64;
        let oxygen_atoms = self.get_oxygen_atoms();
        let distances    = Structure::get_distances(oxygen_atoms);
        let a = 0.37122; 
        let b = 0.34280;
        for distance in distances {
            lj_oo_energy +=  -(a/distance).powi(6) + (b/distance).powi(12)
        }
        lj_oo_energy
    }

    /// Implementation electrostatic contribution to the energy: 
    ///
    ///                    [A atom] [B atom]
    ///                                    Kc * q_{i} * q_{j}
    /// energy_electro =      ∑        ∑   -------------------         [kcal / mol]
    ///                       i        j           r
    ///
    /// Kc           = 332.1    [angstrom * kcal / ( mol * e^{2} )]
    /// Kc*  = 10*Kc = 3321     [   nm    * kcal / ( mol * e^{2} )]
    /// Kc** = 42*Kc = 3321*4.2 [   nm    * kJ   / ( mol * e^{2} ) ]
    /// q_{i/j} in [e]
    /// r in angstrom???
    ///
    fn get_electrostatic_contrib(&self) -> f64 {
        let mut electro_energy = 0f64;
        let Kc = 4.2*3321f64; // nm * kJ / (mol * e^{2})

        for (idx_i, atom_i) in self.atoms.iter().enumerate() {
            for (idx_j, atom_j) in self.atoms.iter().enumerate() {
                if idx_i < idx_j {
                    electro_energy += Kc * atom_i.get_charge() * atom_j.get_charge();
                }
            }
        }

        electro_energy
    }

    /// Calculating SPC energy
    ///
    /// SPC_energy = energy_LJ [kJ/mol] + energy_electro [kJ/mol]
    ///
    fn calc_energy(&self) {
        let energy: f64;
        let lj_oo   = self.get_lj_oo_energy();
        let electro = self.get_electrostatic_contrib();

        energy = lj_oo + electro;

        println!("Total Energy[kcal/mol]: {} | hartree {}", energy,    energy  / 2600f64);     // 
        println!("          LJ[kcal/mol]: {} | hartree {}", lj_oo,     lj_oo   / 2600f64);      // 
        println!("     Electro[kcal/mol]: {} | hartree {}\n", electro, electro / 2600f64);  // 

    }
}


#[derive(Debug, PartialEq, Clone)]
struct Atom {
    name:   String,
    charge: f64,
    pos:    Vec<f64>,
    fpos:   Vec<f64>,
}

impl Atom {
    fn new(name: String, charge: f64, pos: Vec<f64>, fpos: Vec<f64>) -> Self{
        Atom {
            name,
            charge,
            pos,
            fpos,
        }
    }
    fn _set_charge(&mut self, q: f64) {
        self.charge = q;
    }

    fn get_charge(&self) -> f64 {
        self.charge.clone()
    }

    fn distance_to(&self, other: &Atom) -> f64 {
        ( (other.pos[0] - self.pos[0]).powi(2) + (other.pos[1] - self.pos[1]).powi(2) + (other.pos[2] - self.pos[2]).powi(2) ).sqrt()
    }

    fn get_name(&self) -> String {
        self.name.clone()
    }
}

impl FromStr for Atom {
    type Err = std::num::ParseIntError;

    fn from_str(string: &str) -> Result<Self, Self::Err> {
        let name: &str;
        let charge: f64;
        let mut pos:  Vec<f64> = Vec::new();
        let mut fpos: Vec<f64> = Vec::new();

        //H 1.42203366 0.459411945 -0.13831059 0.00499080579 -0.00124616951 0.00616623497
        // to instance of 'Atom'
        let vec_from_string: Vec<&str> = string.split_whitespace().collect();

        // ignore next line:
        // <number of atom>
        // Lattice ...
        if vec_from_string.len() < 2 || vec_from_string.len() > 7  {
            let _ = "a".parse::<i32>()?;
        };

        // H
        name = vec_from_string[0];

        match name {
            "H" => charge = 0.4238,
            "O" => charge = -0.8476,
            _   => charge = 0f64,
        }

        for idx in 1..4 {
            let coord: f64 = vec_from_string[idx].parse::<f64>().unwrap();
            pos.push(coord);
        }
        for idx in 4..7 {
            let fcoord: f64 = vec_from_string[idx].parse::<f64>().unwrap();
            fpos.push(fcoord);
        }

        Ok( 
            Atom { 
                name : name.to_string(), 
                charge,
                pos,
                fpos,
        })
    }
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<fs::File>>>
where P: AsRef<Path>, {
    let file = fs::File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn read_file(filename: &Path) -> Result<Vec<Structure>, Box<dyn std::error::Error + 'static>>{
    println!("Provided:\n [File XYZ] {}", filename.display());

    let mut structures: Vec<Structure> = Vec::new();
    let mut structure: Structure = Structure::new();

    match read_lines(filename) {
     // Consumes the iterator, returns an (Optional) String
        Ok(lines) => {
            for line in lines {
                if let Ok(line) = line {
                    //println!("{:?}", line);
                    match Atom::from_str(&line) {
                        Ok(atom) => {
                            //println!("Atom name: {}, charge: {}, pos: {:?}, fpos: {:?}", atom.name, atom.charge, atom.pos, atom.fpos);
                            structure.add_atom(atom);
                        }
                        Err(_) => {
                            //println!("{} is not a valid atom!", line);
                            let vec_from_string: Vec<&str> = line.split_whitespace().collect();
                            if vec_from_string.len() != 1 {
                                // Lattice ...
                                structure.set_lat_info(line.clone());
                            }

                            if !structure.is_empty() {
                                structure.update_num_atoms();
                                structures.push(structure.clone());
                                structure.reset();
                            }
                        }
                    }
                }
            }
            // last 
            if !structure.is_empty() {
                structure.update_num_atoms();
                structures.push(structure.clone());
                structure.reset();
            }
            //Ok(())
            Ok(structures)
        },
        Err(e) => Err(Box::new(e)),
    }

}


fn help() {
    println!("Usage: ./spce pos.xyz");
}

fn main() {
    // 
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        help();
        process::exit(1);
    }
    let filename = PathBuf::from(&args[1]);
    println!("SPC/E model calculation");
    match read_file(filename.as_path()) {
        Ok(structures) => {
            //println!("{:?}", structures);
            for structure in structures {
                structure.calc_energy();
            }
        },
        Err(e)         => println!("Problem {}", e)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_distance() {
        let atom1 = Atom::new("H1".to_string(), 1.0, vec![1.0, 1.0, 1.0], vec![1.0, 1.0, 1.0]);
        let atom2 = Atom::new("H2".to_string(), 2.0, vec![0.0, 0.0, 0.0], vec![1.0, 1.0, 1.0]);
        let atom3 = Atom::new("H3".to_string(), 3.0, vec![6.0, 4.0, 5.0], vec![1.0, 1.0, 1.0]);

        let distance1 = atom1.distance_to(&atom2);
        let distance2 = atom1.distance_to(&atom3);

        assert!( (distance1 - 3f64.sqrt()) < 1e-10  );
        assert!( (distance2 - 50f64.sqrt()) < 1e-10  );

        assert_eq!("H1", atom1.get_name());
        assert_eq!("H2", atom2.get_name());
        assert_eq!("H3", atom3.get_name());


        assert_eq!(1.0f64, atom1.get_charge());
        assert_eq!(2.0f64, atom2.get_charge());
        assert_eq!(3.0f64, atom3.get_charge());

    }


    #[test]
    fn structure_distances() {
        let mut str1 = Structure::new();

        let atom1 = Atom::new("O".to_string(), 1.0, vec![1.0, 1.0, 1.0], vec![1.0, 1.0, 1.0]);
        let atom2 = Atom::new("O".to_string(), 1.0, vec![0.0, 0.0, 0.0], vec![1.0, 1.0, 1.0]);
        let atom3 = Atom::new("H".to_string(), 1.0, vec![6.0, 4.0, 5.0], vec![1.0, 1.0, 1.0]);
        let atom4 = Atom::new("H".to_string(), 1.0, vec![1.0, 2.0, 3.0], vec![1.0, 1.0, 1.0]);

        str1.add_atom(atom1);
        str1.add_atom(atom2);
        str1.add_atom(atom3);
        str1.add_atom(atom4);

        let oxygen_atoms = str1.get_oxygen_atoms();
        assert_eq!(2, oxygen_atoms.clone().len() );
        assert_eq!(1, Structure::get_distances(oxygen_atoms.clone()).len() );
        assert_eq!(3f64.sqrt(), Structure::get_distances(oxygen_atoms.clone())[0]);

    }
}
