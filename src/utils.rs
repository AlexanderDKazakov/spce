///! Utils for use
///!
///  

use std::{
    fs,
    io::{self, BufRead},
    path::Path,
    error::Error,
};

mod representations;
use representations::{
    Structure,
    Molecule,
    Atom
};
use std::str::FromStr;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<fs::File>>>
where P: AsRef<Path>, {
    let file = fs::File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn read_file(filename: &Path) -> Result<Vec<Structure>, Box<dyn std::error::Error + 'static>>{
    println!("Provided:\n [File XYZ] {}", filename.display());
    let mut structures: Vec<Structure> = Vec::new();

    let mut structure: Structure = Structure::new();
    let mut molecule: Molecule   = Molecule::new();

    match read_lines(filename) {
        Ok(lines) => {
            for line in lines {
                if let Ok(line) = line {
                    //println!("{:?}", line);
                    match Atom::from_str(&line) {
                        Err(_) => {
                            //println!("{} is not a valid atom!", line);
                            let info: Vec<&str> = line.split_whitespace().collect();
                            if info.len() != 1 {
                                if structure.has_lat_info() {
                                    structures.push(structure.clone());
                                    structure.reset();
                                }
                                structure.set_lat_info(line.clone());
                            }
                        }
                        Ok(atom) => {
                            // println!("Atom name: {}, charge: {}, pos: {:?}, fpos: {:?}", atom.name, atom.charge, atom.pos, atom.fpos);
                            molecule.add_atom(atom);
                            if molecule.len() == 3 {
                                molecule.update_self_energy_with_model("SPC");
                                // println!("Molecule name: {}, {:?}, energy: {}", molecule.name, molecule.atoms, molecule.energy);
                                structure.add_molecule(molecule.clone());
                                molecule.reset();
                            }
                        }
                    }
                }
            }
            // last 
            if !structure.is_empty() {
                structures.push(structure.clone());
                structure.reset();
            }
            
            Ok(structures)
        },
        Err(e) => Err(Box::new(e)),
    }
}

pub fn help() {
    println!("Usage: ./spce pos.xyz");
}

pub fn run(filename: &Path) -> Result<(), Box<dyn Error>> {
    println!("SPC/E model calculation");
    println!("File provided: {}", filename.display());

    match read_file(filename) {
        Ok(structures) => {
            //println!("{:?}", structures);
            for structure in structures {
                structure.calc_energy();
            }
            Ok(())
        },
        Err(e) => Err(e)
    }
}


