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
    let mut structures: Vec<Structure> = Vec::new();
    let mut structure:  Structure      = Structure::new();

    let mut atoms:      Vec<Atom>      = Vec::new();
    let mut molecules:  Vec<Molecule>  = Vec::new();
    let mut molecule:   Molecule       = Molecule::new();

    let mut idx: usize = 0;
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
                            //idx = 0;
                        }
                        Ok(atom) => {
                            //println!("Atom name: {}, pos: {:?}", atom.name, atom.pos);
                            match atom.get_name() {

                                "H" => {
                                    if !molecule.has_natoms(2, "H") {
                                        molecule.add_atom(atom);
                                    } else {
                                        molecules.push(molecule.clone());
                                        molecule.reset();
                                        molecule.add_atom(atom);
                                    }
                                },
                                "O" => {
                                    // last 2 H
                                    if molecule.has_natoms(2, "H") {
                                        molecules.push(molecule.clone());
                                        molecule.reset();
                                    } 
                                    match molecules.get_mut(idx) {
                                        Some(molecule) => {
                                            molecule.add_atom(atom);
                                            if molecule.len() == 3 {
                                                structure.add_molecule(molecule.clone());
                                            }
                                        },
                                        None => panic!("Strange, but there is more O than molecules...")
                                    }
                                    idx += 1;  // molecules 
                                },
                                _ => panic!("Not implemented case :( ")
                            }
                        }
                    }
                }
            }
        },
        Err(e) => println!("Problem: {}", e),
    }

    // last 
    if !structure.is_empty() {
        structures.push(structure.clone());
        structure.reset();
    }
    
    Ok(structures)


}

pub fn help() {
    println!("Usage: ./spce extenden_xyz_file.xyz");
}

pub fn run(filename: &Path) -> Result<(), Box<dyn Error>> {
    eprintln!("### SPC/E model calculation ###");
    eprintln!("File provided: {}\n", filename.display());

    match read_file(filename) {
        Ok(structures) => {
            //println!("{:?}", structures);
            for mut structure in structures {
                structure.calc_energy();
                structure.make_exyz();
            }
            Ok(())
        },
        Err(e) => Err(e)
    }
}


#[cfg(test)]
mod tests {
    use super::*;


    fn test_read_file() {

    }

}


