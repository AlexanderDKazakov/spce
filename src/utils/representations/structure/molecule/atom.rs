use std::str::FromStr;

#[derive(Debug, PartialEq, Clone)]
pub struct Atom {
    pub name:   String,
    pub charge: f64,
    pub pos:    Vec<f64>,
    pub fpos:   Vec<f64>,
}

impl Atom {
    pub fn new(name: String, charge: f64, pos: Vec<f64>, fpos: Vec<f64>) -> Self {
        Atom {
            name,
            charge,
            pos,
            fpos,
        }
    }
    pub fn _set_charge(&mut self, q: f64) {
        self.charge = q;
    }

    pub fn get_charge(&self) -> f64 {
        self.charge.clone()
    }

    pub fn distance_to(&self, other: &Atom) -> f64 {
        ( (other.pos[0] - self.pos[0]).powi(2) + (other.pos[1] - self.pos[1]).powi(2) + (other.pos[2] - self.pos[2]).powi(2) ).sqrt()
    }

    /// Get model distance 
    ///  model             TIPS     SPC    TIP3P    SPC/E
    /// NOTE: r(OH), Å    0.9572    1.0    0.9572    1.0 
    ///       α(HOH) deg  104.52   109.47  104.52   109.47 
    ///       r(HH), Å    1.5816   1.6330  1.5816   1.6330
    ///
    pub fn model_distance_to(&self, model: &str, other: &Atom) -> Option<f64> {
        match (model, self.get_name(), other.get_name()) {
            ("SPC", "O", "H")    => Some(1.0000f64),
            ("SPC", "H", "O")    => Some(1.0000f64),
            ("SPC", "H", "H")    => Some(1.6330f64),
            ("SPC", "O", "O")    => None,
            ("SPC/E", "O", "H")  => Some(1.0000f64),
            ("SPC/E", "H", "O")  => Some(1.0000f64),
            ("SPC/E", "H", "H")  => Some(1.6330f64),
            ("SPC/E", "O", "O")  => None,
            (  _,     _,   _)  => None,
        }

    }

    pub fn get_name(&self) -> &str {
        self.name.as_str()
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
            }
        )
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
}