mod lib;

use std::process;
use std::env;
use std::path::Path;

fn main() {

    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        lib::help();
        process::exit(1);
    }
    let filename = Path::new(&args[1]);
    
    if let Err(e) = lib::run(filename) {
        eprintln!("Application error: {}", e);
        process::exit(1);
    }
}
