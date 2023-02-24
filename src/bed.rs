use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use gfa_reader::{NGfa};
use log::{info};




#[derive(Debug)]
/// Whole BED file
///
/// For better post-processing it is already sorted by chromosome (firts column)
pub struct BedFile {
    pub size: usize,
    pub data: HashMap<String, Vec<BedEntry>>
}


#[derive(Debug)]
/// Bed entry
/// chrom (or scaffold) name are stored in the BedFile
pub struct BedEntry {
    pub start: u32,
    pub end: u32,
    pub tag: String
}




impl BedFile {

    pub fn new() -> Self{
        Self{
            size: 0,
            data: HashMap::new(),
        }
    }


    /// Reading a bed file
    /// Format description: https://en.wikipedia.org/wiki/BED_(file_format)
    /// Fast and dirty
    pub fn read_bed(filename: &str) -> Self{
        info!("Reading BED file");
        let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");
        let reader = BufReader::new(file);

        // Fasta_entry -> Vec<Bed_Entry>
        let mut result = HashMap::new();
        for line in reader.lines() {
            let line_data = line.unwrap();
            //let p: Vec<&str> = l.split("\t").collect();
            let col: Vec<&str> = line_data.split("\t").collect();
            //let ko3: HashSet<String> = p2.nth(4).unwrap().split(del).map(|s| s.to_string()).collect();
            // If you are not empty
            if col.len() != 4{
                eprintln!("ERROR")
            }
            let tag = col[3].to_string();
            result.entry(col[0].to_string()).or_insert(Vec::new()).push(BedEntry { start: col[1].parse().unwrap(), end: col[2].parse().unwrap(), tag: tag })
        }



        Self{
            size: get_size(&result),
            data: result,

        }
    }


}





/// Get the total size of the HashMap
pub fn get_size(data: &HashMap<String, Vec<BedEntry>>) -> usize{
    let mut count = 0;
    for (_name, entry) in data.iter(){
        count += entry.len();
    }
    return count;
}


/// Node to feature
pub struct Node2Feature{
    pub data: HashMap<u32, HashSet<String>>
}

impl Node2Feature{

    /// From graph for empty data structure
    pub fn new(graph: &NGfa) -> Self{
        let mut k: HashMap<u32, HashSet<String>> = HashMap::new();
        for (id, _node) in graph.nodes.iter(){
            k.insert(id.clone(), HashSet::new());
        }
        Self{
            data: k,
        }
    }
}






