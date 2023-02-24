
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;
use std::process;
use clap::{App, AppSettings, Arg};
use gfa_reader::{NGfa, NNode, NPath};
use log::{error, info};
use crate::bed::{BedFile, Node2Feature};



pub mod bed;

fn main() {
    let matches = App::new("gfa_annotate")
        .version("0.1.0")
        .author("Sebastian V")
        .about("Overlap annotation and genome graphs")
        .setting(AppSettings::ArgRequiredElseHelp)
        .help_heading("Input")
        .arg(Arg::new("gfa")
            .short('g')
            .required(true)
            .about("input gfa")
            .takes_value(true))
        .arg(Arg::new("bed")
            .short('b')
            .about("bed file")
            .takes_value(true)
            .required(true))
        .arg(Arg::new("output")
            .short('o')
            .required(true)
            .takes_value(true)
            .about("Output file"))
        .get_matches();


    let gfa = matches.value_of("gfa").unwrap();
    let bed = matches.value_of("bed").unwrap();


    if !Path::new(gfa).exists(){
        error!("No gfa file");
        process::exit(0x0100);

    }
    let bedfile;
    if !Path::new(bed).exists(){
        error!("No bed file");

        process::exit(0x0100);
    } else {// Bed file
    info!("Read the gff/bed file");
        bedfile = BedFile::read_bed(bed);
    }


    // Running the graph
    info!("Read the gfa file");
    let mut graph = NGfa::new();
    graph.from_file_direct(gfa);
    let gfa2pos_btree = node2pos(&graph);
    // Bed file
    info!("Read the gff/bed file");

    // For each genome
    let u = bed_intersection(& graph, &bedfile, &gfa2pos_btree);
    writer_v2(u, matches.value_of("output").unwrap());

}

pub fn node2length(nodes: &HashMap<u32, NNode>) -> HashMap<u32, usize>{
    let node_length: HashMap<u32, usize> = nodes.iter().map(|x| (x.0.clone(), x.1.len)).collect();
    return node_length
}

/// Intersecting the bed file and with positions in the graph
///
/// # Arguments:
/// * 'graph': NGfa data structure
/// * 'path2pos': {genome_id -> {pos -> node_id}}
///
/// # Output
/// - 'node2data'
///     - {u32 -> {u32 -> u32
pub fn bed_intersection<'a>(graph: &'a NGfa, bed: &'a BedFile, path2pos: &'a HashMap<String, BTreeMap<u32, usize>>) ->  Vec<(&'a String, u32, u32, String, &'a [u32], &'a [bool])>{

    //let mut k: HashMap<&'a u32, Vec<BTreeMap<String, String>>> = HashMap::new();
    let mut result = Vec::new();

    for (name, data)  in bed.data.iter(){
        let mut ff=  0;

        for x in graph.paths.iter().enumerate(){
            if &x.1.name == name{
                ff = x.0;
            }
        }
        let ff = &graph.paths[ff];
        if path2pos.contains_key(name){

            let index = path2pos.get(name).unwrap();
            for entry in data{
                let interval: Vec<_> = index.range(entry.start+1..entry.end).collect();

                // This is the node coming after the end of the entry - this might also be covered
                let bigger = index.range(entry.end..).next().unwrap();


                // This is the specific case of you are within a node
                if interval.len() == 0{
                    result.push((name, entry.start, entry.end, entry.tag.clone(), &ff.nodes[*bigger.1..*bigger.1+1], &ff.dir[*bigger.1..*bigger.1+1]));


                } else {

                    result.push((name, entry.start, entry.end, entry.tag.clone(), &ff.nodes[*interval.first().unwrap().1..*bigger.1+1], &ff.dir[*interval.first().unwrap().1..*bigger.1+1]));



                }
            }
        }

    }
    result
}

/// Writing output
///
/// **Comment**
/// Tabular output format
/// Node Tag1 Tag2 Tag3
///
/// **Arguments**
/// * index: Index structure for column name
/// * data: Containing node_id + tags
pub fn writer_v2(data: Vec<(&String, u32, u32, String, &[u32], &[bool])>, output: &str){
    let f = File::create(output).expect("Unable to create file");
    let mut f = BufWriter::new(f);

    for x in data.iter(){
        write!(f, "{}\t{}\t{}\t{}\t{:?}\t{:?}\n", x.0, x.1, x.2, x.3, x.4, x.5).expect("Can't write file");
    }

}




/// Position to node for each genome in the graph
///
/// # Arguments
/// * 'graph' - A NGfa data structure
///
/// # Output
/// {Genome_id ->  BtreeMap(position -> node)}
///
pub fn node2pos(graph: &NGfa) -> HashMap<String, BTreeMap<u32, usize>>{
    let mut result: HashMap<String, BTreeMap<u32, usize>> = HashMap::new();

    for path in graph.paths.iter(){
        let mut btree = BTreeMap::new();
        let mut position = 0;
        for (i, node) in path.nodes.iter().enumerate(){
            // {"End"-position of the node -> node_id}
            btree.insert(position + graph.nodes.get(node).unwrap().len as u32, i);
            // Update position
            position += graph.nodes.get(node).unwrap().len as u32
        }
        // Add btree to corresponding path
        btree.insert(position+1, 0);
        result.insert(path.name.clone(), btree);
    }

    return result
}
