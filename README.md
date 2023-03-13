# biosurfer
"Surf" the biological network, from genome to transcriptome to proteome and back to gain insights into human disease biology.

## Installation
Clone the project repository and create a new conda environment if needed. Then, from the top directory of the project, run `pip install --editable .`

The dependency on [`graph-tool`](https://graph-tool.skewed.de/) currently requires a separate installation step. Run `conda install -c conda-forge graph-tool`.

## Bone proteogenomics analysis ([GitHub](https://github.com/aa9gj/Bone_proteogenomics_manuscript))

### Database generation 

Generate SQLite databases for the list of full-length isoforms derived from the [LRP pipeline](https://github.com/sheynkman-lab/Long-Read-Proteogenomics).

This includes two sets of GTF annotation files:
    1. Full-length isoforms with ORFs that *align* to Leafcutter junction cluster
    2. Full-length isoforms with ORFs that *do not align* to Leafcutter junction cluster
    
Link to processed data can be found in [here](https://github.com/aa9gj/Bone_proteogenomics_manuscript#data-availability).

### NMD and protein truncation analysis

Conducted ORF-related analysis to identify sQTLs that affect full-length isoforms that potentially undergo nonsense mediated decay (NMD) and may result in new expression. 
    

    
    
