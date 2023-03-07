import pandas as pd
import shutil

"""
        Initializes the LeafCutterJunctionCluster object.
        
        Args:
        chrom (str): The chromosome where the cluster is located.
        start (int): The start position of the cluster.
        end (int): The end position of the cluster.
        strand (str): The strand of the cluster ('+' or '-').
        name (str): The name of the cluster.
        junctions (list): A list of junctions that belong to the cluster.
"""

# Open the input file and a temporary file for writing
with open("CDS_corrected.gtf") as infile, open("temp.gtf", "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            # Write comment lines to the output file
            outfile.write(line)
        else:
            # Split the line by tabs
            fields = line.strip().split("\t")
            
            # Get the attributes column and split by semicolon
            attributes = fields[8].split(";")
            
            # Extract gene_name and transcript_id attributes
            gene_name = transcript_id = None
            cds_type = "sequence_feature"
            for attr in attributes:
                if "gene_name" in attr:
                    gene_name = attr.strip().split(" ")[1].strip('"')
                elif "transcript_id" in attr:
                    transcript_id = attr.strip().split(" ")[1].strip('"')
                elif "CDS_type" in attr:
                    cds_type = attr.strip().split(" ")[1].strip('"')
            
            # Update the attributes column
            new_attrs = f"gene_id \"{gene_name}\"; transcript_id \"{gene_name}|{transcript_id}|1\";"
            
            # Write the updated line to the output file
            outfile.write("\t".join([fields[0], fields[1], cds_type, fields[3], fields[4], fields[5], fields[6], fields[7], new_attrs]) + "\n")

# Replace the input file with the contents of the temporary file
shutil.move("temp.gtf", "CDS_corrected.gtf")

with open('new_gtf_file_with_CDS_reformatted.gtf', 'r') as infile, open('new_gtf_file_with_CDS_reformatted_2.gtf', 'w') as outfile:
    lines_seen = set()  # To keep track of lines already written to the output file
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)  # Write header lines as is
        else:
            line_parts = line.strip().split('\t')
            feature = line_parts[2]
            attributes = line_parts[-1]
            # Create a unique key for the current line by combining chromosome, feature type, start, and end positions
            key = (line_parts[0], feature, line_parts[3], line_parts[4])
            if key not in lines_seen:
                lines_seen.add(key)
                # Write the line to the output file without redundancies
                outfile.write(f"{line_parts[0]}\t{line_parts[1]}\t{feature}\t{line_parts[3]}\t{line_parts[4]}\t{line_parts[5]}\t{line_parts[6]}\t{line_parts[7]}\t{attributes}\n")

