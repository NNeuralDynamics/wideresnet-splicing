import os
import requests
import csv

# Define the input and output paths
input_file = "encode-shRNA-metadata.tsv"  # Input TSV file path
output_file = "experiment_controls.tsv"  # Output TSV file path
base_dir = "/scratch/runyan.m/controls/bams"  # Base directory for BAM files

# Function to fetch JSON from the ENCODE API
def fetch_encode_data(resource):
    base_url = "https://www.encodeproject.org/{}?format=json"
    headers = {"accept": "application/json"}
    response = requests.get(base_url.format(resource), headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Failed to fetch data for {resource}. HTTP {response.status_code}")
        return None

# Process the input file and fetch possible controls
def process_file(input_path, output_path, base_dir):
    processed_experiments = set()  # Keep track of processed experiment accessions
    unique_controls = set()  # Keep track of unique possible controls

    with open(input_path, "r") as infile, open(output_path, "w", newline="") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")
        
        # Write the header row
        writer.writerow([
            "Experiment Accession",
            "Possible Control Accessions",
            "Control BAM File 1",
            "Control BAM File 2",
            "Knockdown BAM File 1",
            "Knockdown BAM File 2",
            "Biosample Term Name",
            "Library Strand Specific"
        ])
        
        # Iterate over the rows in the input file
        for row in reader:
            exp_accession = row["Experiment accession"]
            
            # Skip if this experiment accession has already been processed
            if exp_accession in processed_experiments:
                continue
            
            print(f"Processing {exp_accession}...")
            
            # Extract Biosample Term Name and Library Strand Specific from the row
            biosample_term_name = row.get("Biosample term name", "Unknown")
            library_strand_specific = row.get("Library strand specific", "Unknown")
            
            # Fetch experiment data
            exp_data = fetch_encode_data(f"experiments/{exp_accession}")
            if exp_data:
                # Extract possible controls and remove the 'experiments/' prefix
                possible_controls = [
                    control["@id"].replace("experiments/", "").strip("/")
                    for control in exp_data.get("possible_controls", [])
                ]
                unique_controls.update(possible_controls)  # Add to unique controls set
                
                # Find BAM file paths
                control_bams = find_bam_files(exp_accession, base_dir, "control")
                knockdown_bams = find_bam_files(exp_accession, base_dir, "KD")
                
                # Write the row to the output file
                writer.writerow([
                    exp_accession, 
                    ", ".join(possible_controls) if possible_controls else "No controls found",
                    control_bams[0] if len(control_bams) > 0 else "None",
                    control_bams[1] if len(control_bams) > 1 else "None",
                    knockdown_bams[0] if len(knockdown_bams) > 0 else "None",
                    knockdown_bams[1] if len(knockdown_bams) > 1 else "None",
                    biosample_term_name,
                    library_strand_specific
                ])
            else:
                writer.writerow([exp_accession, "No controls found", "None", "None", "None", "None", biosample_term_name, library_strand_specific])
            
            # Mark this experiment as processed
            processed_experiments.add(exp_accession)

    # Print the number of unique possible controls
    print(f"Number of unique possible controls: {len(unique_controls)}")

# Function to find BAM files for a specific type (control or KD)
def find_bam_files(exp_accession, base_dir, sample_type):
    bam_files = []
    exp_dir = os.path.join(base_dir, exp_accession, "star")
    
    # Check if the experiment directory exists
    if os.path.exists(exp_dir):
        for rbp in os.listdir(exp_dir):  # Iterate over RBPs
            rbp_dir = os.path.join(exp_dir, rbp, sample_type)
            if os.path.exists(rbp_dir):
                for file in os.listdir(rbp_dir):
                    if file.endswith(".bam"):
                        bam_files.append(os.path.join(rbp_dir, file))
    return bam_files

# Run the script
if __name__ == "__main__":
    if not os.path.exists(input_file):
        print(f"Input file '{input_file}' not found.")
    else:
        process_file(input_file, output_file, base_dir)
        print(f"Output written to '{output_file}'")
