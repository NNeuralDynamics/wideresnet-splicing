import os
import subprocess
import logging
import argparse
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv

def load_config(config_path):
    """
    Load configuration from a JSON file.
    """
    with open(config_path, 'r') as f:
        return json.load(f)


def create_directory(path):
    """
    Safely create a directory if it does not already exist.
    """
    try:
        os.makedirs(path, exist_ok=True)
        logging.info(f"Created directory: {path}")
    except Exception as e:
        logging.warning(f"Directory {path} already exists or cannot be created: {e}")


def extract_sample_name(bam_file):
    """
    Extract the sample name from the BAM file name.
    """
    return os.path.basename(bam_file).replace('.bam', '')


def create_spliser_junctions(bam_file, config):
    """
    Create junction files using regtools.
    """
    sample_name = extract_sample_name(bam_file)
    junc_file = os.path.join(config["junction_files_dir"], f"{sample_name}.spliser.junc")

    if os.path.exists(junc_file):
        logging.info(f"Junction file {junc_file} already exists, skipping creation.")
        return junc_file

    logging.info(f"Creating junction file {junc_file} from {bam_file}...")
    command_junc = ["regtools", "junctions", "extract"] + config["regtools_args"] + [bam_file, "-o", junc_file]
    subprocess.run(command_junc, check=True)
    logging.info(f"Junction file created at {junc_file}")
    return junc_file


def run_spliser_process(bam_file, junc_file, config):
    """
    Run SpliSER process on a BAM file.
    """
    sample_name = extract_sample_name(bam_file)
    output_tsv = os.path.join(config["process_output_dir"], sample_name)

    # Skip if the TSV file already exists
    if os.path.exists(f"{output_tsv}.SpliSER.tsv"):
        logging.info(f"SpliSER TSV {output_tsv}.SpliSER.tsv already exists, skipping processing.")
        return f"{output_tsv}.SpliSER.tsv"

    logging.info(f"Running SpliSER process on {bam_file}...")
    command_process = ["python", config["spliser_script"], "process"] + [
        "-B", bam_file,
        "-b", junc_file,
        "-o", output_tsv,
        "-A", config["gtf_file"]
    ] + config["spliser_process_args"]
    subprocess.run(command_process, check=True)
    logging.info(f"SpliSER process completed, output saved to {output_tsv}.SpliSER.tsv")
    return output_tsv



def process_bam_file(bam_file, config, samples_file):
    """
    Process a single BAM file: create junctions and run SpliSER process.
    """
    junc_file = create_spliser_junctions(bam_file, config)
    output_tsv = run_spliser_process(bam_file, junc_file, config)
    sample_name = extract_sample_name(bam_file)

    if junc_file and output_tsv:
        logging.info(f"Appending {bam_file} to the samples file.")
        with open(samples_file, 'a') as f:
            f.write(f"{sample_name}\t{output_tsv}.SpliSER.tsv\t{bam_file}\n")
        logging.info(f"Successfully appended {bam_file} to the samples file.")


def run_spliser_combine(config):
    """
    Run the SpliSER combine step.
    """
    create_directory(config["output_combined"])
    logging.info(f"Running SpliSER combine...")
    command_combine = ["python", config["spliser_script"], "combine"] + [
        "-S", config["samples_file"],
        "-o", config["output_combined"]
    ] + config["spliser_combine_args"]
    subprocess.run(command_combine, check=True)
    logging.info(f"SpliSER combine finished, results are in {config['output_combined']}/combined.tsv")


def run_spliser_output(config):
    """
    Run the SpliSER output step.
    """
    create_directory(config["output_diffspliser"])
    combined_file = os.path.join(config["output_combined"], "combined.tsv")
    logging.info(f"Running SpliSER output...")
    command_output = ["python", config["spliser_script"], "output"] + [
        "-S", config["samples_file"],
        "-C", combined_file,
        "-t", "DiffSpliSER",
        "-o", config["output_diffspliser"]
    ] + config["spliser_output_args"]
    subprocess.run(command_output, check=True)
    logging.info(f"SpliSER output finished, results are in {config['output_diffspliser']}")


def get_bam_files_from_tsv(tsv_file):
    """
    Get full paths for BAM files listed in the TSV file.
    """
    bam_files = []
    with open(tsv_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row:
                bam_files.extend(row)
    return bam_files


def main():
    parser = argparse.ArgumentParser(description="Run SpliSER pipeline with optional steps.")
    parser.add_argument('--config', type=str, required=True, help="Path to the configuration JSON file.")
    args = parser.parse_args()

    config = load_config(args.config)

    # Create main working directories
    create_directory(config["spliser_dir"])
    create_directory(config["junction_files_dir"])  # Create once
    create_directory(config["process_output_dir"])  # Create once

    log_file = os.path.join(config["spliser_dir"], "spliser_pipeline.log")
    logging.basicConfig(
        filename=log_file,
        filemode='a',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO,
        force=True
    )
    logging.info("SpliSER pipeline started")

    # Get BAM files to process
    bam_files_to_process = get_bam_files_from_tsv(config["bam_file_list"])

    # Process BAM files
    if not config["skip_processing"]:
        create_directory(os.path.dirname(config["samples_file"]))
        with ProcessPoolExecutor(max_workers=config["n_threads"]) as executor:
            futures = [
                executor.submit(process_bam_file, bam, config, config["samples_file"])
                for bam in bam_files_to_process
            ]
            for future in as_completed(futures):
                future.result()

    # Run combine and output steps if required
    if not config["only_process"]:
        run_spliser_combine(config)
        run_spliser_output(config)

    logging.info("SpliSER pipeline finished")


if __name__ == "__main__":
    main()


# python run_spliser.py --config spliser_config.json