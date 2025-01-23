import os
import requests
import json
import yaml


def get(resource):

    base_url = 'https://www.encodeproject.org/{}?format=json'
    headers = {'accept': 'application/json'}
    response = requests.get(base_url.format(resource), headers=headers)
    return response.json()


def get_rbp(exp_acc):

    response = get(os.path.join('experiments', exp_acc))
    return response.get('target', {}).get('label', 'UnknownRBP')


def get_file_accessions(exp_acc, out_dir='input_files'):

    # Create the output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)

    rbp = get_rbp(exp_acc)
    response = get(os.path.join('experiments', exp_acc))
    
    with open(os.path.join(out_dir, '..', exp_acc+'_configuration.json'), 'w') as f:
        json.dump(response, f, indent=2)
    with open(os.path.join(out_dir, '..', exp_acc+'_configuration.yaml'), 'w') as f:
        yaml.dump(response, f, sort_keys=False)

    # Function to download a file
    def download_file(accession, file_name, sample_type):
        url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.fastq.gz"
        dir_path = os.path.join(out_dir, rbp, sample_type)
        file_path = os.path.join(dir_path, file_name)

        # Create directories if they don't exist
        os.makedirs(dir_path, exist_ok=True)

        # Downloading the file
        print(f"Downloading {file_name} to {file_path}...")
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(file_path, 'wb') as f:
                for chunk in response.iter_content(1024):
                    f.write(chunk)
        else:
            print(f"Failed to download {file_name}. Status code: {response.status_code}")

    # Organize files by replicate and direction, including type (KD or control)
    files_by_replicate = {}
    pair_by_replicate = {}
    for file in response['files']:
        if file['file_type'] == 'fastq' and 'replicate' in file:
            replicate_no = file['replicate']['biological_replicate_number']
            paired_end = file.get('paired_end', '')
            direction = 'forward' if paired_end == '1' else 'reverse'
            file_type = 'KD' if 'control' not in file['dataset'] else 'control'

            if replicate_no not in files_by_replicate:
                files_by_replicate[replicate_no] = {'KD': {'forward': [], 'reverse': []},
                                                    'control': {'forward': [], 'reverse': []}}
                pair_by_replicate[replicate_no] = {'KD': {'forward': [], 'reverse': []},
                                                'control': {'forward': [], 'reverse': []}}
            files_by_replicate[replicate_no][file_type][direction].append(file['accession'])
            pair_by_replicate[replicate_no][file_type][direction].append(file['paired_with'].split('/')[2])

    # Process control samples
    for ctrl_acc in response.get('possible_controls', []):
        ctrl_exp = get(ctrl_acc['@id'])
        for file in ctrl_exp['files']:
            if file['file_type'] == 'fastq' and 'replicate' in file:
                replicate_no = file['replicate']['biological_replicate_number']
                tech_replicate_no = file['replicate']['technical_replicate_number']
                paired_end = file.get('paired_end', '')
                direction = 'forward' if paired_end == '1' else 'reverse'
                file_type = 'control'
                if replicate_no not in files_by_replicate:
                    files_by_replicate[replicate_no] = {'KD': {'forward': [], 'reverse': []},
                                                        'control': {'forward': [], 'reverse': []}}
                    pair_by_replicate[replicate_no] = {'KD': {'forward': [], 'reverse': []},
                                                'control': {'forward': [], 'reverse': []}}
                files_by_replicate[replicate_no][file_type][direction].append(file['accession'])
                pair_by_replicate[replicate_no][file_type][direction].append(file['paired_with'].split('/')[2])
                
    # make sure the files are stacked in the right order if there are multiple fastq files        
    for replicate in files_by_replicate.keys():
        for file_type in files_by_replicate[replicate].keys():
            files_by_replicate[replicate][file_type]['reverse'] = pair_by_replicate[replicate][file_type]['forward']

    # Pairing forward and reverse files and triggering download
    for replicate_no, types in files_by_replicate.items():
        for type_tag, directions in types.items():
            forward_files = directions['forward']
            reverse_files = directions['reverse']
            pair_counter = 1  # Initialize a counter for each replicate

            for fwd_acc, rev_acc in zip(forward_files, reverse_files):
                fwd_file_name = f"{rbp}_{type_tag}_{replicate_no}_forward_{pair_counter}.fastq.gz"
                rev_file_name = f"{rbp}_{type_tag}_{replicate_no}_reverse_{pair_counter}.fastq.gz"
                download_file(fwd_acc, fwd_file_name, type_tag)
                download_file(rev_acc, rev_file_name, type_tag)
                pair_counter += 1  # Increment the counter for the next pair
