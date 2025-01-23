import sys
import argparse
import pandas as pd
import gzip
import re
import os

def extract_attribute(attr_string, attr_name):
    pattern = f'{attr_name} "([^"]+)"'
    match = re.search(pattern, attr_string)
    return match.group(1) if match else None

def main(gtf_file, output_file, gene_transcript_output_file):
    paralog_file = "./paralogs_GRCh38.txt"  # Hardcoded local paralog file

    print("Reading GTF file...")
    with gzip.open(gtf_file, 'rt') as f:
        gtf_data = pd.read_csv(
            f, sep='\t', comment='#', header=None, 
            names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
        )

    print("Original columns: ", gtf_data.columns.tolist())

    # extract attributes and filter for protein-coding genes
    gtf_data_parsed = gtf_data.copy()
    for attr in ['gene_id', 'transcript_id', 'gene_name', 'gene_type', 'tag']:
        gtf_data_parsed[attr] = gtf_data_parsed['attributes'].apply(lambda x: extract_attribute(x, attr))

    protein_coding_genes = gtf_data_parsed[(gtf_data_parsed['gene_type'] == 'protein_coding') & (gtf_data_parsed['feature'].isin(['transcript', 'exon']))]
    unique_protein_coding_genes = protein_coding_genes['gene_id'].nunique()
    print(f"Number of unique protein-coding genes: {unique_protein_coding_genes}")
    print(f"Filtered protein-coding entries: {protein_coding_genes.shape[0]}")

    # prioritize transcripts
    transcripts = protein_coding_genes[protein_coding_genes['feature'] == 'transcript'].copy()
    transcripts['priority'] = transcripts['tag'].apply(
        lambda x: 1 if 'basic' in str(x) else
                  2 if 'MANE_Select' in str(x) else
                  3 if 'Ensembl_canonical' in str(x) else
                  4 if 'MANE_Plus_Clinical' in str(x) else 5
    )
    
    transcripts['selected_tag'] = transcripts['priority'].map({
        1: 'basic', 
        2: 'MANE_Select', 
        3: 'Ensembl_canonical', 
        4: 'MANE_Plus_Clinical', 
        5: 'Length fallback'
    })

    transcripts = transcripts.sort_values(['gene_id', 'priority', 'end', 'start'], ascending=[True, True, False, True])
    primary_transcripts = transcripts.drop_duplicates(subset='gene_id').reset_index(drop=True)
    print(f"Selected primary transcripts: {primary_transcripts.shape[0]}")

    selection_counts = primary_transcripts['selected_tag'].value_counts()
    print("Selection counts by tag used:")
    for tag, count in selection_counts.items():
        print(f"{tag}: {count}")


    exons = protein_coding_genes[protein_coding_genes['feature'] == 'exon']
    exons = exons[exons['transcript_id'].isin(primary_transcripts['transcript_id'])]
    
    exon_sites = exons.groupby('transcript_id').agg({
        'start': lambda x: sorted(x),
        'end': lambda x: sorted(x)
    }).reset_index()

    single_exon_transcripts_count = 0  

    def process_exons(row, count_dict):
        start_sites = row['start']
        end_sites = row['end']
        if len(start_sites) > 1:
            start_sites = start_sites[1:]  # Remove smallest start
            end_sites = end_sites[:-1]     # Remove largest end
        # Check if there's only one exon start and end left, and they match the transcript start and end
        if len(start_sites) == 1 and len(end_sites) == 1:
            if start_sites[0] == row['transcript_start'] and end_sites[0] == row['transcript_end']:
                count_dict['single_exon_count'] += 1
                return None, None  # Signal to exclude this transcript
        return ','.join(map(str, start_sites)) + ',', ','.join(map(str, end_sites)) + ','

    # add transcript start and end columns to exon_sites for matching
    primary_transcripts = primary_transcripts.rename(columns={'start': 'transcript_start', 'end': 'transcript_end'})
    exon_sites = exon_sites.merge(primary_transcripts[['transcript_id', 'transcript_start', 'transcript_end']], on='transcript_id')

    # dictionary to keep track of the single exon count
    count_dict = {'single_exon_count': 0}

    exon_sites[['exon_start_sites', 'exon_end_sites']] = exon_sites.apply(lambda row: process_exons(row, count_dict), axis=1, result_type="expand")

    # remove transcripts where exon processing returned None (single-exon cases matching start and end)
    exon_sites = exon_sites.dropna(subset=['exon_start_sites', 'exon_end_sites'])

    primary_transcripts_with_exons = primary_transcripts.merge(
        exon_sites[['transcript_id', 'exon_start_sites', 'exon_end_sites']], on='transcript_id', how='inner'
    )

    final_df = primary_transcripts_with_exons[['gene_id', 'seqname', 'strand', 'transcript_start', 'transcript_end', 'exon_start_sites', 'exon_end_sites', 'gene_name']].copy()

    print("Loading paralog data...")
    paralog_data = pd.read_csv(paralog_file, sep='\t')

    # strip version suffix from 'gene_id' in final_df and 'Gene stable ID version' in paralog_data for matching
    final_df['gene_id_no_version'] = final_df['gene_id'].str.split('.').str[0]
    paralog_data['Gene stable ID no version'] = paralog_data['Gene stable ID version'].str.split('.').str[0]

    # assign paralog status based on whether the gene_id without version matches in the paralog data
    final_df['paralog_status'] = final_df['gene_id_no_version'].isin(
        paralog_data.loc[paralog_data['Human paralogue gene stable ID'].notna(), 'Gene stable ID no version']
    ).astype(int)


    num_with_paralog = final_df['paralog_status'].sum()
    num_without_paralog = len(final_df) - num_with_paralog
    print(f"Number of genes with a paralog: {num_with_paralog}")
    print(f"Number of genes without a paralog: {num_without_paralog}")

    # drop column 'gene_id_no_version' before saving
    final_df_ordered = final_df[['gene_id', 'paralog_status', 'seqname', 'strand', 'transcript_start', 'transcript_end', 'exon_end_sites', 'exon_start_sites']]
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # count where the first exon start equals the transcription start
    matching_start_count = final_df_ordered.apply(
        lambda row: int(row['exon_start_sites'].split(',')[0]) == row['transcript_start'], axis=1
    ).sum()
    print(f"Number of transcripts where the first exon start equals the transcription start: {matching_start_count}")
    
    # count where the greatest exon end equals the transcription end
    matching_end_count = final_df_ordered.apply(
        lambda row: int(row['exon_end_sites'].split(',')[-2]) == row['transcript_end'], axis=1  # Use -2 to account for trailing comma
    ).sum()
    print(f"Number of transcripts where the largest exon end equals the transcription end: {matching_end_count}")
    
    #  count of single-exon transcripts that were excluded
    print(f"Number of single-exon transcripts matching transcription start and end that were excluded: {count_dict['single_exon_count']}")
    
    #  include only chromosomes with 'chr' prefix followed by 1-22, X, or Y
    valid_chromosomes = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    initial_row_count = len(final_df_ordered)
    final_df_ordered = final_df_ordered[final_df_ordered['seqname'].isin(valid_chromosomes)]
    excluded_rows_count = initial_row_count - len(final_df_ordered)
    print(f"Number of rows excluded based on chromosome filter: {excluded_rows_count}")
    

    total_rows_to_save = len(final_df_ordered)
    print(f"\nTotal number of rows being saved to the output file: {total_rows_to_save}")
    

    final_df_ordered.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"File saved as {output_file}")

    #save gene-to-transcript mapping file
    gene_transcript_mapping = primary_transcripts[['gene_id', 'transcript_id', 'gene_name', 'selected_tag']].copy()
    gene_transcript_mapping.to_csv(gene_transcript_output_file, sep='\t', index=False, header=True)
    print(f"Gene-to-transcript mapping file saved as {gene_transcript_output_file}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GTF file and check paralog status using a local paralog file.")
    parser.add_argument('-g', '--gtf_file', type=str, required=True, help='Path to the input GTF file (gzipped).')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to the output file.')
    parser.add_argument('-m', '--gene_transcript_output_file', type=str, required=True, help='Path to save the gene-to-transcript mapping file.')
    args = parser.parse_args()

    main(args.gtf_file, args.output_file, args.gene_transcript_output_file)



# usage:
# python make_splice_table.py -g ../../GRCh38/gencode.v29.primary_assembly.annotation.gtf.gz -o ../../splice_tables/GRCh38_v29_splice_table.txt -m gene_2_tx.tsv


