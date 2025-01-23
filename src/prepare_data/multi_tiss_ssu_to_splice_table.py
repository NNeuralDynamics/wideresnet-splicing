import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor

# Function to process each row of the splice table for a specific tissue
def process_row(row, aggregated_sse, index, tissue_suffix):
    chrom = row['Region']
    strand = row['Strand']
    tx_start = row['tx_start']
    tx_end = row['tx_end']

    # Filter aggregated_sse DataFrame for matching chromosome, strand, and within transcript boundaries
    filtered_sse = aggregated_sse[(aggregated_sse['Region'] == chrom) &
                                  (aggregated_sse['Strand'] == strand) &
                                  (aggregated_sse['Site'] >= tx_start) &
                                  (aggregated_sse['Site'] <= tx_end)]

    # Separate exon_end and exon_start annotations
    exon_end_sse = filtered_sse[filtered_sse['Annotation'] == 'exon_end']
    exon_start_sse = filtered_sse[filtered_sse['Annotation'] == 'exon_start']

    # Create comma-separated strings for splice sites and SSE values
    exon_end_positions = ','.join(map(str, exon_end_sse['Site'])) + ',' if not exon_end_sse.empty else ''
    exon_end_values = ','.join(map(str, exon_end_sse['Average_SSE'])) + ',' if not exon_end_sse.empty else ''
    exon_start_positions = ','.join(map(str, exon_start_sse['Site'])) + ',' if not exon_start_sse.empty else ''
    exon_start_values = ','.join(map(str, exon_start_sse['Average_SSE'])) + ',' if not exon_start_sse.empty else ''

    # Log progress every 500 rows
    if index % 500 == 0:
        print(f"Processed {index} rows for tissue {tissue_suffix}...")

    return row

# Function to add missing positions from one tissue to the other, setting SSE to 0 for missing values
def synchronize_tissues_positions(row):
    def synchronize_positions(pos1, sse1, pos2, sse2):
        # Remove trailing commas and convert comma-separated strings to lists, handling empty strings gracefully
        pos1_list = list(map(int, pos1.rstrip(',').split(','))) if pos1 else []
        pos2_list = list(map(int, pos2.rstrip(',').split(','))) if pos2 else []
    
        sse1_list = list(map(float, sse1.rstrip(',').split(','))) if sse1 else []
        sse2_list = list(map(float, sse2.rstrip(',').split(','))) if sse2 else []
    
        # Ensure all lists are synchronized to the combined position set
        combined_positions = sorted(set(pos1_list + pos2_list))
        
        synced_sse1 = []
        synced_sse2 = []
    
        for pos in combined_positions:
            # SSE for tissue 1
            if pos in pos1_list:
                synced_sse1.append(sse1_list[pos1_list.index(pos)])
            else:
                synced_sse1.append(0.0)
    
            # SSE for tissue 2
            if pos in pos2_list:
                synced_sse2.append(sse2_list[pos2_list.index(pos)])
            else:
                synced_sse2.append(0.0)
    
        # Only join the lists if there are actual values to avoid trailing commas
        pos_synced = ','.join(map(str, combined_positions)) if combined_positions else ''
        sse1_synced = ','.join(map(str, synced_sse1)) if synced_sse1 else ''
        sse2_synced = ','.join(map(str, synced_sse2)) if synced_sse2 else ''
    
        return pos_synced, sse1_synced, sse2_synced


    # Synchronize exon_end positions and SSE values
    row['exon_end_positions'], row['exon_end_sse_t1'], row['exon_end_sse_t2'] = synchronize_positions(
        row['exon_end_positions_t1'], row['exon_end_sse_t1'], row['exon_end_positions_t2'], row['exon_end_sse_t2']
    )

    # Synchronize exon_start positions and SSE values
    row['exon_start_positions'], row['exon_start_sse_t1'], row['exon_start_sse_t2'] = synchronize_positions(
        row['exon_start_positions_t1'], row['exon_start_sse_t1'], row['exon_start_positions_t2'], row['exon_start_sse_t2']
    )

    # Drop the intermediate tissue-specific position columns after synchronization
    row.drop(['exon_end_positions_t1', 'exon_end_positions_t2', 
              'exon_start_positions_t1', 'exon_start_positions_t2'], inplace=True)

    return row

# Main function
def main(splice_table_path, aggregated_sse_path_t1, aggregated_sse_path_t2, output_path, annotated_only):
    print("Reading aggregated SSE data for Tissue 1...")
    aggregated_sse_t1 = pd.read_csv(aggregated_sse_path_t1, sep='\t')

    print("Reading aggregated SSE data for Tissue 2...")
    aggregated_sse_t2 = pd.read_csv(aggregated_sse_path_t2, sep='\t')

    print("Reading splice table data...")
    splice_table = pd.read_csv(splice_table_path, sep='\t', header=None)

    splice_table.columns = [
        'Transcript_ID', 'Paralog_status', 'Region', 'Strand',
        'tx_start', 'tx_end', 'Exon_ends', 'Exon_starts'
    ]

    print("Processing rows in parallel for Tissue 1...")
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(process_row, row, aggregated_sse_t1, idx, 't1') 
            for idx, row in splice_table.iterrows()
        ]
        updated_rows_t1 = [future.result() for future in futures]

    print("Processing rows in parallel for Tissue 2...")
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(process_row, row, aggregated_sse_t2, idx, 't2') 
            for idx, row in pd.DataFrame(updated_rows_t1).iterrows()
        ]
        updated_rows_t2 = [future.result() for future in futures]

    updated_splice_table = pd.DataFrame(updated_rows_t2)

    print("Synchronizing exon positions between tissues...")
    updated_splice_table = updated_splice_table.apply(synchronize_tissues_positions, axis=1)

    if annotated_only:
        print("Filtering the splice table to annotated sites...")
        updated_splice_table = filter_to_annotated_sites(updated_splice_table)

    print("Filtering out rows with empty new columns...")
    filtered_splice_table = updated_splice_table[~(
        (updated_splice_table['exon_end_positions'] == '') |
        (updated_splice_table['exon_end_sse_t1'] == '') |
        (updated_splice_table['exon_start_positions'] == '') |
        (updated_splice_table['exon_start_sse_t1'] == '') |
        (updated_splice_table['exon_end_sse_t2'] == '') |
        (updated_splice_table['exon_start_sse_t2'] == '')
    )]

    print("Saving the updated splice table to the output path...")
    filtered_splice_table.to_csv(output_path, sep='\t', index=False, header=False)
    print(f"Updated splice table saved successfully to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine SSE data from two tissues with splice table and analyze results.")
    parser.add_argument("-s", "--splice-table", help="Path to the splice table file", required=True)
    parser.add_argument("-a1", "--aggregated-sse-t1", help="Path to the aggregated SSE file for Tissue 1", required=True)
    parser.add_argument("-a2", "--aggregated-sse-t2", help="Path to the aggregated SSE file for Tissue 2", required=True)
    parser.add_argument("-o", "--output", help="Path to save the output splice table", required=True)
    parser.add_argument("--annotated-only", action="store_true", help="Only keep exon positions and SSE values that match annotated positions in the splice table")

    args = parser.parse_args()
    main(args.splice_table, args.aggregated_sse_t1, args.aggregated_sse_t2, args.output, args.annotated_only)
