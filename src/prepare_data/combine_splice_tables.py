import pandas as pd
import argparse

def load_and_prepare(filepath):
    """Load dataset and prepare columns for combining."""
    df = pd.read_csv(filepath, sep='\t', header=None, names=[
        'Transcript_ID', 'Paralog_status', 'Region', 'Strand', 'tx_start', 'tx_end', 
        'Exon_ends', 'Exon_starts', 'exon_end_positions', 'exon_end_sse', 
        'exon_start_positions', 'exon_start_sse'
    ])
    df['exon_end_positions'] = df['exon_end_positions'].fillna('').apply(lambda x: [int(pos) for pos in str(x).split(',') if pos])
    df['exon_start_positions'] = df['exon_start_positions'].fillna('').apply(lambda x: [int(pos) for pos in str(x).split(',') if pos])
    df['exon_end_sse'] = df['exon_end_sse'].fillna('').apply(lambda x: [float(sse) for sse in str(x).split(',') if sse])
    df['exon_start_sse'] = df['exon_start_sse'].fillna('').apply(lambda x: [float(sse) for sse in str(x).split(',') if sse])
    return df

def merge_positions_and_sse(positions1, sse1, positions2, sse2):
    """Merge exon positions and SSE values, preserving order and counting matches and unique positions."""
    combined_positions = []
    combined_sse1 = []
    combined_sse2 = []
    match_count = 0
    unique1_count = 0
    unique2_count = 0

    # Add positions and SSE from the first dataset
    for pos, sse in zip(positions1, sse1):
        if pos not in combined_positions:
            combined_positions.append(pos)
            combined_sse1.append(sse)
            combined_sse2.append(sse2[positions2.index(pos)] if pos in positions2 else 0)  # Get matching SSE or 0
            if pos in positions2:
                match_count += 1
            else:
                unique1_count += 1

    # Add positions and SSE from the second dataset that weren't already added
    for pos, sse in zip(positions2, sse2):
        if pos not in combined_positions:
            combined_positions.append(pos)
            combined_sse1.append(0)
            combined_sse2.append(sse)
            unique2_count += 1

    return combined_positions, combined_sse1, combined_sse2, match_count, unique1_count, unique2_count

def combine_datasets(df1, df2):
    combined_df = pd.merge(
        df1, df2,
        on=['Transcript_ID', 'Region', 'Strand', 'tx_start', 'tx_end'],
        how='outer',
        suffixes=('_1', '_2')
    )

    combined_data = []
    total_end_matches = 0
    total_start_matches = 0
    unique_end1 = 0
    unique_end2 = 0
    unique_start1 = 0
    unique_start2 = 0

    for _, row in combined_df.iterrows():
        combined_ends, dataset1_end_sse, dataset2_end_sse, end_match_count, end_unique1, end_unique2 = merge_positions_and_sse(
            row['exon_end_positions_1'] if isinstance(row['exon_end_positions_1'], list) else [],
            row['exon_end_sse_1'] if isinstance(row['exon_end_sse_1'], list) else [],
            row['exon_end_positions_2'] if isinstance(row['exon_end_positions_2'], list) else [],
            row['exon_end_sse_2'] if isinstance(row['exon_end_sse_2'], list) else []
        )
        
        combined_starts, dataset1_start_sse, dataset2_start_sse, start_match_count, start_unique1, start_unique2 = merge_positions_and_sse(
            row['exon_start_positions_1'] if isinstance(row['exon_start_positions_1'], list) else [],
            row['exon_start_sse_1'] if isinstance(row['exon_start_sse_1'], list) else [],
            row['exon_start_positions_2'] if isinstance(row['exon_start_positions_2'], list) else [],
            row['exon_start_sse_2'] if isinstance(row['exon_start_sse_2'], list) else []
        )

        total_end_matches += end_match_count
        total_start_matches += start_match_count
        unique_end1 += end_unique1
        unique_end2 += end_unique2
        unique_start1 += start_unique1
        unique_start2 += start_unique2

        combined_data.append({
            'Transcript_ID': row['Transcript_ID'],
            # Ensure 'Paralog_status' is cast to int if not NaN
            'Paralog_status': int(row['Paralog_status_1']) if pd.notnull(row['Paralog_status_1']) else int(row['Paralog_status_2']),
            'Region': row['Region'],
            'Strand': row['Strand'],
            'tx_start': row['tx_start'],
            'tx_end': row['tx_end'],
            'Annotated_Exon_ends': (row['Exon_ends_1'] if pd.notnull(row['Exon_ends_1']) else row['Exon_ends_2']) + ',',
            'Annotated_Exon_starts': (row['Exon_starts_1'] if pd.notnull(row['Exon_starts_1']) else row['Exon_starts_2']) + ',',
            'Combined_Exon_ends': ','.join(map(str, combined_ends)) + ',',
            'Combined_Exon_starts': ','.join(map(str, combined_starts)) + ',',
            'Dataset1_Exon_end_sse': ','.join(map(str, dataset1_end_sse)) + ',',
            'Dataset1_Exon_start_sse': ','.join(map(str, dataset1_start_sse)) + ',',
            'Dataset2_Exon_end_sse': ','.join(map(str, dataset2_end_sse)) + ',',
            'Dataset2_Exon_start_sse': ','.join(map(str, dataset2_start_sse)) + ',',
        })

    print(f"\nTotal exon end matches: {total_end_matches}")
    print(f"Total exon start matches: {total_start_matches}")
    
    print(f"\nTotal unique exon end sites in dataset 1: {unique_end1}")
    print(f"Total unique exon end sites in dataset 2: {unique_end2}")
    
    print(f"\nTotal unique exon start sites in dataset 1: {unique_start1}")
    print(f"Total unique exon start sites in dataset 2: {unique_start2}\n")

    return pd.DataFrame(combined_data)



def main(splice_table1_path, splice_table2_path, output_path):
    print("Loading datasets...")
    df1 = load_and_prepare(splice_table1_path)
    df2 = load_and_prepare(splice_table2_path)
    
    print("Combining datasets...")
    combined_df = combine_datasets(df1, df2)
    
    print("Saving the combined dataset...")
    combined_df.to_csv(output_path, sep='\t', index=False, header=False)
    print(f"Combined splice table saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine exon positions and SSE values from two splice tables.")
    parser.add_argument("-s1", "--splice-table1", help="Path to the first splice table", required=True)
    parser.add_argument("-s2", "--splice-table2", help="Path to the second splice table", required=True)
    parser.add_argument("-o", "--output", help="Path to save the combined splice table", required=True)
    args = parser.parse_args()
    main(args.splice_table1, args.splice_table2, args.output)

   

# usage:
'''
python combine_splice_tables.py \
    -s1 ../../splice_tables/k562_filt_rna_seq_splice_table.txt \
    -s2 ../../splice_tables/hepg2_filt_rna_seq_splice_table.txt \
    -o ../../splice_tables/k562_and_hepg2_filt_rna_seq_splice_table.txt

'''

