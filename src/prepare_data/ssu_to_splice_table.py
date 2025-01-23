from tqdm import tqdm
import pandas as pd
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed


def preprocess_aggregated_sse(aggregated_sse):
    """Index aggregated_sse by chromosome and strand for faster lookup."""
    sse_dict = {}
    for (chrom, strand), group in aggregated_sse.groupby(['Region', 'Strand']):
        sse_dict[(chrom, strand)] = group
    return sse_dict


def process_row(row, sse_dict, index):
    chrom = row['Region']
    strand = row['Strand']
    tx_start = int(row['tx_start'])
    tx_end = int(row['tx_end'])


    if (chrom, strand) in sse_dict:
        filtered_sse = sse_dict[(chrom, strand)]
        filtered_sse = filtered_sse[(filtered_sse['Site'] >= tx_start) & (filtered_sse['Site'] <= tx_end)]
    else:
        return row  


    exon_end_sse = filtered_sse[filtered_sse['Annotation'] == 'exon_end']
    exon_start_sse = filtered_sse[filtered_sse['Annotation'] == 'exon_start']

    row['exon_end_positions'] = ','.join(map(str, exon_end_sse['Site'])) + ',' if not exon_end_sse.empty else ''
    row['exon_end_sse'] = ','.join(map(str, exon_end_sse['Average_SSE'])) + ',' if not exon_end_sse.empty else ''
    row['exon_start_positions'] = ','.join(map(str, exon_start_sse['Site'])) + ',' if not exon_start_sse.empty else ''
    row['exon_start_sse'] = ','.join(map(str, exon_start_sse['Average_SSE'])) + ',' if not exon_start_sse.empty else ''


    return row

def add_one_to_positions(positions):
    values = str(positions).split(',')
    values = [str(int(x) + 1) for x in values if x.isdigit()]
    return ','.join(values) + ',' if values else ''

def filter_to_annotated_sites(splice_table):
    def process_positions_and_sse(exon_positions, annotated_positions, sse_values):
        exon_positions_list = exon_positions.split(',')
        annotated_positions_list = annotated_positions.split(',')
        sse_values_list = sse_values.split(',')

        filtered_positions = []
        filtered_sse = []

        for pos, sse in zip(exon_positions_list, sse_values_list):
            if pos.strip().isdigit() and pos in annotated_positions_list:
                filtered_positions.append(pos)
                filtered_sse.append(sse)

        return ','.join(filtered_positions) + (',' if filtered_positions else ''), ','.join(filtered_sse) + (',' if filtered_sse else '')

    splice_table[['exon_end_positions', 'exon_end_sse']] = splice_table.apply(
        lambda row: process_positions_and_sse(row['exon_end_positions'], row['Exon_ends'], row['exon_end_sse']), axis=1, result_type='expand'
    )
    splice_table[['exon_start_positions', 'exon_start_sse']] = splice_table.apply(
        lambda row: process_positions_and_sse(row['exon_start_positions'], row['Exon_starts'], row['exon_start_sse']), axis=1, result_type='expand'
    )

    return splice_table

def main(splice_table_path, aggregated_sse_path, output_path, annotated_only):
    print("Reading aggregated SSE data...")
    aggregated_sse = pd.read_csv(aggregated_sse_path, sep='\t')

    print("Reading splice table data...")
    splice_table = pd.read_csv(splice_table_path, sep='\t', header=None, names=[
        'Transcript_ID', 'Paralog_status', 'Region', 'Strand',
        'tx_start', 'tx_end', 'Exon_ends', 'Exon_starts'
    ])
    splice_table[['exon_end_positions', 'exon_end_sse', 'exon_start_positions', 'exon_start_sse']] = ''

    print("Preprocessing aggregated SSE data...")
    sse_dict = preprocess_aggregated_sse(aggregated_sse)

    total_rows = len(splice_table)
    print(f"Total rows to process: {total_rows}")

    print("Processing rows in parallel...")

    with ProcessPoolExecutor() as executor:
        print(f"Number of cores being used: {executor._max_workers}")
        futures = {
            executor.submit(process_row, row, sse_dict, idx): idx 
            for idx, row in splice_table.iterrows()
        }
        with tqdm(total=total_rows, desc="Processing rows") as pbar:
            for future in as_completed(futures):
                row_index = futures[future]
                try:
                    splice_table.iloc[row_index] = future.result()
                except Exception as exc:
                    print(f"Row {row_index} generated an exception: {exc}")
                # Update progress bar
                pbar.update(1)

    print("Adding 1 to each value in exon_start_positions...")
    splice_table['exon_start_positions'] = splice_table['exon_start_positions'].apply(add_one_to_positions)

    if annotated_only:
        print("Filtering the splice table to annotated sites...")
        splice_table = filter_to_annotated_sites(splice_table)

    filtered_splice_table = splice_table[~(
        (splice_table['exon_end_positions'] == '') |
        (splice_table['exon_start_positions'] == '')
    )]

    print("Saving the updated splice table to the output path...")
    filtered_splice_table.to_csv(output_path, sep='\t', index=False, header=False)
    print(f"Updated splice table saved successfully to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine SSE data with splice table and analyze results.")
    parser.add_argument("-s", "--splice-table", help="Path to the splice table file", required=True)
    parser.add_argument("-a", "--aggregated-sse", help="Path to the aggregated SSE file", required=True)
    parser.add_argument("-o", "--output", help="Path to save the output splice table", required=True)
    parser.add_argument("--annotated-only", action="store_true", help="Only keep exon positions and SSE values that match gencode annotated positions in the splice table")
    args = parser.parse_args()
    main(args.splice_table, args.aggregated_sse, args.output, args.annotated_only)

# usage:
#  python ssu_to_splice_table.py -s ../../splice_tables/GRCh38_v29_splice_table.txt -a ./hepg2_filt_aggregated_ssu.tsv -o ../../splice_tables/hepg2_filt_rna_seq_splice_table.txt
#  python ssu_to_splice_table.py -s ../../splice_tables/GRCh38_v29_splice_table.txt -a ./hepg2_filt_aggregated_ssu.tsv -o ../../splice_tables/hepg2_filt_annotated_splice_table.txt --annotated-only
#
# python ssu_to_splice_table.py -s ../../splice_tables/GRCh38_v29_splice_table.txt -a ./k562_filt_aggregated_ssu.tsv -o ../../splice_tables/k562_filt_rna_seq_splice_table.txt
# python ssu_to_splice_table.py -s ../../splice_tables/GRCh38_v29_splice_table.txt -a ./k562_filt_aggregated_ssu.tsv -o ../../splice_tables/k562_filt_annotated_splice_table.txt --annotated-only








