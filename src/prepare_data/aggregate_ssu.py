import argparse
import os
import pandas as pd
import time
from multiprocessing import Pool

def process_file(file_path, min_alpha_count, log_file):
    site_data = {}
    skipped_sites = []
    alpha_count_exclusions = 0  # Counter for sites excluded due to min_alpha_count

    try:
        df = pd.read_csv(file_path, sep='\t')
        expected_columns = ['Region', 'Site', 'Strand', 'SSE', 'alpha_count', 'beta1_count', 'beta2Simple_count', 'Partners', 'Competitors']
        if not all(col in df.columns for col in expected_columns):
            print(f"Error: One or more expected columns are missing in file {file_path}")
            return site_data, skipped_sites, alpha_count_exclusions

        if df.empty:
            print(f"Warning: File {file_path} is empty, skipping...")
            return site_data, skipped_sites, alpha_count_exclusions

        for _, row in df.iterrows():
            region, site, strand = row['Region'], row['Site'], row['Strand']
            alpha_count, beta1_count, beta2Simple_count = row['alpha_count'], row['beta1_count'], row['beta2Simple_count']
            sse = row['SSE']
            
            if alpha_count < min_alpha_count:
                alpha_count_exclusions += 1
                continue  # Skip if alpha_count is below the threshold

            if sse == 0:
                with open(log_file, 'a') as log:
                    log.write(f"{region}\t{site}\t{strand}\tAlpha: {alpha_count}\tBeta1: {beta1_count}\tBeta2Simple: {beta2Simple_count}\n")

            partners = [int(p) for p in eval(row['Partners']).keys()] if isinstance(row['Partners'], str) else []
            if not partners:
                continue

            upstream_partners = [p for p in partners if p < site]
            downstream_partners = [p for p in partners if p > site]

            if upstream_partners and downstream_partners:
                skipped_sites.append((region, site, strand, file_path, "both"))
                continue

            annotation = "exon_end" if downstream_partners else "exon_start" if upstream_partners else None
            site_key = (region, site, strand)
            site_data[site_key] = {"SSE": sse, "Annotation": annotation}

        return site_data, skipped_sites, alpha_count_exclusions
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return site_data, skipped_sites, alpha_count_exclusions

def process_file_with_alpha(args):
    file_path, min_alpha_count, log_file = args
    return process_file(file_path, min_alpha_count, log_file)

def aggregate_sse(input_dir, output_file, skipped_output_file, log_file, n_threads, min_alpha_count, min_file_count):
    start_time = time.time()
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".SpliSER.tsv")]
    total_files = len(files)
    print(f"Total number of files to process: {total_files}")

    with Pool(processes=n_threads) as pool:
        results = []
        skipped_sites_all = []
        total_alpha_count_exclusions = 0  # Total sites excluded due to min_alpha_count

        args = [(f, min_alpha_count, log_file) for f in files]
        for idx, (site_result, skipped_result, alpha_exclusions) in enumerate(pool.imap_unordered(process_file_with_alpha, args), 1):
            if site_result:
                results.append(site_result)
            if skipped_result:
                skipped_sites_all.extend(skipped_result)
            total_alpha_count_exclusions += alpha_exclusions
            print(f"Processed {idx}/{total_files} files")

    # Track occurrences of each site with required alpha reads across files
    site_occurrences = {}
    site_data = {}
    for file_site_data in results:
        for site_key, site_info in file_site_data.items():
            if site_key not in site_data:
                site_data[site_key] = {"SSE": [], "Annotation": site_info["Annotation"]}
                site_occurrences[site_key] = 0
            site_data[site_key]["SSE"].append(site_info["SSE"])
            site_occurrences[site_key] += 1

    # Filter to retain only sites that appear in min_file_count or more files
    filtered_site_data = {}
    min_file_count_exclusions = 0  # Counter for sites excluded due to min_file_count
    for site_key, site_info in site_data.items():
        if site_occurrences[site_key] >= min_file_count:
            filtered_site_data[site_key] = site_info
        else:
            min_file_count_exclusions += 1

    # Calculate average SSE for each unique site
    aggregated_data = []
    for site_key, site_info in filtered_site_data.items():
        sse_values = site_info["SSE"]
        if sse_values:
            average_sse = sum(sse_values) / len(sse_values)
            aggregated_data.append((*site_key, average_sse, site_info["Annotation"]))

    if aggregated_data:
        aggregated_df = pd.DataFrame(aggregated_data, columns=['Region', 'Site', 'Strand', 'Average_SSE', 'Annotation'])
        aggregated_df.to_csv(output_file, sep='\t', index=False)
        print(f"Aggregated SSE data saved to {output_file}")

        total_sites = len(aggregated_df)
        sse_zero_count = len(aggregated_df[aggregated_df['Average_SSE'] == 0])
        sse_0_0_1_count = len(aggregated_df[(aggregated_df['Average_SSE'] > 0) & (aggregated_df['Average_SSE'] <= 0.1)])
        sse_0_1_0_9_count = len(aggregated_df[(aggregated_df['Average_SSE'] > 0.1) & (aggregated_df['Average_SSE'] < 0.9)])
        sse_0_9_1_count = len(aggregated_df[(aggregated_df['Average_SSE'] >= 0.9) & (aggregated_df['Average_SSE'] <= 1)])

        print(f"Sites with Average_SSE == 0: {sse_zero_count} ({(sse_zero_count / total_sites) * 100:.2f}%)")
        print(f"Sites with 0 < Average_SSE <= 0.1: {sse_0_0_1_count} ({(sse_0_0_1_count / total_sites) * 100:.2f}%)")
        print(f"Sites with 0.1 < Average_SSE < 0.9: {sse_0_1_0_9_count} ({(sse_0_1_0_9_count / total_sites) * 100:.2f}%)")
        print(f"Sites with 0.9 <= Average_SSE <= 1: {sse_0_9_1_count} ({(sse_0_9_1_count / total_sites) * 100:.2f}%)")
    else:
        print("No data to write to output file.")

    if skipped_sites_all:
        skipped_df = pd.DataFrame(skipped_sites_all, columns=['Region', 'Site', 'Strand', 'File', 'Reason'])
        skipped_df.to_csv(skipped_output_file, sep='\t', index=False)
        print(f"Skipped sites data saved to {skipped_output_file}")
    else:
        print("No skipped sites to write to output file.")

    print(f"Total sites filtered out due to min_alpha_count: {total_alpha_count_exclusions}")
    print(f"Total sites filtered out due to min_file_count: {min_file_count_exclusions}")
    print(f"Total time taken: {time.time() - start_time:.2f} seconds")

def main():
    parser = argparse.ArgumentParser(description="Aggregate SSE values from multiple SpliSER output files.")
    parser.add_argument('--input_dir', type=str, required=True, help="Path to the directory containing SpliSER output files.")
    parser.add_argument('--output_file', type=str, required=True, help="Path to the output file for aggregated SSE values.")
    parser.add_argument('--skipped_output_file', type=str, required=True, help="Path to the output file for skipped sites.")
    parser.add_argument('--log_file', type=str, required=True, help="Path to the log file for SSE=0 entries.")
    parser.add_argument('--n_threads', type=int, default=4, help="Number of threads for parallel processing (default: 4)")
    parser.add_argument('--min_alpha_count', type=int, default=5, help="Minimum alpha count to include a site (default: 5)")
    parser.add_argument('--min_file_count', type=int, default=2, help="Minimum files with required alpha count for a site (default: 2)")
    args = parser.parse_args()

    aggregate_sse(args.input_dir, args.output_file, args.skipped_output_file, args.log_file, args.n_threads, args.min_alpha_count, args.min_file_count)

if __name__ == "__main__":
    main()


# usage:
# python aggregate_sse.py --input_dir ../../k562 --output_file k562_filt_aggregated_ssu.tsv --skipped_output_file k562_up_and_down_partners.tsv --log_file zero_sse_log.txt --n_threads 56 --min_alpha_count 5 --min_file_count 2

# python aggregate_sse.py --input_dir ../../hepg2 --output_file hepg2_filt_aggregated_ssu.tsv --skipped_output_file hepg2_up_and_down_partners.tsv --n_threads 56 --min_alpha_count 5 --min_file_count 2 --log_file zero_sse_log.tx



