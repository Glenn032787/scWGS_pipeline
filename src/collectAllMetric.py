import pandas as pd

def main():
    file_paths = snakemake.input
    merged_df = pd.DataFrame()
    for file_path in file_paths:
        # Read each TSV file into a DataFrame
        metric = pd.read_csv(file_path, sep='\t')
        # Check if "cell_id" column exists in the DataFrame
        if 'cell_id' not in metric.columns:
            continue
        
        # Merge metrics
        if 'cell_id' not in merged_df.columns:
            merged_df = metric
        else:
            merged_df = pd.merge(merged_df, metric, on='cell_id', how='left')
            
    merged_df.to_csv(snakemake.output[0], sep="\t", index = False) 

if __name__ == "__main__":
    main()