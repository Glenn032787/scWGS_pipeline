import pandas as pd
import os
import pandas as pd

def main():
    # Import flagstat data to get total mapped reads and total reads 
    flagstat = pd.read_csv(snakemake.input["flagstat"], sep='\t', comment='#')
    flagstat = flagstat[["cell_id", "total_mapped_reads", "total_reads"]]

    # Import insert size metric to get average read length
    insertSize = pd.read_csv(snakemake.input["insertSize"], sep='\t', comment='#')
    insertSize = insertSize[["cell_id", "mean_insert_size"]]

    # Get genome length from reference genome fasta index
    fai = pd.read_csv(snakemake.input["genomeIndex"], sep = "\t", header = None)
    genomeLength = fai[1].sum()

    # Calculate actual and expected coverage
    # (numReads * readLength)/genomeLength
    result_df = pd.merge(flagstat, insertSize, on='cell_id', how='left')
    result_df["actual_coverage"] = (result_df["total_mapped_reads"] * result_df["mean_insert_size"])/genomeLength
    result_df["expected_coverage"] = (result_df["total_reads"] * result_df["mean_insert_size"])/genomeLength
    coverage_metric = result_df[["cell_id", "actual_coverage", "expected_coverage"]]

    # Save coverage metrics
    coverage_metric.to_csv(snakemake.output[0], sep = "\t", index = False)

if __name__ == "__main__":
    main()