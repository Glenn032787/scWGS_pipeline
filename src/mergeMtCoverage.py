import pandas as pd

def parse_mt_coverage_bed(cell_id, mt_coverage_file):
    df = pd.read_csv(mt_coverage_file, sep = '\t', names = ['chr', 'start', 'end', 'mt_coverage_depth', 'n_bases_covered', 'len_coverage', 'fraction_covered'])
    df.insert(0, 'cell_id', cell_id)
    df = df[['cell_id', 'mt_coverage_depth']]
    return df

def collate_mt_coverage(mt_coverage_dict):
    df_collate = pd.DataFrame()
    for cell_id in mt_coverage_dict.keys():
        try:
            df = parse_mt_coverage_bed(cell_id, mt_coverage_dict[cell_id])
            df_collate = pd.concat([df_collate, df], axis=0)
        except pd.errors.EmptyDataError:
            print(f"{cell_id} has an empty mt coverage file")
    return df_collate
    

def get_input(cell_ids, sample_id):
    paths = {}
    for cell_id in cell_ids:
        fpaths = f"output/{sample_id}/{cell_id}/0_qc/mtDNA/{cell_id}.mtCov.bed"
        paths[cell_id] = fpaths
    return paths

def main():
    input_paths = get_input(snakemake.params["cell_id_lst"], snakemake.wildcards["sample_id"])
    metrics = collate_mt_coverage(input_paths)
    metrics.to_csv(snakemake.output[0], sep="\t", index = False) 

if __name__ == "__main__":
    main()