import pandas as pd
import os

def compute_coverage_breadth(df_hist_col):
    '''Compute the coverage breadth from the WGS metrics historgram.'''
    genome_territory = sum(df_hist_col['count'])
    empty_territory = df_hist_col['count'][0]

    covered_count = genome_territory - empty_territory
    coverage_breadth = covered_count / genome_territory

    return coverage_breadth

def parse_wgs_metrics(cell_id, wgs_metrics_file):
    '''Parse the Picard WGS metrics file.'''
    df = pd.read_csv(wgs_metrics_file, sep='\t', comment='#')

    df_metrics = df.iloc[[0]]
    df_metrics.insert(0, 'cell_id', cell_id)
    df_metrics.columns = [x.lower() for x in df_metrics.columns]

    df_hist_col = df.iloc[2:, 0:2].astype(int).reset_index(drop=True)
    df_hist_col.columns = ['coverage', 'count']

    coverage_breadth = compute_coverage_breadth(df_hist_col)
    df_metrics.insert(2, 'coverage_breadth', coverage_breadth)

    df_hist = df_hist_col.set_index('coverage')
    df_hist.index.name = None
    df_hist = df_hist.T
    df_hist = df_hist.reset_index(drop=True)
    df_hist.insert(0, 'cell_id', cell_id)

    return df_metrics, df_hist


def collate_wgs_metrics(wgs_metrics_dict):
    '''Collate the Picard WGS metrics output from multiple cells.'''
    df_metrics_collate = pd.DataFrame()
    df_hist_collate = pd.DataFrame()

    for cell_id in wgs_metrics_dict.keys():
        df_metrics, df_hist = parse_wgs_metrics(cell_id, wgs_metrics_dict[cell_id])
        df_metrics_collate = pd.concat([df_metrics_collate, df_metrics], axis=0)
        df_hist_collate = pd.concat([df_hist_collate, df_hist], axis=0)
    return df_metrics_collate

def get_input(cell_ids, sample_id):
    paths = {}
    for cell_id in cell_ids:
        fpaths = f"output/{sample_id}/{cell_id}/0_qc/WgsMetrics/{cell_id}.wgsMetric.txt"
        paths[cell_id] = fpaths
    return paths

def main():
    input_paths = get_input(snakemake.params["cell_id_lst"], snakemake.wildcards["sample_id"])
    metric = collate_wgs_metrics(input_paths)
    metric = metric.loc[:, ["cell_id", "mean_coverage", "coverage_breadth"]]
    metric.rename(columns={"mean_coverage": "coverage_depth"}, inplace=True)
    metric.to_csv(snakemake.output[0], sep="\t", index = False) 

if __name__ == "__main__":
    main()