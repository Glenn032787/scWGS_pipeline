import pandas as pd
import os

def parse_gc_bias_metrics(cell_id, gc_bias_metrics_file):
    """Parse the Picard GC bias metrics file."""
    df = pd.read_csv(gc_bias_metrics_file, sep='\t', comment='#')

    df.columns = [x.lower() for x in df.columns]

    df.insert(0, 'cell_id', cell_id)

    return df


def collate_gc_bias_metrics(gc_bias_metrics_dict):
    '''Collate the Picard GC bias metrics output from multiple cells.'''
    df_metrics_collate = pd.DataFrame()
    for cell_id in gc_bias_metrics_dict.keys():
        df = pd.read_csv(gc_bias_metrics_dict[cell_id], sep='\t', comment='#')
        metric = parse_gc_bias_metrics(cell_id, gc_bias_metrics_dict[cell_id])
        metric = metric.loc[:, ['windows', 'normalized_coverage',  'gc']]
        metric = metric.set_index('gc')
        metric = metric.rename(columns={'normalized_coverage': cell_id, 'windows': 'reference'})
        metric = metric.T
        metric['cell_id'] = metric.index
        metric.columns = metric.columns.astype(str)
        metric = metric.replace('?', 0)
        df_metrics_collate = pd.concat([df_metrics_collate, metric], axis=0)
    return df_metrics_collate.drop_duplicates()


def get_input(cell_ids, sample_id):
    paths = {}
    for cell_id in cell_ids:
        fpaths = f"output/{sample_id}/{cell_id}/0_qc/gcBias/{cell_id}.gcBias_metric.txt"
        paths[cell_id] = fpaths
    return paths

def main():
    input_paths = get_input(snakemake.params["cell_id_lst"], snakemake.wildcards["sample_id"])
    metric = collate_gc_bias_metrics(input_paths)
    metric.to_csv(snakemake.output[0], sep=",", index = False) 

if __name__ == "__main__":
    main()