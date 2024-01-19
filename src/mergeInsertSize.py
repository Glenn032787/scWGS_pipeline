import pandas as pd
import os

def _is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def _find_section_start_end_lines(file_path, section_header):
    lines = []
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            if line.startswith(section_header):
                lines.append(i)
            elif (len(lines) == 1) and (line.strip() == ''):
                lines.append(i) 
    return lines


def parse_insert_metrics(cell_id, insert_metrics_file):
    '''Parse the Picard duplicate metrics file.'''
    metrics_positions = _find_section_start_end_lines(insert_metrics_file, 'MEDIAN_INSERT_SIZE')
    metric_nrows = (metrics_positions[1]) - (metrics_positions[0])
    metrics_df = pd.read_csv(insert_metrics_file, sep='\t', skiprows=metrics_positions[0], nrows=metric_nrows-1)
    metrics_df.insert(0, 'cell_id', cell_id)
    metrics_df.columns = [x.lower() for x in metrics_df.columns]
    
    histogram_positions = _find_section_start_end_lines(insert_metrics_file, 'insert_size')
    hist_nrows = histogram_positions[1] - histogram_positions[0]
    hist_df = pd.read_csv(insert_metrics_file, sep='\t', skiprows=histogram_positions[0], nrows=hist_nrows)
    hist_df = hist_df.astype(int).reset_index(drop=True)
    hist_df.columns = ['insert_size', 'read_count']
    hist_df = hist_df.set_index('insert_size')
    hist_df.index.name = None
    hist_df = hist_df.T
    hist_df = hist_df.reset_index(drop=True)
    hist_df.insert(0, 'cell_id', cell_id)
    
    return metrics_df, hist_df


def collate_insert_metrics(insert_metrics_dict):
    df_metrics_collate = pd.DataFrame()
    df_hist_collate = pd.DataFrame()
    for cell_id in insert_metrics_dict.keys():
        if _is_non_zero_file(insert_metrics_dict[cell_id]):
            df_metrics, df_hist = parse_insert_metrics(cell_id, insert_metrics_dict[cell_id])
            df_metrics_collate = pd.concat([df_metrics_collate, df_metrics], axis=0)
            df_hist_collate = pd.concat([df_hist_collate, df_hist], axis=0)
    return df_metrics_collate, df_hist_collate


def get_input(cell_ids, sample_id):
    paths = {}
    for cell_id in cell_ids:
        insertSizePath = f"output/{sample_id}/{cell_id}/0_qc/insertsize/{cell_id}.insertSize.txt"
        paths[cell_id] = insertSizePath
    return paths

def main():
    insert_metric_paths = get_input(snakemake.params["cell_id_lst"], snakemake.wildcards["sample_id"])
    metrics, hist = collate_insert_metrics(insert_metric_paths)
    metrics = metrics.loc[:, ["cell_id", "mean_insert_size", "median_insert_size", "standard_deviation"]]
    metrics.rename(columns={'standard_deviation': 'standard_deviation_insert_size'}, inplace=True)
    metrics.to_csv(snakemake.output["metric"], sep="\t", index = False) 
    hist.to_csv(snakemake.output["hist"], sep = "\t", index = False)

if __name__ == "__main__":
    main()