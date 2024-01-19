import pandas as pd
import numpy as np
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


def parse_duplicate_metrics(cell_id, duplicate_metrics_file):
    '''Parse the Picard duplicate metrics file.'''
    df = pd.read_csv(duplicate_metrics_file, sep='\t', comment='#', nrows=1)
    df.columns = [x.lower() for x in df.columns.values]
    df = df.rename(columns={'percent_duplication': 'duplicate_fraction'})
    df.insert(0, 'cell_id', cell_id)
    return df

   
def collate_duplicate_metrics(duplicate_metrics_dict):
    df_collate = pd.DataFrame()

    for cell_id in duplicate_metrics_dict.keys():
        df = parse_duplicate_metrics(cell_id, duplicate_metrics_dict[cell_id])
        df_collate = pd.concat([df_collate, df], axis=0)     
    return df_collate 

def get_input(cell_ids, sample_id):
    paths = {}
    for cell_id in cell_ids:
        fpaths = f"output/{sample_id}/{cell_id}/0_qc/markDup/{cell_id}.markDup_metric.txt"
        paths[cell_id] = fpaths
    return paths

def main():
    input_paths = get_input(snakemake.params["cell_id_lst"], snakemake.wildcards["sample_id"])
    metrics = collate_duplicate_metrics(input_paths)
    
    renameCol = {"cell_id": "cell_id", 
                 "unpaired_reads_examined": "unpaired_mapped_reads", 
                 "read_pairs_examined": "paired_mapped_reads", 
                 "unpaired_read_duplicates": "unpaired_duplicate_reads", 
                 "read_pair_duplicates": "paired_duplicate_reads", 
                 "unmapped_reads": "unmapped_reads", 
                "estimated_library_size": "estimated_library_size"}
    metrics = metrics.loc[:, list(renameCol.keys()) + ["read_pair_optical_duplicates"]]
    metrics.rename(columns=renameCol, inplace=True)

    def percent_duplicate(row):
        numMapped = row["unpaired_mapped_reads"] + (row["paired_mapped_reads"] * 2)
        numDupplicate = row["unpaired_duplicate_reads"] + ((row["paired_duplicate_reads"] + row["read_pair_optical_duplicates"]) * 2) 
        return np.where(numMapped != 0, numDupplicate / numMapped, 0)
    
    metrics["percent_duplicate_reads"] = metrics.apply(percent_duplicate, axis=1)
    metrics.drop('read_pair_optical_duplicates', axis=1, inplace=True)
    metrics.to_csv(snakemake.output[0], sep="\t", index = False) 

if __name__ == "__main__":
    main()