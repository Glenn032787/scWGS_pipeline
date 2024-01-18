
import pandas as pd

def getRG(wildcards, metadata_path):
    metadata = pd.read_table(metadata_path)
    metadata_col = list(metadata.columns.values)
    cell_id = wildcards.cell_id

    assert "sample_id" in metadata_col
    assert "chip_id" in metadata_col 
    assert "flowcell_id" in metadata_col
    assert "lane_number" in metadata_col
    assert "machine" in metadata_col
    assert "cell_id" in metadata_col

    curr_metadata = metadata[metadata['cell_id'] == cell_id]

    assert curr_metadata.shape[0] == 1

    sample_id = curr_metadata['sample_id'].iloc[0]
    chip_id = curr_metadata['chip_id'].iloc[0]
    flowcell_id = curr_metadata['flowcell_id'].iloc[0]
    lane_number = curr_metadata['lane_number'].iloc[0]
    machine = curr_metadata['machine'].iloc[0]

    RG_ID = sample_id + "_" + chip_id + "_" + flowcell_id + "_" + str(lane_number)
    RG_LB = chip_id
    RG_PL = machine
    RG_PU = flowcell_id + "_" + str(lane_number) + "_" + sample_id
    RG_SM = sample_id

    RG = f"@RG\\tID:{RG_ID}\\tLB:{RG_LB}\\tPL:{RG_PL}\\tPU:{RG_PU}\\tSM:{RG_SM}"

    return RG