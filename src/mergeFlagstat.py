import pandas as pd

def parse_flagstat(cell_id, flagstat_file):
    "Parse samtools flagstat metric output, code taken from modrian utilities"
    df = pd.read_csv(
        flagstat_file,
        sep=r'\s\+\s0\s',
        header=None,
        names=['value', 'type'],
        engine = "python"
    )

    tot_reads = df[df['type'] == 'in total (QC-passed reads + QC-failed reads)']['value']
    tot_mpd_reads = df[
        (df['type'].str.contains('mapped') == True) &
        (df['type'].str.contains('primary mapped') == False) &
        (df['type'].str.contains('mate mapped') == False)
        ]
    tot_dup_reads = df[df['type'] == 'duplicates']['value']
    tot_prop_paired = df[df['type'].str.contains('properly paired')]

    assert len(tot_reads) == 1
    assert len(tot_mpd_reads) == 1
    assert len(tot_dup_reads) == 1
    assert len(tot_prop_paired) == 1

    tot_reads = tot_reads.iloc[0]
    tot_mpd_reads = tot_mpd_reads['value'].iloc[0]
    tot_dup_reads = tot_dup_reads.iloc[0]
    tot_prop_paired = tot_prop_paired['value'].iloc[0]

    outdata = {
        'cell_id': cell_id,
        'total_reads': tot_reads,
        'total_mapped_reads': tot_mpd_reads,
        'total_duplicate_reads': tot_dup_reads,
        'total_properly_paired': tot_prop_paired
    }

    flagstat = pd.DataFrame.from_dict(outdata, orient='index').T
    return flagstat

def collate_flagstat_metrics(flagstat_metrics_dict):
    df_collate = pd.DataFrame()

    for cell_id in flagstat_metrics_dict.keys():
        df = parse_flagstat(cell_id, flagstat_metrics_dict[cell_id])
        df_collate = pd.concat([df_collate, df], axis=0)     
    return df_collate 

def get_input(cell_ids, sample_id):
    paths = {}
    for cell_id in cell_ids:
        fpaths = f"output/{sample_id}/{cell_id}/0_qc/flagstat/{cell_id}.flagstat.txt"
        paths[cell_id] = fpaths
    return paths

def main():
    input_paths = get_input(snakemake.params["cell_id_lst"], snakemake.wildcards["sample_id"])
    metrics = collate_flagstat_metrics(input_paths)
    metrics.to_csv(snakemake.output[0], sep="\t", index = False) 

if __name__ == "__main__":
    main()