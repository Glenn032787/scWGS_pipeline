import pandas as pd
import sys
pd.options.mode.chained_assignment = None  # default='warn'

CONTROLS = ["gDNA", "NCC", "NTC", "control"]

def addStatus(cell_id, sample_id, metadata_path, bam_file, fastq_path_1, fastq_path_2, reference = "human"):
    #######
    # Check if sample is control
    ######
    metadata = pd.read_table(metadata_path)
    curr_metadata = metadata[metadata['cell_id'] == cell_id]

    if curr_metadata["cell_condition"].iloc[0] in CONTROLS:
        file_path = f'output/{sample_id}/control_sample.txt'
        with open(file_path, 'a') as file:
            file.write(bam_file + '\n')
        sys.exit()

    #######
    # Check if sample has contamination 
    ######
    singleContamination = []
    multiHitContamination = []

    for fastq in [fastq_path_1, fastq_path_2]:
        fastq_screen = pd.read_table(fastq, skiprows=1)
        fastq_screen = fastq_screen.drop(fastq_screen.index[-1])

        # Find proportion of reads that hits contaimination genome only 
        contamination = fastq_screen[fastq_screen['Genome'] != reference]
        prop_hit = (contamination['#One_hit_one_genome'] + contamination['#Multiple_hits_one_genome'])/contamination['#Reads_processed']
        contamination['prop_hit'] = prop_hit
        sumOneGenomeContamination = contamination['prop_hit'].sum()
        singleContamination.append(sumOneGenomeContamination)

        # Find proportion of reads that multihit in reference genome (human) 
        multiHit = fastq_screen.filter(['Genome','#Reads_processed', "#One_hit_multiple_genomes", "Multiple_hits_multiple_genomes"])
        multiHit = multiHit[multiHit['Genome'] == 'human']
        prop_multihit = (multiHit['#One_hit_multiple_genomes'] + contamination['Multiple_hits_multiple_genomes'])/contamination['#Reads_processed']
        multiHit["prop_multihit"] = prop_multihit
        multihitHuman = multiHit["prop_multihit"].iloc[0]
        multiHitContamination.append(multihitHuman)

    # Get average between fq1 and fq2
    singleContaimination = sum(singleContamination)/len(singleContamination)
    multiHitContamination = sum(multiHitContamination)/len(multiHitContamination)

    # Contamination  when 5% of reads map exclusively to non-human genome 
    # Contamination if 25% reads map to human genome but also map with other genome  
    if singleContaimination > 0.05 or multiHitContamination > 0.25:
        file_path = f'output/{sample_id}/contamination_sample.txt'
        with open(file_path, 'a') as file:
            file.write(bam_file + '\n')
    else:
        file_path = f'output/{sample_id}/experimental_sample.txt'
        with open(file_path, 'a') as file:
            file.write(bam_file + '\n')



with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    addStatus(snakemake.wildcards["cell_id"],
          snakemake.wildcards["sample_id"],
          snakemake.config["metadata"], 
          snakemake.input[0], 
          snakemake.input[1], 
          snakemake.input[2], 
          snakemake.config["organism"])