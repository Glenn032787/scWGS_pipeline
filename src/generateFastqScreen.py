import os

def create_fastq_screen_config(config, output_file):
    fastq_config = {'database': dict(), 'aligner_paths': {"bowtie2": "bowtie2"}}
    for organism in config:
        fastq_config['database'][organism] = {"bowtie2": config[organism]}

    with open(output_file, 'w') as fout:
        for label, indexes in fastq_config['database'].items():
            for aligner, index in indexes.items():
                fout.write('\t'.join([
                    'DATABASE', label, index, aligner.upper()]) + '\n')
        for aligner, path in fastq_config['aligner_paths'].items():
            fout.write('\t'.join([aligner.upper(), path]) + '\n')


directory = os.path.dirname(snakemake.output[0])
if not os.path.exists(directory):
    os.makedirs(directory)
    
create_fastq_screen_config(snakemake.config["fastq_screen"], snakemake.output[0])