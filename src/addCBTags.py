import pysam

def add_cb_tag(bam_file, output_bam, cell_barcode):
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
            # Iterate through each read in the input BAM file
            for read in infile:
                # Add the CB tag to each read
                read.set_tag("CB", cell_barcode, value_type="Z", replace=False)
                # Write the modified read to the output BAM file
                outfile.write(read)

add_cb_tag(snakemake.input[0], snakemake.output[0], snakemake.wildcards["cell_id"])