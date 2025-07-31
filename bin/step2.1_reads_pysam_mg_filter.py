import pysam

bam_in = pysam.AlignmentFile("run30_chr16.bam", "rb")
bam_out = pysam.AlignmentFile("run30_chr16.mg95.bam", "wb", template=bam_in)

threshold = 95.0  # e.g., minimum mg score

for read in bam_in:
    if read.has_tag('mg'):
        mg_val = read.get_tag('mg')
        if mg_val >= threshold:
            bam_out.write(read)
    else:
        # keep reads without mg tag, or skip them? Your choice
        bam_out.write(read)

bam_in.close()
bam_out.close()
