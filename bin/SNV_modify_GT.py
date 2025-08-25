#!/usr/bin/env python3

from cyvcf2 import VCF, Writer
import sys
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Adjust SNV genotypes within a deletion region.")
parser.add_argument("--sv_vcf", required=True, help="Input SV VCF file (bgzipped and indexed).")
parser.add_argument("--snv_vcf", required=True, help="Input SNV VCF file (bgzipped and indexed).")
parser.add_argument("--output_vcf", required=True, help="Output VCF file with adjusted genotypes.")
args = parser.parse_args()

# Region to check for overlap
query_chrom = "chr16"
query_start = 89735549
query_end = 89818647

# Step 1: Parse SV VCF to find overlapping deletion
sv_vcf = VCF(args.sv_vcf)
region_start = None
region_end = None

for variant in sv_vcf(f"{query_chrom}:{query_start}-{query_end}"):
    if variant.is_sv and variant.INFO.get("SVTYPE") == "DEL":
        sv_start = variant.POS
        sv_end = int(variant.INFO.get("END", 0))
        if sv_start <= query_end and sv_end >= query_start:
            region_start = sv_start
            region_end = sv_end
            print(f"Found deletion: {query_chrom}:{region_start}-{region_end}")
            break

sv_vcf.close()

if region_start is None or region_end is None:
    print("No overlapping deletion found in SV VCF.")
    sys.exit(1)

# Step 2: Modify SNV VCF based on detected region
snv_vcf = VCF(args.snv_vcf)
out = open(args.output_vcf, "w")

# Write header
for line in snv_vcf.raw_header.strip().split("\n"):
    out.write(line + "\n")

# Process variants
for variant in snv_vcf:
    chrom = variant.CHROM
    pos = variant.POS
    fields = str(variant).strip().split("\t")

    if chrom == query_chrom and region_start <= pos <= region_end:
        format_keys = fields[8].split(":")
        if "GT" in format_keys:
            gt_index = format_keys.index("GT")
            for i in range(9, len(fields)):
                sample_data = fields[i].split(":")
                gt = sample_data[gt_index]
                if gt in ("0/0", "0|0", "1/1", "1|1"):
                    sample_data[gt_index] = "0/1"
                    fields[i] = ":".join(sample_data)

    out.write("\t".join(fields) + "\n")

out.close()
snv_vcf.close()
