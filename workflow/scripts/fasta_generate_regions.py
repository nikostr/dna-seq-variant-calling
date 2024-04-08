from math import ceil
import os


def generate_regions(fasta_index_file, size, chunks=False, chromosomes=None, bed_files=None):

    if not fasta_index_file.endswith(".fai"):
        fasta_index_file = fasta_index_file + ".fai"

    with open(fasta_index_file, "r") as fasta_index:
        for line in fasta_index:
            fields = line.strip().split("\t")
            chrom_name = fields[0]
            chrom_length = int(fields[1])
            if chromosomes is not None and chrom_name not in chromosomes:
                continue
            region_start = 0
            if chunks is True:
                region_size = ceil(chrom_length / size)  # have to make sure this works
            else:
                region_size = size
            while region_start < chrom_length:
                region_end = region_start + region_size
                if region_end > chrom_length:
                    region_end = chrom_length
                start = str(region_start)
                end = str(region_end)
                if bed_files is not None:
                    region = str(ceil(region_end / region_size))
                    file_path = f"{bed_files}.{chrom_name}.region.{region}.bed"
                    os.makedirs(os.path.dirname(file_path), exist_ok=True)
                    with open(file_path, "w") as f:
                        f.write("\t".join([chrom_name, start, end]))
                else:
                    print(f"{chrom_name}:{start}-{end}")
                region_start = region_end

generate_regions(snakemake.input.ref_idx[0], snakemake.params.chunks, chunks=True, bed_files="results/freebayes/regions/genome")
