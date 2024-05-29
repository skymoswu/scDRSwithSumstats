configfile: "scdrs.yaml"
import re
import pandas as pd

rule all:
    input:
        expand(
            "{fullscore}/{trait}.full_score.gz",
            fullscore = config["fullscoredir"],
            trait = config["traits"].keys())

rule prepare_sumstats:
    input:
        config["traits"].values()
    output:
        expand("sumstats/{trait}.tsv.gz", trait = config["traits"].keys())
    run:
        for k, v in config["traits"].items():
            shell("cp {v} sumstats/{k}.tsv.gz")

rule generate_magma_input:
    input:
        "sumstats/{trait}.tsv.gz"
    output:
        "magma_files/{trait}.loc.tsv",
        "magma_files/{trait}.p.txt"
    run:
        for f in input:
            print(f)
            sumstat = pd.read_table(f, sep = "\t", header = 0)
            sumstat.loc[:, ["SNP", "CHR", "BP"]].to_csv(f"magma_files/{wildcards.trait}.loc.tsv", sep = "\t", header = False, index = None)
            sumstat.loc[:, ["SNP", "P", "N"]].to_csv(f"magma_files/{wildcards.trait}.p.txt", sep = "\t", index = None)

rule generate_gene_locs:
    input:
        "magma_files/{trait}.loc.tsv"
    output:
        "magma_files/{trait}.genes.annot"
    shell:
        "magma --annotate --snp-loc {input} --gene-loc {config[gene_loc]} --out magma_files/{wildcards.trait}"

rule generate_gene_out:
    input:
        "magma_files/{trait}.p.txt",
        "magma_files/{trait}.genes.annot"
    output:
        "magma_files/{trait}.genes.out"
    shell:
        "magma --bfile {config[bfile]} --pval {input[0]} ncol=N --gene-annot {input[1]} --out magma_files/{wildcards.trait}"

rule generate_z_file:
    input:
        "magma_files/{trait}.genes.out"
    output:
        "magma_files/{trait}_z.tsv"
    shell:
        "Rscript utils/convert_p_stat.R {input} {output}"

rule munge_gs:
    input:
        "magma_files/{trait}_z.tsv"
    output:
        "gs/{trait}.gs"
    shell:
        "scdrs munge-gs --out-file {output} --zscore-file {input} --weight zscore --n-max 1000"

rule merge_gs:
    input:
        expand("gs/{trait}.gs", trait = config["traits"].keys())
    output:
        "temp/full_traits.gs"
    run:
        with open('temp/full_traits.gs', "a") as summary_file:
            summary_file.write("TRAIT	GENESET\n")
            for f in input:
                stem = re.match("gs/(.*)\.gs", f).group(1)
                summary_file.write(open(f).readlines()[1].replace("TRAIT", stem))

rule run_scdrs:
    input:
        "temp/full_traits.gs",
        config["adata"]
    output:
        expand("{fullscore}/{trait}.full_score.gz", fullscore = config["fullscoredir"], trait = config["traits"].keys())
    shell:
        '''
        scdrs compute-score \
        --h5ad-file {input[1]} \
        --h5ad-species human \
        --gs-file {input[0]} \
        --gs-species human \
        --flag-filter-data True \
        --flag-raw-count True \
        --flag-return-ctrl-raw-score False \
        --flag-return-ctrl-norm-score True \
        --out-folder {config[fullscoredir]}/
        '''