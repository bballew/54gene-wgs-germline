# from snakemake.utils import validate
# import pandas as pd


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"


##### load config and sample sheets #####


# configfile: "config/config.yaml"


# validate(config, schema="../schemas/config.schema.yaml")

# samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
# samples.index.names = ["sample_id"]
# validate(samples, schema="../schemas/samples.schema.yaml")


rule get_resources:
    """Programmatically retrieve Broad resources from an AWS s3 bucket.
    May want to allow a switch to their GCP bucket, depending on
    which cloud provider we're using.
    """
    output:
        "resources/Homo_sapiens_assembly38.fasta",
        "resources/Homo_sapiens_assembly38.dict",
        "resources/Homo_sapiens_assembly38.fasta.64.alt",
        "resources/Homo_sapiens_assembly38.fasta.64.amb",
        "resources/Homo_sapiens_assembly38.fasta.64.ann",
        "resources/Homo_sapiens_assembly38.fasta.64.bwt",
        "resources/Homo_sapiens_assembly38.fasta.64.pac",
        "resources/Homo_sapiens_assembly38.fasta.64.sa",
        "resources/Homo_sapiens_assembly38.fasta.fai",
        "resources/Homo_sapiens_assembly38.known_indels.vcf.gz",
        "resources/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi",
        "resources/Homo_sapiens_assembly38.dbsnp138.vcf",
        "resources/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
    benchmark:
        "results/performance_benchmarks/get_resources/benchmarks.tsv"
    conda:
        "../envs/aws.yaml"
    shell:
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dict resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf resources/ --no-sign-request && "
        "aws s3 cp s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx resources/ --no-sign-request"
