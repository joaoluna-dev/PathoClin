import os
import glob

#obtenção das amostras existentes no diretório de input
#obtém todos os arquivos vcf, sendo .vcf ou .vcf.gz
raw_files = glob.glob("data/raw/*.vcf") + glob.glob("data/raw/*.vcf.gz")

#itera sobre os arquivos obtidos anteriormente, e remove a extensão
SAMPLES = []
for f in raw_files:
    base_name = os.path.basename(f)
    if base_name.endswith(".vcf.gz"):
        SAMPLES.append(base_name.replace(".vcf.gz", ""))
    elif base_name.endswith(".vcf"):
        SAMPLES.append(base_name.replace(".vcf", ""))

#remove duplicatas, importante caso haja um arquivo1.vcf e um arquivo1.vcf.gz na mesma pasta
SAMPLES = list(set(SAMPLES))

#verifica se o arquivo de input é um .vcf ou um .vcf.gz, e passa para a o input da regra filter_vcf
def get_vcf_input(wildcards):
    gz_file = f"data/raw/{wildcards.sample}.vcf.gz"
    if os.path.exists(gz_file):
        return gz_file
    return f"data/raw/{wildcards.sample}.vcf"

#regra all: define o arquivo final para a interrupção do pipeline
rule all:
    input:
        expand("data/temp/{sample}_report.json",sample=SAMPLES)

#regra filter_vcf: passa o vcf bruto para a filtragem pelo filter.py
rule filter_vcf:
    input:
        vcf = get_vcf_input
    params:
        min_depth = 20,
        min_qual = 30,
        min_qd = 2.0,
        min_mq = 40.0,
        max_fs_snp = 60.0,
        max_fs_indel = 200.0,
        max_sor = 3.0,
        min_mq_rank_sum = -12.5,
        min_read_pos_rank_sum_snp = -8.0,
        min_read_pos_rank_sum_indel = -20.0
    output:
        output_file = "data/temp/{sample}_filtered.vcf"
    script:
        "scripts/filter.py"

#regra annotate_vcf: faz a anotação do vcf utilizando o ANNOVAR e o InterVar
rule annotate_vcf:
    input:
        "data/temp/{sample}_filtered.vcf"
    params:
        output_base = "data/temp/{sample}"
    output:
        "data/temp/{sample}.norm.vcf",
        "data/temp/{sample}.hg38_multianno.txt",
        "data/temp/{sample}.hg38_multianno.txt.intervar"
    shell:
        "./scripts/annotate.sh {input} {params.output_base}"

#regra parse_annotations: faz o parsing das anotações para criar o json anotado
rule parse_annotations:
    input:
        norm_vcf = "data/temp/{sample}.norm.vcf",
        annovar_file = "data/temp/{sample}.hg38_multianno.txt",
        intervar_file = "data/temp/{sample}.hg38_multianno.txt.intervar"
    output:
        output_json = "data/temp/{sample}_report.json"
    script:
        "scripts/parser.py"

#fim do pipeline








