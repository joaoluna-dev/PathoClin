# importação dos módulos
from cyvcf2 import VCF, Writer

# referências científicas:
#AUWERA, G. A. et al. From FastQ Data to High‐Confidence Variant Calls:
# The Genome Analysis Toolkit Best Practices Pipeline. Current Protocols in Bioinformatics,
# out. 2013. v. 43, n. 1. doi.org/10.1002/0471250953.bi1110s43
#
# ROY, S. et al. Standards and Guidelines for Validating Next-Generation Sequencing Bioinformatics Pipelines:
# A Joint Recommendation of the Association for Molecular Pathology and the College of American Pathologists.
# The Journal of molecular diagnostics: JMD, 1 jan. 2018. v. 20, n. 1, p. 4–27.
# Disponível em: <https://www.ncbi.nlm.nih.gov/pubmed/29154853>. doi.org/10.1016/j.jmoldx.2017.11.003.


print("==============================================================================")
print("GenoLaudo - Iniciando filtragem de vcf...")
print(
    "Referências para para a filtragem: \n"
    "AUWERA, G. A. et al. From FastQ Data to High‐Confidence Variant Calls: \n "
    "The Genome Analysis Toolkit Best Practices Pipeline. Current Protocols in Bioinformatics, out. 2013. v. 43, n. 1. \n"
    "\n ROY, S. et al. Standards and Guidelines for Validating Next-Generation Sequencing Bioinformatics Pipelines: \n "
    "A Joint Recommendation of the Association for Molecular Pathology and the College of American Pathologists. \n "
    "The Journal of molecular diagnostics: JMD, 1 jan. 2018. v. 20, n. 1, p. 4–27. \n "
    "Disponível em: <https://www.ncbi.nlm.nih.gov/pubmed/29154853>. "
)
print("==============================================================================")
# arquivos de entrada e saída
input_file = snakemake.input.vcf
output_file = snakemake.output.output_file

# parâmetros
min_depth = snakemake.params.min_depth
min_qual = snakemake.params.min_qual
min_qd = snakemake.params.min_qd
min_mq = snakemake.params.min_mq
max_fs_snp = snakemake.params.max_fs_snp
max_fs_indel = snakemake.params.max_fs_indel
max_sor = snakemake.params.max_sor
min_mq_rank_sum = snakemake.params.min_mq_rank_sum
min_read_pos_rank_sum_snp = snakemake.params.min_read_pos_rank_sum_snp
min_read_pos_rank_sum_indel = snakemake.params.min_read_pos_rank_sum_indel

# criação do iterável a partir do vcf de entrada, e o escritor que criará o vcf de saída
print("GenoLaudo - Lendo arquivo de entrada...")
vcf_file = VCF(input_file)

print("GenoLaudo - Criando arquivo de saida...")
writer = w = Writer(output_file, vcf_file)

# contadores para indicar o usuário a quantidade de variantes filtradas pelo sistema
total = 0
passed = 0
filtered = 0

print("GenoLaudo - Iniciando filtragem das variantes...")
for variant in vcf_file:
    total += 1
    # hard filters: se a condição não está presente, a variante é filtrada do resultado final
    # FILTER
    if variant.FILTER and variant.FILTER != "PASS" and variant.FILTER != ".":
        filtered += 1
        continue

    # DP
    try:
        format_dp = variant.format("DP")
    except KeyError:
        format_dp = None

    if format_dp is not None:
        # se temos dados de DP da amostra, cada uma quem tem DP baixo é individualmente mascarada
        gt = variant.genotypes
        valid_samples = 0
        for i, dp_val in enumerate(format_dp):
            if dp_val[0] < min_depth:
                gt[i] = [-1] * variant.ploidy + [
                    False
                ]  # transforma uma amostra com DP ruim no genótipo missing ./.
            else:
                valid_samples += 1

        # se todos os genótipos forem convertidos em missing, a variante é filtrada para remoção
        if valid_samples == 0:
            filtered += 1
            continue
        variant.genotypes = gt  # salva a substituição

    # caso format(DP) esteja vazio, usa como callback o campo INFO
    elif variant.INFO.get("DP") is None or variant.INFO.get("DP") < min_depth:
        filtered += 1
        continue

    # QUAL
    if variant.QUAL is None or variant.QUAL < min_qual:
        filtered += 1
        continue

    if variant.is_indel or variant.is_deletion:
        # soft filters: se a condição não estiver presente, a variante passa
        # QD
        if (
            variant.INFO.get("QD") is not None
            and variant.INFO.get("QD") < min_qd
        ):
            filtered += 1
            continue

        # MQ
        if (
            variant.INFO.get("MQ") is not None
            and variant.INFO.get("MQ") < min_mq
        ):
            filtered += 1
            continue

        # FS
        if (
            variant.INFO.get("FS") is not None
            and variant.INFO.get("FS") > max_fs_indel
        ):
            filtered += 1
            continue

        # SOR
        if (
            variant.INFO.get("SOR") is not None
            and variant.INFO.get("SOR") > max_sor
        ):
            filtered += 1
            continue

        # MQRankSum
        if (
            variant.INFO.get("MQRankSum") is not None
            and variant.INFO.get("MQRankSum") < min_mq_rank_sum
        ):
            filtered += 1
            continue

        # ReadPosRankSum
        if (
            variant.INFO.get("ReadPosRankSum") is not None
            and variant.INFO.get("ReadPosRankSum")
            < min_read_pos_rank_sum_indel
        ):
            filtered += 1
            continue

        # se a variante passar em todas as verificações, ela é escrita no vcf de saída usando o writer
        writer.write_record(variant)
        passed += 1
    elif variant.is_snp:
        # soft filters: se a condição não estiver presente, a variante passa
        # QD
        if (
            variant.INFO.get("QD") is not None
            and variant.INFO.get("QD") < min_qd
        ):
            filtered += 1
            continue

        # MQ
        if (
            variant.INFO.get("MQ") is not None
            and variant.INFO.get("MQ") < min_mq
        ):
            filtered += 1
            continue

        # FS
        if (
            variant.INFO.get("FS") is not None
            and variant.INFO.get("FS") > max_fs_snp
        ):
            filtered += 1
            continue

        # SOR
        if (
            variant.INFO.get("SOR") is not None
            and variant.INFO.get("SOR") > max_sor
        ):
            filtered += 1
            continue

        # MQRankSum
        if (
            variant.INFO.get("MQRankSum") is not None
            and variant.INFO.get("MQRankSum") < min_mq_rank_sum
        ):
            filtered += 1
            continue

        # ReadPosRankSum
        if (
            variant.INFO.get("ReadPosRankSum") is not None
            and variant.INFO.get("ReadPosRankSum") < min_read_pos_rank_sum_snp
        ):
            filtered += 1
            continue

        # se a variante passar em todas as verificações, ela é escrita no vcf de saída usando o writer
        writer.write_record(variant)
        passed += 1
    else:
        writer.write_record(variant)
        passed += 1

print("=========================================================================================================")
print("GenoLaudo - Filtragem finalizada com sucesso!")
print(f"GenoLaudo - Total de variantes processadas: {total}")
print(f"GenoLaudo - Número de variantes que passaram na verificação: {passed} ({(passed / total) * 100:.2f}%)")
print(f"GenoLaudo - Número de variantes reprovadas na verificação {filtered} ({(filtered / total) * 100:.2f}%)")
print("=========================================================================================================")

writer.close()
vcf_file.close()
