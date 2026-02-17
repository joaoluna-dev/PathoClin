# script de filtragem do arquivo .vcf
# referências:
# Van der Auwera et al. (2013) - GATK Best Practices. 10.1002/0471250953.bi1110s43
# Roy et al. (2018) - AMP/CAP Guidelines. 10.1016/j.jmoldx.2017.11.003.

# importação dos módulos
from cyvcf2 import VCF, Writer
import sys

print("GenoLaudo - Iniciando filtragem de vcf...")
print("Referências para os valores definidos para a filtragem: \n"
      "Van der Auwera et al. (2013) - GATK Best Practices. DOI: doi.org/10.1002/0471250953.bi1110s43 \n"
      "Roy et al. (2018) - AMP/CAP Guidelines. DOI: doi.org/10.1016/j.jmoldx.2017.11.003")
# arquivos de entrada e saída
# input_file = snakemake.input.vcf
input_file = sys.argv[1]

# output_file = snakemake.output.output_file
output_file = sys.argv[2]

# parameters
# min_depth = snakemake.params.min_depth
min_depth = 20

# min_qual = snakemake.params.min_qual
min_qual = 30

# min_qd = snakemake.params.min_qd
min_qd = 2.0

# min_mq = snakemake.params.min_mq
min_mq = 40.0

# max_fs_snp = snakemake.params.max_fs_snp
max_fs_snp = 60.0

# max_fs_indel = snakemake.params.max_fs_indel
max_fs_indel = 200.0

# max_sor = snakemake.params.max_sor
max_sor = 3.0

# min_mq_rank_sum = snakemake.params.min_mq_rank_sum
min_mq_rank_sum = -12.5

# min_read_pos_rank_sum_snp = snakemake.params.min_read_pos_rank_sum_snp
min_read_pos_rank_sum_snp = -8.0

# min_read_pos_rank_sum_indel = snakemake.params.min_read_pos_rank_sum_indel
min_read_pos_rank_sum_indel = -20.0

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
    format_dp = variant.format('DP')

    if format_dp is not None:
        # se temos dados de DP da amostra, cada uma quem tem DP baixo é individualmente mascarada
        gt = variant.genotypes
        valid_samples = 0
        for i, dp_val in enumerate(format_dp):
            if dp_val[0] < min_depth:
                gt[i] = [-1] * variant.ploidy + [False]  #transforma uma amostra com DP ruim no genótipo missing ./.
            else:
                valid_samples += 1

        # se todos os genótipos forem convertidos em missing, a variante é filtrada para remoção
        if valid_samples == 0:
            filtered += 1
            continue
        variant.genotypes = gt # salva a substituição

    # caso format(DP) esteja vazio, usa como callback o campo INFO
    elif variant.INFO.get('DP') is None or variant.INFO.get('DP') < min_depth:
        filtered += 1
        continue

    # QUAL
    if variant.QUAL is None or variant.QUAL < min_qual:
        filtered += 1
        continue

    if variant.is_indel or variant.is_deletion:
        # soft filters: se a condição não estiver presente, a variante passa
        # QD
        if variant.INFO.get('QD') is not None and variant.INFO.get('QD') < min_qd:
            filtered += 1
            continue

        # MQ
        if variant.INFO.get('MQ') is not None and variant.INFO.get('MQ') < min_mq:
            filtered += 1
            continue

        # FS
        if variant.INFO.get('FS') is not None and variant.INFO.get('FS') > max_fs_indel:
            filtered += 1
            continue

        # SOR
        if variant.INFO.get('SOR') is not None and variant.INFO.get('SOR') > max_sor:
            filtered += 1
            continue

        # MQRankSum
        if variant.INFO.get('MQRankSum') is not None and variant.INFO.get('MQRankSum') < min_mq_rank_sum:
            filtered += 1
            continue

        # ReadPosRankSum
        if variant.INFO.get('ReadPosRankSum') is not None and variant.INFO.get(
                'ReadPosRankSum') < min_read_pos_rank_sum_indel:
            filtered += 1
            continue

        # se a variante passar em todas as verificações, ela é escrita no vcf de saída usando o writer
        writer.write_record(variant)
        passed += 1
    elif variant.is_snp:
        # soft filters: se a condição não estiver presente, a variante passa
        # QD
        if variant.INFO.get('QD') is not None and variant.INFO.get('QD') < min_qd:
            filtered += 1
            continue

        # MQ
        if variant.INFO.get('MQ') is not None and variant.INFO.get('MQ') < min_mq:
            filtered += 1
            continue

        # FS
        if variant.INFO.get('FS') is not None and variant.INFO.get('FS') > max_fs_snp:
            filtered += 1
            continue

        # SOR
        if variant.INFO.get('SOR') is not None and variant.INFO.get('SOR') > max_sor:
            filtered += 1
            continue

        # MQRankSum
        if variant.INFO.get('MQRankSum') is not None and variant.INFO.get('MQRankSum') < min_mq_rank_sum:
            filtered += 1
            continue

        # ReadPosRankSum
        if variant.INFO.get('ReadPosRankSum') is not None and variant.INFO.get(
                'ReadPosRankSum') < min_read_pos_rank_sum_snp:
            filtered += 1
            continue

        # se a variante passar em todas as verificações, ela é escrita no vcf de saída usando o writer
        writer.write_record(variant)
        passed += 1
    else:
        writer.write_record(variant)
        passed += 1

print("GenoLaudo - Filtragem finalizada com sucesso!")
print(f"GenoLaudo - Total de variantes processadas: {total}")
print(f"GenoLaudo - Número de variantes que passaram na verificação: {passed} ({(passed / total) * 100:.2f}%)")
print(f"GenoLaudo - Número de variantes reprovadas na verificação {filtered} ({(filtered / total) * 100:.2f}%)")

writer.close()
vcf_file.close()