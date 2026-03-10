# script de parsing dos dados de anotação

import sys
import os
import json
import gc
import pandas as pd
from cyvcf2 import VCF


# referências científicas
# RICHARDS, S. et al. Standards and guidelines for the interpretation of sequence variants:
# a joint consensus recommendation of the American College of Medical Genetics and Genomics and the
# Association for Molecular Pathology. Genetics in medicine :
# official journal of the American College of Medical Genetics, 2015. v. 17, n. 5, p. 405–24.


def parse_aachange(aachange):
    """
    Divide o AAChange em Transcript e Variant (c. / p.) segundo a ACMG.
    Param:aachange: recebe a notação HGVS do intervar
    return: string formatada em transcript e variant
    """

    # verifica se o campo é nulo ou não, caso sim, atribui o valor unknown para Transcript e Variant
    if pd.isna(aachange) or aachange == "." or aachange == "":
        return "Unknown", "Unknown"

    # splitting da string de aachange para obter a notação HGVS
    parts = str(aachange).split(",")[0].split(":")
    if len(parts) >= 5:
        transcript = parts[1]
        # Monta a string "c.XXX (p.XXX)" ou apenas "c.XXX"
        variant_hgvs = (
            f"{parts[3]} ({parts[4]})" if len(parts) > 4 else f"{parts[3]}"
        )
        return transcript, variant_hgvs

    return "Unknown", str(aachange)


def parse_disease(clndn_str):
    """
    processa a coluna CLNDN do ClinVar para obtenção da doença/fenótipo.
    param clndn_str: recebe a string do clinvar (CLNDN) do arquivo do annovar
    return: string formatada das doenças (separadas por " | ")
    """
    # verifica se o valor é válido ou nulo
    if pd.isna(clndn_str) or clndn_str == ".":
        return "Not provided"

    # recebe a coluna CLNDN do clinvar, e divide por |, para o caso de mais de uma doença
    parts = str(clndn_str).split("|")
    valid_diseases = []

    # itera sobre os elementos divididos pelo |
    for p in parts:
        clean_name = p.replace("_", " ").strip()  # remove caracteres _
        # filtra dados que não são relevantes
        if clean_name.lower() not in ["not provided", "not specified", ""]:
            # adiciona as doenças válidas à lista valid_diseases
            valid_diseases.append(clean_name)

    # caso nenhuma doença válida, retorna "not provided"
    if not valid_diseases:
        return "Not provided"

    # reúne os nomes de doenças válidas, ordena os nomes e remove duplicatas
    return " | ".join(sorted(list(set(valid_diseases))))


def load_orpha_mapping(intervar_dir):
    """
    carrega a tradução de IDs do OMIM para nomes de doenças e herança usando a base do orphanet.
    callback para tentar completar o nome de doenças e herança do arquivo .json final
    param intervar_dir: recebe a localização do intervar no projeto
    return: dicionário de mapeamento com os OMIM IDs, fenótipos e padrão de herança
    """
    # cria o diciońário de mapeamento
    mapping = {}

    # obtém o caminho do arquivo do orphanet (orpha.txt)
    filepath = os.path.join(intervar_dir, "intervardb", "orpha.txt")

    # caso o arquivo exista, processa o mesmo
    if os.path.exists(filepath):
        # ignora erros de encoding caso o arquivo tenha caracteres especiais
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            next(f, None)  # pula a linha de cabeçalho usando o next
            for line in f:
                # remove caracteres de espaço e separa a linha original
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    info = parts[1].split("|")
                    if len(info) >= 6:
                        syndrome = info[1]  # obtém o nome do fenótipo
                        inheritance = (
                            info[3] if info[3] else "Unknown"
                        )  # obtém a herança
                        omims = info[5].split()

                        for omim in omims:
                            # salva síndrome e herança num mini-dicionário para cada OMIM
                            mapping[omim] = {
                                "syndrome": syndrome,
                                "inheritance": inheritance,
                            }
    else:
        print(f"[AVISO] Arquivo de tradução não encontrado em: {filepath}")

    return mapping


def parse_omim_disease(omim_str, orpha_mapping):
    """
    extrai e traduz o ID do fenótipo OMIM e a herança genética usando o Orphanet.
    param omim_str: recebe o ID OMIM a partir do parsing do arquivo intervar
    param orpha_mapping: dicionário contendo o mapeamento OMIM ID: Fenótipo/Herança
    return: o fenótipo correspondente, a herança correspondente limpa
    """

    # verifica se a linha obtida é um valor nulo ou válido
    if pd.isna(omim_str) or str(omim_str).strip() in [".", ""]:
        return "Not provided", "Unknown"

    # obtém os OMIM ids o processamento
    ids = str(omim_str).replace("_", " ").strip(";").split(";")
    names = []
    inheritances = []

    # itera sobre os OMIM ids, caso o ID seja nulo, pula a execução
    # caso o ID esteja no dicionário de mapeamento, obtém a(as) síndrome correspondente
    for mim_id in ids:
        mim_id = mim_id.strip()
        if not mim_id:
            continue

        if mim_id in orpha_mapping:
            names.append(orpha_mapping[mim_id]["syndrome"])

            # tenta capturar a herança da doença, e limpa tags HTML e caracteres de preenchimento indesejados
            inh = (
                str(orpha_mapping[mim_id]["inheritance"])
                .replace("<br>", " ")
                .replace("&nbsp;", " ")
            )
            # verifica se a herança é valida
            if inh.strip() and inh.strip() not in [
                "Unknown",
                "Not applicable",
                "-",
            ]:
                inheritances.append(inh.strip())
        else:
            names.append(f"OMIM ID: {mim_id}")

    # une as doenças correspondentes ao ID e as heranças das doenças
    final_disease = " | ".join(sorted(list(set(names))))
    final_inheritance = (
        " | ".join(sorted(list(set(inheritances))))
        if inheritances
        else "Unknown"
    )

    return final_disease, final_inheritance


def parse_intervar(iv_str):
    """
    extrai a classificação final e a lista de evidências (ex: PVS1, PM2)
    param iv_str: recebe o elemento iterado da coluna de classificação do intervar
    return: a classificação do intervar (ex: Benign) e a(s) evidência(s) (ex: AB1)
    """
    # verifica se a string é válida ou inválida
    if pd.isna(iv_str) or iv_str == ".":
        return "Not classified", ""

    # ex: InterVar: Benign PVS1=0 PS=[0, 0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0, 0] PP=[0, 0, 1, 0, 0, 0] BA1=1 BS=[1, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]
    # remove o "InterVar", e remove espaços
    clean = str(iv_str).replace("InterVar: ", "").strip()
    # separa os elementos
    parts = clean.split()

    classification_words = []
    evidence_list = []

    # obtém as evidências da string processada anteriormente
    hit_evidence = False  # A nossa trava lógica

    for p in parts:
        if "=" in p:
            # obtém a evidência, e avalia o seu valor
            hit_evidence = True
            # caso exista evidência, adiciona ela a lista de evidências
            if "=1" in p:
                evidence_list.append(p.split("=")[0])
        elif not hit_evidence:
            # guarda palavras se ainda não tiver chegado na seção de evidências
            classification_words.append(p)

    # reúne as strings de evidência e classificação
    classification = " ".join(classification_words).capitalize()
    evidence_str = ", ".join(evidence_list)

    return classification, evidence_str


def main(vcf_path, annovar_path, intervar_path, output_json_path):
    print(
        "=============================================================================="
    )
    print("GenoLaudo - Iniciando Parsing das anotações...")
    print(
        "Referências para os campos do parsing: \n "
        "RICHARDS, S. et al. Standards and guidelines for the interpretation of sequence variants: \n "
        "a joint consensus recommendation of the American College of Medical Genetics and Genomics and the \n "
        "Association for Molecular Pathology. Genetics in medicine : official journal of the American College of \n "
        "Medical Genetics, 2015. v. 17, n. 5, p. 405–24. "
    )
    print("==============================================================================")

    # obtenção dinâmica do orpha.txt
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    real_intervar_dir = os.path.join(project_root, "data", "intervar")

    # criação do dicionário com o mapeamento OMIM ID: fenótipo a partir do orpha.txt
    print("[0/4] Carregando dicionário Orphanet (OMIM to Syndrome)...")
    orpha_mapping = load_orpha_mapping(real_intervar_dir)
    print(f"      -> {len(orpha_mapping)} IDs médicos mapeados.")

    # =========================================================================
    # 1. LER O VCF: OBTENÇÃO DO HGVS E GENÓTIPO DA VARIANTE PARA CADA AMOSTRA
    # =========================================================================
    # leitura do vcf com o cyvcf2
    print(f"[1/4] Lendo VCF: {vcf_path}")
    vcf = VCF(vcf_path)
    sample_names = vcf.samples  # obtenção dos nomes das amostras no vcf.norm

    # cria a base do json com a estrutura {sample: {dados}...}
    samples_dict = {sample: {} for sample in sample_names}
    # set para remover duplicatas de coodernadas genômicas
    all_valid_keys = set()

    # iterando sobre o VCF para obter a notação HGVS das variantes
    for variant in vcf:
        chrom = variant.CHROM.replace("chr", "")
        pos = int(variant.POS)
        ref = variant.REF
        alt = variant.ALT[0] if variant.ALT else "."

        # Emulação da matemática do ANNOVAR para geração da chave de junção
        annovar_start = pos
        annovar_ref = ref
        annovar_alt = alt

        # Aplica a quebra apenas para Indels (tamanhos de REF e ALT diferentes)
        if len(annovar_ref) != len(annovar_alt):
            # Verifica se a primeira base é idêntica (âncora do VCF)
            if annovar_ref and annovar_alt and annovar_ref[0] == annovar_alt[0]:
                annovar_ref = annovar_ref[1:]  # Remove a âncora de REF
                annovar_alt = annovar_alt[1:]  # Remove a âncora de ALT

                if len(annovar_ref) == 0:
                    annovar_ref = "-"
                    # Inserções puras: o ANNOVAR mantém o Start igual à posição da âncora
                elif len(annovar_alt) == 0:
                    annovar_alt = "-"
                    # Deleções puras: o ANNOVAR soma 1 (a variante começa na base seguinte)
                    annovar_start += 1
                else:
                    # Substituições em bloco (delins): soma 1 também
                    annovar_start += 1

        # A chave de cruzamento usa estritamente os dados normalizados sem a âncora
        key = f"{chrom}_{annovar_start}_{annovar_ref}_{annovar_alt}"

        # itera sobre as amostras, se o genótipo for missing, a variante não é registrada
        for idx, sample_name in enumerate(sample_names):
            gt = variant.genotypes[idx]
            a, b = gt[0], gt[1]  # genótipo diploide obtido

            # pula genótipos missing ou 0/0 (referência)
            if a == -1 or b == -1 or (a == 0 and b == 0):
                continue

            # define a zigosidade da amostra
            zygosity = "Homozygous" if a == b else "Heterozygous"

            # a amostra é inicializada no dicionário, mantendo a coordenada visual do VCF original
            samples_dict[sample_name][key] = {
                "Coordenada Genômica": f"chr{chrom}:{pos}:{ref}:{alt}",
                "Gene": "Unknown",
                "Transcrito": "Unknown",
                "Localização": "Unknown",
                "Variante": "Unknown",
                "Zigosidade": zygosity,
                "Classificação": "Not classified",
                "Evidência ACMG": "",
                "Doença": "Not provided",
                "Herança": "Unknown",
                "Frequência ABraOM": 0,
                "Origem Parental": "-",
            }

            # armazena a chave modificada para o cruzamento nos Passos 2 e 3
            all_valid_keys.add(key)

    print(f"      -> Total de amostras: {len(sample_names)}")
    print(f"      -> {len(all_valid_keys)} variantes válidas retidas.")

    gc.collect()  # limpeza de resíduos na memória

    # =========================================================================
    # 2. LER O ANNOVAR: OBTENÇÃO DA DOENÇA, FREQUÊNCIA, REVEL SCORE E CADD
    # =========================================================================
    # leitura do annovar, e processamento do header para obter as colunas
    print(f"[2/4] Lendo anotações ANNOVAR: {annovar_path}")
    with open(annovar_path, "r") as f:
        header = f.readline().strip().split("\t")

    # identificando as colunas presentes no ANNOVAR
    # next(): itera sobre o elemento até obter a coluna correspondente, caso o elemento não exista, retorna None
    clndn_col = next((c for c in header if c == "CLNDN"), None)
    abraom_col = next((c for c in header if c == "abraom_freq"), None)
    clnsig_col = next((c for c in header if c == "CLNSIG"), None)  # Nova detecção da coluna do ClinVar

    # colunas padrão do ANNOVAR
    # adiciona a cols_anno as colunas obtidas anteriormente caso elas existam
    cols_anno = [
        "Chr",
        "Start",
        "Ref",
        "Alt",
        "Gene.refGene",
        "Func.refGene",
        "ExonicFunc.refGene",
        "AAChange.refGene",
    ]
    if clndn_col:
        cols_anno.append(clndn_col)
    if abraom_col:
        cols_anno.append(abraom_col)
    if clnsig_col:
        cols_anno.append(clnsig_col)  # Adicionada para extração

    # itera sobre as linhas do arquivo do ANNOVAR, retornando chunks de 100000 linhas, para evitar estouro da memória
    # retorna apenas os valores de colunas em cols_anno
    chunk_iterator = pd.read_csv(
        annovar_path, sep="\t", usecols=cols_anno, dtype=str, chunksize=100000
    )

    # itera sobre a lista de chunks geradas e obtém o HGVS das variantes
    # verifica se ela está no set de variantes com genótipo válido
    for chunk in chunk_iterator:
        chunk["Key"] = (
                chunk["Chr"].str.replace("chr", "")
                + "_"
                + chunk["Start"]
                + "_"
                + chunk["Ref"]
                + "_"
                + chunk["Alt"]
        )
        chunk = chunk[chunk["Key"].isin(all_valid_keys)]

        # criação de arrays com os valores obtidos das colunas de cada chunk
        keys_arr = chunk["Key"].values
        genes_arr = chunk["Gene.refGene"].fillna("Unknown").values
        funcs_arr = chunk["Func.refGene"].fillna("").values
        exonics_arr = chunk["ExonicFunc.refGene"].fillna("").values
        aachanges_arr = chunk["AAChange.refGene"].fillna("").values

        # verifica se os valores das colunas armazenados nos arrays são válidos
        # substitui valores nulos por "."
        clndn_arr = (
            chunk[clndn_col].fillna("").values
            if clndn_col
            else [""] * len(keys_arr)
        )
        abraom_arr = (
            chunk[abraom_col].fillna(".").values
            if abraom_col
            else ["."] * len(keys_arr)
        )
        clnsig_arr = (
            chunk[clnsig_col].fillna("").values
            if clnsig_col
            else [""] * len(keys_arr)
        )  # Array com os dados do ClinVar

        # iteração sobre os arrays do ANNOVAR, utilizando o zip para iteração paralela
        for (
                k,
                gene,
                func,
                exonic,
                aachange,
                clndn_val,
                abraom_val,
                clnsig_val,  # Variável iteradora do ClinVar
        ) in zip(
            keys_arr,
            genes_arr,
            funcs_arr,
            exonics_arr,
            aachanges_arr,
            clndn_arr,
            abraom_arr,
            clnsig_arr,  # Inserida na iteração
        ):
            # obtenção de valores e saneamento
            location = (
                f"{func} ({exonic})" if exonic and exonic != "." else func
            )  # obtenção do location
            transcript, variant_hgvs = parse_aachange(
                aachange
            )  # obtenção do valor transcript
            disease = (
                parse_disease(clndn_val) if clndn_col else "Not provided"
            )  # obtenção da doença

            # extração e formatação da classificação do ClinVar
            clinvar_sig = ""
            if clnsig_col and str(clnsig_val).strip() not in [".", ""]:
                # extrai a primeira classificação (caso haja múltiplas) e formata adequadamente
                raw_sig = str(clnsig_val).split(",")[0].replace("_", " ").capitalize()
                # filtra classificações inúteis ou não-conclusivas
                if raw_sig.lower() not in ["not provided", "conflicting interpretations of pathogenicity", ""]:
                    clinvar_sig = raw_sig

            # formatação da frequência do ABraOM como numérico (float/int) para garantir a filtragem
            try:
                abraom_float = (
                    float(abraom_val) if abraom_val not in [".", ""] else 0
                )
            except ValueError:
                abraom_float = 0

            # atualizando o dicionário com os dados da variante
            # utiliza a HGVS para identificar a variante correta
            for sample_name in sample_names:
                if k in samples_dict[sample_name]:
                    samples_dict[sample_name][k]["Gene"] = gene
                    samples_dict[sample_name][k]["Localização"] = location
                    samples_dict[sample_name][k]["Transcrito"] = transcript
                    samples_dict[sample_name][k]["Variante"] = variant_hgvs
                    samples_dict[sample_name][k]["Doença"] = disease
                    samples_dict[sample_name][k]["Frequência ABraOM"] = abraom_float

                    # Prioridade Clínica: se houver classificação válida do ClinVar, ela é registrada imediatamente
                    if clinvar_sig:
                        samples_dict[sample_name][k]["Classificação"] = clinvar_sig

        # remoção dos arrays armazenados
        del chunk, keys_arr, genes_arr, funcs_arr, exonics_arr, aachanges_arr
        del clndn_arr, abraom_arr, clnsig_arr  # Memória do array do ClinVar limpa
        gc.collect()  # limpeza de resíduos na memória

    # ===========================================================================
    # 3. LER INTERVAR: OBTER CLASSIFICAÇÃO INTERVAR E FENÓTIPO OMIM
    # ===========================================================================
    print(
        f"[3/4] Lendo classificações InterVar e cruzando Phenotype_MIM: {intervar_path}"
    )
    # obtendo header do intervar
    with open(intervar_path, "r") as f:
        iv_header = f.readline().strip("\n").split("\t")

    # obtenção das colunas do intervar
    chr_col = iv_header[0]
    intervar_col = next(
        (c for c in iv_header if "InterVar" in c and "Evidence" in c), None
    )
    phenotype_col = next((c for c in iv_header if "Phenotype_MIM" in c), None)

    if not intervar_col:
        print("[ERRO FATAL] Coluna do InterVar não encontrada no cabeçalho!")
        sys.exit(1)

    # definição das colunas de interesse do intervar
    # adiciona as colunas obtidas anteriormente caso elas existam
    iv_cols = [chr_col, "Start", "Ref", "Alt", intervar_col]
    if phenotype_col:
        iv_cols.append(phenotype_col)

    # itera sobre as linhas do arquivo do intervar, retornando chunks de 100000 linhas, para evitar estouro da memória
    # retorna apenas os valores de colunas em iv_cols
    chunk_iterator_iv = pd.read_csv(
        intervar_path, sep="\t", usecols=iv_cols, dtype=str, chunksize=100000
    )

    # itera sobre a lista de chunks geradas e obtém o HGVS da variante nos dados do intervar
    # verifica se ela está no set de variantes com genótipo válido para as amostras
    for iv_chunk in chunk_iterator_iv:
        iv_chunk["Key"] = (
                iv_chunk[chr_col].str.replace("chr", "")
                + "_"
                + iv_chunk["Start"]
                + "_"
                + iv_chunk["Ref"]
                + "_"
                + iv_chunk["Alt"]
        )
        iv_chunk = iv_chunk[iv_chunk["Key"].isin(all_valid_keys)]

        # cria os arrays com os valores das colunas obtidas durante a iteração
        # substitui valores nulos por "."
        keys_arr_iv = iv_chunk["Key"].values
        classes_arr = iv_chunk[intervar_col].fillna("").values
        pheno_arr = (
            iv_chunk[phenotype_col].fillna("").values
            if phenotype_col
            else [""] * len(keys_arr_iv)
        )

        # iteração sobre os arrays do intervar, utilizando o zip para iteração paralela
        for k, raw_class, pheno_val in zip(
                keys_arr_iv, classes_arr, pheno_arr
        ):
            classification, acmg_evidence = parse_intervar(
                raw_class
            )  # obtenção da classificação e evidência ACMG

            # repassa o mapping do orphanet para tradução visual e extração de herança limpa
            omim_disease, omim_inheritance = parse_omim_disease(
                pheno_val, orpha_mapping
            )

            # atualiza o dicionário com os valores obtidos da iteração
            # utiliza a HGVS para identificar a variante correta
            for sample_name in sample_names:
                if k in samples_dict[sample_name]:

                    # Override Secundário: Utiliza a classificação do InterVar APENAS se o ClinVar for omisso (Not classified) no Passo 2
                    if samples_dict[sample_name][k]["Classificação"] == "Not classified":
                        samples_dict[sample_name][k]["Classificação"] = classification

                    # As evidências calculadas pelo InterVar continuam sendo mantidas para suplementar o laudo
                    samples_dict[sample_name][k]["Evidência ACMG"] = (
                        acmg_evidence
                    )

                    # preenche a herança do orphanet se estiver disponível
                    if omim_inheritance != "Unknown":
                        samples_dict[sample_name][k]["Herança"] = (
                            omim_inheritance
                        )

                    # fallback do OMIM, para tentar assimilar um valor
                    current_disease = samples_dict[sample_name][k].get(
                        "Doença", "Not provided"
                    )
                    if (
                            current_disease == "Not provided"
                            or current_disease == "Unknown"
                    ):
                        if omim_disease != "Not provided":
                            samples_dict[sample_name][k]["Doença"] = (
                                omim_disease
                            )
        # remoção dos arrays da memória
        del iv_chunk, keys_arr_iv, classes_arr, pheno_arr
        gc.collect()  # limpeza de resíduos na memória

    # =========================================================================
    # PASSO 4: EXPORTAR JSON
    # =========================================================================
    print("[4/4] Montando JSON final e gravando no disco...")
    # criação da lista que conterá os dicionários das amostras
    final_output = []
    # itera sobre as amostras e suas variantes válidas
    for sample_name, variants_dict in samples_dict.items():
        # cria uma entrada válida para a amostra no json
        final_output.append(
            {f"{sample_name}_vcf": list(variants_dict.values())}
        )
    # cria o arquivo json
    os.makedirs(os.path.dirname(output_json_path), exist_ok=True)
    with open(output_json_path, "w", encoding="utf-8") as f:
        json.dump(final_output, f, indent=4, ensure_ascii=False)

    print("==============================================================================")
    print(f"GenoLaudo - CONCLUÍDO! O laudo formatado está em: {output_json_path}")
    print("==============================================================================")


if __name__ == "__main__":
    # chamada da função main, que inicia o parsing
    main(
        snakemake.input.norm_vcf,
        snakemake.input.annovar_file,
        snakemake.input.intervar_file,
        snakemake.output.output_json,
    )