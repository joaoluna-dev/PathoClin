import sys
import os
import json
import gc
import pandas as pd
from cyvcf2 import VCF


def parse_aachange(aachange):
    """
    Divide o AAChange em Transcript e Variant (c. / p.) segundo a ACMG.
    Param:aachange: recebe a notação HGVS do intervar
    return: string formatada
    """
    if pd.isna(aachange) or aachange == '.' or aachange == '':
        return "Unknown", "Unknown"

    parts = str(aachange).split(',')[0].split(':')
    if len(parts) >= 5:
        transcript = parts[1]
        # Monta a string "c.XXX (p.XXX)" ou apenas "c.XXX"
        variant_hgvs = f"{parts[3]} ({parts[4]})" if len(parts) > 4 else f"{parts[3]}"
        return transcript, variant_hgvs

    return "Unknown", str(aachange)


def parse_disease(clinvar_str):
    if pd.isna(clinvar_str) or clinvar_str == '.':
        return "Not provided"
    parts = str(clinvar_str).split('|_|')
    if len(parts) >= 3:
        return parts[2].replace('_', ' ')
    return str(clinvar_str).replace('_', ' ')


def load_orpha_mapping(intervar_dir):
    """Carrega a tradução de IDs do OMIM para Nomes de Doenças usando a base do Orphanet."""
    mapping = {}
    filepath = os.path.join(intervar_dir, "intervardb", "orpha.txt")

    if os.path.exists(filepath):
        # Ignora erros de encoding caso o arquivo tenha caracteres especiais
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            next(f, None)  # Pula a linha de cabeçalho
            for line in f:
                parts = line.strip().split('\t')
                # Exemplo da linha: OrphaID \t OrphaID|SyndromeName|Prev|Inheritance|Age|OMIM
                if len(parts) >= 2:
                    info = parts[1].split('|')
                    if len(info) >= 6:
                        syndrome = info[1]
                        omims = info[5].split()  # Pode haver mais de um OMIM por doença
                        for omim in omims:
                            mapping[omim] = syndrome
    else:
        print(f"[AVISO] Arquivo de tradução não encontrado em: {filepath}")

    return mapping


def parse_omim_disease(omim_str, orpha_mapping):
    """Extrai e traduz a string de fenótipo usando o Orphanet."""
    if pd.isna(omim_str) or str(omim_str).strip() in ['.', '']:
        return "Not provided"

    ids = str(omim_str).replace('_', ' ').strip(';').split(';')
    names = []

    for mim_id in ids:
        mim_id = mim_id.strip()
        if not mim_id:
            continue

        if mim_id in orpha_mapping:
            names.append(orpha_mapping[mim_id])
        else:
            names.append(f"OMIM ID: {mim_id}")

    # Usa set para evitar duplicatas, caso 2 códigos apontem para a mesma síndrome
    return " | ".join(sorted(list(set(names))))


def parse_intervar(iv_str):
    """Extrai a Classificação Final e a lista de Evidências (ex: PVS1, PM2)."""
    if pd.isna(iv_str) or iv_str == '.':
        return "Not classified", ""

    clean = str(iv_str).replace("InterVar: ", "").strip()
    parts = clean.split()

    classification_words = []
    evidence_list = []

    hit_evidence = False  # A nossa trava lógica

    for p in parts:
        if '=' in p:
            hit_evidence = True  # A trava desce ao achar a primeira regra
            if '=1' in p:
                evidence_list.append(p.split('=')[0])
        elif not hit_evidence:
            # Só guarda palavras se ainda não tivermos chegado na seção de evidências
            classification_words.append(p)

    classification = " ".join(classification_words).capitalize()
    evidence_str = ", ".join(evidence_list)

    return classification, evidence_str


def main(vcf_path, annovar_path, intervar_path, output_json_path):
    print("==============================================================================")
    print("GenoLaudo - Iniciando Parsing das anotações...")
    print("==============================================================================")

    # Resolve o caminho do projeto dinamicamente para achar o orpha.txt
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    real_intervar_dir = os.path.join(project_root, "data", "intervar")

    print("[0/4] Carregando dicionário Orphanet (OMIM to Syndrome)...")
    orpha_mapping = load_orpha_mapping(real_intervar_dir)
    print(f"      -> {len(orpha_mapping)} IDs médicos mapeados.")

    # =========================================================================
    # PASSO 1: LER O VCF E CRIAR O DICIONÁRIO EM MEMÓRIA
    # =========================================================================
    print(f"[1/4] Lendo VCF: {vcf_path}")
    vcf = VCF(vcf_path)
    sample_names = vcf.samples

    samples_dict = {sample: {} for sample in sample_names}
    all_valid_keys = set()

    for variant in vcf:
        chrom = variant.CHROM.replace('chr', '')
        pos = str(variant.POS)
        ref = variant.REF
        alt = variant.ALT[0]
        key = f"{chrom}_{pos}_{ref}_{alt}"

        for idx, sample_name in enumerate(sample_names):
            gt = variant.genotypes[idx]
            a, b = gt[0], gt[1]

            # pula Missing ou Homo Reference
            if a == -1 or b == -1 or (a == 0 and b == 0):
                continue

            zygosity = "Homozygous" if a == b else "Heterozygous"

            samples_dict[sample_name][key] = {
                "Genomic_Coordinate": f"chr{chrom}:{pos}:{ref}:{alt}",
                "Gene": "Unknown",
                "Transcript": "Unknown",
                "Location": "Unknown",
                "Variant": "Unknown",
                "Zygosity": zygosity,
                "Classification": "Not classified",
                "ACMG_Evidence": "",
                "Disease": "Not provided",
                "Inheritance": "Unknown",
                "Parental_Origin": "Unknown"
            }
            all_valid_keys.add(key)

    print(f"      -> Total de amostras: {len(sample_names)}")
    print(f"      -> {len(all_valid_keys)} variantes válidas retidas.")

    gc.collect()  # limpeza de resíduos na memória

    # =========================================================================
    # PASSO 2: LER ANNOVAR
    # =========================================================================
    print(f"[2/4] Lendo anotações ANNOVAR: {annovar_path}")
    with open(annovar_path, 'r') as f:
        header = f.readline().strip().split('\t')

    clinvar_col = next((c for c in header if 'clinvar' in c.lower()), None)
    cols_anno = ['Chr', 'Start', 'Ref', 'Alt', 'Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 'AAChange.refGene']
    if clinvar_col: cols_anno.append(clinvar_col)

    chunk_iterator = pd.read_csv(annovar_path, sep='\t', usecols=cols_anno, dtype=str, chunksize=100000)

    for chunk in chunk_iterator:
        chunk['Key'] = chunk['Chr'].str.replace('chr', '') + "_" + chunk['Start'] + "_" + chunk['Ref'] + "_" + chunk[
            'Alt']
        chunk = chunk[chunk['Key'].isin(all_valid_keys)]

        # criação de arrays com os valores obtidos das colunas de cada chunk
        keys_arr = chunk['Key'].values
        genes_arr = chunk['Gene.refGene'].fillna('Unknown').values
        funcs_arr = chunk['Func.refGene'].fillna('').values
        exonics_arr = chunk['ExonicFunc.refGene'].fillna('').values
        aachanges_arr = chunk['AAChange.refGene'].fillna('').values
        clinvars_arr = chunk[clinvar_col].fillna('').values if clinvar_col else [''] * len(keys_arr)

        # iteração sobre as linhas do arquivo ANNOVAR, utilizando o zip para iteração paralela de cada array
        for k, gene, func, exonic, aachange, clin_val in zip(keys_arr, genes_arr, funcs_arr, exonics_arr, aachanges_arr,
                                                             clinvars_arr):
            location = f"{func} ({exonic})" if exonic and exonic != '.' else func
            transcript, variant_hgvs = parse_aachange(aachange)
            disease = parse_disease(clin_val) if clinvar_col else "Not provided"

            # atualizando o dicionário com os dados da variante
            for sample_name in sample_names:
                if k in samples_dict[sample_name]:
                    samples_dict[sample_name][k]["Gene"] = gene
                    samples_dict[sample_name][k]["Location"] = location
                    samples_dict[sample_name][k]["Transcript"] = transcript
                    samples_dict[sample_name][k]["Variant"] = variant_hgvs
                    samples_dict[sample_name][k]["Disease"] = disease

        del chunk, keys_arr, genes_arr, funcs_arr, exonics_arr, aachanges_arr, clinvars_arr
        gc.collect()  # limpeza de resíduos na memória

    # ===========================================================================
    # PASSO 3: LER INTERVAR E APLICAR NOME DA DOENÇA USANDO ORPHANET DO INTERVAR
    # ===========================================================================
    print(f"[3/4] Lendo classificações InterVar e cruzando Phenotype_MIM: {intervar_path}")
    with open(intervar_path, 'r') as f:
        iv_header = f.readline().strip('\n').split('\t')

    chr_col = iv_header[0]
    intervar_col = next((c for c in iv_header if 'InterVar' in c and 'Evidence' in c), None)
    phenotype_col = next((c for c in iv_header if 'Phenotype_MIM' in c), None)

    if not intervar_col:
        print("[ERRO FATAL] Coluna do InterVar não encontrada no cabeçalho!")
        sys.exit(1)

    iv_cols = [chr_col, 'Start', 'Ref', 'Alt', intervar_col]
    if phenotype_col:
        iv_cols.append(phenotype_col)

    chunk_iterator_iv = pd.read_csv(intervar_path, sep='\t', usecols=iv_cols, dtype=str, chunksize=100000)

    for iv_chunk in chunk_iterator_iv:
        iv_chunk['Key'] = iv_chunk[chr_col].str.replace('chr', '') + "_" + iv_chunk['Start'] + "_" + iv_chunk[
            'Ref'] + "_" + iv_chunk['Alt']
        iv_chunk = iv_chunk[iv_chunk['Key'].isin(all_valid_keys)]

        keys_arr_iv = iv_chunk['Key'].values
        classes_arr = iv_chunk[intervar_col].fillna('').values
        pheno_arr = iv_chunk[phenotype_col].fillna('').values if phenotype_col else [''] * len(keys_arr_iv)

        for k, raw_class, pheno_val in zip(keys_arr_iv, classes_arr, pheno_arr):
            classification, acmg_evidence = parse_intervar(raw_class)

            # Repassa o mapping do orphanet para tradução visual
            omim_disease = parse_omim_disease(pheno_val, orpha_mapping)

            for sample_name in sample_names:
                if k in samples_dict[sample_name]:
                    samples_dict[sample_name][k]["Classification"] = classification
                    samples_dict[sample_name][k]["ACMG_Evidence"] = acmg_evidence

                    # Sistema de resgate (Fallback) do OMIM
                    current_disease = samples_dict[sample_name][k].get("Disease", "Not provided")
                    if current_disease == "Not provided" or current_disease == "Unknown":
                        if omim_disease != "Not provided":
                            samples_dict[sample_name][k]["Disease"] = omim_disease

        del iv_chunk, keys_arr_iv, classes_arr, pheno_arr
        gc.collect()

    # =========================================================================
    # PASSO 4: EXPORTAR JSON
    # =========================================================================
    print(f"[4/4] Montando JSON final e gravando no disco...")

    final_output = []
    for sample_name, variants_dict in samples_dict.items():
        final_output.append({
            f"{sample_name}_vcf": list(variants_dict.values())
        })

    os.makedirs(os.path.dirname(output_json_path), exist_ok=True)
    with open(output_json_path, 'w', encoding='utf-8') as f:
        json.dump(final_output, f, indent=4, ensure_ascii=False)

    print(f"==============================================================================")
    print(f"GenoLaudo - CONCLUÍDO! O laudo formatado está em: {output_json_path}")
    print(f"==============================================================================")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Uso: python parser.py <input.norm.vcf> <input.multianno.txt> <input.intervar> <output.json>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])