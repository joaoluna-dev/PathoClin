import sys
import os
import json
import gc
import pandas as pd
from cyvcf2 import VCF


def parse_aachange(aachange):
    if pd.isna(aachange) or aachange == '.' or aachange == '':
        return "Unknown"
    parts = str(aachange).split(',')[0].split(':')
    if len(parts) >= 5:
        return f"{parts[1]}:{parts[3]} ({parts[4]})" if len(parts) > 4 else f"{parts[1]}:{parts[3]}"
    return str(aachange)


def parse_disease(clinvar_str):
    if pd.isna(clinvar_str) or clinvar_str == '.':
        return "Not provided"
    parts = str(clinvar_str).split('|_|')
    if len(parts) >= 3:
        return parts[2].replace('_', ' ')
    return str(clinvar_str).replace('_', ' ')


def parse_intervar(iv_str):
    if pd.isna(iv_str) or iv_str == '.':
        return "Not classified"
    clean = str(iv_str).replace("InterVar: ", "").strip()
    parts = clean.split()
    classification = []
    for p in parts:
        if '=' in p or '[' in p:
            break
        classification.append(p)
    return " ".join(classification).capitalize()


def main(vcf_path, annovar_path, intervar_path, output_json_path):
    print("==============================================================================")
    print("GenoLaudo - Iniciando Parsing das variantes anotadas pelo ANNOVAR e InterVar")
    print("==============================================================================")

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

            #pula Missing ou Homo Reference
            if a == -1 or b == -1 or (a == 0 and b == 0):
                continue

            zygosity = "Homozygous" if a == b else "Heterozygous"

            samples_dict[sample_name][key] = {
                "Variant": f"chr{chrom}:{pos}:{ref}:{alt}",
                "Gene": "Unknown",
                "Transcript_Variant": "Unknown",
                "Location": "Unknown",
                "Zygosity": zygosity,
                "Classification": "Not classified",
                "Disease": "Not provided",
                "Inheritance": "Unknown",
                "Parental_Origin": "Unknown"
            }
            all_valid_keys.add(key)

    print(f"      -> Total de amostras: {len(sample_names)}")
    print(f"      -> {len(all_valid_keys)} variantes válidas retidas.")

    gc.collect()  #limpeza de resíduos na memória RAM

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

        #extração em massa via vetores
        keys_arr = chunk['Key'].values
        genes_arr = chunk['Gene.refGene'].fillna('Unknown').values
        funcs_arr = chunk['Func.refGene'].fillna('').values
        exonics_arr = chunk['ExonicFunc.refGene'].fillna('').values
        aachanges_arr = chunk['AAChange.refGene'].fillna('').values
        clinvars_arr = chunk[clinvar_col].fillna('').values if clinvar_col else [''] * len(keys_arr)

        #iteração sobre as linhas do arquivo ANNOVAR, utilizando o zip para iteração paralela
        for k, gene, func, exonic, aachange, clin_val in zip(keys_arr, genes_arr, funcs_arr, exonics_arr, aachanges_arr,
                                                             clinvars_arr):
            location = f"{func} ({exonic})" if exonic and exonic != '.' else func
            transcript = parse_aachange(aachange)
            disease = parse_disease(clin_val) if clinvar_col else "Not provided"

            for sample_name in sample_names:
                if k in samples_dict[sample_name]:
                    samples_dict[sample_name][k]["Gene"] = gene
                    samples_dict[sample_name][k]["Location"] = location
                    samples_dict[sample_name][k]["Transcript_Variant"] = transcript
                    samples_dict[sample_name][k]["Disease"] = disease

        del chunk, keys_arr, genes_arr, funcs_arr, exonics_arr, aachanges_arr, clinvars_arr
        gc.collect()  #limpeza de resíduos na memória RAM

    # =========================================================================
    # PASSO 3: LER INTERVAR
    # =========================================================================
    print(f"[3/4] Lendo classificações InterVar: {intervar_path}")
    with open(intervar_path, 'r') as f:
        iv_header = f.readline().strip('\n').split('\t')

    chr_col = iv_header[0]
    intervar_col = next((c for c in iv_header if 'InterVar' in c and 'Evidence' in c), None)

    if not intervar_col:
        print("[ERRO FATAL] Coluna do InterVar não encontrada no cabeçalho!")
        sys.exit(1)

    iv_cols = [chr_col, 'Start', 'Ref', 'Alt', intervar_col]
    chunk_iterator_iv = pd.read_csv(intervar_path, sep='\t', usecols=iv_cols, dtype=str, chunksize=100000)

    for iv_chunk in chunk_iterator_iv:
        iv_chunk['Key'] = iv_chunk[chr_col].str.replace('chr', '') + "_" + iv_chunk['Start'] + "_" + iv_chunk[
            'Ref'] + "_" + iv_chunk['Alt']
        iv_chunk = iv_chunk[iv_chunk['Key'].isin(all_valid_keys)]

        keys_arr_iv = iv_chunk['Key'].values
        classes_arr = iv_chunk[intervar_col].fillna('').values

        for k, raw_class in zip(keys_arr_iv, classes_arr):
            classification = parse_intervar(raw_class)

            for sample_name in sample_names:
                if k in samples_dict[sample_name]:
                    samples_dict[sample_name][k]["Classification"] = classification

        del iv_chunk, keys_arr_iv, classes_arr
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