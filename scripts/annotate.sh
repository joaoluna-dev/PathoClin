#!/bin/bash
set -e

# ==============================================================================
# CONFIGURAÇÃO DE CAMINHOS
# ==============================================================================
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANNOVAR_DIR="$PROJECT_ROOT/data/annovar"
INTERVAR_DIR="$PROJECT_ROOT/data/intervar"
TABLE_ANNOVAR="$ANNOVAR_DIR/table_annovar.pl"
CONVERT_ANNOVAR="$ANNOVAR_DIR/convert2annovar.pl"
DB_DIR="$ANNOVAR_DIR/humandb/dbs"
REF_GENOME="$PROJECT_ROOT/data/genome/hg38.fa"

if [ "$#" -ne 2 ]; then
    echo "Uso: $0 <input_vcf> <output_basename>"
    exit 1
fi

INPUT_VCF=$1
OUTPUT_BASE=$2

if [[ -e "$(dirname "$OUTPUT_BASE")" ]]; then
  echo "O diretório $OUTPUT_BASE já existe. Os arquivos serão sobrescritos."
  :
fi

mkdir -p "$(dirname "$OUTPUT_BASE")"

echo "=============================================================================="
echo "GenoLaudo - Iniciando Anotação das variantes filtradas em $INPUT_VCF..."
echo "=============================================================================="

# ==============================================================================
# 1. VERIFICAÇÃO DOS BANCOS OBRIGATÓRIOS
# ==============================================================================
PROTOCOLS=""
OPERATIONS=""

get_db_filename() {
    ls "$DB_DIR"/hg38_*"$1"*.txt 2>/dev/null | head -n 1
}

# RefGene
if [ -z "$(get_db_filename "refGene")" ]; then echo "[ERRO] refGene não encontrado."; exit 1; fi
PROTOCOLS="refGene"
OPERATIONS="g"

# ClinVar
CLINVAR_FILE=$(get_db_filename "clinvar")
if [ -z "$CLINVAR_FILE" ]; then echo "[ERRO] ClinVar não encontrado."; exit 1; fi
CLINVAR_NAME=$(basename "$CLINVAR_FILE" | sed 's/hg38_//;s/.txt//')
PROTOCOLS="$PROTOCOLS,$CLINVAR_NAME"
OPERATIONS="$OPERATIONS,f"

# ABraOM
ABRAOM_FILE=$(get_db_filename "abraom")
if [ -z "$ABRAOM_FILE" ]; then echo "[ERRO] ABraOM não encontrado."; exit 1; fi
ABRAOM_NAME=$(basename "$ABRAOM_FILE" | sed 's/hg38_//;s/.txt//')
PROTOCOLS="$PROTOCOLS,$ABRAOM_NAME"
OPERATIONS="$OPERATIONS,f"

# dbNSFP
DBNSFP_FILE=$(get_db_filename "dbnsfp")
if [ -z "$DBNSFP_FILE" ]; then echo "[ERRO] dbNSFP não encontrado."; exit 1; fi
DBNSFP_NAME=$(basename "$DBNSFP_FILE" | sed 's/hg38_//;s/.txt//')
PROTOCOLS="$PROTOCOLS,$DBNSFP_NAME"
OPERATIONS="$OPERATIONS,f"

# ==============================================================================
# 2. NORMALIZAÇÃO E CONVERSÃO PARA AVINPUT
# ==============================================================================
echo "GenoLaudo - [1/4] Normalizando e convertendo para formato nativo ANNOVAR..."
NORM_VCF="${OUTPUT_BASE}.norm.vcf"
AV_INPUT="${OUTPUT_BASE}.avinput"

bcftools norm -m -any -f "$REF_GENOME" "$INPUT_VCF" > "$NORM_VCF"

# Converte para o formato seguro (impede o erro de reconstrução do VCF)
perl "$CONVERT_ANNOVAR" -format vcf4 "$NORM_VCF" > "$AV_INPUT" 2>/dev/null

# ==============================================================================
# 3. ANOTAÇÃO DAS VARIANTES
# ==============================================================================
echo "GenoLaudo - [2/4] Executando ANNOVAR..."

# Nota: Sem a flag -vcfinput. O output será um arquivo tabular .txt perfeito.
perl "$TABLE_ANNOVAR" "$AV_INPUT" "$DB_DIR" \
    -buildver hg38 \
    -out "$OUTPUT_BASE" \
    -remove \
    -protocol "$PROTOCOLS" \
    -operation "$OPERATIONS" \
    -nastring .

RAW_TXT="${OUTPUT_BASE}.hg38_multianno.txt"

if [ ! -f "$RAW_TXT" ]; then
    echo "  [ERRO] Falha na anotação do ANNOVAR."
    exit 1
fi

# ==============================================================================
# 4. FILTRAGEM DE BENIGNAS NO ARQUIVO TXT
# ==============================================================================
echo "GenoLaudo - [3/4] Removendo variantes Benignas do relatório..."

CLEAN_TXT="${OUTPUT_BASE}.clean_multianno.txt"

# Lê o arquivo tabular gerado e remove linhas que contêm "Benign" ou "Likely_benign"
awk -F'\t' '{
    if (NR==1) { print $0; next } # Mantém o cabeçalho
    if ($0 ~ /\tBenign/ || $0 ~ /\tLikely_benign/ || $0 ~ /Benign;/ || $0 ~ /Likely_benign;/) { next }
    print $0
}' "$RAW_TXT" > "$CLEAN_TXT"

# Substitui o original pelo limpo para que o InterVar leia apenas as suspeitas
mv "$CLEAN_TXT" "$RAW_TXT"

# ==============================================================================
# 5. CLASSIFICAÇÃO DAS VARIANTES USANDO O INTERVAR
# ==============================================================================
echo "GenoLaudo - [4/4] Executando InterVar (ACMG)..."

# 1. Cria links simbólicos temporários para satisfazer a checagem rígida do InterVar
ln -sf "$ANNOVAR_DIR/annotate_variation.pl" ./annotate_variation.pl
ln -sf "$ANNOVAR_DIR/convert2annovar.pl" ./convert2annovar.pl

# 2. Localiza o script do InterVar
INTERVAR_EXEC=""
if [ -f "$INTERVAR_DIR/InterVar.py" ]; then INTERVAR_EXEC="$INTERVAR_DIR/InterVar.py"; fi
if [ -f "$INTERVAR_DIR/Intervar.py" ]; then INTERVAR_EXEC="$INTERVAR_DIR/Intervar.py"; fi

if [ -n "$INTERVAR_EXEC" ]; then
    # Adicionamos a flag --skip_annovar
    # Isso força o InterVar a ler diretamente o arquivo .txt que limpamos no passo 4
    python3 "$INTERVAR_EXEC" \
        -b hg38 \
        -i "$AV_INPUT" \
        -d "$DB_DIR" \
        --table_annovar "$TABLE_ANNOVAR" \
        -o "$OUTPUT_BASE" \
        -t "$INTERVAR_DIR/intervardb" \
        --input_type AVinput \
        --skip_annovar
else
    echo "  [AVISO] Intervar.py não encontrado."
fi

# 3. Limpeza dos links temporários
echo "  [AVISO] Removendo arquivos temporários..."
rm -f ./annotate_variation.pl ./convert2annovar.pl
rm "$RAW_TXT.grl_p" "$RAW_TXT" "$AV_INPUT"

echo "=============================================================================="
echo "Anotação concluída."