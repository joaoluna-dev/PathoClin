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
GENOME_DIR="$PROJECT_ROOT/data/genome"

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

# RefGene (BUSCA DINÂMICA)
REFGENE_FILE=$(get_db_filename "refGene")
if [ -z "$REFGENE_FILE" ]; then echo "[ERRO] refGene não encontrado."; exit 1; fi
REFGENE_NAME=$(basename "$REFGENE_FILE" | sed 's/hg38_//;s/.txt//')
PROTOCOLS="$REFGENE_NAME"
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
echo "GenoLaudo - [1/3] Preparando genoma, normalizando e convertendo para ANNOVAR..."

REF_GENOME=$(ls "$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fasta 2>/dev/null | head -n 1)

if [ -z "$REF_GENOME" ]; then
    echo "[ERRO] Nenhum genoma de referência (.fa ou .fasta) encontrado em $GENOME_DIR."
    exit 1
fi

echo "  -> Genoma de referência localizado: $(basename "$REF_GENOME")"

if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "  -> Índice .fai não encontrado. Gerando índice dinamicamente com samtools..."
    samtools faidx "$REF_GENOME"
    echo "  -> Índice criado com sucesso!"
fi

NORM_VCF="${OUTPUT_BASE}.norm.vcf"
AV_INPUT="${OUTPUT_BASE}.avinput"

bcftools norm -m -any -f "$REF_GENOME" "$INPUT_VCF" > "$NORM_VCF"
perl "$CONVERT_ANNOVAR" -format vcf4 "$NORM_VCF" > "$AV_INPUT"

# ==============================================================================
# 3. ANOTAÇÃO DAS VARIANTES
# ==============================================================================
echo "GenoLaudo - [2/3] Executando ANNOVAR..."

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
# 4. CORREÇÃO DE CABEÇALHO E INJEÇÃO DE COLUNAS PARA INTERVAR
# ==============================================================================
echo "GenoLaudo - Adaptando multianno.txt para compatibilidade com InterVar..."

# Ajusta nomenclatura do refGene se necessário
if [ "$REFGENE_NAME" != "refGene" ]; then
    sed -i "1s/$REFGENE_NAME/refGene/g" "$RAW_TXT"
fi

# Injeção de colunas de frequência fantasmas para evitar o erro BA1 (Cromossomo 7 interpretado como freq)
# Esta implementação via Python inline garante que o InterVar leia freq=0 em vez de colunas deslocadas.
python3 - <<EOF
import sys
import os

file_path = "$RAW_TXT"
temp_path = file_path + ".tmp"

with open(file_path, 'r') as f_in, open(temp_path, 'w') as f_out:
    header = f_in.readline().strip().split('\t')

    # Colunas que o InterVar busca via hardcoding para calcular critérios populacionais
    missing_cols = ['esp6500siv2_all', '1000g2015aug_all', 'gnomAD_genome_ALL', 'gnomAD_exome_ALL', 'AF']

    # Identifica quais colunas realmente precisam ser adicionadas
    to_add = [c for c in missing_cols if c not in header]

    # Escreve o novo cabeçalho
    new_header = header + to_add
    f_out.write('\t'.join(new_header) + '\n')

    # Escreve as linhas injetando "0" nas colunas fantasmas
    for line in f_in:
        f_out.write(line.strip() + '\t' + '\t'.join(['0'] * len(to_add)) + '\n')

os.replace(temp_path, file_path)
EOF

echo "  -> Colunas populacionais injetadas com sucesso (freq=0)."

# ==============================================================================
# 5. CLASSIFICAÇÃO DAS VARIANTES USANDO O INTERVAR
# ==============================================================================
echo "GenoLaudo - [3/3] Executando InterVar (ACMG)..."

ln -sf "$ANNOVAR_DIR/annotate_variation.pl" ./annotate_variation.pl
ln -sf "$ANNOVAR_DIR/convert2annovar.pl" ./convert2annovar.pl

INTERVAR_EXEC=""
if [ -f "$INTERVAR_DIR/InterVar.py" ]; then INTERVAR_EXEC="$INTERVAR_DIR/InterVar.py"; fi
if [ -f "$INTERVAR_DIR/Intervar.py" ]; then INTERVAR_EXEC="$INTERVAR_DIR/Intervar.py"; fi

if [ -n "$INTERVAR_EXEC" ]; then
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

echo "  [AVISO] Removendo arquivos temporários..."
rm -f ./annotate_variation.pl ./convert2annovar.pl
rm -f "$RAW_TXT.grl_p" "$AV_INPUT"

echo "=============================================================================="
echo "Anotação concluída."