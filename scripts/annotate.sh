#!/bin/bash
set -euo
set -x

# ==============================================================================
# CONFIGURAÇÃO DE CAMINHOS DE ARQUIVOS A SEREM LIDOS E CRIADOS
# ==============================================================================
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANNOVAR_DIR="$PROJECT_ROOT/data/annovar"
INTERVAR_DIR="$PROJECT_ROOT/data/intervar"
TABLE_ANNOVAR="$ANNOVAR_DIR/table_annovar.pl"
CONVERT_ANNOVAR="$ANNOVAR_DIR/convert2annovar.pl"
DB_DIR="$ANNOVAR_DIR/humandb/dbs"
GENOME_DIR="$PROJECT_ROOT/data/genome"

# Verifica se os parâmetros foram passados corretamente

if [ "$#" -ne 2 ]; then
    echo "Uso: $0 <input_vcf> <output_basename>"
    exit 1
fi

INPUT_VCF=$1
OUTPUT_BASE=$2 # O nome do arquivo, que sera utilizado para nomear os outros arquivos gerados

# Aviso ao usuário, de que em caso de repetição, os arquivos serão sobrescritos

if [[ -e "$(dirname "$OUTPUT_BASE")" ]]; then
  echo "O diretório $OUTPUT_BASE já existe. Os arquivos serão sobrescritos."
  :
fi

mkdir -p "$(dirname "$OUTPUT_BASE")"

echo "=============================================================================="
echo "PathoClin - Iniciando Anotação das variantes filtradas em $INPUT_VCF..."
echo "=============================================================================="

# ==============================================================================
# 1. VERIFICAÇÃO DOS BANCOS DE DADOS ANNOVAR
# ==============================================================================
# String vazia para o comando do annovar
PROTOCOLS=""
OPERATIONS=""

#Função para identificar os bancos de dados dentro da pasta humandb/dbs, verificando o arquivo txt correspondente
get_db_filename() {
    ls "$DB_DIR"/hg38_*"$1"*.txt 2>/dev/null | head -n 1
}

# Busca dos bancos a partir da chamada da função anterior
# Ao identificar o banco correspondente, a string de protocols e operations é atualizada com as informacoes correspondentes

# RefGene
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

# ===================================================================================================
# 2. OBTENÇÃO DO GENOMA DE REFERENCIA, NORMALIZAÇÃO DO VCF FILTRADO E CONVERSÃO PARA FORMATO AVINPUT
# ===================================================================================================
echo "PathoClin - [1/3] Verificando mismatches, preparando genoma, normalizando e convertendo para ANNOVAR..."

# Obtenção dos genomas presentes na pasta genomes
REF_GENOME=$(ls "$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fasta 2>/dev/null | head -n 1)

if [ -z "$REF_GENOME" ]; then
    echo "[ERRO] Nenhum genoma de referência (.fa ou .fasta) encontrado em $GENOME_DIR."
    exit 1
fi

echo "  -> Genoma de referência localizado: $(basename "$REF_GENOME")"

#Verifica se o arquivo .fai existe, caso não, ele é criado automaticamente
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "  -> Índice .fai não encontrado. Gerando índice dinamicamente com samtools..."
    samtools faidx "$REF_GENOME"
    echo "  -> Índice criado com sucesso!"
fi

#Realizando a normalização
echo "  -> Iniciando normalização..."
NORM_VCF="${OUTPUT_BASE}.norm.vcf"
bcftools norm -m -any --check-ref w -f "$REF_GENOME" "$INPUT_VCF" > "$NORM_VCF"
echo "  -> VCF normalizado com sucesso!"

# Criação do AVINPUT a partir do VCF filtrado
echo "  -> Iniciando conversão para formato AVINPUT..."
AV_INPUT="${OUTPUT_BASE}.avinput"
perl "$CONVERT_ANNOVAR" -format vcf4 "$NORM_VCF" > "$AV_INPUT"
echo "  -> Conversão para formato AVINPUT concluída."

# ==============================================================================
# 3. ANOTAÇÃO DAS VARIANTES UTILIZANDO O ANNOVAR
# ==============================================================================
echo "PathoClin - [2/3] Executando anotação ANNOVAR..."

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
echo "  -> Adaptando multianno.txt para compatibilidade com InterVar..."

# Ajusta nomenclatura do refGene se necessário para evitar erros do intervar
if [ "$REFGENE_NAME" != "refGene" ]; then
    sed -i "1s/$REFGENE_NAME/refGene/g" "$RAW_TXT"
fi

# Injeção de colunas de frequência para evitar erros durante a leitura da anotação 
# Caso as colunas de frequencia não utilizadas permaneçam vazias, ocorre um erro de leitura
# este erro torna todas variantes benigas (criterio BA1 e BS1 sempre)
# O código python abaixo obtém as colunas de frequencia não preenchidas, injeta 0 em cada uma, e retorna um novo arquivo

python3 - <<EOF
import sys
import os

file_path = "$RAW_TXT"
temp_path = file_path + ".tmp"

with open(file_path, 'r') as f_in, open(temp_path, 'w') as f_out:
    header = f_in.readline().strip().split('\t')

    # Lista das colunas de frequências lidas pelo InterVar
    missing_cols = [
        'esp6500siv2_all', '1000g2015aug_all', 'gnomAD_genome_ALL', 'gnomAD_exome_ALL',
        'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr',
        'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'Kaviar_AF', 'CG46'
    ]

    # Identifica quais colunas precisam ter um 0 adicionado
    to_add = [c for c in missing_cols if c not in header]

    # Escreve o novo cabeçalho
    new_header = header + to_add
    f_out.write('\t'.join(new_header) + '\n')

    # Escreve as linhas injetando 0 nas colunas vazias
    for line in f_in:
        f_out.write(line.strip() + '\t' + '\t'.join(['0'] * len(to_add)) + '\n')

os.replace(temp_path, file_path)
EOF

echo "  -> Colunas populacionais injetadas com sucesso (freq=0)."

# ==============================================================================
# 5. CLASSIFICAÇÃO DAS VARIANTES USANDO O INTERVAR
# ==============================================================================
echo "PathoClin - [3/3] Executando InterVar..."

#Cria atalhos para os scripts annovar na raiz do projeto, para evitar erros de dependencia do InterVar, que utiliza scripts annovar
#Teoricamente, nós pulamos a etapa de anotação annovar, que ja foi realizada, mas esta funcionalidade evita erros em caso de verificacoes de dependências
ln -sf "$ANNOVAR_DIR/annotate_variation.pl" ./annotate_variation.pl
ln -sf "$ANNOVAR_DIR/convert2annovar.pl" ./convert2annovar.pl

# Armazena o caminho para o script intervar.py
INTERVAR_EXEC=""

# Busca o script do intervar utilizando letras minusculas e maiusculas (Var e var)
# Evita erros de localização por problemas de case sensitive
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
echo "=============================================================================="
