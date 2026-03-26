#!/bin/bash

HOST="0.0.0.0"
DATA_CLEANING="$1" #ON ou OFF para limpeza de dados temporários
#definição de cores para o splash screen
CYAN='\033[1;36m'
NC='\033[0m' #serve para voltar o terminal para a cor normal

echo -e "${CYAN}"
cat << "EOF"
  █████████                             █████                             █████
 ███░░░░░███                           ░░███                             ░░███
███     ░░░   ██████  ████████   ██████ ░███       ██████   █████ ████ ███████  ██████
███          ███░░███░░███░░███ ███░░███░███      ░░░░░███ ░░███ ░███ ███░░███ ███░░███
███    █████░███████  ░███ ░███░███ ░███░███       ███████  ░███ ░███░███ ░███░███ ░███
░███  ░░███ ░███░░░   ░███ ░███░███ ░███░███      ███░░███  ░███ ░███░███ ░███░███ ░███
 ░░█████████ ░░██████ ████ ████░░██████ █████████░░████████ ░░████████░░████████░░██████
  ░░░░░░░░░   ░░░░░░ ░░░░ ░░░░  ░░░░░░ ░░░░░░░░░  ░░░░░░░░   ░░░░░░░░  ░░░░░░░░  ░░░░░░
EOF
echo "Bem vindo ao GenoLaudo!"
echo "Versão: 0.1.0"
echo "Criado por: João Gabriel Barbosa de Luna"
echo "Referências bibliográficas em REFERENCES.md"
echo -e "${NC}"

#execução do snakemake para processar as variantes
echo -e "${CYAN}"
echo "Executando pipeline snakemake..."
echo -e "${NC}"
snakemake --cores all

#inicialização do servidor do streamlit app.py, com a interface visual
echo -e "${CYAN}"
echo "Executando interface visual..."
echo -e "${NC}"
# streamlit run app.py --server.address $HOST

#remoção dos arquivos temporários, para manter a segurança dos dados
if [[ "$DATA_CLEANING" == "ON" ]]; then
  echo -e "${CYAN}"
  echo "Iniciando limpeza dos dados temporários..."
  echo -e "${NC}"
  rm data/temp/*.intervar
  rm data/temp/*.vcf
  rm data/temp/*.txt
  rm data/temp/*.json
elif [[ "$DATA_CLEANING" == "OFF" ]]; then
  :
else
  echo "Valor inválido em DATA_CLEANING: $DATA_CLEANING"
  exit 1
fi




