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

#remoção dos arquivos temporários, para manter a segurança dos dados
if [[ "$DATA_CLEANING" == "ON" ]]; then
  ./clean_data.sh
elif [[ "$DATA_CLEANING" == "OFF" ]]; then
  :
else
  echo "Valor inválido em DATA_CLEANING: $DATA_CLEANING"
  exit 1
fi




