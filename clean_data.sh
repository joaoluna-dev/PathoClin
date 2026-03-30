#!/bin/bash

target=$1

CYAN='\033[1;36m'
NC='\033[0m' #serve para voltar o terminal para a cor normal


if [[ "$target" == "temp_files" ]]; then
  echo -e "${CYAN}"
  echo "Iniciando limpeza dos dados temporários..."
  echo -e "${NC}"
  rm data/temp/*.intervar
  rm data/temp/*.vcf
  rm data/temp/*.txt
  rm data/temp/*.json
elif [[ "$target" == "vcfs" ]]; then
  echo -e "${CYAN}"
  echo "Iniciando limpeza dos vcfs..."
  echo -e "${NC}"
  rm data/raw/*
fi
