#!/bin/bash

echo "Bem vindo ao assistente de instalação do GenoLaudo!"
echo "São necessários no mínimo 70Gb de espaço livre de armazenamento para continuar"

#verifica se o usuário possui o armazenamento necessário para a execução
read -rp "Deseja prosseguir? (s/n): " selection
if [[ $selection == "s" ]]; then
  :
elif [[ $selection == "n" ]]; then
  exit 1
fi

echo "É necessário um link de download do ANNOVAR para prosseguir."
echo "Para mais informações, acesse https://annovar.openbioinformatics.org/en/latest/user-guide/download/"

#obtendo link de download do ANNOVAR
read -rp "Insira o seu link de download do ANNOVAR: " annovar_link

#baixando o ANNOVAR
curl "$annovar_link" -o annovar.latest.tar.gz

#extraindo o ANNOVAR
tar -xvzf annovar.latest.tar.gz

#movendo o ANNOVAR para o diretório /data
mkdir data
mv annovar data/

#sanity check do ANNOVAR
if [[ -e data/annovar ]]; then
  :
else
  echo "O diretório do annovar não foi corretamente instalado. Interrompendo execução."
  exit 1
fi

mkdir data/annovar/humandb/dbs

#movendo os bancos de dados padrões do annovar para o diretório dbs, que na montagem receberá os bancos de dados do usuário
mv data/annovar/humandb/hg38_refGeneWithVer.txt data/annovar/humandb/dbs/
mv data/annovar/humandb/hg38_refGeneWIthVerMrna.fa data/annovar/humandb/dbs/

#sanity check das databases padrão do annovar
if [[ "$(ls -A data/annovar/humandb/dbs | wc -l)" -eq 0 ]]; then
  echo "Os bancos de dados padrões do annovar não estão na pasta annovar/humandb/dbs. Interrompendo execução."
  exit 1
else
  :
fi

# --baixando e instalando o intervar--
echo "Iniciando a instalação do InterVar..."
#criando diretório
mkdir data/intervar

#baixando o mim2gene.txt
echo "Baixando arquivo OMIM..."
curl https://www.omim.org/static/omim/data/mim2gene.txt -o mim2gene.txt

#baixando o intervar
echo "Baixando InterVar..."
wget https://github.com/WGLab/InterVar/archive/master.zip -O InterVar.zip

#extraindo o intervar
unzip InterVar.zip

#movendo os arquivos do intervar e removendo arquivos residuais
echo "Instalando InterVar..."
mv InterVar-master/* data/intervar/

#passando o arquivo mim2gene do OMIM para o diretório /intervar
mv mim2gene.txt data/intervar/intervardb/

#removendo arquivos temporários
rm InterVar.zip
rm annovar.latest.tar.gz
rm -rf InterVar-master/

#criando estrutura de pastas para o projeto
echo "Criando estrutura de pastas para o projeto"
mkdir data/temp
mkdir data/raw
mkdir data/genome
mkdir data/output

#baixando arquivos de suporte para o pipeline
echo "A instalação das ferramentas auxiliares está concluida."
read -rp "Deseja baixar os genomas de referência (hg19 e hg38)? (s/n): " selection_genome
if [[ $selection_genome == "s" ]]; then
  #baixando arquivos de genoma de referência
  #hg19
  curl -C - -L -o hg19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
  gunzip hg19.fa.gz

  #hg38
  curl -C - -L -o hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip hg38.fa.gz

  #movendo os arquivos para um diretório de interesse do usuário

  read -rp "Indique o local onde os arquivos devem ser armazenados (ex: /home/user/data/): " genome_place

  mv hg38.fa "$genome_place"
  mv hg19.fa "$genome_place"

  echo "Dados movidos com sucesso."

elif [[ $selection_genome == "n" ]]; then
  :
fi

#baixando arquivos de bancos de dados do annovar
read -rp "Deseja baixar os bancos de dados ANNOVAR (ABraOM, ClinVar e dbNSFP) mais recentes? (s/n): " selection_dbs
if [[ $selection_dbs == "s" ]]; then
  #verifica se o perl está instalado no sistema, caso contrário, instala
  if ! command -v perl &> /dev/null; then
    sudo apt update
    sudo apt install perl -y
  fi
  #baixando arquivos de banco de dados hg38
  echo "Baixando bancos de dados para a versão hg38..."
  perl data/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar abraom .
  perl data/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20250721 .
  perl data/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42c .

  #baixando arquivos de banco de dados hg19
  echo "Baixando bancos de dados para hg19..."
  perl data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar abraom .
  perl data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20250721 .
  perl data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp42c .

  echo "Bancos de dados baixados com sucesso!"

  #movendo os arquivos baixados para o diretório de interesse
  read -rp "Indique o local onde os arquivos devem ser armazenados (ex: /home/user/data/): " data_place

  mv hg38_abraom* "$data_place"
  mv hg38_clinvar* "$data_place"
  mv hg38_dbnsfpc* "$data_place"
  mv hg19_abraom* "$data_place"
  mv hg19_clinvar* "$data_place"
  mv hg19_dbnsfpc* "$data_place"

  echo "Dados movidos com sucesso."
  echo "Todos os arquivos baixados devem ser montados em conjunto com o container que executará a aplicação."

elif [[ $selection_dbs == "n" ]]; then
  echo "Finalizando instalação..."
  exit 1
fi

echo "Finalizando instalação"




