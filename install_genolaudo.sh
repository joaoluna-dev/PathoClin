#!/bin/bash

echo "Bem vindo ao assistente de instalação do GenoLaudo!"
echo "São necessários no mínimo 70Gb de espaço livre de armazenamento para continuar"

#verifica se o usuário possui o armazenamento necessário para a execução
read -p "Deseja prosseguir? (s/n)" selection
if [[ $selection == "s" ]]; then
  :
elif [[ $selection == "n" ]]; then
  exit 1
fi

echo "É necessário um link de download do ANNOVAR para prosseguir."
echo "Para mais informações, acesse https://annovar.openbioinformatics.org/en/latest/user-guide/download/"

#obtendo link de download do ANNOVAR
read -p "Insira o seu link de download do ANNOVAR: " annovar_link

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
mv hg38_refGeneWithVer.txt dbs/
mv hg38_refGeneWIthVerMrna.fa dbs/

#sanity check das databases padrão do annovar
if [[ "$(ls -A data/annovar/humandb/dbs | wc -l)" -eq 0 ]]; then
  :
else
  echo "Os bancos de dados padrões do annovar não estão na pasta annovar/humandb/dbs. Interrompendo execução."
  exit 1
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
rm -rf InterVar-master/

#criando estrutura de pastas para o projeto
echo "Criando estrutura de pastas para o projeto"
mkdir data/temp
mkdir data/raw
mkdir data/genome
mkdir data/output