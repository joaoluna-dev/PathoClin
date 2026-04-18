# PathoClin 📄🧬 

![Python](https://img.shields.io/badge/Python-3.10-blue)
![Streamlit](https://img.shields.io/badge/Streamlit-1.54-red)
![Weasyprint](https://img.shields.io/badge/weasyprint-67.0-red)
![Snakemake](https://img.shields.io/badge/Snakemake-7.32.4-teal)
![Docker](https://img.shields.io/badge/Docker-Container-2496ED)
![Apptainer](https://img.shields.io/badge/Apptainer-Container-2496ED)
![bcftools](https://img.shields.io/badge/bcftools-red)
![samtools](https://img.shields.io/badge/samtools-red)
![ANNOVAR](https://img.shields.io/badge/ANNOVAR-green)
![Intervar](https://img.shields.io/badge/Intervar-green)


Repositório oficial do projeto PathoClin, um pipeline opensource 100% automatizado pra a anotação de variantes em arquivos .vcf e geração de laudos genéticos 
padronizados conforme as diretrizes laboratoriais brasileiras da ANVISA.

## Destaques

- Processamento automatizado de arquivos .vcf fornecidos pelo usuário.
- Filtragem de variantes com base nas diretrizes da ACMG e Broad Institute.
- Anotação das variantes utilizando o ANNOVAR e o InterVar, utilizando bancos de dados ABraOm, ClinVar e dbnsfp.
- Processamento de dados da anotação e criação de laudo automatizado em PDF.
- Interface visual de fácil utilização para preenchimento de dados clínicos de cada amostra.
- Conteinerização para maior versatilidade e reprodutibilidade.
- Sistema completamente automatizado. Foque apenas no que importa.
- Instalação simples e modular. Utilize as versões de bancos de dados de seu interesse.

## Visão geral do projeto

O PathoClin é um pipeline bioinformático automatizado, desenvolvido em Python, dedicado à anotação e interpretação de 
variantes genômicas humanas a partir de arquivos VCF. Projetado para processar dados de fluxos padrão de variant calling. 
O software foca na triagem de mutações em genes de relevância médica. O seu principal objetivo é assessorar a tomada de 
decisão clínica de forma ágil e abrangente, reduzindo a complexidade operacional na rotina de diagnóstico molecular 
através de uma interface de fácil preenchimento de informações clínicas. Como diferencial tecnológico e de privacidade, 
o sistema realiza a anotação das variantes empregando bancos de dados executados localmente por meio do ANNOVAR e Intervar, 
eliminando a necessidade de submissão de dados genéticos, garantindo a segurança e sigilo às informações sensíveis dos pacientes, 
além de acelerar substancialmente o tempo de processamento analítico.

### Algorítmo
1. Entrada de Dados (Input): A inicialização do sistema requer três componentes essenciais, o diretório com os bancos de dados do ANNOVAR, o arquivo FASTA/FA do genoma de referência (na mesma versão de montagem utilizada na etapa de chamada de variantes) e o arquivo VCF bruto da amostra. A entrada destes arquivos no conteiner ocorre através de mapeamento de volumes (bind mount).
2. Filtragem de Variantes (Filtering): O VCF bruto é submetido a um processo de controle de qualidade. As variantes são avaliadas e filtradas com base nas métricas do arquivo, seguindo as diretrizes técnicas do Broad Institute (GATK) e as recomendações de boas práticas da ACMG/AMP. As variantes aprovadas são exportadas para um arquivo VCF temporário.
3. Normalização e Anotação Genômica: As variantes filtradas são primeiramente normalizadas (alinhadas à esquerda e parsimoniosas) para garantir uma representação padronizada. Em seguida, ocorre a anotação. As variantes são enriquecidas com dados clínicos (ClinVar), frequências alélicas da população brasileira [ABraOM](https://abraom.ib.usp.br/) e escores de predição computacional de impacto patogênico (dbNSFP). O produto desta etapa é salvo nos formatos nativos das ferramentas, incluindo arquivos .intervar, .multianno.txt e o .vcf normalizado.
4. Processamento e Estruturação de Dados (Parsing): Os arquivos gerados na etapa de anotação são processados para extrair os dados exigidos pelas diretrizes da ACMG para a elaboração de laudos (Coordenada Genômica, Gene, Transcrito, Nomenclatura da Variante HGVS, Zigosidade, Classificação de Patogenicidade e Doença Associada). O sistema possui suporte a VCFs multiamostra, segregando automaticamente as informações por paciente e serializando os dados processados de forma estruturada em um arquivo JSON.
5. Integração Fenotípica e Emissão de Laudo: A interface gráfica interativa é executada automaticamente. O ambiente disponibiliza formulários para o preenchimento dos dados clínicos e de identificação de cada amostra. Após a confirmação, o painel cruza os dados clínicos com o JSON estruturado, aplicando as informações em um template HTML. Finalmente, o sistema escreve esse template em um documento PDF, liberado para download direto na interface.

### Fluxograma das análises

![Fluxograma Pathoclin.png](images/Fluxograma%20Pathoclin.png)

### Dependências 

- snakemake
- cyvcf2
- bcftools
- ANNOVAR
- InterVar
- pandas
- streamlit
- weasyprint

### Equipe

João Gabriel Barbosa de Luna <br>
João Luiz de Lemos Padilha Pitta, PhD. <br>
Túlio de Lima Campos, PhD.

### Sobre mim

Biomédico, Bioinformata e Desenvolvedor Python. Trabalha com ciência de dados aplicada à genômica, transcriptômica, epigenômica e proteômica. Esta ferramenta é desenvolvida por profissionais para profissionais, para facilitar a tomada de decisões em genômica médica, democratizando o acesso à este tipo de ferramenta para centros de pesquisa e pequenos laboratórios.

### Informação dos autores

![Github](https://img.shields.io/badge/github-gray) <br>
[GitHub João Gabriel](https://github.com/joaoluna-dev) <br>
[GitHub OncogenSUS]() <br>

![Linkedin](https://img.shields.io/badge/linkedin-blue) <br>
[Linkedin João Gabriel](https://www.linkedin.com/in/jo%C3%A3o-gabriel-a3b28130a/) <br>
[Linkedin João Pitta]() <br>
[Linkedin Túlio Campos]() <br>

![Lattes](https://img.shields.io/badge/lattes-white) <br>
[Lattes João Gabriel](https://lattes.cnpq.br/3598905984017635) <br>
[Lattes João Pitta]() <br>
[Lattes Túlio Campos]() <br>

## Instalação

### Pré-requisitos

- [Docker](https://docs.docker.com/get-docker/) ou [Apptainer](https://apptainer.org/docs/admin/main/installation.html#) instalados.

### Configuração

1.  **Clone o repositório:**
    ```bash
    git clone https://github.com/joaoluna-dev/GenoLaudo.git
    cd GenoLaudo
    ```

2. **Realizar a instalação das ferramentas auxiliares**
   Execute o script abaixo para iniciar o download e instalação do ANNOVAR e do Intervar, e a criação da estrutura de pastas dos dados gerados
    ```bash
    chmod +x install_pathoclin.sh
    ./install_pathoclin.sh
    ```
> **ANNOVAR e Bancos de dados**: para completar a instalação, você precisa de um link para download do ANNOVAR, que pode ser obtido [aqui](https://www.openbioinformatics.org/annovar/annovar_download_form.php). Para obter os bancos de dados ClinVar, ABraOM e dbSNFP nas versões de interesse, clique [aqui](https://annovar.openbioinformatics.org/en/latest/user-guide/download/).

3. **Criar a imagem:**
   Execute um dos comandos abaixo para iniciar a criação da imagem:

- Docker
    ```bash
    docker build -t genolaudo .
    ```

- Apptainer
    ```bash
    apptainer build ./PathoClin.sif ./PathoClin.def
    ```

> Este procedimento pode levar um tempo.

## Utilização

1. **Criando o contêiner e iniciando as análises:**
    Execute um dos comandos abaixo para executar a imagem criada:

- Docker:

    ```bash
    docker run -d --name genolaudo -p 8501:8501 -v /caminho/local/para/diretorio/bancos_annovar:/app/data/annovar/humandb/dbs:ro -v /caminho/local/para/os/genomas_de_referencia:/app/data/genome:ro genolaudo:latest
    ```
    `-d`: executa o container em modo detached (desacoplado), liberando o terminal para uso <br>
    `--name genolaudo`: define o nome do container <br>
    `-p`: abre a porta da interface visual, para que ela possa ser executada <br>
    `-v`: faz a montagem dos arquivos necessários para a execução (o ou os VCF, o diretório com os bancos de dados do ANNOVAR, e o genoma de referência da amostra). Os arquivos são montados em somente leitura (:ro), evitando alterações nos dados originais, e a execução de scripts maliciosos injetados nos arquivos <br>
    `genolaudo:latest`: indica qual imagem deve ser utilizada <br>

- Apptainer:

    ```bash
    nohup apptainer run --bind /caminho/local/para/diretorio/bancos_annovar:/app/data/annovar/humandb/dbs:ro --bind /caminho/local/para/os/genomas_de_referencia:/app/data/genome:ro genolaudo.sif > genolaudo.log 2>&1 &
    ```
    `nohup ... &`: executa o container em segundo plano, liberando o terminal para uso e salvando os registros de execução no arquivo `genolaudo.log` <br>
    `apptainer run`: inicia a execução do container acionando o script principal da imagem <br>
    `--bind`: faz a montagem dos arquivos necessários para a execução (o VCF, o diretório com os bancos de dados do ANNOVAR, e o genoma de referência da amostra). Os arquivos são montados em somente leitura (:ro), evitando alterações nos dados originais e a execução de scripts maliciosos injetados nos arquivos <br>
    `genolaudo.sif`: indica qual imagem compilada do Apptainer deve ser utilizada <br>

2. **Acessando a interface visual:**
   Aguarde até que o terminal indique que o servidor do Streamlit está em execução e, em seguida, abra o seu navegador web:

   - 📊 Interface Gráfica: http://localhost:8501

> **Acesso em Rede Local (LAN)**: A aplicação pode ser acessada a partir de qualquer dispositivo conectado à mesma rede local. Para isso, basta utilizar o endereço IP da máquina *host* (onde o container está rodando) seguido da porta exposta:
  http://<ENDERECO_IP_DO_HOST>:8501


3. **Preenchimento de dados e geração do laudo:**
   Utilize a interface gráfica interativa para preencher as informações clínicas do paciente e compilar o resultado final:

   - **Informações Clínicas:** Navegue pelos formulários disponibilizados na interface para preencher os dados de identificação do paciente e as informações fenotípicas/clínicas referentes à amostra.
     - Identificação do paciente: Nome, Registro, Sexo, Data de nascimento, Nome social e Nome da mãe
     - Dados da amostra e pedido médico: Nome do médico solicitante, Exame solicitado (Gene, Painel ou Exoma), Data da coleta e Tipo de amostra
     - Dados clínicos: Método analítico, Genoma de base, Informações clínicas adicionais
     - Dados institucionais e assinatura: Nome do laboratório, CNES, Endereço e contato, Responsável técnico, Responsável pela liberação do laudo
   - **Revisão e Integração:** Após o preenchimento, confirme os dados. O sistema irá cruzar automaticamente as suas requisições médicas com as variantes genéticas que foram previamente filtradas e anotadas pelo *pipeline* para a amostra correspondente, retornando apenas as variantes não-benignas do laudo, com frequência < 0.01.
   - **Download do Documento:** O sistema aplicará os dados integrados em um *template* HTML que compõe o laudo genômico final, escrevendo-o em um arquivo PDF. O arquivo PDF resultante ficará imediatamente disponível para *download* na própria interface.


4. **Encerrando o container:**
   Você pode interromper a execução do *software* a qualquer momento rodando o seguinte comando no seu terminal:
   
   ```bash
   docker stop genolaudo
   ```

> Após a finalização da execução da ferramenta, o sistema limpará todos os arquivos temporários, para garantir a segurança dos dados. Certifique-se de que os laudos foram preenchidos e processados corretamente!

## Próximas implementações

1. **Integração de Inteligência Artificial (AI Overview):**
    - Implementação de um modelo de linguagem (LLM) executado localmente, operando em conjunto com a arquitetura RAG (*Retrieval-Augmented Generation*) conectada à base de dados do GeneReviews. O objetivo é gerar resumos clínicos automatizados, contextualizados e seguros sobre as síndromes e patologias identificadas nas variantes do paciente.

## Arquivos gerados

Durante a execução do *pipeline*, o GenoLaudo gera arquivos para a análise das variantes. Estes dados são gerenciados automaticamente:

* `_filtered.vcf`: Arquivo VCF temporário que contém as variantes aprovadas na etapa de filtragem do pipeline
* `.norm.vcf`: Arquivo VCF temporário que contém as variantes aprovadas na etapa de filtragem de qualidade e que foram estruturalmente normalizadas pelo BCFTools (alinhadas à esquerda e parsimoniosas). Garante a padronização das variantes antes da etapa de anotação.
* `.multianno.txt`: Relatório tabular gerado pelo ANNOVAR. Contém a anotação de cada variante, cruzando as informações com bancos de dados utilizados (ClinVar), frequência populacional (ABraOM) e predições *in silico* de impacto patogênico (dbNSFP).
* `.intervar`: Arquivo gerado pelo algoritmo InterVar. Apresenta o cálculo de patogenicidade automatizado para cada variante identificada, seguindo os 28 critérios estabelecidos pelas diretrizes da ACMG.
* `._report.json`: Arquivo estruturado final, gerado pelo script de processamento do pipeline (parser.py). Ele une os dados de anotação e patogenicidade apenas das variantes aprovadas, segregando as informações por amostra. É este arquivo que alimenta a interface visual e o laudo em PDF. A sua estrutura obedece ao seguinte formato:

```json
[
    {
        "nome_da_amostra_no_vcf": [
            {
                "Coordenada Genômica": "ex: chr7:117559590:ATCT:A",
                "Gene": "ex: CFTR",
                "Transcrito": "ex: NM_000492.4",
                "Localização": "ex: exonic (nonframeshift deletion)",
                "Variante": "ex: c.1520_1522del (p.F508del)",
                "Zigosidade": "ex: Homozygous",
                "Classificação": "ex: Pathogenic",
                "Evidência ACMG": "ex: Classificação obtida da base de dados do Clinvar",
                "Doença": "ex: Cystic fibrosis...",
                "Herança": "ex: Autosomal dominant...",
                "Frequência ABraOM": 0.002463
            }, 
            {
                "dados da variante 2": "..."
            }
        ]
    }
]
```

## Feedback e Contribuição

Contribuições são bem-vindas! Por favor, reporte bugs através do [GitHub Issues](https://github.com/seu-usuario/genolaudo/issues), incluindo os logs do Docker e os passos para reprodução.

**Como contribuir:**
1. Faça um Fork do repositório e crie uma branch de funcionalidade (feature branch).
2. Certifique-se de que o seu código passa nos padrões básicos de formatação (ferramentas de Lint).
3. Submeta um Pull Request direcionado à branch `main`.

**Contato:** João Gabriel ([joaogabrieldeluna@gmail.com](mailto:joaogabrieldeluna@gmail.com))
