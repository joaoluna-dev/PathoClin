import streamlit as st
import json
from weasyprint import HTML
import pandas as pd
import pathlib
import os
import datetime

#referências científicas
#MACARTHUR, D. G. et al. Guidelines for investigating causality of sequence variants in human disease.
# Nature, abr. 2014. v. 508, n. 7497, p. 469–476. Acesso em: 4 nov. 2019.
#
#PEDERSEN, B. S. et al. Effective variant filtering and expected candidate variant yield in studies of
# rare human disease. npj Genomic Medicine, 15 jul. 2021. v. 6, n. 1. doi.org/10.1038/s41525-021-00227-3
#
#RICHARDS, S. et al. Standards and guidelines for the interpretation of sequence variants:
# a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association
# for Molecular Pathology. Genetics in medicine : official journal of the American College of Medical Genetics,
# v. 17, n. 5, p. 405–24, 2015.


if 'info_message' not in st.session_state:
    st.session_state.info_message = None


#funções principais

def filter_data(variants_df, exame, genes):
    #(RICHARDS et al., 2015), (MACARTHUR et al., 2014), (PEDERSEN et al., 2021)
    variants_df = variants_df[variants_df["Frequência ABraOM"].astype(float) <= 0.01]
    variants_df = variants_df.drop('Frequência ABraOM', axis=1)

    acmg_classification = ["pathogenic", "likely pathogenic", "uncertain significance"]
    variants_df = variants_df[variants_df["Classificação"].str.lower().isin(acmg_classification)]

    if (exame == "Painel Genético" or exame == "Gene") and not genes:
        st.session_state.info_message = "Nenhum gene foi selecionado para o exame {exame}"
        return None

    elif exame == "Gene" and genes:
        genes = genes.strip()
        filtered_variants_df = variants_df[variants_df["Gene"].str.contains(rf'\b{genes}\b', na=False, regex=True)]
        return filtered_variants_df

    elif exame == "Painel Genético" and genes:
        parsed_genes = genes.split()
        padrao_regex = rf'\b(?:{"|".join(parsed_genes)})\b'
        filtered_variants_df = variants_df[variants_df["Gene"].str.contains(padrao_regex, na=False, regex=True)]
        return filtered_variants_df

    elif exame == "Exoma":
        padrao_regex = r'\b(?:exonic|splicing)\b'
        filtered_variants_df = variants_df[variants_df["Localização"].str.contains(padrao_regex, na=False, regex=True)]
        return filtered_variants_df

    return None


def create_html(filtered_dataframe):
    in_house_str = "Sim" if metodologia_in_house else "Não"
    restricao_str = "Sim" if material_com_restricao else "Não"
    material_final = tipo_material_outro if tipo_material == "Outro" else tipo_material

    nome_social_html = f"<strong>Nome Social:</strong> {nome_social}<br>" if nome_social else ""
    alvos_html = f"<strong>Alvos (Genes):</strong> {targets}<br>" if targets is not None else ""
    hora_coleta_html = f"<strong>Horário da Coleta:</strong> {horario_coleta}<br>" if horario_coleta else ""
    info_clinicas_html = f"<strong>Informações Clínicas:</strong> {info_clinicas}<br>" if info_clinicas else ""

    if filtered_dataframe is not None and not filtered_dataframe.empty:
        tabela_variantes_html = filtered_dataframe.to_html(index=False, justify='center', border=0,
                                                           classes='tabela-resultados')
    else:
        tabela_variantes_html = "<p style='text-align:center;'><strong>Nenhuma variante encontrada para os critérios selecionados.</strong></p>"

    html_template = f"""
    <!DOCTYPE html>
    <html lang="pt-BR">
    <head>
        <meta charset="UTF-8">
        <title>Laudo Genômico - {nome_paciente}</title>
        <style>
            @page {{
                size: A4;
                margin: 20mm 15mm 20mm 15mm;
            }}
            body {{
                font-family: Arial, Helvetica, sans-serif;
                font-size: 11pt;
                color: #2c3e50;
                line-height: 1.4;
            }}
            .header {{
                text-align: center;
                border-bottom: 2px solid #2c3e50;
                padding-bottom: 10px;
                margin-bottom: 20px;
            }}
            .header h1 {{ margin: 0; font-size: 18pt; }}
            .header p {{ margin: 3px 0; font-size: 10pt; color: #7f8c8d; }}

            .section {{
                margin-bottom: 20px;
                page-break-inside: avoid;
            }}
            .section-title {{
                background-color: #ecf0f1;
                padding: 5px;
                font-size: 12pt;
                font-weight: bold;
                border-left: 4px solid #3498db;
                margin-bottom: 10px;
                text-transform: uppercase;
            }}
            .grid-container {{
                display: table;
                width: 100%;
            }}
            .grid-item {{
                display: table-cell;
                width: 50%;
                padding-right: 10px;
            }}

            .tabela-resultados {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                font-size: 8pt;
                text-align: center;
                table-layout: fixed;
                word-wrap: break-word;
            }}
            .tabela-resultados th, .tabela-resultados td {{
                border: 1px solid #bdc3c7;
                padding: 6px;
                overflow-wrap: break-word;
            }}
            .tabela-resultados th {{
                background-color: #34495e;
                color: white;
                font-weight: bold;
            }}
            .tabela-resultados tr:nth-child(even) {{ background-color: #f9f9f9; }}

            .footer {{
                margin-top: 40px;
                border-top: 1px solid #bdc3c7;
                padding-top: 10px;
                font-size: 9pt;
                page-break-inside: avoid;
            }}
            .signatures {{
                display: table;
                width: 100%;
                text-align: center;
                margin-top: 30px;
            }}
            .signature-box {{
                display: table-cell;
                width: 50%;
                padding: 0 20px;
            }}
            .signature-line {{
                border-top: 1px solid #2c3e50;
                margin: 0 auto 5px auto;
                width: 80%;
            }}
            .acmg-note {{
                text-align: justify;
                font-size: 8pt;
                color: #7f8c8d;
                margin-top: 15px;
            }}
        </style>
    </head>
    <body>

        <div class="header">
            <h1>{nome_servico}</h1>
            <p>CNES: {cnes_servico} | Contato: {endereco_telefone}</p>
        </div>

        <div class="section">
            <div class="section-title">Dados do Paciente</div>
            <div class="grid-container">
                <div class="grid-item">
                    <strong>Nome:</strong> {nome_paciente}<br>
                    {nome_social_html}
                    <strong>Registro (RG/CPF/Prontuário):</strong> {rg_identificacao}<br>
                    <strong>Data de Nascimento:</strong> {data_nascimento}<br>
                </div>
                <div class="grid-item">
                    <strong>Sexo Biológico:</strong> {sexo_biologico}<br>
                    <strong>Nome da Mãe:</strong> {nome_mae}<br>
                </div>
            </div>
        </div>

        <div class="section">
            <div class="section-title">Dados da Amostra e Solicitação</div>
            <div class="grid-container">
                <div class="grid-item">
                    <strong>Médico Solicitante:</strong> {nome_solicitante}<br>
                    <strong>Exame:</strong> {nome_exame}<br>
                    {alvos_html}
                    <strong>Material Aceito com Restrição:</strong> {restricao_str}<br>
                </div>
                <div class="grid-item">
                    <strong>Data da Coleta:</strong> {data_coleta}<br>
                    {hora_coleta_html}
                    <strong>Tipo de Material Biológico:</strong> {material_final}<br>
                </div>
            </div>
        </div>

        <div class="section">
            <div class="section-title">Informações Analíticas e Clínicas</div>
            <strong>Método Analítico:</strong> {metodo_analitico}<br>
            <strong>Valores de Referência / Genoma Base:</strong> {valores_referencia}<br>
            <strong>Metodologia Própria (in house):</strong> {in_house_str}<br>
            {info_clinicas_html}
        </div>

        <div class="section">
            <div class="section-title" style="text-align: center; border-left: none; background-color: #34495e; color: white;">
                Resultados da Análise Genômica
            </div>
            {tabela_variantes_html}
        </div>

        <div class="section">
            <div class="section-title">Limitações Técnicas e Interpretação</div>
            <p style="text-align: justify; font-size: 10pt;">{limitacoes_teste}</p>
        </div>

        <div class="footer">
            <p><strong>Data de Emissão do Laudo:</strong> {data_emissao}</p>

            <div class="signatures">
                <div class="signature-box">
                    <div class="signature-line"></div>
                    <strong>Responsável Técnico (RT)</strong><br>
                    {dados_rt}
                </div>
                <div class="signature-box">
                    <div class="signature-line"></div>
                    <strong>Profissional Responsável</strong><br>
                    {dados_assinatura}
                </div>
            </div>

            <div class="acmg-note">
                * A classificação de patogenicidade das variantes reportadas neste laudo segue as diretrizes 
                conjuntas do American College of Medical Genetics and Genomics (ACMG) e da Association for 
                Molecular Pathology (AMP). Os resultados devem ser correlacionados com os achados clínicos, 
                histórico familiar e exames complementares pelo médico assistente.
            </div>
        </div>

    </body>
    </html>
    """

    return html_template


def create_pdf(html):
    HTML(string=html).write_pdf(pdf_file)
    if os.path.exists(pdf_file) and os.path.getsize(pdf_file) > 0:
        st.session_state.info_message = f"Arquivo PDF criado com sucesso em {pdf_file}"
    else:
        st.session_state.info_message = "Arquivo pdf não criado. Tente novamente."


def execute():
    try:
        filtered_data = filter_data(variants_dataframe, nome_exame, targets)
        html_string = create_html(filtered_data)
        create_pdf(html_string)
    except Exception as e:
        st.session_state.info_message = f"Ocorreu um erro durante a criação do laudo: {e}"


#header do site
st.set_page_config(
	page_title='GenoLaudo: Bioinformática Clínica',
	page_icon='🧬',
	layout='wide'
)
st.logo("logo/logo1.png", size="large")

#título principal
st.title("Genolaudo")
st.text("João Gabriel, 2026")

#caminho absoluto do json com as anotações (estático)
try:
    app_path = pathlib.Path(__file__).resolve()
    root = app_path.parent
    json_file_path = root / "data" / "temp"
    output_path = root / "data" / "output"
except Exception as e:
    st.error(f"Ocorreu um erro durante a leitura dos arquivos de anotação: {e}")
#lista com os arquivos de laudos
files = []

#obtendo os arquivos json com o laudo de informações
for file in os.listdir(json_file_path):
    if file.endswith(".json"):
        files.append(file)

counter = 0

#lendo os arquivos .json obtidos
for json_file in files:
    json_folder = root / output_path / json_file.strip("_report.json")
    if not os.path.exists(json_folder):
        os.mkdir(json_folder)
    with st.expander(f"Arquivo {json_file}", expanded=False):
        with open(f'{json_file_path}/{json_file}', 'r') as file:
            data = json.load(file)
            for amostra in data:
                for sample, variants in amostra.items():
                    pdf_file = json_folder / f"{sample}.pdf"
                    variants_dataframe = pd.DataFrame(variants)
                    counter += 1
                    with st.form(f"sample_{counter}_form"):
                        st.header(f"Formulário da amostra {sample}")
                        st.subheader("1. Dados de Identificação do Paciente")
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            nome_paciente = st.text_input("Nome completo do paciente*")
                            rg_identificacao = st.text_input("Nº de Registro (RG/CPF/Prontuário)*")
                        with col2:
                            data_nascimento = st.date_input("Data de nascimento*", min_value=datetime.date(1900, 1, 1))
                            sexo_biologico = st.selectbox("Sexo biológico*",
                                                          ["Feminino", "Masculino", "Intersexo", "Não informado"])
                        with col3:
                            nome_social = st.text_input("Nome social (se aplicável)")
                            nome_mae = st.text_input("Nome da mãe")

                        st.markdown("---")

                        st.subheader("2. Dados da Amostra e do Pedido Médico")
                        col4, col5, col6 = st.columns(3)
                        with col4:
                            nome_solicitante = st.text_input("Nome do médico solicitante*")
                            nome_exame = st.selectbox("Nome do exame*",
                                                      options=("Painel Genético",
                                                               "Exoma",
                                                               "Gene"),
                                                      index=0)
                        targets = st.text_area(label="Insira os alvos do exame. Caso tenha selecionado 'Gene', "
                                                          "Escreva o nome do Gene(s) de interesse. Caso tenha selecionado "
                                                          "'Painel Genético', insira o nome dos genes, separados por espaço. "
                                                          "Caso tenha selecionado 'Exoma', deixe em branco.*").upper()
                        if not targets:
                            targets = None
                        with col5:
                            data_coleta = st.date_input("Data da coleta*")
                            horario_coleta = st.time_input("Horário da coleta")
                        with col6:
                            tipo_material = st.selectbox("Tipo de material biológico*",
                                                         ["Sangue periférico", "Saliva", "Swab bucal", "Tecido fresco",
                                                          "FFPE", "Outro"])
                        tipo_material_outro = st.text_area("Especifique o material, caso tenha selecionado 'Outro'*")

                        st.markdown("---")

                        st.subheader("3. Dados Analíticos e Clínicos")
                        col7, col8 = st.columns(2)
                        with col7:
                            metodo_analitico = st.text_input("Método analítico*",
                                                             value="Sequenciamento de Nova Geração (NGS)")
                            valores_referencia = st.text_input("Valores de referência / Genoma Base*", value="GRCh38")
                        with col8:
                            info_clinicas = st.text_area(
                                "Informações clínicas adicionais (Uso de medicamentos, histórico familiar, etc.)")

                        limitacoes_teste = st.text_area("Limitações técnicas e dados para interpretação*", height=100)

                        col9, col10 = st.columns(2)
                        with col9:
                            metodologia_in_house = st.checkbox("Exame realizado por Metodologia Própria (in house)")
                        with col10:
                            material_com_restricao = st.checkbox("Material biológico aceito com restrição")
                            if material_com_restricao:
                                motivo_restricao = st.text_input("Especifique a restrição (ex: hemólise, volume baixo)")

                        st.markdown("---")

                        st.subheader("4. Dados Institucionais e Assinatura (Preenchimento Manual)")
                        col11, col12 = st.columns(2)
                        with col11:
                            nome_servico = st.text_input("Nome do Serviço / Laboratório*")
                            cnes_servico = st.text_input("Número do CNES*")
                            endereco_telefone = st.text_area("Endereço e Telefone de contato*")
                        with col12:
                            dados_rt = st.text_input("Nome e Registro do Responsável Técnico (RT)*")
                            dados_assinatura = st.text_input("Nome e Registro do profissional que assina o laudo*")
                            data_emissao = st.date_input("Data de emissão do laudo*", value=datetime.date.today())

                        st.text("* Campos obrigatórios.")
                        # Botão de submissão do formulário
                        submit_button = st.form_submit_button(label=f"Salvar Laudo: {sample}")

                    if submit_button:
                        execute()
                        st.info(st.session_state.info_message)
                        if "Arquivo PDF criado com sucesso" in st.session_state.info_message:
                            with open(pdf_file, "rb") as report_file:
                                st.download_button(
                                    label="Baixar laudo",
                                    data=report_file,
                                    file_name=f"{nome_paciente}.pdf",
                                    mime="application/pdf",
                                    icon="📄"
                                )






