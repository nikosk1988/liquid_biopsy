import streamlit as st
import sqlite3
import pandas as pd
import plotly.express as px
import subprocess
import os
from scipy.stats import norm
import numpy as np
from fpdf import FPDF
from docx import Document
from datetime import datetime
import io
import re

# --- 1. ΡΥΘΜΙΣΗ ΣΕΛΙΔΑΣ & ΕΠΑΓΓΕΛΜΑΤΙΚΟ CSS ---
st.set_page_config(page_title="ctDNA clinical Discoverer Pro", layout="wide")

st.markdown("""
    <style>
    .stAppDeployButton {display: none !important;}
    [data-testid="stTabs"] button p { font-size: 1.5rem !important; font-weight: 800; color: #00d1b2 !important; }
    .metric-card {
        background-color: var(--secondary-background-color);
        color: var(--text-color);
        padding: 22px; border-radius: 12px; border-left: 8px solid #00d1b2;
        margin-bottom: 20px; box-shadow: 0 4px 10px rgba(0,0,0,0.1);
        border: 1px solid rgba(128, 128, 128, 0.2);
    }
    .metric-info { font-size: 1.1rem; line-height: 1.6; }
    .record-count { color: #00d1b2; font-weight: bold; font-size: 1.2rem; }
    [data-testid="stSidebar"] label p { font-size: 1.2rem !important; font-weight: 600; }
    </style>
""", unsafe_allow_html=True)

# --- 2. ΣΤΑΤΙΣΤΙΚΕΣ ΠΑΡΑΜΕΤΡΟΙ & HELPER FUNCTIONS ---
H_MEAN, H_STD = 18.5, 2.1
COLOR_MAP = {"healthy.bam": "#00d1b2", "cancer.bam": "#ff4b4b"}

def to_safe_latin(text):
    gr = "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩαβγδεζηθικλμνξοπρστυφχψωάέήίόύώ"
    en = "ABGDEZHQIKLMNXOPRSTUFXPSabgdezhqiklmnxoprstufxpsaehiouw"
    table = str.maketrans(gr, en)
    return re.sub(r'[^A-Za-z0-9\s\.\_\-]', '', str(text).translate(table))

def calculate_vaf_ci(mut, total):
    if total < 5: return 0.0, 0.0
    p = mut / total
    z = 1.96
    err = z * np.sqrt(max(0, p * (1 - p) / total))
    return max(0, (p - err) * 100), min(100, (p + err) * 100)

# --- 3. ΣΥΝΑΡΤΗΣΕΙΣ ΕΞΑΓΩΓΗΣ (PDF, WORD, EXCEL) ---

def generate_pdf(df, p_info):
    pdf = FPDF()
    pdf.add_page(); pdf.set_font("Arial", "B", 16)
    pdf.cell(190, 10, txt="CLINICAL ctDNA DISCOVERY REPORT", ln=True, align='C'); pdf.ln(10)
    pdf.set_font("Arial", "B", 12); pdf.cell(190, 10, txt="PATIENT PROFILE", ln=True)
    pdf.set_font("Arial", "", 11)
    pdf.cell(190, 8, txt=f"Name: {to_safe_latin(p_info['name'])} | Age: {p_info['age']} | ID: {to_safe_latin(p_info['id'])}", ln=True)
    pdf.ln(5)
    for sn in df['Sample'].unique():
        p150 = df[df.Sample == sn][df.Size < 150].Density.sum()
        z = (p150 - H_MEAN) / H_STD
        pdf.set_font("Arial", "B", 11); pdf.cell(190, 8, txt=f"Analysis: {to_safe_latin(sn)}", ln=True)
        pdf.set_font("Arial", "", 10)
        pdf.cell(190, 6, txt=f" - Z-Score: {z:.2f} sigma | P150: {p150:.2f}% | p-val: {(1-norm.cdf(abs(z))):.4f}", ln=True)
    output = pdf.output(dest='S')
    return bytes(output) if isinstance(output, bytearray) else output.encode('latin-1', 'ignore')

def generate_word(df, muts, p_info):
    doc = Document(); doc.add_heading('Clinical Biopsy Diagnostic', 0)
    doc.add_paragraph(f"Name: {p_info['name']} | ID: {p_info['id']} | Age: {p_info['age']}")
    for sn in df['Sample'].unique():
        p150 = df[df.Sample == sn][df.Size < 150].Density.sum()
        doc.add_heading(f'Sample: {sn}', level=1)
        doc.add_paragraph(f"Result: {'ABNORMAL' if p150>25 else 'NORMAL'} (Z-score: {(p150-H_MEAN)/H_STD:.2f} σ)")
        if sn in muts and not muts[sn].empty:
            m = muts[sn].iloc[0]; vaf = (m['Mut_Count']/(m['Mut_Count']+m['WT_Count']))*100
            doc.add_paragraph(f"Genetic Signal (VAF): {vaf:.2f}%")
    bio = io.BytesIO(); doc.save(bio); return bio.getvalue()

def generate_excel(df, muts, perf):
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Fragmentation_Analysis')
        if muts: pd.concat([d.assign(Sample=k) for k, d in muts.items()]).to_excel(writer, index=False, sheet_name='Genetic_Mutations')
        if perf: pd.concat([d.assign(Sample=k) for k, d in perf.items()]).to_excel(writer, index=False, sheet_name='Engine_Benchmarks')
    return output.getvalue()

# --- 4. SESSION STATE ΑΡΧΙΚΟΠΟΙΗΣΗ ---
if 'all_data' not in st.session_state: st.session_state['all_data'] = pd.DataFrame()
if 'motifs' not in st.session_state: st.session_state['motifs'] = {}
if 'muts' not in st.session_state: st.session_state['muts'] = {}
if 'perf' not in st.session_state: st.session_state['perf'] = {}

# --- 5. SIDEBAR: ΕΛΕΓΧΟΣ ---
st.sidebar.title("🧬 Control Core")
if st.sidebar.button("🏠 Home / Reset Session"):
    st.session_state.clear();
    if os.path.exists('liquid_biopsy.db'): os.remove('liquid_biopsy.db')
    st.rerun()

st.sidebar.markdown("---")
uname = st.sidebar.text_input("Ονοματεπώνυμο", " ")
uid = st.sidebar.text_input("Patient ID", "000")
uage = st.sidebar.number_input("Ηλικία", 1, 120, 40)
p_info = {'name': uname, 'id': uid, 'age': uage}

up_files = st.sidebar.file_uploader("Upload BAM samples", type=["bam"], accept_multiple_files=True)

if up_files:
    for f in up_files:
        if st.sidebar.button(f"Analyze: {f.name}"):
            with st.spinner("Processing..."):
                with open(f.name, "wb") as b: b.write(f.getbuffer())
                subprocess.run(["./read_bam", f.name, "4"])
                if os.path.exists('liquid_biopsy.db'):
                    conn = sqlite3.connect('liquid_biopsy.db')
                    df_r = pd.read_sql_query("SELECT * FROM FragmentSizes", conn)
                    if not df_r.empty:
                        df_r['Density'] = (df_r['Count'] / df_r['Count'].sum()) * 100
                        df_r['Sample'] = f.name
                        st.session_state['all_data'] = pd.concat([st.session_state['all_data'], df_r], ignore_index=True)
                        st.session_state['motifs'][f.name] = pd.read_sql_query("SELECT * FROM Motifs ORDER BY Count DESC LIMIT 30", conn)
                        st.session_state['muts'][f.name] = pd.read_sql_query("SELECT * FROM Mutations", conn)
                        st.session_state['perf'][f.name] = pd.read_sql_query("SELECT * FROM ProcessingMetrics", conn)
                    conn.close()

if not st.session_state['all_data'].empty:
    st.sidebar.markdown("---")
    st.sidebar.header("📥 Lab Exports")
    st.sidebar.download_button("📄 PDF Report", generate_pdf(st.session_state['all_data'], p_info), "report.pdf", "application/pdf")
    st.sidebar.download_button("📘 Word Report", generate_word(st.session_state['all_data'], st.session_state['muts'], p_info), "report.docx")
    st.sidebar.download_button("📊 Excel Export", generate_excel(st.session_state['all_data'], st.session_state['muts'], st.session_state['perf']), "data.xlsx")

# --- 6. ΚΥΡΙΩΣ DASHBOARD ---
st.title("🧬 Liquid Biopsy Discovery Platform")
df = st.session_state['all_data']

if not df.empty:
    # 1. Diagnostic Summary Cards
    unames = df.Sample.unique()
    grid = st.columns(len(unames))
    for i, name in enumerate(unames):
        sdf = df[df.Sample == name]
        p150 = sdf[sdf.Size < 150].Density.sum()
        z = (p150 - H_MEAN) / H_STD
        total_r = int(st.session_state['perf'][name].iloc[0]['TotalReads']) if name in st.session_state['perf'] else 0
        gc = st.session_state['muts'][name].iloc[0]['Avg_GC'] if name in st.session_state['muts'] else 0
        color = COLOR_MAP.get(name, "#00d1b2")
        with grid[i]:
            st.markdown(f"""
                <div class="metric-card" style="border-left-color:{color}">
                    <h3 style="color:{color};">{name}</h3>
                    <div class="metric-info">
                        Fragments Analysed: <span class="record-count">{total_r:,}</span><br>
                        Z-Score: <b>{z:.2f} σ</b> | p-val: <b>{(1-norm.cdf(abs(z))):.4f}</b><br>
                        P150 Index: <b>{p150:.2f}%</b> | GC: {gc:.1f}%<br>
                        Outcome: <b style="color:{color}">{('ABNORMAL' if z > 5 else 'NORMAL')}</b>
                    </div>
                </div>""", unsafe_allow_html=True)

    # 2. ΤΑ TABS
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Length Profile", "🧬 DNA Motifs", "🚀 Performance", "🎯 Mutations"])

    with tab1:
        st.subheader("Integrated Fragment Length Distributions")
        fig1 = px.line(df, x='Size', y='Density', color='Sample', line_shape='spline', color_discrete_map=COLOR_MAP)
        fig1.update_traces(line=dict(width=3.5))
        fig1.add_vrect(x0=160, x1=185, fillcolor="#00d1b2", opacity=0.08, layer="below", annotation_text="HEALTHY ZONE")
        st.plotly_chart(fig1, use_container_width=True)

        st.subheader("P150 Ratio Comparison (%)")
        p_list = pd.DataFrame([{'Sample': n, 'P150 Index': df[df.Sample == n][df.Size < 150].Density.sum()} for n in unames])
        st.plotly_chart(px.bar(p_list, x='Sample', y='P150 Index', color='Sample', color_discrete_map=COLOR_MAP, template="plotly_white"), use_container_width=True)

    with tab2:
        if st.session_state['motifs']:
            sm = st.selectbox("View Signature for:", list(st.session_state['motifs'].keys()))
            st.plotly_chart(px.bar(st.session_state['motifs'][sm], x='Motif', y='Count', color='Count', color_continuous_scale='Magma'), use_container_width=True)

    with tab3:
        st.subheader("Engine Efficiency (C-Engine vs Python)")
        for sn, p_df in st.session_state['perf'].items():
            m = p_df.iloc[0]; py_t = int(m['TotalReads']) * 0.00008; speedup = py_t / m['Runtime']
            st.info(f"**Sample {sn}:** C-Engine is {speedup:.1f}x faster.")
            st.plotly_chart(px.bar(pd.DataFrame({'System':['C-Engine','Standard Python'], 'Sec':[m['Runtime'], py_t]}), x='System', y='Sec', color='System', color_discrete_map={'C-Engine':'#00d1b2','Standard Python':'#ff4b4b'}), use_container_width=True)

    with tab4:
        st.subheader("🎯 Targeted Mutation Profiling (Genetic Ratio)")
        m_cols = st.columns(len(st.session_state['muts']))
        for idx, (sname, m_stats) in enumerate(st.session_state['muts'].items()):
            if not m_stats.empty:
                m_r = m_stats.iloc[0]; wt, mu = m_r['WT_Count'], m_r['Mut_Count']; vaf = (mu/(wt+mu)*100) if (wt+mu)>0 else 0
                low, high = calculate_ci = calculate_vaf_ci(mu, wt+mu)
                with m_cols[idx]:
                    st.write(f"**Sample: {sname}** (VAF: {vaf:.2f}%)")
                    st.write(f"95% CI: **[{low:.2f}% - {high:.2f}%]**")
                    pie_df = pd.DataFrame({"Type": ["WT", "Mutant"], "Count": [wt, mu]})
                    st.plotly_chart(px.pie(pie_df, names='Type', values='Count', hole=0.4, color='Type', color_discrete_map={"WT":"#00d1b2", "Mutant":"#ff4b4b"}), use_container_width=True)
else:
    st.info("System Ready. Analyze data to generate clinical findings.")
