"""Interactive dashboard for exploring DESeq2 results.  Run: streamlit run app.py"""
import numpy as np, pandas as pd, plotly.express as px, streamlit as st

st.set_page_config(page_title="RNA-seq Explorer", layout="wide")
st.title("RNA-seq Differential Expression Explorer")

@st.cache_data
def load(path):
    return pd.read_csv(path, index_col=0)

up = st.sidebar.file_uploader("Upload DESeq2 results CSV (optional)", type="csv")
try:
    res = pd.read_csv(up, index_col=0) if up else load("deseq2_results.csv")
except FileNotFoundError:
    st.error("Run `python rnaseq_pipeline.py` first, or upload a results CSV."); st.stop()

padj_cut = st.sidebar.slider("Adjusted p-value cutoff", 0.001, 0.2, 0.05, 0.001)
lfc_cut = st.sidebar.slider("|log2FC| cutoff", 0.0, 4.0, 1.0, 0.1)

d = res.dropna(subset=["padj"]).copy()
d["nlog10padj"] = -np.log10(d["padj"].clip(lower=1e-300))
d["status"] = np.where((d["padj"] < padj_cut) & (d["log2FoldChange"] > lfc_cut), "Up",
              np.where((d["padj"] < padj_cut) & (d["log2FoldChange"] < -lfc_cut), "Down", "NS"))

c1, c2, c3 = st.columns(3)
c1.metric("Genes tested", len(d)); c2.metric("Upregulated", (d.status == "Up").sum())
c3.metric("Downregulated", (d.status == "Down").sum())

tab1, tab2, tab3 = st.tabs(["Volcano plot", "MA plot", "Results table"])
cmap = {"Up": "#C44E52", "Down": "#4C72B0", "NS": "#BBBBBB"}
with tab1:
    d2 = d.reset_index().rename(columns={d.index.name or "index": "gene"})
    fig = px.scatter(d2, x="log2FoldChange", y="nlog10padj", color="status",
                     color_discrete_map=cmap, hover_name="gene", opacity=.7)
    st.plotly_chart(fig, width="stretch")
with tab2:
    fig = px.scatter(d2, x="baseMean", y="log2FoldChange", color="status",
                     color_discrete_map=cmap, hover_name="gene", log_x=True, opacity=.7)
    st.plotly_chart(fig, width="stretch")
with tab3:
    q = st.text_input("Search gene")
    t = d[d.status != "NS"]
    if q: t = d[d.index.str.contains(q, case=False)]
    st.dataframe(t.drop(columns="nlog10padj"), width="stretch")
    st.download_button("Download table", t.to_csv().encode(), "filtered_results.csv")
