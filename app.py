"""Interactive dashboard for exploring DESeq2 results.  Run: streamlit run app.py"""
import os
import numpy as np, pandas as pd, plotly.express as px, streamlit as st

st.set_page_config(page_title="RNA-seq Explorer", layout="wide")
st.title("RNA-seq Differential Expression Explorer")

def find(name):
    """Look for a file in the app folder or in a results/ subfolder."""
    for p in (name, os.path.join("results", name)):
        if os.path.exists(p):
            return p
    return None

@st.cache_data
def load(path):
    return pd.read_csv(path, index_col=0)

up = st.sidebar.file_uploader("Upload DESeq2 results CSV (optional)", type="csv")
res_path = find("deseq2_results.csv")
if up:
    res = pd.read_csv(up, index_col=0)
elif res_path:
    res = load(res_path)
else:
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

tab1, tab2, tab3, tab4 = st.tabs(["Volcano plot", "MA plot", "Heatmap", "Results table"])
cmap = {"Up": "#C44E52", "Down": "#4C72B0", "NS": "#BBBBBB"}
d2 = d.reset_index().rename(columns={d.index.name or "index": "gene"})

with tab1:
    fig = px.scatter(d2, x="log2FoldChange", y="nlog10padj", color="status",
                     color_discrete_map=cmap, hover_name="gene", opacity=.7)
    st.plotly_chart(fig, width="stretch")
with tab2:
    fig = px.scatter(d2, x="baseMean", y="log2FoldChange", color="status",
                     color_discrete_map=cmap, hover_name="gene", log_x=True, opacity=.7)
    st.plotly_chart(fig, width="stretch")
with tab3:
    norm_path = find("normalized_counts.csv")
    if not norm_path:
        st.info("Heatmap needs `normalized_counts.csv` in the same folder as the app.")
    else:
        norm = load(norm_path)
        n = st.slider("Number of top genes", 10, 100, 30, 5)
        top = d.sort_values("padj").head(n).index.intersection(norm.index)
        z = np.log2(norm.loc[top] + 1)
        z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1) + 1e-9, axis=0)
        z = z[sorted(z.columns)]
        z = z.loc[d.loc[top].sort_values("log2FoldChange").index]
        fig = px.imshow(z, aspect="auto", color_continuous_scale="RdBu_r",
                        color_continuous_midpoint=0, labels=dict(color="z-score"))
        fig.update_layout(height=max(450, 18 * len(z)))
        st.plotly_chart(fig, width="stretch")
        st.caption("Top genes by adjusted p-value, row z-scored log2 normalized counts.")
with tab4:
    q = st.text_input("Search gene")
    t = d[d.status != "NS"]
    if q: t = d[d.index.str.contains(q, case=False)]
    st.dataframe(t.drop(columns="nlog10padj"), width="stretch")
    st.download_button("Download table", t.to_csv().encode(), "filtered_results.csv")
