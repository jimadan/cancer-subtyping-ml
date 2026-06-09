import pandas as pd
import numpy as np
import gseapy as gp
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from pathlib import Path
from analysis.pathway_analysis import assign_pathway_group


def differential_expression(expr, labels, cluster):

    g1 = expr[labels == cluster]
    g2 = expr[labels != cluster]

    out = []
    eps = 1e-9

    for gene in expr.columns:

        try:
            _, p = ttest_ind(
                g1[gene],
                g2[gene],
                equal_var=False,
                nan_policy="omit"
            )

            logfc = np.log2(
                (g1[gene].mean() + eps) /
                (g2[gene].mean() + eps)
            )

            out.append({
                "gene": gene,
                "logFC": logfc,
                "pval": p
            })

        except Exception:
            continue

    df = pd.DataFrame(out)

    if df.empty:
        return None

    df["padj"] = multipletests(
        df["pval"],
        method="fdr_bh"
    )[1]

    df["score"] = (
        df["logFC"]
        * -np.log10(df["padj"] + 1e-300)
    )

    return df.sort_values("score", ascending=False)


def run_gsea_pipeline(expr, clusters):

    Path("results/gsea").mkdir(parents=True, exist_ok=True)

    labels = clusters["cluster"]
    results_all = []

    gene_sets = [
        "MSigDB_Hallmark_2020",
        "KEGG_2021_Human",
        "Reactome_2022",
        "GO_Biological_Process_2023"
    ]

    for cl in sorted(labels.unique()):

        print(f"[GSEA] cluster {cl}")

        de = differential_expression(expr, labels, cl)
        if de is None:
            continue

        ranking = (
            de[["gene", "score"]]
            .assign(gene=lambda x: x["gene"].astype(str).str.upper())
            .drop_duplicates("gene")
            .set_index("gene")["score"]
        )

        for gs in gene_sets:

            try:
                pre = gp.prerank(
                    rnk=ranking,
                    gene_sets=gs,
                    permutation_num=1000,
                    min_size=5,
                    max_size=2000,
                    seed=42,
                    outdir=f"results/gsea/{gs}/cluster_{cl}",
                    verbose=False
                )

                res = pre.res2d
                if res is None or res.empty:
                    continue

                res = res.loc[res["FDR q-val"] < 0.25].copy()

                if res.empty:
                    continue

                res["cluster"] = cl
                res["gene_set"] = gs
                res["pathway_group"] = res["Term"].apply(assign_pathway_group)

                results_all.append(res)

            except Exception as e:
                print(f"[ERROR] cluster {cl} {gs}: {e}")

    if not results_all:
        print("[WARNING] No significant GSEA results found")
        final = pd.DataFrame(
            columns=[
                "Name",
                "Term",
                "ES",
                "NES",
                "NOM p-val",
                "FDR q-val",
                "FWER p-val",
                "Tag %",
                "Gene %",
                "Lead_genes",
                "cluster",
                "gene_set",
                "pathway_group",
            ]
        )
        final.to_csv("results/gsea/all_gsea_results.csv", index=False)
        return final

    final = pd.concat(results_all, ignore_index=True)
    final.to_csv("results/gsea/all_gsea_results.csv", index=False)

    return final
