import pandas as pd

FILE = "human_genetics_gwas_associations copy.xlsx"
df = pd.read_excel(FILE)

print("=== BETA EXAMPLES ===")
print(df["beta"].astype(str).value_counts().head(15))

print("\n=== P-VALUE ZERO CHECK ===")
print("Reported zeros:", (df["p_value"] == 0).sum())
print(df.loc[df["p_value"] == 0,
             ["queried_gene", "p_value", "p_value_formatted",
              "pvalue_mantissa", "pvalue_exponent"]].head(10))

print("\n=== SIGNIFICANCE BY GENE ===")
print(pd.crosstab(
    df["queried_gene"],
    df["genome_wide_significant"],
    margins=True
))

print("\n=== ODDS-RATIO EXAMPLES ===")
print(df.loc[df["or_value"].notna(),
             ["queried_gene", "reported_trait", "or_value",
              "or_per_copy_num", "ci_lower", "ci_upper"]].head(10))

print("\n=== EFFECT-SIZE RAW VALUES ===")
for col in ["beta", "or_value", "or_per_copy_num",
            "ci_lower", "ci_upper", "risk_frequency"]:
    print(f"\n{col}:")
    print(df[col].drop_duplicates().head(10).to_list())
