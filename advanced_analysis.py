
from pathlib import Path
import re

import pandas as pd
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
DATA = BASE / "human_genetics_gwas_associations copy.xlsx"
OUT = BASE / "results"
OUT.mkdir(exist_ok=True)

df = pd.read_excel(DATA)

# 1. Parse beta text without treating unavailable values as zero.
beta_text = df["beta"].astype(str)
df["beta_numeric"] = pd.to_numeric(
    beta_text.str.extract(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")[0],
    errors="coerce"
)
df["beta_direction"] = beta_text.str.extract(
    r"\b(increase|decrease)\b", flags=re.IGNORECASE
)[0].str.lower()

# A decrease is represented as a negative signed effect.
df.loc[df["beta_direction"] == "decrease", "beta_numeric"] *= -1

print("=== BETA PARSING ===")
print("Numeric beta values:", df["beta_numeric"].notna().sum())
print(df[["beta", "beta_numeric", "beta_direction"]].head(12).to_string(index=False))

# 2. Reconstruct p-values using the reported mantissa and exponent.
# Log10 representation remains usable even when p-values underflow to zero.
df["log10_p"] = (
    df["pvalue_mantissa"].astype(float).map(lambda x: __import__("math").log10(x))
    + df["pvalue_exponent"].astype(float)
)

print("\n=== EXTREME P-VALUES ===")
print(
    df.loc[df["p_value"].eq(0),
           ["queried_gene", "p_value_formatted", "log10_p"]]
    .to_string(index=False)
)

# 3. Summarise significance labels by queried gene.
sig = pd.crosstab(df["queried_gene"], df["genome_wide_significant"])
sig["total"] = sig.sum(axis=1)
sig["significant_percent"] = (
    sig.get("Yes", 0) / sig["total"] * 100
)
sig.to_csv(OUT / "significance_by_gene.csv")

print("\n=== SIGNIFICANCE BY GENE ===")
print(sig.to_string())

ax = sig["significant_percent"].sort_values().plot(kind="barh")
ax.set_xlim(0, 100)
ax.set_xlabel("Records labelled genome-wide significant (%)")
ax.set_ylabel("Queried gene")
ax.set_title("Significance labels in this collected dataset")
plt.tight_layout()
plt.savefig(OUT / "significance_by_gene.png", dpi=200)
plt.close()

# 4. Check confidence intervals against odds ratios.
or_value = pd.to_numeric(df["or_value"], errors="coerce")
lower = pd.to_numeric(df["ci_lower"], errors="coerce")
upper = pd.to_numeric(df["ci_upper"], errors="coerce")

has_all = or_value.notna() & lower.notna() & upper.notna()
df["or_ci_contains_estimate"] = (
    (lower <= or_value) & (or_value <= upper)
).where(has_all)

print("\n=== ODDS-RATIO / CI CHECK ===")
print("Rows with OR and both limits:", int(has_all.sum()))
print(
    "Rows where CI contains OR:",
    int(df.loc[has_all, "or_ci_contains_estimate"].sum())
)

print("\nExamples requiring source verification:")

invalid_ci = has_all & df["or_ci_contains_estimate"].eq(False)

if invalid_ci.any():
    print(
        df.loc[
            invalid_ci,
            [
                "queried_gene",
                "reported_trait",
                "or_value",
                "ci_lower",
                "ci_upper",
            ],
        ].head(15).to_string(index=False)
    )
else:
    print(
        "No inconsistencies found among records with complete "
        "odds ratios and confidence intervals."
    )



# 5. Export the cleaned, derived fields for further inspection.
df.to_csv(OUT / "derived_gwas_data.csv", index=False)
print(f"\nOutputs saved in {OUT}")
