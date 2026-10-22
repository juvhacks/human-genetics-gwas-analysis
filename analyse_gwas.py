
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "human_genetics_gwas_associations copy.xlsx"
OUTPUT_DIR = BASE_DIR / "results"
OUTPUT_DIR.mkdir(exist_ok=True)

# Load data
df = pd.read_excel(DATA_FILE)

# Remove columns that contain no useful data
df = df.dropna(axis=1, how="all")

print(f"Dataset: {df.shape[0]} rows, {df.shape[1]} usable columns")
print(f"Unique association IDs: {df['association_id'].nunique()}")

# 1. Records per queried gene
gene_counts = df["queried_gene"].value_counts(dropna=False)
print("\nRecords per queried gene:")
print(gene_counts.to_string())

gene_counts.rename_axis("queried_gene").reset_index(
    name="record_count"
).to_csv(OUTPUT_DIR / "records_per_gene.csv", index=False)

ax = gene_counts.sort_values().plot(kind="barh")
ax.set_title("GWAS records per queried gene")
ax.set_xlabel("Number of records")
ax.set_ylabel("Queried gene")
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "records_per_gene.png", dpi=200)
plt.close()

# 2. Genome-wide significance labels
print("\nGenome-wide significance labels:")
print(df["genome_wide_significant"].value_counts(dropna=False).to_string())

# 3. P-value summary
pvals = pd.to_numeric(df["p_value"], errors="coerce")
valid_pvals = pvals[pvals.between(0, 1)]

print(f"\nValid p-values: {len(valid_pvals)} / {len(df)}")
print(valid_pvals.describe().to_string())

# 4. Effect-size availability
for col in ["beta", "or_value", "or_per_copy_num", "ci_lower", "ci_upper"]:
    numeric = pd.to_numeric(df[col], errors="coerce")
    print(f"{col}: {numeric.notna().sum()} numeric values")

# 5. Missingness report
missing = df.isna().sum().sort_values(ascending=False)
missing.rename("missing_count").to_csv(OUTPUT_DIR / "missing_values.csv")

print(f"\nResults saved to: {OUTPUT_DIR}")
