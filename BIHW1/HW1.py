import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# Load dataset
UTA_ID = "1002170533"
filename = f"{UTA_ID}.csv"
df = pd.read_csv(filename)

# Ensure correct column names
expected_columns = {"SNP", "Case_C", "Case_T", "Control_C", "Control_T"}
if not expected_columns.issubset(df.columns):
    print(f"Error: Missing expected columns. Found columns: {df.columns}")
    exit()

# Check for missing values and fill with zeros
df.fillna(0, inplace=True)

# Perform Fisher's Exact Test for each SNP
significance_threshold = 5e-8
num_tests = len(df)
bonferroni_threshold = significance_threshold / num_tests

results = []
for index, row in df.iterrows():
    table = [[row["Case_C"], row["Case_T"]], [row["Control_C"], row["Control_T"]]]
    odds_ratio, p = stats.fisher_exact(table, alternative='two-sided')
    significant = p < significance_threshold
    bonferroni_significant = p < bonferroni_threshold
    results.append([row["SNP"], p, significant, bonferroni_significant])

# Convert results to DataFrame and save
results_df = pd.DataFrame(results, columns=["SNP", "P-value", "Significant", "Bonferroni-Significant"])
results_df.to_csv("results.csv", index=False)
print("Results saved in results.csv")

# Generate Manhattan Plot
plt.figure(figsize=(10, 5))
plt.scatter(range(len(results_df)), -np.log10(results_df["P-value"]), color="blue", alpha=0.7)
plt.axhline(y=-np.log10(significance_threshold), color="red", linestyle="dashed", label="Significance (5e-8)")
plt.axhline(y=-np.log10(bonferroni_threshold), color="green", linestyle="dashed", label="Bonferroni Threshold")
plt.xlabel("SNP Index")
plt.ylabel("-log10(P-value)")
plt.title("Manhattan Plot")
plt.legend()
plt.savefig("manhattan_plot.png")
plt.show()
print("Manhattan plot saved.")

# Substitution Matrix for Custom Alignment
characters = "abcdefghijklmnopqrstuvwxyz "
name = "pritdesai"
pangram = "thequickbrownfoxjumpsoverthelazydog"

S = np.zeros((len(characters), len(characters)), dtype=int)
char_to_index = {char: idx for idx, char in enumerate(characters)}

for i, c1 in enumerate(characters):
    for j, c2 in enumerate(characters):
        if c1 == c2:
            S[i][j] = 2  # Match
        elif c1 in name and c2 in name:
            S[i][j] = 1  # Semi-Match
        else:
            S[i][j] = -1  # Mismatch

# Save S matrix
with open(f"{UTA_ID}_S.txt", "w") as f:
    f.write("  " + " ".join(characters) + "\n")
    for i, char in enumerate(characters):
        f.write(char + " " + " ".join(map(str, S[i])) + "\n")

# Smith-Waterman Algorithm
def smith_waterman(seq1, seq2, S, gap_penalty=-1):
    m, n = len(seq1), len(seq2)
    D = np.zeros((m+1, n+1), dtype=int)
    max_score, max_pos = 0, (0, 0)

    for i in range(1, m+1):
        for j in range(1, n+1):
            match = D[i-1][j-1] + S[char_to_index[seq1[i-1]], char_to_index[seq2[j-1]]]
            delete = D[i-1][j] + gap_penalty
            insert = D[i][j-1] + gap_penalty
            D[i][j] = max(0, match, delete, insert)
            if D[i][j] > max_score:
                max_score, max_pos = D[i][j], (i, j)

    return D

D = smith_waterman(name, pangram, S)
with open(f"{UTA_ID}_D.txt", "w") as f:
    for row in D:
        f.write(" ".join(map(str, row)) + "\n")
print("Alignment matrix saved.")

# Needleman-Wunsch Algorithm
def needleman_wunsch(seq1, seq2, S, gap_penalty=-1):
    m, n = len(seq1), len(seq2)
    D = np.zeros((m+1, n+1), dtype=int)
    for i in range(m+1):
        D[i][0] = i * gap_penalty
    for j in range(n+1):
        D[0][j] = j * gap_penalty

    for i in range(1, m+1):
        for j in range(1, n+1):
            match = D[i-1][j-1] + S[char_to_index[seq1[i-1]], char_to_index[seq2[j-1]]]
            delete = D[i-1][j] + gap_penalty
            insert = D[i][j-1] + gap_penalty
            D[i][j] = max(match, delete, insert)

    return D

D_global = needleman_wunsch(name, pangram, S)
with open(f"{UTA_ID}_NW.txt", "w") as f:
    for row in D_global:
        f.write(" ".join(map(str, row)) + "\n")
print("Global alignment matrix saved.")
