import csv

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

class Solution:
    def custom_alignment(self):
        name = "pritdesai"  # Replace with your concatenated first and last name
        alphabet = "abcdefghijklmnopqrstuvwxyz"
        semi_matches = set(name)

        S = {}
        for char1 in alphabet:
            S[char1] = {}
            for char2 in alphabet:
                if char1 == char2:
                    S[char1][char2] = 2  # Match
                elif char2 in semi_matches:
                    S[char1][char2] = 1  # Semi-match
                else:
                    S[char1][char2] = -1  # Mismatch

        return S, name

    def local_alignment(self, sequence_A: str, sequence_B: str, substitution: dict, gap: int) -> [tuple]:
        m, n = len(sequence_A), len(sequence_B)
        D = [[0] * (n + 1) for _ in range(m + 1)]

        max_score = 0
        max_pos = (0, 0)

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = D[i - 1][j - 1] + substitution[sequence_A[i - 1]][sequence_B[j - 1]]
                delete = D[i - 1][j] + gap
                insert = D[i][j - 1] + gap
                D[i][j] = max(0, match, delete, insert)

                if D[i][j] > max_score:
                    max_score = D[i][j]
                    max_pos = (i, j)

        align_A, align_B = "", ""
        i, j = max_pos

        while D[i][j] > 0:
            if D[i][j] == D[i - 1][j - 1] + substitution[sequence_A[i - 1]][sequence_B[j - 1]]:
                align_A = sequence_A[i - 1] + align_A
                align_B = sequence_B[j - 1] + align_B
                i -= 1
                j -= 1
            elif D[i][j] == D[i - 1][j] + gap:
                align_A = sequence_A[i - 1] + align_A
                align_B = "-" + align_B
                i -= 1
            else:
                align_A = "-" + align_A
                align_B = sequence_B[j - 1] + align_B
                j -= 1

        return [(align_A, align_B)]


# Read the CSV file
data = []
with open('1002170533.csv', 'r') as file:
    csv_reader = csv.reader(file)
    next(csv_reader)  # Skip header
    for row in csv_reader:
        data.append(row)

# Perform Fisher's exact test
p_values = []
effective_p_value = 5e-8
bonferroni_corrected_p_value = effective_p_value / len(data)

results = []
for row in data:
    snp, case_c, case_t, control_c, control_t = row
    _, p_value = stats.fisher_exact([[int(case_c), int(case_t)], [int(control_c), int(control_t)]])
    p_values.append(p_value)

    is_significant_original = p_value < effective_p_value
    is_significant_corrected = p_value < bonferroni_corrected_p_value

    results.append([snp, p_value, is_significant_original, is_significant_corrected])

# Save results to CSV
with open('results.csv', 'w', newline='') as file:
    csv_writer = csv.writer(file)
    csv_writer.writerow(['SNP', 'p-value', 'Significant (Original)', 'Significant (Corrected)'])
    csv_writer.writerows(results)


# Create an instance of the Solution class
solution = Solution()

# Generate custom substitution matrix and name
S, name = solution.custom_alignment()

# Perform local alignment
result = solution.local_alignment(name, "thequickbrownfoxjumpsoverthelazydog", S, -2)


# Manhattan Plot
plt.figure(figsize=(12, 6))
plot_x = range(1, len(p_values) + 1)
plot_y = [-np.log10(p) for p in p_values]

plt.axhline(y=-np.log10(effective_p_value), color='r', linestyle='--', label=f'Original Threshold (5e-8)')
plt.axhline(y=-np.log10(bonferroni_corrected_p_value), color='b', linestyle='--', label=f'Bonferroni Correction ({bonferroni_corrected_p_value:.2e})')
plt.scatter(plot_x, plot_y, s=10, alpha=0.6, edgecolors='k')
plt.xlabel("SNP Index")
plt.ylabel("-log10(p-values)")
plt.title("Manhattan Plot for GWAS")
plt.legend()
plt.tight_layout()
plt.savefig('manhattan_plot.png')
plt.close()


# Print the result
print("Alignment result:", result)

# Pretty print the substitution matrix S
S_matrix = [['-'] + list(S.keys())] + [[key] + [S[key][k] for k in S[key]] for key in S]


# Save S matrix to file
np.savetxt('1002170533_S.txt', S_matrix, fmt='%s', delimiter='\t')

# Pretty print the alignment matrix D
sequence_A = name
sequence_B = "thequickbrownfoxjumpsoverthelazydog"
m, n = len(sequence_A), len(sequence_B)
D = [[0] * (n + 1) for _ in range(m + 1)]

for i in range(1, m + 1):
    for j in range(1, n + 1):
        match = D[i - 1][j - 1] + S[sequence_A[i - 1]][sequence_B[j - 1]]
        delete = D[i - 1][j] - 2
        insert = D[i][j - 1] - 2
        D[i][j] = max(0, match, delete, insert)

D_matrix = [['-'] + list(sequence_B)] + [[sequence_A[i - 1]] + D[i][1:] for i in range(1, m + 1)]


# Save D matrix to file
np.savetxt('1002170533_D.txt', D_matrix, fmt='%s', delimiter='\t')
