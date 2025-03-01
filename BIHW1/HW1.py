import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import csv

# Constants
effective_p_value = 5 * (10 ** -8)
num_tests = 1000  # Number of SNPs based on the dataset description
bonferroni_corrected_p_value = effective_p_value / num_tests

# Read SNP Data
with open("1002170533.csv", "r") as mydata:
    dataByLine = csv.reader(mydata, delimiter=";")
    snps = []
    next(dataByLine)  # Skip header
    for line in dataByLine:
        r = line[0].split(",")
        snp_name = r[0]
        snp_counts = list(map(int, r[1:5]))
        snps.append([snp_name] + snp_counts)

# Output file
with open('results.csv', 'w', newline='') as outputFile:
    outputHook = csv.writer(outputFile)
    outputHook.writerow(["SNP", "p_value", "Significant (Original Threshold)", "Significant (Bonferroni Correction)"])

    plot_x = []
    plot_y = []

    significant_original = 0
    significant_corrected = 0

    for idx, snp in enumerate(snps, start=1):
        snp_name, c_case, t_case, c_control, t_control = snp
        table = np.array([[c_case, t_case], [c_control, t_control]])
        _, p_value = fisher_exact(table, alternative="two-sided")

        is_significant_original = p_value <= effective_p_value
        is_significant_corrected = p_value <= bonferroni_corrected_p_value

        if is_significant_original:
            significant_original += 1
        if is_significant_corrected:
            significant_corrected += 1

        outputHook.writerow([snp_name, p_value, is_significant_original, is_significant_corrected])

        plot_x.append(idx)
        plot_y.append(-np.log10(p_value))

# Manhattan Plot
plt.figure(figsize=(12, 6))
plt.axhline(y=-np.log10(effective_p_value), color='r', linestyle='--', label=f'Original Threshold (5e-8)')
plt.axhline(y=-np.log10(bonferroni_corrected_p_value), color='b', linestyle='--', label=f'Bonferroni Correction ({bonferroni_corrected_p_value:.2e})')
plt.scatter(plot_x, plot_y, s=10, alpha=0.6, edgecolors='k')
plt.xlabel("SNP Index")
plt.ylabel("-log10(p-values)")
plt.title("Manhattan Plot for GWAS")
plt.legend()
plt.tight_layout()
plt.show()

print(f"Significant SNPs (Original Threshold): {significant_original}")
print(f"Significant SNPs (Bonferroni Correction): {significant_corrected}")