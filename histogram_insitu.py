import argparse
import pandas as pd
import matplotlib.pyplot as plt
import argcomplete
import time
import csv 
import os
# ===============================
# Parser les arguments de la ligne de commande
# ===============================
parser = argparse.ArgumentParser(description="Génère un graphique à partir d'un fichier CSV de sismos.")
parser.add_argument('--input', type=str, required=True, help="Chemin vers le fichier CSV d'entrée")

argcomplete.autocomplete(parser)
args = parser.parse_args()

input_file = args.input
output_file = args.output
bench = {}

# ===============================
# Lecture du CSV
# ===============================
start = time.time()
data = pd.read_csv(input_file)
bench["read_csv"] = time.time() - start


if "p" not in data.columns:
    raise ValueError("La colonne 'p' n'existe pas dans le CSV.")

p_values = data["p"].astype(float)

pmin = p_values.min()
pmax = p_values.max()

NBINS = 10

if pmax == pmin:
    print(f"Toutes les valeurs sont identiques : p = {pmin}")
    print(f"Bin 0 [{pmin} , {pmax}] : {len(p_values)} points")
else:
    bin_width = (pmax - pmin) / NBINS

    # Création manuelle des bins comme en C++
    hist = [0] * NBINS

    for p in p_values:
        bin_index = int((p - pmin) / bin_width)
        if bin_index == NBINS:  # Cas limite, identique à ton code C++
            bin_index = NBINS - 1
        hist[bin_index] += 1

    # Affichage du résultat
    for i in range(NBINS):
        bmin = pmin + i * bin_width
        bmax = bmin + bin_width
        print(f"Bin {i} [{bmin} , {bmax}] : {hist[i]} points")
