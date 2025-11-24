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
parser.add_argument('--output', type=str, required=False, help="Chemin pour sauvegarder le graphique (ex: sortie.png)")
parser.add_argument('--benchmark', type=str, required=False, help="Chemin vers le fichiers CSV pour sauvegarder les benchmarks")

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

# Vérification des colonnes
if 'timestep' not in data.columns or 'pressure' not in data.columns:
    raise ValueError("Le fichier doit contenir les colonnes 'timestep' et 'pressure'.")

# ===============================
# Création du graphique (simple ligne)
# ===============================
start = time.time()
plt.figure(figsize=(10,6))
plt.plot(data['timestep'], data['pressure'], linestyle='-', color='blue')  # <- uniquement ligne
plt.xlabel("Timestep")
plt.ylabel("Pressure")
plt.title("Evolution de la pression en fonction du temps")
plt.grid(True)
plt.tight_layout()

bench["plot_creation"] = time.time() - start
# ===============================
# Sauvegarde ou affichage
# ===============================
if output_file:
    start = time.time()
    plt.savefig(output_file)
    bench["save_plot"] = time.time() - start
    print(f"Graphique sauvegardé dans {output_file}")
else:
    start = time.time()
    plt.show()
    bench["show_plot"] = time.time() - start

# ===============================
# Sauvegarde du benchmark csv
# ===============================
if args.benchmark:
    write_header = not os.path.exists(args.benchmark)
    with open(args.benchmark, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=bench.keys())
        if write_header:
            writer.writeheader()
        writer.writerow(bench)

    print(f"Benchmark écrit dans {args.benchmark}")


