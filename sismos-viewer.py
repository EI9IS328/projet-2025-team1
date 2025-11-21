import argparse
import pandas as pd
import matplotlib.pyplot as plt
import argcomplete

# ===============================
# Parser les arguments de la ligne de commande
# ===============================
parser = argparse.ArgumentParser(description="Génère un graphique à partir d'un fichier CSV de sismos.")
parser.add_argument('--input', type=str, required=True, help="Chemin vers le fichier CSV d'entrée")
parser.add_argument('--output', type=str, required=False, help="Chemin pour sauvegarder le graphique (ex: sortie.png)")

argcomplete.autocomplete(parser)
args = parser.parse_args()

input_file = args.input
output_file = args.output

# ===============================
# Lecture du CSV
# ===============================
data = pd.read_csv(input_file)

# Vérification des colonnes
if 'timestep' not in data.columns or 'pressure' not in data.columns:
    raise ValueError("Le fichier doit contenir les colonnes 'timestep' et 'pressure'.")

# ===============================
# Création du graphique (simple ligne)
# ===============================
plt.figure(figsize=(10,6))
plt.plot(data['timestep'], data['pressure'], linestyle='-', color='blue')  # <- uniquement ligne
plt.xlabel("Timestep")
plt.ylabel("Pressure")
plt.title("Evolution de la pression en fonction du temps")
plt.grid(True)
plt.tight_layout()

# ===============================
# Sauvegarde ou affichage
# ===============================
if output_file:
    plt.savefig(output_file)
    print(f"Graphique sauvegardé dans {output_file}")
else:
    plt.show()


