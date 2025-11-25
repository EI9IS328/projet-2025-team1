#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import numpy as np
import pandas as pd
from plotnine import *
import argparse
import argcomplete
import time
import csv
import os

#todo avec une touche passer direct au snapshot suivant/précédent

def main():
    # Configuration des arguments de la ligne de commande
    parser = argparse.ArgumentParser(description='Génère un graphique à partir d\'un fichier CSV de snapshot.')
    parser.add_argument('--input', type=str, required=True, help='Chemin vers le fichier CSV d\'entrée')
    parser.add_argument('--output', type=str, required=True, help='Chemin pour sauvegarder le graphique')
    parser.add_argument('--benchmark', type=str, required=False, help='Chemin du CSV benchmark à générer')

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    bench = {}


    # Lecture du fichier CSV avec gestion des erreurs
    start = time.time()
    try:
        df = pd.read_csv(args.input, delimiter=',')
    except FileNotFoundError:
        print(f"Erreur : Le fichier d'entrée '{args.input}' n'a pas été trouvé.")
        return
    except Exception as e:
        print(f"Une erreur est survenue lors de la lecture du fichier CSV : {e}")
        return
    bench["read_csv"] = time.time() - start

    # Vérification que les colonnes nécessaires (x, y, z, p) existent
    required_columns = ['x', 'y', 'z', 'p']
    if not all(col in df.columns for col in required_columns):
        print(f"Erreur : Le fichier CSV doit contenir les colonnes suivantes : {required_columns}")
        return

    ##############################################################################
    start = time.time()
    # 1. Déterminer la plage min/max de l'axe Z
    z_min = df['z'].min()
    z_max = df['z'].max()
    z_range = z_max - z_min
    bench["z_calcul"] = time.time() - start

    print(f"Plage de Z détectée : min={z_min:.4f}, max={z_max:.4f}")

    # 2. Définir les pourcentages et calculer les coordonnées Z cibles
    percentages = [0.25, 0.50, 0.75, 1.0]
    # Pour la dernière coupe, on prend le max directement pour éviter les problèmes de précision
    z_targets = [z_min + p * z_range for p in percentages[:-1]] + [z_max]

    start = time.time()
    # 3. Pour chaque cible, trouver la valeur Z existante la plus proche dans les données
    unique_z_values = df['z'].unique()
    z_slices_actual = []
    for target in z_targets:
        closest_z_index = np.abs(unique_z_values - target).argmin()
        z_slices_actual.append(unique_z_values[closest_z_index])

    # S'assurer que les valeurs sont uniques et triées
    z_slices_actual = sorted(list(set(z_slices_actual)))
    bench["z_slice"] = time.time() - start

    print(f"Coordonnées Z cibles pour {percentages}: {[f'{z:.4f}' for z in z_targets]}")
    print(f"Tranches Z les plus proches trouvées dans les données : {[f'{z:.4f}' for z in z_slices_actual]}")

    start = time.time()
    # 4. Filtrer le dataframe pour ne conserver que les données de ces tranches Z
    df_slices = df[df['z'].isin(z_slices_actual)].copy()

    # Conversion de la colonne 'z' en type Catégoriel pour un meilleur affichage avec facet_wrap
    df_slices['z'] = pd.Categorical(df_slices['z'])
    bench["df_filter"] = time.time() - start

    if df_slices.empty:
        print("Avertissement : Aucune donnée trouvée pour les tranches Z calculées. Vérifiez vos données.")
        return

    start = time.time()
    # 5. (Optionnel) Calculer le min et le max de la pression sur les coupes pour information
    p_min_slices = df_slices['p'].min()
    p_max_slices = df_slices['p'].max()
    bench["pressure_range_s"] = time.time() - start

    print(f"Échelle de pression détectée sur les coupes (min/max) : [{p_min_slices:.4f}, {p_max_slices:.4f}]")

    start = time.time()
    # 6. Créer le graphique avec plotnine
    p = (
            ggplot(df_slices, aes(x='x', y='y', fill='p'))
            + geom_raster()
            + facet_wrap('~z', labeller='label_both', ncol=2)
            + scale_fill_cmap('viridis')
            + coord_fixed()
            + theme_minimal()
            + labs(
        title="Coupes de Pression à différentes profondeurs (Z)",
        subtitle=f"Coupes aux niveaux Z les plus proches de 25%, 50%, 75% et 100%",
        x="Axe X",
        y="Axe Y",
        fill="Pression"
    )
    )
    bench["plot_creation"] = time.time() - start 
    # Sauvegarder le graphique avec gestion des erreurs

    start = time.time()
    try:
        p.save(args.output, width=12, height=10, units='in', dpi=300, verbose=False)
        print(f"Graphique sauvegardé dans {args.output}")
    except Exception as e:
        print(f"Une erreur est survenue lors de la sauvegarde du graphique : {e}")
    bench["plot_save"] = time.time() - start

    # Sauvegarder le benchmark dans un fichier csv
    if args.benchmark:
        write_header = not os.path.exists(args.benchmark)

        with open(args.benchmark, 'a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=bench.keys())
            if write_header:
                writer.writeheader()
            writer.writerow(bench)

        print(f"Benchmark sauvegardé dans {args.benchmark}")
if __name__ == "__main__":
    main()