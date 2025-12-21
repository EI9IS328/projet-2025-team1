#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import numpy as np
import pandas as pd
import argparse
import argcomplete
import os
import glob
import re
import matplotlib.pyplot as plt

class SnapshotViewer:
    def __init__(self, files):
        self.files = files
        self.index = 0
        self.fig, self.axes = plt.subplots(2, 2, figsize=(12, 10))
        self.axes = self.axes.flatten() # Pour itérer facilement sur les 4 subplots

        # Connecter les touches
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)

        # Charger la première image
        self.update_plot()

    def get_snapshot_time(self, filename):
        match = re.search(r'snapshot_(\d+)\.csv', filename)
        return int(match.group(1)) if match else -1

    def load_data(self, filepath):
        print(f"--- [{self.index+1}/{len(self.files)}] Chargement : {os.path.basename(filepath)} ---")
        try:
            # Lecture optimisée
            df = pd.read_csv(filepath, delimiter=',', usecols=['x', 'y', 'z', 'p'])
        except Exception as e:
            print(f"Erreur : {e}")
            return None, None

        if df.empty: return None, None

        # Calcul des tranches (exactement comme ton script original)
        z_min, z_max = df['z'].min(), df['z'].max()
        z_range = z_max - z_min
        percentages = [0.25, 0.50, 0.75, 1.0]
        z_targets = [z_min + p * z_range for p in percentages[:-1]] + [z_max]

        unique_z = df['z'].unique()
        z_slices_actual = []
        for target in z_targets:
            idx = (np.abs(unique_z - target)).argmin()
            z_slices_actual.append(unique_z[idx])

        z_slices_actual = sorted(list(set(z_slices_actual)))

        df_slices = df[df['z'].isin(z_slices_actual)].copy()
        return df_slices, z_slices_actual

    def update_plot(self):
        filepath = self.files[self.index]
        df_slices, z_vals = self.load_data(filepath)

        if df_slices is None:
            return

        time_ms = self.get_snapshot_time(filepath)
        self.fig.suptitle(f"Snapshot: {os.path.basename(filepath)} ({time_ms}ms)\nFlèches G/D pour naviguer", fontsize=16)

        # Nettoyer les axes
        for ax in self.axes:
            ax.clear()
            ax.set_aspect('equal')

        # Trouver min/max global pour la couleur (pour que ce soit cohérent sur les 4 graphs)
        vmin, vmax = df_slices['p'].min(), df_slices['p'].max()

        # Dessiner chaque tranche
        for i, z_val in enumerate(z_vals):
            if i >= 4: break # Max 4 subplots

            ax = self.axes[i]
            data = df_slices[df_slices['z'] == z_val]

            # Utilisation de scatter avec marqueurs carrés pour imiter geom_raster
            sc = ax.scatter(data['x'], data['y'], c=data['p'], cmap='viridis',
                            marker='s', s=10, vmin=vmin, vmax=vmax) # s=taille des points

            ax.set_title(f"Z = {z_val:.2f}")
            ax.set_xlabel("X")
            ax.set_ylabel("Y")

        # Cacher les axes vides si moins de 4 tranches
        for j in range(len(z_vals), 4):
            self.axes[j].axis('off')

        # Rafraîchir
        self.fig.canvas.draw()

    def on_key(self, event):
        if event.key == 'right':
            self.index = (self.index + 1) % len(self.files)
            self.update_plot()
        elif event.key == 'left':
            self.index = (self.index - 1) % len(self.files)
            self.update_plot()
        elif event.key in ['escape', 'q']:
            plt.close('all')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True)
    # Arguments legacy ignorés
    parser.add_argument('--output', type=str, required=False)
    parser.add_argument('--benchmark', type=str, required=False)
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # Trouver les fichiers
    input_path = os.path.abspath(args.input)
    input_dir = os.path.dirname(input_path)
    files = glob.glob(os.path.join(input_dir, "snapshot_*.csv"))

    # Tri numérique
    files.sort(key=lambda f: int(re.search(r'snapshot_(\d+)\.csv', f).group(1)) if re.search(r'snapshot_(\d+)\.csv', f) else -1)

    if not files:
        print("Aucun fichier trouvé.")
        return

    # Lancer le viewer
    try:
        start_index = files.index(input_path)
    except ValueError:
        start_index = 0

    print("=== Viewer Matplotlib ===")
    print("Chargement...")

    viewer = SnapshotViewer(files)
    viewer.index = start_index
    # Un dernier update pour se caler sur le bon index de départ
    if start_index != 0:
        viewer.update_plot()

    plt.show()

if __name__ == "__main__":
    main()