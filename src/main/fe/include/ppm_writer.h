#ifndef PPM_WRITER_H
#define PPM_WRITER_H

#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

/**
 * @brief Classe pour générer des images PPM in-situ
 * 
 * Permet de convertir des slices 2D de champs de pression en images PPM
 * avec différentes palettes de couleurs (colormaps).
 */
class PPMWriter {
public:
  enum class Colormap {
    JET,       // Bleu -> Cyan -> Vert -> Jaune -> Rouge
    VIRIDIS,   // Violet -> Bleu -> Vert -> Jaune
    COOLWARM,  // Bleu -> Blanc -> Rouge
    GRAYSCALE  // Noir -> Blanc
  };

  enum class SlicePlane {
    XY,  // Plan horizontal (z constant)
    XZ,  // Plan vertical (y constant)
    YZ   // Plan vertical (x constant)
  };

  /**
   * @brief Écrit une image PPM binaire (P6) à partir d'une slice 2D
   * 
   * @param filename Nom du fichier de sortie
   * @param data Données 2D (ligne majeure)
   * @param width Largeur de l'image
   * @param height Hauteur de l'image
   * @param colormap Palette de couleurs à utiliser
   * @param vmin Valeur minimale pour normalisation (auto si NaN)
   * @param vmax Valeur maximale pour normalisation (auto si NaN)
   */
  static void writeSlice(const std::string& filename,
                         const std::vector<float>& data,
                         int width, int height,
                         Colormap colormap = Colormap::VIRIDIS,
                         float vmin = NAN, float vmax = NAN) {
    
    // Calcul automatique des bornes si non spécifiées
    if (std::isnan(vmin) || std::isnan(vmax)) {
      auto minmax = std::minmax_element(data.begin(), data.end());
      if (std::isnan(vmin)) vmin = *minmax.first;
      if (std::isnan(vmax)) vmax = *minmax.second;
    }

    // Éviter division par zéro
    if (vmax == vmin) vmax = vmin + 1.0f;

    std::ofstream out(filename, std::ios::binary);
    if (!out) {
      throw std::runtime_error("Cannot open PPM file: " + filename);
    }

    // En-tête PPM P6 (binaire)
    out << "P6\n" << width << " " << height << "\n255\n";

    // Écriture des pixels
    for (int i = 0; i < width * height; ++i) {
      float value = data[i];
      
      // Normalisation [vmin, vmax] -> [0, 1]
      float normalized = (value - vmin) / (vmax - vmin);
      normalized = std::max(0.0f, std::min(1.0f, normalized));

      // Application de la colormap
      auto rgb = applyColormap(normalized, colormap);
      
      out.put(static_cast<char>(rgb[0]));
      out.put(static_cast<char>(rgb[1]));
      out.put(static_cast<char>(rgb[2]));
    }

    out.close();
  }

  /**
   * @brief Extrait une slice 2D d'un champ 3D
   * 
   * @param field3D Champ 3D complet
   * @param nx, ny, nz Dimensions du champ
   * @param plane Plan de coupe
   * @param sliceIndex Position de la coupe (index dans la dimension perpendiculaire)
   * @return std::vector<float> Données de la slice 2D
   */
  static std::vector<float> extractSlice(const std::vector<float>& field3D,
                                         int nx, int ny, int nz,
                                         SlicePlane plane, int sliceIndex) {
    std::vector<float> slice;

    switch (plane) {
      case SlicePlane::XY: {
        // Plan z = sliceIndex
        slice.reserve(nx * ny);
        for (int j = 0; j < ny; ++j) {
          for (int i = 0; i < nx; ++i) {
            int idx = i + j * nx + sliceIndex * nx * ny;
            slice.push_back(field3D[idx]);
          }
        }
        break;
      }
      case SlicePlane::XZ: {
        // Plan y = sliceIndex
        slice.reserve(nx * nz);
        for (int k = 0; k < nz; ++k) {
          for (int i = 0; i < nx; ++i) {
            int idx = i + sliceIndex * nx + k * nx * ny;
            slice.push_back(field3D[idx]);
          }
        }
        break;
      }
      case SlicePlane::YZ: {
        // Plan x = sliceIndex
        slice.reserve(ny * nz);
        for (int k = 0; k < nz; ++k) {
          for (int j = 0; j < ny; ++j) {
            int idx = sliceIndex + j * nx + k * nx * ny;
            slice.push_back(field3D[idx]);
          }
        }
        break;
      }
    }

    return slice;
  }

private:
  /**
   * @brief Applique une colormap à une valeur normalisée [0,1]
   * @return RGB values in [0, 255]
   */
  static std::array<uint8_t, 3> applyColormap(float t, Colormap cmap) {
    // Clamp to [0, 1]
    t = std::max(0.0f, std::min(1.0f, t));

    switch (cmap) {
      case Colormap::JET:
        return jet(t);
      case Colormap::VIRIDIS:
        return viridis(t);
      case Colormap::COOLWARM:
        return coolwarm(t);
      case Colormap::GRAYSCALE:
        return grayscale(t);
      default:
        return {0, 0, 0};
    }
  }

  // Implémentations des colormaps

  static std::array<uint8_t, 3> jet(float t) {
    float r = std::max(0.0f, std::min(1.0f, 1.5f - 4.0f * std::abs(t - 0.75f)));
    float g = std::max(0.0f, std::min(1.0f, 1.5f - 4.0f * std::abs(t - 0.5f)));
    float b = std::max(0.0f, std::min(1.0f, 1.5f - 4.0f * std::abs(t - 0.25f)));
    return {static_cast<uint8_t>(r * 255),
            static_cast<uint8_t>(g * 255),
            static_cast<uint8_t>(b * 255)};
  }

  static std::array<uint8_t, 3> viridis(float t) {
    // Approximation simplifiée de Viridis
    float r = 0.267004f + t * (0.993248f - 0.267004f);
    float g = 0.004874f + t * (0.906157f - 0.004874f);
    float b = 0.329415f + t * (0.143936f - 0.329415f);
    
    if (t < 0.5f) {
      r = 0.267004f + 2.0f * t * (0.127568f - 0.267004f);
      g = 0.004874f + 2.0f * t * (0.566949f - 0.004874f);
      b = 0.329415f + 2.0f * t * (0.550556f - 0.329415f);
    } else {
      float t2 = (t - 0.5f) * 2.0f;
      r = 0.127568f + t2 * (0.993248f - 0.127568f);
      g = 0.566949f + t2 * (0.906157f - 0.566949f);
      b = 0.550556f + t2 * (0.143936f - 0.550556f);
    }
    
    return {static_cast<uint8_t>(r * 255),
            static_cast<uint8_t>(g * 255),
            static_cast<uint8_t>(b * 255)};
  }

  static std::array<uint8_t, 3> coolwarm(float t) {
    // Bleu -> Blanc -> Rouge
    float r, g, b;
    if (t < 0.5f) {
      // Bleu -> Blanc
      float t2 = t * 2.0f;
      r = 0.23f + t2 * (1.0f - 0.23f);
      g = 0.30f + t2 * (1.0f - 0.30f);
      b = 0.75f + t2 * (1.0f - 0.75f);
    } else {
      // Blanc -> Rouge
      float t2 = (t - 0.5f) * 2.0f;
      r = 1.0f;
      g = 1.0f - t2 * (1.0f - 0.34f);
      b = 1.0f - t2 * (1.0f - 0.29f);
    }
    return {static_cast<uint8_t>(r * 255),
            static_cast<uint8_t>(g * 255),
            static_cast<uint8_t>(b * 255)};
  }

  static std::array<uint8_t, 3> grayscale(float t) {
    uint8_t gray = static_cast<uint8_t>(t * 255);
    return {gray, gray, gray};
  }
};

#endif // PPM_WRITER_H