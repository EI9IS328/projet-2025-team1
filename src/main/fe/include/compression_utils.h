#ifndef COMPRESSION_UTILS_H
#define COMPRESSION_UTILS_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include <fstream>
#include <iostream>

/**
 * @brief Utilitaires pour la compression de données avec/sans perte
 */
class CompressionUtils {
public:
  
  /**
   * @brief Métadonnées pour la quantification
   */
  struct QuantizationMetadata {
    float vmin;
    float vmax;
    int nbits;  // 8, 16, 32
  };

  /**
   * @brief Quantifie un tableau de float en entiers sur N bits
   * 
   * @param data Données en float
   * @param nbits Nombre de bits (8, 16 ou 32)
   * @param global_min Min pour normalisation (auto si NaN)
   * @param global_max Max pour normalisation (auto si NaN)
   * @return std::pair<std::vector<uint32_t>, QuantizationMetadata>
   */
  static std::pair<std::vector<uint32_t>, QuantizationMetadata>
  quantize(const std::vector<float>& data, int nbits,
           float global_min = NAN, float global_max = NAN) {
    
    if (nbits != 8 && nbits != 16 && nbits != 32) {
      throw std::invalid_argument("nbits must be 8, 16, or 32");
    }

    // Calcul des bornes
    float vmin = global_min;
    float vmax = global_max;
    
    if (std::isnan(vmin) || std::isnan(vmax)) {
      vmin = std::numeric_limits<float>::max();
      vmax = std::numeric_limits<float>::lowest();
      for (float v : data) {
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
      }
    }

    // Éviter division par zéro
    if (vmax == vmin) vmax = vmin + 1.0f;

    // Nombre de niveaux
    uint32_t max_val = (1u << nbits) - 1;

    std::vector<uint32_t> quantized;
    quantized.reserve(data.size());

    for (float v : data) {
      // Normalisation [vmin, vmax] -> [0, 1]
      float normalized = (v - vmin) / (vmax - vmin);
      normalized = std::max(0.0f, std::min(1.0f, normalized));
      
      // Quantification
      uint32_t q = static_cast<uint32_t>(normalized * max_val + 0.5f);
      quantized.push_back(q);
    }

    QuantizationMetadata meta{vmin, vmax, nbits};
    return {quantized, meta};
  }

  /**
   * @brief Dequantifie les données
   */
  static std::vector<float> dequantize(
      const std::vector<uint32_t>& quantized,
      const QuantizationMetadata& meta) {
    
    uint32_t max_val = (1u << meta.nbits) - 1;
    
    std::vector<float> data;
    data.reserve(quantized.size());

    for (uint32_t q : quantized) {
      float normalized = static_cast<float>(q) / max_val;
      float v = normalized * (meta.vmax - meta.vmin) + meta.vmin;
      data.push_back(v);
    }

    return data;
  }

  /**
   * @brief Calcule l'erreur RMSE entre données originales et reconstruites
   */
  static float computeRMSE(const std::vector<float>& original,
                           const std::vector<float>& reconstructed) {
    if (original.size() != reconstructed.size()) {
      throw std::invalid_argument("Size mismatch in RMSE computation");
    }

    double sum_sq_error = 0.0;
    for (size_t i = 0; i < original.size(); ++i) {
      double diff = original[i] - reconstructed[i];
      sum_sq_error += diff * diff;
    }

    return std::sqrt(sum_sq_error / original.size());
  }

  /**
   * @brief Calcule l'erreur relative maximale
   */
  static float computeMaxRelativeError(const std::vector<float>& original,
                                       const std::vector<float>& reconstructed) {
    if (original.size() != reconstructed.size()) {
      throw std::invalid_argument("Size mismatch in error computation");
    }

    float max_rel_err = 0.0f;
    for (size_t i = 0; i < original.size(); ++i) {
      if (std::abs(original[i]) > 1e-10f) {
        float rel_err = std::abs(original[i] - reconstructed[i]) / std::abs(original[i]);
        max_rel_err = std::max(max_rel_err, rel_err);
      }
    }

    return max_rel_err;
  }

  /**
   * @brief Run-Length Encoding (RLE)
   * 
   * Encode les séquences de valeurs identiques.
   * Format: [valeur, count, valeur, count, ...]
   */
  static std::vector<uint32_t> rleEncode(const std::vector<uint32_t>& data) {
    if (data.empty()) return {};

    std::vector<uint32_t> encoded;
    encoded.reserve(data.size() / 2);  // Estimation

    uint32_t current_value = data[0];
    uint32_t count = 1;

    for (size_t i = 1; i < data.size(); ++i) {
      if (data[i] == current_value && count < 0xFFFFFFFF) {
        count++;
      } else {
        encoded.push_back(current_value);
        encoded.push_back(count);
        current_value = data[i];
        count = 1;
      }
    }

    // Dernière séquence
    encoded.push_back(current_value);
    encoded.push_back(count);

    return encoded;
  }

  /**
   * @brief Décode RLE
   */
  static std::vector<uint32_t> rleDecode(const std::vector<uint32_t>& encoded,
                                         size_t expected_size = 0) {
    std::vector<uint32_t> decoded;
    if (expected_size > 0) {
      decoded.reserve(expected_size);
    }

    for (size_t i = 0; i < encoded.size(); i += 2) {
      uint32_t value = encoded[i];
      uint32_t count = encoded[i + 1];
      
      for (uint32_t j = 0; j < count; ++j) {
        decoded.push_back(value);
      }
    }

    return decoded;
  }

  /**
   * @brief Calcule le taux de compression RLE
   */
  static float rleCompressionRatio(size_t original_size, size_t rle_size) {
    return static_cast<float>(original_size) / rle_size;
  }

  /**
   * @brief Sauvegarde données quantifiées avec métadonnées
   */
  static void saveQuantized(const std::string& filename,
                            const std::vector<uint32_t>& data,
                            const QuantizationMetadata& meta) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
      throw std::runtime_error("Cannot open file: " + filename);
    }

    // Écriture métadonnées
    out.write(reinterpret_cast<const char*>(&meta.vmin), sizeof(float));
    out.write(reinterpret_cast<const char*>(&meta.vmax), sizeof(float));
    out.write(reinterpret_cast<const char*>(&meta.nbits), sizeof(int));

    // Taille des données
    size_t size = data.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    // Données quantifiées
    if (meta.nbits == 8) {
      for (uint32_t v : data) {
        uint8_t val = static_cast<uint8_t>(v);
        out.write(reinterpret_cast<const char*>(&val), sizeof(uint8_t));
      }
    } else if (meta.nbits == 16) {
      for (uint32_t v : data) {
        uint16_t val = static_cast<uint16_t>(v);
        out.write(reinterpret_cast<const char*>(&val), sizeof(uint16_t));
      }
    } else {  // 32 bits
      out.write(reinterpret_cast<const char*>(data.data()), 
                data.size() * sizeof(uint32_t));
    }

    out.close();
  }

  /**
   * @brief Sauvegarde données avec RLE
   */
  static void saveRLE(const std::string& filename,
                     const std::vector<uint32_t>& rle_data,
                     const QuantizationMetadata& meta,
                     size_t original_size) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
      throw std::runtime_error("Cannot open file: " + filename);
    }

    // Métadonnées
    out.write(reinterpret_cast<const char*>(&meta.vmin), sizeof(float));
    out.write(reinterpret_cast<const char*>(&meta.vmax), sizeof(float));
    out.write(reinterpret_cast<const char*>(&meta.nbits), sizeof(int));
    out.write(reinterpret_cast<const char*>(&original_size), sizeof(size_t));

    // Taille RLE
    size_t rle_size = rle_data.size();
    out.write(reinterpret_cast<const char*>(&rle_size), sizeof(size_t));

    // Données RLE
    out.write(reinterpret_cast<const char*>(rle_data.data()), 
              rle_data.size() * sizeof(uint32_t));

    out.close();
  }

  /**
   * @brief Statistiques de compression
   */
  struct CompressionStats {
    size_t original_bytes;
    size_t compressed_bytes;
    float compression_ratio;
    float rmse;
    float max_rel_error;
    double compression_time_ms;
    double decompression_time_ms;
  };

  /**
   * @brief Affiche les statistiques
   */
  static void printStats(const CompressionStats& stats, 
                        const std::string& method_name) {
    std::cout << "\n=== Compression Stats: " << method_name << " ===" << std::endl;
    std::cout << "Original size:     " << stats.original_bytes / 1024.0 << " KB" << std::endl;
    std::cout << "Compressed size:   " << stats.compressed_bytes / 1024.0 << " KB" << std::endl;
    std::cout << "Compression ratio: " << stats.compression_ratio << "x" << std::endl;
    std::cout << "RMSE:              " << stats.rmse << std::endl;
    std::cout << "Max relative err:  " << stats.max_rel_error * 100 << "%" << std::endl;
    std::cout << "Compression time:  " << stats.compression_time_ms << " ms" << std::endl;
    std::cout << "Decompression time:" << stats.decompression_time_ms << " ms" << std::endl;
  }
};

#endif // COMPRESSION_UTILS_H