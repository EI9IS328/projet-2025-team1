#pragma once

#include <cxxopts.hpp>
#include <stdexcept>
#include <string>

class SemProxyOptions
{
 public:
  // Defaults
  int order = 2;
  int ex = 50, ey = 50, ez = 50;
  float lx = 2000.f, ly = 2000.f, lz = 2000.f;
  float srcx = 1010.f, srcy = 1010.f, srcz = 1010.f;
  float rcvx = 1410.f, rcvy = 1010.f, rcvz = 1010.f;
  std::string implem = "makutu";  // makutu|shiva
  std::string method = "sem";     // sem|dg
  std::string mesh = "cartesian";
  float dt = 0.001;
  float timemax = 1.5;
  bool autodt = false;
  // sponge boundaries parameters
  float boundaries_size = 0;
  bool surface_sponge = false;
  float taper_delta = 0.015;
  // Boolean to tell if the model is charged on nodes or on element
  bool isModelOnNodes = false;
  bool isElastic = false;

  // snapshot
  bool saveSnapshot = false;
  string snapFolder = "";
  int snapInterval = 50;

  // sismo
  string sismoFile = "";
  string sismoFolder = "";

  // save perf
  string perfFile = "";

  // insitu
  string insituFolder = "";
  bool insituHistogram = false;
  int insituInterval = 50;

  // PPM 
  bool savePPM = false;
  string ppmFolder = "";
  int ppmInterval = 50;
  string ppmPlane = "xy";        // xy, xz, yz
  int ppmSliceIndex = 0;         // 0 = milieu automatique
  string ppmColormap = "viridis"; // viridis, jet, coolwarm, grayscale

  // Compression 
  bool useCompression = false;
  int compressionBits = 16;      // 8, 16, 32
  bool useRLE = false;
  string compressionStatsFile = "";

  void validate() const
  {
    if (order < 1) throw std::runtime_error("order must be >= 1");
    if (ex <= 0 || ey <= 0 || ez <= 0)
      throw std::runtime_error("ex/ey/ez must be > 0");
    if (lx <= 0 || ly <= 0 || lz <= 0)
      throw std::runtime_error("lx/ly/lz must be > 0");
    
    // Validation PPM
    if (ppmPlane != "xy" && ppmPlane != "xz" && ppmPlane != "yz")
      throw std::runtime_error("ppm-plane must be: xy, xz, or yz");
    
    if (ppmColormap != "viridis" && ppmColormap != "jet" && 
        ppmColormap != "coolwarm" && ppmColormap != "grayscale")
      throw std::runtime_error("ppm-colormap must be: viridis, jet, coolwarm, or grayscale");
    
    // Validation compression
    if (compressionBits != 8 && compressionBits != 16 && compressionBits != 32)
      throw std::runtime_error("compression-bits must be: 8, 16, or 32");
  }

  // Bind CLI flags to this instance (no --help here)
  static void bind_cli(cxxopts::Options& opts, SemProxyOptions& o)
  {
    opts.add_options()
        ("o,order", "Order of approximation", cxxopts::value<int>(o.order))
        ("ex", "Number of elements on X (Cartesian mesh)", cxxopts::value<int>(o.ex))
        ("ey", "Number of elements on Y (Cartesian mesh)", cxxopts::value<int>(o.ey))
        ("ez", "Number of elements on Z (Cartesian mesh)", cxxopts::value<int>(o.ez))
        ("lx", "Domain size X (Cartesian)", cxxopts::value<float>(o.lx))
        ("ly", "Domain size Y (Cartesian)", cxxopts::value<float>(o.ly))
        ("lz", "Domain size Z (Cartesian)", cxxopts::value<float>(o.lz))
        ("implem", "Implementation: makutu|shiva", cxxopts::value<std::string>(o.implem))
        ("method", "Method: sem|dg", cxxopts::value<std::string>(o.method))
        ("mesh", "Mesh: cartesian|ucartesian", cxxopts::value<std::string>(o.mesh))
        ("dt", "Time step selection in s (default = 0.001s)", cxxopts::value<float>(o.dt))
        ("timemax", "Duration of the simulation in s (default = 1.5s)", cxxopts::value<float>(o.timemax))
        ("auto-dt", "Select automatique dt via CFL equation.", cxxopts::value<bool>(o.autodt))
        ("boundaries-size", "Size of absorbing boundaries (meters)", cxxopts::value<float>(o.boundaries_size))
        ("sponge-surface", "Considere the surface's nodes as non sponge nodes", cxxopts::value<bool>(o.surface_sponge))
        ("taper-delta", "Taper delta for sponge boundaries value", cxxopts::value<float>(o.taper_delta))
        ("is-model-on-nodes", "Boolean to tell if the model is charged on nodes (true) or on element (false)", cxxopts::value<bool>(o.isModelOnNodes))
        ("is-elastic", "Elastic simulation", cxxopts::value<bool>(o.isElastic))

        // Snapshots
        ("save-snapshot", "Save snapshots during the simulation", cxxopts::value<bool>(o.saveSnapshot))
        ("snap-folder", "Folder to save snapshots (default = current folder)", cxxopts::value<string>(o.snapFolder))
        ("snap-interval", "Step interval between snapshots save (default = 50)", cxxopts::value<int>(o.snapInterval))

        // Sismogrammes
        ("sismo-file", "Sismo file for source definition", cxxopts::value<string>(o.sismoFile))
        ("sismo-folder", "Folder where sismo file is located", cxxopts::value<string>(o.sismoFolder))

        // Performance
        ("perf-file", "File to save performance", cxxopts::value<string>(o.perfFile))

        // In-situ histogram
        ("insitu-histogram", "Enable in-situ histogram processing", cxxopts::value<bool>(o.insituHistogram))
        ("insitu-interval", "Interval between in-situ processing", cxxopts::value<int>(o.insituInterval))
        ("insitu-folder", "Folder to save in-situ results", cxxopts::value<string>(o.insituFolder))

        // PPM Visualization 
        ("save-ppm", "Enable PPM visualization (in-situ image generation)", cxxopts::value<bool>(o.savePPM))
        ("ppm-folder", "Folder to save PPM images (default = current folder)", cxxopts::value<string>(o.ppmFolder))
        ("ppm-interval", "Step interval between PPM images (default = 50)", cxxopts::value<int>(o.ppmInterval))
        ("ppm-plane", "Slice plane: xy, xz, or yz (default = xy)", cxxopts::value<string>(o.ppmPlane))
        ("ppm-slice-index", "Position of the slice (0 = automatic center)", cxxopts::value<int>(o.ppmSliceIndex))
        ("ppm-colormap", "Colormap: viridis, jet, coolwarm, grayscale (default = viridis)", cxxopts::value<string>(o.ppmColormap))

        // Compression 
        ("use-compression", "Enable data compression for snapshots", cxxopts::value<bool>(o.useCompression))
        ("compression-bits", "Quantization bits: 8, 16, or 32 (default = 16)", cxxopts::value<int>(o.compressionBits))
        ("use-rle", "Enable Run-Length Encoding (RLE) compression", cxxopts::value<bool>(o.useRLE))
        ("compression-stats", "File to save compression statistics", cxxopts::value<string>(o.compressionStatsFile))
    ;
  }
};