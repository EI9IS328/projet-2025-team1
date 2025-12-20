//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "sem_proxy.h"

#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <sem_solver_acoustic.h>
#include <source_and_receiver_utils.h>
#include <stdio.h>

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

#include "compression_utils.h"
#include "ppm_writer.h"

using namespace SourceAndReceiverUtils;

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  const int order = opt.order;
  nb_elements_[0] = opt.ex;
  nb_elements_[1] = opt.ey;
  nb_elements_[2] = opt.ez;
  nb_nodes_[0] = opt.ex * order + 1;
  nb_nodes_[1] = opt.ey * order + 1;
  nb_nodes_[2] = opt.ez * order + 1;

  const float spongex = opt.boundaries_size;
  const float spongey = opt.boundaries_size;
  const float spongez = opt.boundaries_size;
  const std::array<float, 3> sponge_size = {spongex, spongey, spongez};
  src_coord_[0] = opt.srcx;
  src_coord_[1] = opt.srcy;
  src_coord_[2] = opt.srcz;

  domain_size_[0] = opt.lx;
  domain_size_[1] = opt.ly;
  domain_size_[2] = opt.lz;

  rcv_coord_[0] = opt.rcvx;
  rcv_coord_[1] = opt.rcvy;
  rcv_coord_[2] = opt.rcvz;

  bool isModelOnNodes = opt.isModelOnNodes;
  isElastic_ = opt.isElastic;
  cout << boolalpha;
  bool isElastic = isElastic_;

  // snapshots
  save_snapshot = opt.saveSnapshot;
  snapFolder = opt.snapFolder;
  snapInterval = opt.snapInterval;

  savePPM = opt.savePPM;
  ppmFolder = opt.ppmFolder;
  ppmInterval = opt.ppmInterval;
  ppmSliceIndex = opt.ppmSliceIndex;

  if (opt.ppmPlane == "xy")
    ppmPlane = PPMWriter::SlicePlane::XY;
  else if (opt.ppmPlane == "xz")
    ppmPlane = PPMWriter::SlicePlane::XZ;
  else if (opt.ppmPlane == "yz")
    ppmPlane = PPMWriter::SlicePlane::YZ;
  else
    ppmPlane = PPMWriter::SlicePlane::XY;  // default

  // Conversion string -> enum pour la colormap
  if (opt.ppmColormap == "viridis")
    ppmColormap = PPMWriter::Colormap::VIRIDIS;
  else if (opt.ppmColormap == "jet")
    ppmColormap = PPMWriter::Colormap::JET;
  else if (opt.ppmColormap == "coolwarm")
    ppmColormap = PPMWriter::Colormap::COOLWARM;
  else if (opt.ppmColormap == "grayscale")
    ppmColormap = PPMWriter::Colormap::GRAYSCALE;
  else
    ppmColormap = PPMWriter::Colormap::VIRIDIS;

  // sismos
  sismosFile = opt.sismoFile;
  sismosFolder = opt.sismoFolder;
  sismosPoints.clear();
  sismosNodeIndex.clear();

  // perf
  perfFile = opt.perfFile;

  // histogramme
  insituHistogram = opt.insituHistogram;
  insituInterval = opt.insituInterval;
  insituFolder = opt.insituFolders;

  // slice insitu
  sliceSnapshot = opt.sliceSnapshot;
  sliceInterval = opt.sliceInterval;
  axe = opt.axe;
  values = opt.value;
  sliceFolder = opt.sliceFolder;

  const SolverFactory::methodType methodType = getMethod(opt.method);
  const SolverFactory::implemType implemType = getImplem(opt.implem);
  const SolverFactory::meshType meshType = getMesh(opt.mesh);
  const SolverFactory::modelLocationType modelLocation =
      isModelOnNodes ? SolverFactory::modelLocationType::OnNodes
                     : SolverFactory::modelLocationType::OnElements;
  const SolverFactory::physicType physicType =
      SolverFactory::physicType::Acoustic;

  float lx = domain_size_[0];
  float ly = domain_size_[1];
  float lz = domain_size_[2];
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  if (meshType == SolverFactory::Struct)
  {
    switch (order)
    {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      default:
        throw std::runtime_error(
            "Order other than 1 2 3 is not supported (semproxy)");
    }
  }
  else if (meshType == SolverFactory::Unstruct)
  {
    model::CartesianParams<float, int> param(order, ex, ey, ez, lx, ly, lz,
                                             isModelOnNodes);
    model::CartesianUnstructBuilder<float, int> builder(param);
    m_mesh = builder.getModel();
  }
  else
  {
    throw std::runtime_error("Incorrect mesh type (SEMproxy ctor.)");
  }

  // time parameters
  if (opt.autodt)
  {
    float cfl_factor = (order == 2) ? 0.5 : 0.7;
    dt_ = find_cfl_dt(cfl_factor);
  }
  else
  {
    dt_ = opt.dt;
  }
  timemax_ = opt.timemax;
  num_sample_ = timemax_ / dt_;

  m_solver = SolverFactory::createSolver(methodType, implemType, meshType,
                                         modelLocation, physicType, order);
  m_solver->computeFEInit(*m_mesh, sponge_size, opt.surface_sponge,
                          opt.taper_delta);

  initFiniteElem();

  std::cout << "Number of node is " << m_mesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << m_mesh->getNumberOfElements()
            << std::endl;
  std::cout << "Launching the Method " << opt.method << ", the implementation "
            << opt.implem << " and the mesh is " << opt.mesh << std::endl;
  std::cout << "Model is on " << (isModelOnNodes ? "nodes" : "elements")
            << std::endl;
  std::cout << "Physics type is " << (isElastic ? "elastic" : "acoustic")
            << std::endl;
  std::cout << "Order of approximation will be " << order << std::endl;
  std::cout << "Time step is " << dt_ << "s" << std::endl;
  std::cout << "Simulated time is " << timemax_ << "s" << std::endl;
}
void SEMproxy::run()
{

  initSismoPoints();


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////  MAIN LOOP  //////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  SEMsolverDataAcoustic solverData(i1, i2, myRHSTerm, pnGlobal, rhsElement,
                                   rhsWeights);

  for (int indexTimeSample = 0; indexTimeSample < num_sample_;
       indexTimeSample++)
  {
    // ===== CALCUL =====
    startComputeTime = system_clock::now();
    m_solver->computeOneStep(dt_, indexTimeSample, solverData);
    totalComputeTime += system_clock::now() - startComputeTime;

    // ===== SORTIES =====
    startOutputTime = system_clock::now();

    // Output debug (tous les 50 pas)
    if (indexTimeSample % 50 == 0)
    {
      m_solver->outputSolutionValues(indexTimeSample, i1, rhsElement[0],
                                     pnGlobal, "pnGlobal");
    }

    // Save pressure at receiver
    const int order = m_mesh->getOrder();
    float varnp1 = 0.0;
    for (int i = 0; i < order + 1; i++)
    {
      for (int j = 0; j < order + 1; j++)
      {
        for (int k = 0; k < order + 1; k++)
        {
          int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
          int globalNodeOnElement =
              i + j * (order + 1) + k * (order + 1) * (order + 1);
          varnp1 +=
              pnGlobal(nodeIdx, i2) * rhsWeightsRcv(0, globalNodeOnElement);
        }
      }
    }

    pnAtReceiver(0, indexTimeSample) = varnp1;
    swap(i1, i2);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  MAIN LOOP AFTER COMPUTE one step  //////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int time_ms = static_cast<int>((indexTimeSample * dt_ * 1000.0));

    if (save_snapshot && indexTimeSample % snapInterval == 0)
      saveSnapshot(time_ms);

    if (savePPM && indexTimeSample % ppmInterval == 0)
      saveSnapshotPPM(time_ms);


    if (insituHistogram && indexTimeSample % insituInterval == 0)
      saveHistogramInsitu(time_ms);


    if (sliceSnapshot && indexTimeSample % insituInterval == 0)
      saveSliceSnapshot(time_ms, axe, values);

    if (!sismosFile.empty())
      saveSismoPoints(indexTimeSample);


    /////////////////////////////////////////////
    /////////////////////////////////////////////

    // Swap des indices temporels
    auto tmp = solverData.m_i1;
    solverData.m_i1 = solverData.m_i2;
    solverData.m_i2 = tmp;

    totalOutputTime += system_clock::now() - startOutputTime;
  }

  // ===== STATISTIQUES FINALES =====
  float kerneltime_ms = time_point_cast<microseconds>(totalComputeTime)
                            .time_since_epoch()
                            .count();
  float outputtime_ms =
      time_point_cast<microseconds>(totalOutputTime).time_since_epoch().count();

  cout << "------------------------------------------------ " << endl;
  cout << "\n---- Elapsed Kernel Time : " << kerneltime_ms / 1E6 << " seconds."
       << endl;
  cout << "---- Elapsed Output Time : " << outputtime_ms / 1E6 << " seconds."
       << endl;
  cout << "------------------------------------------------ " << endl;

  savePerf(kerneltime_ms, outputtime_ms);
}

void SEMproxy::savePerf(float kerneltime_ms, float outputtime_ms)
{
  if (!perfFile.empty())
  {
    std::ofstream out(perfFile, std::ios::app);
    if (!out.is_open())
    {
      std::cerr << "Error: could not open perf file " << perfFile << std::endl;
      return;
    }

    // if empty file, write header
    out.seekp(0, std::ios::end);
    if (out.tellp() == 0)
    {
      // kernel_time, output_time, total_time, nb_nodes, nb_step, nb_snapshot,
      // save_sismo
      out << "kernel_time_ms,output_time_ms,total_time_ms,nb_nodes,nb_steps,"
             "nb_snapshot,nb_sismo,nb_snapshot_ppm,nb_histo,nb_slice\n";
    }

    float total_time_ms = kerneltime_ms + outputtime_ms;
    int nb_nodes = m_mesh->getNumberOfNodes();

    int nb_snapshot = save_snapshot ? (num_sample_ / snapInterval) + 1 : 0;
    int nb_sismo = sismosPoints.size();
    int nb_snapshot_ppm = savePPM ? (num_sample_ / ppmInterval) + 1 : 0;
    int nb_histo = insituHistogram ? (num_sample_ / insituInterval) + 1 : 0;
    int nb_slice = sliceSnapshot ? (num_sample_ / sliceInterval) + 1 : 0;

    out << std::fixed << std::setprecision(3) << kerneltime_ms << ","
        << outputtime_ms << "," << total_time_ms << "," << nb_nodes << ","
        << num_sample_ << "," << nb_snapshot << "," << nb_sismo << ","
        << nb_snapshot_ppm << "," << nb_histo << "," << nb_slice << "\n";

    out.close();
  }
}

string SEMproxy::getSismoFileName(array<float, 3> point)
{
  std::string sismo_file = "sismo_" + std::to_string(point[0]) + "-" +
                           std::to_string(point[1]) + "-" +
                           std::to_string(point[2]) + ".csv";

  string filename;
  if (sismosFolder.empty())
  {
    filename = sismo_file;
  }
  else
  {
    filename = sismosFolder + "/" + sismo_file;
  }

  return filename;
}

void SEMproxy::initSismoPoints()
{
  if (!sismosFile.empty())
  {
    parsePointSismos(sismosFile);

    sismosNodeIndex.clear();

    for (size_t i = 0; i < sismosPoints.size(); i++)
    {
      const auto& s = sismosPoints[i];
      int index_node_sismos = findClosestNode(s[0], s[1], s[2]);

      sismosNodeIndex.push_back(index_node_sismos);

      string filename = getSismoFileName(s);

      std::ofstream out(filename);
      if (!out.is_open())
      {
        std::cerr << "Error: could not open sismo file " << filename
                  << std::endl;
        continue;
      }

      out << "timestep,pressure\n";
      out.close();
    }
  }
}

void SEMproxy::saveSismoPoints(int timestep)
{
  if (!sismosFile.empty())
  {
    for (size_t i = 0; i < sismosNodeIndex.size(); i++)
    {
      int node = sismosNodeIndex[i];
      const auto& s = sismosPoints[i];

      float pression = pnGlobal(node, i2);

      string filename = getSismoFileName(s);

      std::ofstream out(filename, std::ios::app);
      if (!out.is_open())
      {
        std::cerr << "Error: could not open sismo file " << filename
                  << std::endl;
        continue;
      }

      out << timestep << "," << pression << "\n";
      out.close();
    }
  }
}

void SEMproxy::saveSnapshot(int time_ms)
{
  // Construire le nom du fichier
  string filename;
  if (snapFolder.empty())
  {
    filename = "snapshot_" + std::to_string(time_ms) + ".csv";
  }
  else
  {
    filename = snapFolder + "/snapshot_" + std::to_string(time_ms) + ".csv";
  }
  printf("save in %s\n", filename.c_str());
  std::ofstream out(filename);

  if (!out.is_open())
  {
    std::cerr << "Error: could not open snapshot file " << filename
              << std::endl;
    return;
  }

  out << "x,y,z,p\n";

  int nbNodes = m_mesh->getNumberOfNodes();

  for (int node = 0; node < nbNodes; node++)
  {
    // Récupération des coordonnées du nœud
    float x = m_mesh->nodeCoord(node, 0);
    float y = m_mesh->nodeCoord(node, 1);
    float z = m_mesh->nodeCoord(node, 2);

    float p = pnGlobal(node, i2);

    out << x << "," << y << "," << z << "," << p << "\n";
  }

  out.close();
}

void SEMproxy::saveSliceSnapshot(int time_ms, int axe, float value)
{
  float eps = 1e-6;
  // Construire le nom du fichier
  string filename;
  if (sliceFolder.empty())
  {
    filename = "sliceSnapshot_" + std::to_string(time_ms) + ".csv";
  }
  else
  {
    filename =
        sliceFolder + "/sliceSnapshot_" + std::to_string(time_ms) + ".csv";
  }
  printf("save in %s\n", filename.c_str());
  std::ofstream out(filename);

  if (!out.is_open())
  {
    std::cerr << "Error: could not open snapshot file " << filename
              << std::endl;
    return;
  }

  out << "x,y,z,p\n";

  int nbNodes = m_mesh->getNumberOfNodes();
  for (int node = 0; node < nbNodes; node++)
  {
    // Récupération des coordonnées du nœud
    float coords[3];
    coords[0] = m_mesh->nodeCoord(node, 0);  // x
    coords[1] = m_mesh->nodeCoord(node, 1);  // y
    coords[2] = m_mesh->nodeCoord(node, 2);  // z

    float p = pnGlobal(node, i2);

    if (std::fabs(coords[axe] - value) < eps)
    {
      out << coords[0] << "," << coords[1] << "," << coords[2] << "," << p
          << "\n";
    }
  }

  out.close();
}

void SEMproxy::saveHistogramInsitu(int time_ms)
{
  string filename;
  if (insituFolder.empty())
  {
    filename = "histogram_" + std::to_string(time_ms) + ".csv";
  }
  else
  {
    filename =
        insituFolder + "/histogram_" + std::to_string(time_ms) + ".csv";
  }
  printf("save histogram in %s\n", filename.c_str());

  std::ofstream out(filename);
  if (!out.is_open())
  {
    std::cerr << "Error: could not open histogram file " << filename
              << std::endl;
    // Ne pas faire return ici, continuer la simulation
  }
  else
  {
    int nbNodes = m_mesh->getNumberOfNodes();

    float pmin = std::numeric_limits<float>::max();
    float pmax = std::numeric_limits<float>::lowest();

    for (int node = 0; node < nbNodes; node++)
    {
      float p = pnGlobal(node, i2);
      if (p < pmin) pmin = p;
      if (p > pmax) pmax = p;
    }

    const int NBINS = 10;
    std::vector<int> hist(NBINS, 0);

    float binWidth = (pmax - pmin) / NBINS;

    for (int node = 0; node < nbNodes; node++)
    {
      float p = pnGlobal(node, i2);

      int bin = (int)((p - pmin) / binWidth);
      if (bin == NBINS) bin = NBINS - 1;

      hist[bin]++;
    }

    for (int i = 0; i < NBINS; i++)
    {
      float bmin = pmin + i * binWidth;
      float bmax = bmin + binWidth;
      out << "Bin " << i << " [" << bmin << " , " << bmax
          << "] : " << hist[i] << " points\n";
    }
    out.close();
  }
}

/*TP 3 code*/

void SEMproxy::saveSnapshotPPM(int time_ms)
{
  // Extraire les données du champ de pression
  int nbNodes = m_mesh->getNumberOfNodes();
  std::vector<float> pressure_field;
  pressure_field.reserve(nbNodes);

  for (int node = 0; node < nbNodes; node++)
  {
    pressure_field.push_back(pnGlobal(node, i2));
  }

  // Dimensions de la grille
  int nx = nb_nodes_[0];
  int ny = nb_nodes_[1];
  int nz = nb_nodes_[2];

  // Position de la coupe (par défaut au milieu)
  int sliceIndex = ppmSliceIndex;

  if (sliceIndex == 0)
  {
    switch (ppmPlane)
    {
      case PPMWriter::SlicePlane::XY:
        sliceIndex = nz / 2;
        break;
      case PPMWriter::SlicePlane::XZ:
        sliceIndex = ny / 2;
        break;
      case PPMWriter::SlicePlane::YZ:
        sliceIndex = nx / 2;
        break;
    }
  }

  // Extraction de la slice
  auto slice =
      PPMWriter::extractSlice(pressure_field, nx, ny, nz, ppmPlane, sliceIndex);

  // Dimensions de l'image
  int width, height;
  switch (ppmPlane)
  {
    case PPMWriter::SlicePlane::XY:
      width = nx;
      height = ny;
      break;
    case PPMWriter::SlicePlane::XZ:
      width = nx;
      height = nz;
      break;
    case PPMWriter::SlicePlane::YZ:
      width = ny;
      height = nz;
      break;
  }

  // Nom du fichier
  std::string filename;
  if (ppmFolder.empty())
  {
    filename = "ppm_" + std::to_string(time_ms) + ".ppm";
  }
  else
  {
    filename = ppmFolder + "/ppm_" + std::to_string(time_ms) + ".ppm";
  }

  // Écriture de l'image PPM
  try
  {
    PPMWriter::writeSlice(filename, slice, width, height, ppmColormap);
    std::cout << "✓ PPM saved: " << filename << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cerr << "Error saving PPM: " << e.what() << std::endl;
  }
}

/*TP 3 code end*/

void SEMproxy::parsePointSismos(const std::string& filename)
{
  std::ifstream file(filename);
  if (!file.is_open())
  {
    std::cerr << "Erreur : could not open sismos file " << filename
              << std::endl;
    return;
  }

  sismosPoints.clear();

  std::string line;
  while (std::getline(file, line))
  {
    if (line.empty()) continue;

    float x, y, z;
    char c1, c2;  // for ','

    std::stringstream ss(line);

    if (!(ss >> x >> c1 >> y >> c2 >> z))
    {
      std::cerr << " invalid line in sismos file, line:" << line << std::endl;
      continue;
    }

    if (c1 != ',' || c2 != ',')
    {
      std::cerr << "Format error in sismos file, line:" << line << std::endl;
      continue;
    }

    sismosPoints.push_back({x, y, z});
  }

  file.close();

  std::cout << "✔ " << sismosPoints.size() << " sismos points loaded from "
            << filename << std::endl;
}

int SEMproxy::findClosestNode(float x, float y, float z)
{
  int nNodes = m_mesh->getNumberOfNodes();
  int closestNode = 0;
  float minDist2 = std::numeric_limits<float>::max();

  for (int i = 0; i < nNodes; i++)
  {
    float nx = m_mesh->nodeCoord(i, 0);
    float ny = m_mesh->nodeCoord(i, 1);
    float nz = m_mesh->nodeCoord(i, 2);

    float dx = x - nx;
    float dy = y - ny;
    float dz = z - nz;

    float dist2 = dx * dx + dy * dy + dz * dz;

    if (dist2 < minDist2)
    {
      minDist2 = dist2;
      closestNode = i;
    }
  }
  return closestNode;
}

// Initialize arrays
void SEMproxy::init_arrays()
{
  cout << "Allocate host memory for source and pressure values ..." << endl;

  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(
      myNumberOfRHS, m_mesh->getNumberOfPointsPerElement(), "RHSWeight");
  myRHSTerm = allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTerm");
  pnGlobal =
      allocateArray2D<arrayReal>(m_mesh->getNumberOfNodes(), 2, "pnGlobal");
  pnAtReceiver = allocateArray2D<arrayReal>(1, num_sample_, "pnAtReceiver");
  // Receiver
  rhsElementRcv = allocateVector<vectorInt>(1, "rhsElementRcv");
  rhsWeightsRcv = allocateArray2D<arrayReal>(
      1, m_mesh->getNumberOfPointsPerElement(), "RHSWeightRcv");
}

// Initialize sources
void SEMproxy::init_source()
{
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(1, 3, "RHSLocation");
  // std::cout << "All source are currently are coded on element 50." <<
  // std::endl;
  std::cout << "All source are currently are coded on middle element."
            << std::endl;
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  float lx = domain_size_[0];
  float ly = domain_size_[1];
  float lz = domain_size_[2];

  // Get source element index

  int source_index = floor((src_coord_[0] * ex) / lx) +
                     floor((src_coord_[1] * ey) / ly) * ex +
                     floor((src_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElement[i] = source_index;
  }

  // Get coordinates of the corners of the sourc element
  float cornerCoords[8][3];
  int I = 0;
  int nodes_corner[2] = {0, m_mesh->getOrder()};
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElement[0], i, j, k);
        cornerCoords[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoords[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoords[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  // initialize source term
  vector<float> sourceTerm =
      myUtils.computeSourceTerm(num_sample_, dt_, f0, sourceOrder);
  for (int j = 0; j < num_sample_; j++)
  {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }

  // get element number of source term
  myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myElementSource << endl
       << endl;

  int order = m_mesh->getOrder();

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }

  // Receiver computation
  int receiver_index = floor((rcv_coord_[0] * ex) / lx) +
                       floor((rcv_coord_[1] * ey) / ly) * ex +
                       floor((rcv_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElementRcv[i] = receiver_index;
  }

  // Get coordinates of the corners of the receiver element
  float cornerCoordsRcv[8][3];
  I = 0;
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
        cornerCoordsRcv[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoordsRcv[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoordsRcv[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }
}

SolverFactory::implemType SEMproxy::getImplem(string implemArg)
{
  if (implemArg == "makutu") return SolverFactory::MAKUTU;
  if (implemArg == "shiva") return SolverFactory::SHIVA;

  throw std::invalid_argument(
      "Implentation type does not follow any valid type.");
}

SolverFactory::meshType SEMproxy::getMesh(string meshArg)
{
  if (meshArg == "cartesian") return SolverFactory::Struct;
  if (meshArg == "ucartesian") return SolverFactory::Unstruct;

  std::cout << "Mesh type found is " << meshArg << std::endl;
  throw std::invalid_argument("Mesh type does not follow any valid type.");
}

SolverFactory::methodType SEMproxy::getMethod(string methodArg)
{
  if (methodArg == "sem") return SolverFactory::SEM;
  if (methodArg == "dg") return SolverFactory::DG;

  throw std::invalid_argument("Method type does not follow any valid type.");
}

float SEMproxy::find_cfl_dt(float cfl_factor)
{
  float sqrtDim3 = 1.73;  // to change for 2d
  float min_spacing = m_mesh->getMinSpacing();
  float v_max = m_mesh->getMaxSpeed();

  float dt = cfl_factor * min_spacing / (sqrtDim3 * v_max);

  return dt;
}
