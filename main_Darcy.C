// $Id$
//==============================================================================
//!
//! \file main_Darcy.C
//!
//! \date Mar 27 2015
//!
//! \author Yared Bekele
//!
//! \brief Main program for an isogeometric solver for Darcy flow.
//!
//==============================================================================

#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMDarcy.h"
#include "SIMSolverAdap.h"
#include "Utilities.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "AppCommon.h"
#include "TimeStep.h"


template<class Dim, template<class T> class Solver>
int runSimulator(char* infile)
{
  SIMDarcy<Dim> darcy;
  Solver<SIMDarcy<Dim>> solver(darcy);
  int res = ConfigureSIM(darcy, infile, false);

  if (res)
    return res;

  // HDF5 output
  std::unique_ptr<DataExporter> exporter;

  if (darcy.opt.dumpHDF5(infile))
    exporter.reset(SIM::handleDataOutput(darcy, solver,
                                         darcy.opt.hdf5, false, 1, 1));

  return solver.solveProblem(infile, exporter.get());
}


template<class Dim>
int runSimulator1(char* infile, bool adaptive)
{
  if (adaptive)
    return runSimulator<Dim, SIMSolverAdap>(infile);

  return runSimulator<Dim, SIMSolver>(infile);
}


int main(int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  int  i;
  char ndim = 3;
  char* infile = 0;
  bool adaptive = false;

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      ndim = 2;
    else if (!strcmp(argv[i],"-1D"))
      ndim = 1;
    else if (!strcmp(argv[i],"-adap"))
      adaptive = true;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
              <<" [-free] [-lag|-spec|-LR] [-1D|-2D] [-nGauss <n>]"
              <<"\n       [-vtf <format> [-nviz <nviz>]"
              <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n"
              <<"       [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]\n"
              <<"       [-ignore <p1> <p2> ...] [-fixDup]" << std::endl;
    return 0;
  }

  if (adaptive)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\n >>> IFEM Darcy equation solver <<<"
             <<"\n ====================================\n"
             <<"\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
  IFEM::cout <<"\n\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (ndim == 3)
    return runSimulator1<SIM3D>(infile,adaptive);
  else if (ndim == 2)
    return runSimulator1<SIM2D>(infile,adaptive);
  else
    return runSimulator1<SIM1D>(infile,adaptive);

  return 1;
}
