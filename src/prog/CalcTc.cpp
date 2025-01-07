// Copyright (C) 2023 Lisa Biermann, João Viana
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program calculates the critical temperature for any found coexisting
 * phases
 *
 */

#include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/gravitational_waves/gwUtils/bounce_solution/ActionCalculation.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h> // for atoi, EXIT_FAILURE
#include <string>   // for string, operator<<
#include <utility>  // for pair
#include <vector>   // for vector

using namespace std;
using namespace BSMPT;
using namespace GravitationalWaves;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int firstline{0}, lastline{0};
  double tempstart{0}, tempend{300};
  std::string inputfile, outputfile;
  bool TerminalOutput{false};
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};
  bool UseMultithreading{true};
  bool CheckNLOStability{true};

  CLIOptions(int argc, char *argv[]);
  bool good() const;
};

int main(int argc, char *argv[])
try
{
  const CLIOptions args(argc, argv);
  if (not args.good())
  {
    return EXIT_FAILURE;
  }

  int linecounter = 1;
  std::ifstream infile(args.inputfile);
  if (!infile.good())
  {
    Logger::Write(LoggingLevel::Default,
                  "Input file " + args.inputfile + " not found ");
    return EXIT_FAILURE;
  }

  std::ofstream outfile(args.outputfile);
  if (!outfile.good())
  {
    Logger::Write(LoggingLevel::Default,
                  "Can not create file " + args.outputfile);
    return EXIT_FAILURE;
  }
  std::string linestr;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model);

  auto len_legTemp = (modelPointer->addLegendTemp()).size();

  while (getline(infile, linestr))
  {
    if (linecounter > args.lastline) break;

    if (linecounter == 1)
    {
      outfile << linestr << sep << modelPointer->addLegendCT() << sep
              << modelPointer->addLegendTemp() << std::endl;

      modelPointer->setUseIndexCol(linestr);
    }
    const clock_t begin_time = clock();
    if (linecounter >= args.firstline and linecounter <= args.lastline and
        linecounter != 1)
    {
      if (args.TerminalOutput)
      {
        Logger::Write(LoggingLevel::ProgDetailed,
                      "Currently at line " + std::to_string(linecounter));
      }
      std::pair<std::vector<double>, std::vector<double>> parameters =
          modelPointer->initModel(linestr);
      modelPointer->PerformVCTShift();
      if (args.firstline == args.lastline)
      {
        modelPointer->write();
      }

      Logger::Write(LoggingLevel::GWDetailed,
                    "Currently at line " + std::to_string(linecounter));

      MinimumTracer MinTracer(modelPointer);

      // EW symmetry restoration check
      Logger::Write(
          LoggingLevel::GWDetailed,
          "Check for EW symmetry restoration: " +
              std::to_string(MinTracer.IsThereEWSymmetryRestoration()) + "\n");

      // NLOVEV at T = 0 GeV check
      bool nlostable = modelPointer->CheckNLOVEV(
          MinTracer.ConvertToVEVDim(MinTracer.GetGlobalMinimum(0)));

      if (nlostable or not args.CheckNLOStability)
      {
        // phase tracking
        Logger::Write(LoggingLevel::GWDetailed,
                      "Track phases in between T_low = " +
                          std::to_string(args.tempstart) +
                          " and T_high = " + std::to_string(args.tempend));

        Vacuum vac(args.tempstart, args.tempend, MinTracer, modelPointer);

        Logger::Write(LoggingLevel::GWDetailed,
                      "Found and traced " +
                          std::to_string(vac.PhasesList.size()) +
                          " minima.\n-------------------------------");
        for (auto phase : vac.PhasesList)
        {
          Logger::Write(LoggingLevel::GWDetailed,
                        std::to_string(phase.T_low) + "\t->\t" +
                            std::to_string(phase.T_high));
        }
        Logger::Write(LoggingLevel::GWDetailed,
                      "-------------------------------");

        Logger::Write(LoggingLevel::GWDetailed,
                      "Number of coexisting phase regions : " +
                          std::to_string(vac.CoexPhasesList.size()));

        if (vac.CoexPhasesList.size() == 0)
        {
          outfile << linestr;
          outfile << sep << parameters.second;
          std::vector<int> error_no_coex_phase_found(len_legTemp, -1);
          outfile << sep << error_no_coex_phase_found;
          outfile << std::endl;
        }
        else
        {

          Logger::Write(
              LoggingLevel::GWDetailed,
              "Start the GW Calculation for the strongest EWPT coexPhase.");

          for (auto coex : vac.CoexPhasesList)
          {
            if (coex.num_phases == 1) continue;

            for (std::size_t TruePhaseCandidate = 1;
                 TruePhaseCandidate < coex.num_phases;
                 TruePhaseCandidate++)
            {
              double Tc = coex.CalculateTc(TruePhaseCandidate);
              std::cout << "Tc = " << Tc << std::endl;

              if (Tc > 0)
              {
                // calculate EW PT strengths xi = vEW / T
                double vEW_crit =
                    modelPointer->EWSBVEV(modelPointer->MinimizeOrderVEV(
                        coex.vec_phases[TruePhaseCandidate].Get(Tc).point));

                Logger::Write(LoggingLevel::GWDetailed,
                              "Found Tc = " + std::to_string(Tc) +
                                  "\n vEW_c = " + std::to_string(vEW_crit) +
                                  "\n xic = " + std::to_string(vEW_crit / Tc));

                outfile << linestr;
                outfile << sep << parameters.second;
                outfile << sep << Tc;
                outfile << sep << vEW_crit;
                outfile << sep << vEW_crit / Tc;
                outfile << sep
                        << coex.vec_phases[TruePhaseCandidate].Get(Tc).point;
                outfile << std::endl;
              }
              else // no or only one solution or Tc invalid
              {
                outfile << linestr;
                outfile << sep << parameters.second;
                std::vector<int> error_out_bounce_failed(len_legTemp, -2);
                outfile << sep << error_out_bounce_failed;
                outfile << std::endl;
              }
            }
          }
        }
      }
      else
      {
        Logger::Write(LoggingLevel::GWDetailed, "NLO VEV is unstable!\n");

        std::vector<int> error_out_NLO_unstable(len_legTemp, -3);
        outfile << linestr;
        outfile << sep << parameters.second;
        outfile << sep << error_out_NLO_unstable;
        outfile << std::endl;
      }

      std::cout << "Took\t" << float(clock() - begin_time) / CLOCKS_PER_SEC
                << " seconds.\n";
    }

    linecounter++;
    if (infile.eof()) break;
  }
  outfile.close();
  return EXIT_SUCCESS;
}
catch (int)
{
  return EXIT_SUCCESS;
}
catch (exception &e)
{
  Logger::Write(LoggingLevel::Default, e.what());
  return EXIT_FAILURE;
}

CLIOptions::CLIOptions(int argc, char *argv[])
{

  std::vector<std::string> args;
  for (int i{1}; i < argc; ++i)
    args.push_back(argv[i]);

  if (argc < 6 or args.at(0) == "--help")
  {
    int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
    std::stringstream ss;

    ss << std::boolalpha << "GW calculates the strength of GW signal"
       << std::endl
       << "It is called by " << std::endl
       << "./CalcTc [arguments]" << std::endl
       << "with the following arguments" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--help"
       << "Shows this menu" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--model="
       << "The model you want to investigate" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--input="
       << "The input file in tsv format" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--output="
       << "The output file in tsv format" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--tstart="
       << "The start temperature at which the phases should be evaluated"
       << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--tend="
       << "The end temperature until which the phases should be evaluated"
       << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--firstline="
       << "The first line in the input file to calculate the GW. Expects "
          "line 1 to be a legend."
       << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--lastline="
       << "The last line in the input file to calculate the GW." << std::endl;
    std::string GSLhelp{"--UseGSL="};
    GSLhelp += Minimizer::UseGSLDefault ? "true" : "false";
    ss << std::setw(SizeOfFirstColumn) << std::left << GSLhelp
       << "Use the GSL library to minimize the effective potential"
       << std::endl;
    std::string CMAEShelp{"--UseCMAES="};
    CMAEShelp += Minimizer::UseLibCMAESDefault ? "true" : "false";
    ss << std::setw(SizeOfFirstColumn) << std::left << CMAEShelp
       << "Use the CMAES library to minimize the effective potential"
       << std::endl;
    std::string NLoptHelp{"--UseNLopt="};
    NLoptHelp += Minimizer::UseNLoptDefault ? "true" : "false";
    ss << std::setw(SizeOfFirstColumn) << std::left << NLoptHelp
       << "Use the NLopt library to minimize the effective potential"
       << std::endl;
    ss << std::setw(SizeOfFirstColumn) << std::left
       << "--UseMultithreading = true"
       << "Enables/Disables multi threading for the minimizers" << std::endl;
    ss << std::setw(CheckNLOStability) << std::left
       << "--CheckNLOStability = true"
       << "Enables/Disables check for NLO stability" << std::endl;
    ss << std::setw(SizeOfFirstColumn) << std::left << "--TerminalOutput="
       << "y/n Turns on additional information in the terminal during "
          "the calculation."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    ShowLoggerHelp();
    ShowInputError();
  }

  if (args.size() > 0 and args.at(0) == "--help")
  {
    throw int{0};
  }
  else if (argc < 6)
  {
    throw std::runtime_error("Too few arguments.");
  }

  std::string prefix{"--"};
  bool UsePrefix = StringStartsWith(args.at(0), prefix);
  std::vector<std::string> UnusedArgs;
  if (UsePrefix)
  {
    for (const auto &arg : args)
    {
      auto el = arg;
      std::transform(el.begin(), el.end(), el.begin(), ::tolower);
      if (StringStartsWith(el, "--model="))
      {
        Model =
            BSMPT::ModelID::getModel(el.substr(std::string("--model=").size()));
      }
      else if (StringStartsWith(el, "--input="))
      {
        inputfile = arg.substr(std::string("--input=").size());
      }
      else if (StringStartsWith(el, "--output="))
      {
        outputfile = arg.substr(std::string("--output=").size());
      }
      else if (StringStartsWith(el, "--tstart="))
      {
        tempstart = std::stod(el.substr(std::string("--tstart=").size()));
      }
      else if (StringStartsWith(el, "--tend="))
      {
        tempend = std::stod(el.substr(std::string("--tend=").size()));
      }
      else if (StringStartsWith(el, "--firstline="))
      {
        firstline = std::stoi(el.substr(std::string("--firstline=").size()));
      }
      else if (StringStartsWith(el, "--lastline="))
      {
        lastline = std::stoi(el.substr(std::string("--lastline=").size()));
      }
      else if (StringStartsWith(el, "--terminaloutput="))
      {
        TerminalOutput =
            el.substr(std::string("--terminaloutput=").size()) == "y";
      }
      else if (StringStartsWith(el, "--usegsl="))
      {
        UseGSL = el.substr(std::string("--usegsl=").size()) == "true";
      }
      else if (StringStartsWith(el, "--usecmaes="))
      {
        UseCMAES = el.substr(std::string("--usecmaes=").size()) == "true";
      }
      else if (StringStartsWith(el, "--usenlopt="))
      {
        UseNLopt = el.substr(std::string("--usenlopt=").size()) == "true";
      }
      else if (StringStartsWith(el, "--usemultithreading="))
      {
        UseMultithreading =
            el.substr(std::string("--usemultithreading=").size()) == "true";
      }
      else if (StringStartsWith(el, "--checknlostability="))
      {
        CheckNLOStability =
            el.substr(std::string("--checknlostability=").size()) == "true";
      }
      else
      {
        UnusedArgs.push_back(el);
      }
    }
    WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);
    SetLogger(UnusedArgs);
  }
  else
  {
    Model      = ModelID::getModel(args.at(0));
    inputfile  = args.at(1);
    outputfile = args.at(2);
    firstline  = std::stoi(args.at(3));
    lastline   = std::stoi(args.at(4));
    tempstart  = std::stod(args.at(5));
    tempend    = std::stod(args.at(6));
  }
}

bool CLIOptions::good() const
{
  if (UseGSL and not Minimizer::UseGSLDefault)
  {
    throw std::runtime_error(
        "You set --UseGSL=true but GSL was not found during compilation.");
  }
  if (UseCMAES and not Minimizer::UseLibCMAESDefault)
  {
    throw std::runtime_error("You set --UseCMAES=true but CMAES was not "
                             "found during compilation.");
  }
  if (UseNLopt and not Minimizer::UseNLoptDefault)
  {
    throw std::runtime_error("You set --UseNLopt=true but NLopt was not "
                             "found during compilation.");
  }
  if (WhichMinimizer == 0)
  {
    throw std::runtime_error(
        "You disabled all minimizers. You need at least one.");
  }

  if (Model == ModelID::ModelIDs::NotSet)
  {

    Logger::Write(
        LoggingLevel::Default,
        "Your Model parameter does not match with the implemented Models.");
    ShowInputError();
    return false;
  }
  if (firstline < 1)
  {
    Logger::Write(LoggingLevel::Default, "Start line counting with 1");
    return false;
  }
  if (firstline > lastline)
  {
    Logger::Write(LoggingLevel::Default, "firstline is smaller then lastline ");
    return false;
  }
  if (tempstart > tempend)
  {
    Logger::Write(
        LoggingLevel::Default,
        "Invalid temperature choice. tempend has to be above tempstart.");
    return false;
  }
  return true;
}
