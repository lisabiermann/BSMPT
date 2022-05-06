// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Calculates the electroweak phase transition for a given Inputfile for a given
 * subset of lines in the file and adds it at the end of the line in the format
 * T_c v_c all single vevs. One parameter point per line.
 *
 */

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h>   // for Class_Potential_Origin
#include <BSMPT/models/ClassPotentialR2HDMEFT.h> // for Class_Potential_R2HDMEFT
#include <BSMPT/models/ClassPotentialR2HDMEFTPHI6.h> // for Class_Potential_R2HDMEFTPHI6
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
#include <algorithm> // for copy
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>   // for unique_ptr
#include <stdlib.h> // for EXIT_FAILURE, atoi
#include <string>   // for getline, string
#include <utility>  // for pair
#include <vector>   // for vector
using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int Line{};
  std::string InputFile;
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};

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
  std::ifstream infile(args.InputFile);
  if (!infile.good())
  {
    Logger::Write(LoggingLevel::Default, "Input file not found ");
    return EXIT_FAILURE;
  }

  Logger::Write(LoggingLevel::ProgDetailed, "Found file");

  std::string linestr;
  std::unique_ptr<Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model);

  Logger::Write(LoggingLevel::ProgDetailed, "Created modelpointer ");

  while (getline(infile, linestr))
  {
    if (linecounter > args.Line) break;

    if (linecounter == 1)
    {

      modelPointer->setUseIndexCol(linestr);
    }
    if (linecounter == args.Line and linecounter != 1)
    {
      Logger::Write(LoggingLevel::ProgDetailed, "Found line");
      modelPointer->initModel(linestr);
      modelPointer->write();
      std::vector<double> dummy;
      modelPointer->Debugging(dummy, dummy);
      modelPointer->CheckImplementation(args.WhichMinimizer);
    }
    linecounter++;
    if (infile.eof()) break;
  }
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

  if (argc == 5)
  {
    BSMPT::Models::Class_Potential_R2HDMEFT::Op6_111111 = std::stod(args.at(3));
  }
  else
  {
    BSMPT::Models::Class_Potential_R2HDMEFT::Op6_111111 = 0;
  }

  if (argc == 12)
  {
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_111111 = std::stod(args.at(3));
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_111122 = std::stod(args.at(4));
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_122111 = std::stod(args.at(5));
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_121211 = std::stod(args.at(6));
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_222222 = std::stod(args.at(7));
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_112222 = std::stod(args.at(8));
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_122122 = std::stod(args.at(9));
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_121222 = std::stod(args.at(10));
  }
  else
  {
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_111111 = 0;
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_111122 = 0;
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_122111 = 0;
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_121211 = 0;
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_222222 = 0;
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_112222 = 0;
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_122122 = 0;
    BSMPT::Models::Class_Potential_R2HDMEFTPHI6::Op6_121222 = 0;
  }

  if (argc < 4 or args.at(0) == "--help")
  {
    std::stringstream ss;
    int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
    ss << "Test performs a serious of tests on the given model. "
          "Intended for testing new models."
       << std::endl
       << "It is called either by " << std::endl
       << "./Test model input Line" << std::endl
       << "or with the following arguments" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--help"
       << "Shows this menu" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--model="
       << "The model you want to test" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--input="
       << "The input file in tsv format" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--Line="
       << "The line in the input file with the parameter point used to "
          "check the model."
       << std::endl;
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
    Logger::Write(LoggingLevel::Default, ss.str());
    ShowLoggerHelp();
    ShowInputError();
  }

  if (args.size() > 0 and args.at(0) == "--help")
  {
    throw int{0};
  }
  else if (argc < 4)
  {
    throw std::runtime_error("Too few arguments.");
  }

  const std::string prefix{"--"};
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
        InputFile = arg.substr(std::string("--input=").size());
        Logger::Write(LoggingLevel::ProgDetailed, "Inputfile = " + InputFile);
      }
      else if (StringStartsWith(el, "--line="))
      {
        Line = std::stoi(el.substr(std::string("--line=").size()));
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
    Model     = ModelID::getModel(args.at(0));
    InputFile = args.at(1);
    Line      = std::stoi(args.at(2));
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
    throw std::runtime_error(
        "You set --UseCMAES=true but CMAES was not found during compilation.");
  }
  if (UseNLopt and not Minimizer::UseNLoptDefault)
  {
    throw std::runtime_error(
        "You set --UseNLopt=true but NLopt was not found during compilation.");
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
  if (Line < 1)
  {
    Logger::Write(LoggingLevel::Default, "Start line counting with 1");
    return false;
  }
  return true;
}
