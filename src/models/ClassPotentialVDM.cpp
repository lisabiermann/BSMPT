// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/ClassPotentialVDM.h>
#include <BSMPT/models/SMparam.h> // for C_vev0, C_MassTop, C_g
#include <algorithm>              // for max, copy
#include <iostream>               // for operator<<, endl, basic_o...
#include <memory>                 // for allocator_traits<>::value...
#include <stddef.h>               // for std::size_t

#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

/**
 * @file
 * Template for adding a new model class
 */

namespace BSMPT
{
namespace Models
{
/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();
  */
Class_VDM::Class_VDM(const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model = ModelID::ModelIDs::VDM; // global int constant which will be used to
                                  // tell the program which model is called
  NNeutralHiggs = 4;              // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 2; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 6;  // number of parameters in the tree-Level Lagrangian
  nParCT = 11; // number of parameters in the counterterm potential

  nVEV = 2; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;
  NGauge = 5; // overwrite number of gauge bosons as we have an additional

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  VevOrder[0] = 2;
  VevOrder[1] = 4;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_VDM::~Class_VDM()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_VDM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dmuHsq");
  labels.push_back("dlambdaH");
  labels.push_back("dkappa");
  labels.push_back("dmuSsq");
  labels.push_back("dlambdaS");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  labels.push_back("dT5");
  labels.push_back("dT6");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_VDM::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c"); // Label for the critical temperature
  labels.push_back("v_c"); // Label for the critical vev
  labels.push_back("omega_c/T_c");
  labels.push_back("omega_c");
  labels.push_back(
      "omega_sc"); // Label for v_c/T_c, you could use xi_c also for example
  // out += "Your VEV order"; // Now you have to put the label for your vevs
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 *
 */
std::vector<std::string> Class_VDM::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);
  // here you have to define the particle names in the vector particles

  particles[0] = "G+";
  particles[1] = "G-";
  particles[2] = "G0";
  particles[3] = "H1";
  particles[4] = "H2";
  particles[5] = "H3";

  std::string out = "Tree_";
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = i; j < NHiggs; j++)
    {
      for (std::size_t k = j; k < NHiggs; k++)
      {
        labels.push_back("Tree_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CT_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CW_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
      }
    }
  }

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_VDM::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order";
  labels.push_back("omega");
  labels.push_back("omega_s");
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_VDM::ReadAndSet(const std::string &linestr, std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k < 7; k++)
  {
    // alte zuordnung
    ss >> tmp;

    // if (k == 1)
    //   MH1 = tmp; // MH1
    // else if (k == 2)
    //   MH2 = tmp; // MH2
    // else if (k == 3)
    //   v   = tmp; // v
    // else if (k == 4)
    //   vs  = tmp; // vs
    // else if (k == 5)
    //   alpha = tmp; // alpha
    // else if (k == 6)
    //   MX = tmp; // MX

    if (k == 1)
      MH1 = tmp; //  MH1
    else if (k == 2)
      MH2 = tmp; // MH2
    else if (k == 3)
      MX = tmp; // MX
    else if (k == 4)
      alpha = tmp; // alpha
    else if (k == 5)
      v = tmp; // v
    else if (k == 6)
      C_gX = tmp; // kopplung
  }

  vs = MX / C_gX;

  par[0] = MH1;
  par[1] = MH2;
  par[2] = v;
  par[3] = vs;
  par[4] = alpha;
  par[5] = MX;

  // par[0] =  MH1;
  // par[1] =  MH2;
  // par[2] =  v;
  // par[3] =  MX /gX; //vs
  // par[4] =  alpha;
  // par[5] =  MX;

  set_gen(par); // This you have to call so that everything will be set
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_VDM::set_gen(const std::vector<double> &par)
{

  v       = SMConstants.C_vev0;
  vs      = par[3];
  lambdaH = (par[0] * par[0] * cos(par[4]) * cos(par[4]) +
             par[1] * par[1] * sin(par[4]) * sin(par[4])) /
            (2 * par[2] * par[2]);
  kappa = (par[0] * par[0] - par[1] * par[1]) * sin(par[4]) * cos(par[4]) /
          (par[2] * par[3]);
  lambdaS = (par[1] * par[1] * cos(par[4]) * cos(par[4]) +
             par[0] * par[0] * sin(par[4]) * sin(par[4])) /
            (2 * par[3] * par[3]);
  C_gX = par[5] / par[3];

  muHsq = vs * vs * kappa / 2 + v * v * lambdaH; // tadpole conditions
  muSsq = kappa * v * v / 2 + lambdaS * vs * vs;

  scale = v; // hier evtl vh

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);

  // Here you have to set the vector vevTreeMin. The vector vevTree will then be
  // set by the function MinimizeOrderVEV
  vevTreeMin[0] = v;
  vevTreeMin[1] = vs;

  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_VDM::set_CT_Pot_Par(const std::vector<double> &par)
{

  dmuHsq   = par[0];
  dlambdaH = par[1];
  dkappa   = par[2];
  dmuSsq   = par[3];
  dlambdaS = par[4];
  dT1      = par[5];
  dT2      = par[6];
  dT3      = par[7];
  dT4      = par[8];
  dT5      = par[9];
  dT6      = par[10];

  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;
  Curvature_Higgs_CT_L1[4] = dT5;
  Curvature_Higgs_CT_L1[5] = dT6;

  Curvature_Higgs_CT_L2[0][0] = -dmuHsq;
  Curvature_Higgs_CT_L2[1][1] = -dmuHsq;
  Curvature_Higgs_CT_L2[2][2] = -dmuHsq;
  Curvature_Higgs_CT_L2[3][3] = -dmuHsq;
  Curvature_Higgs_CT_L2[4][4] = -dmuSsq;
  Curvature_Higgs_CT_L2[5][5] = -dmuSsq;

  Curvature_Higgs_CT_L4[0][0][0][0] = 6 * dlambdaH;
  Curvature_Higgs_CT_L4[0][0][1][1] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[0][0][2][2] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[0][0][3][3] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[0][0][4][4] = dkappa;
  Curvature_Higgs_CT_L4[0][0][5][5] = dkappa;
  Curvature_Higgs_CT_L4[0][1][0][1] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[0][1][1][0] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[0][2][0][2] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[0][2][2][0] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[0][3][0][3] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[0][3][3][0] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[0][4][0][4] = dkappa;
  Curvature_Higgs_CT_L4[0][4][4][0] = dkappa;
  Curvature_Higgs_CT_L4[0][5][0][5] = dkappa;
  Curvature_Higgs_CT_L4[0][5][5][0] = dkappa;
  Curvature_Higgs_CT_L4[1][0][0][1] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[1][0][1][0] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[1][1][0][0] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[1][1][1][1] = 6 * dlambdaH;
  Curvature_Higgs_CT_L4[1][1][2][2] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[1][1][3][3] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[1][1][4][4] = dkappa;
  Curvature_Higgs_CT_L4[1][1][5][5] = dkappa;
  Curvature_Higgs_CT_L4[1][2][1][2] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[1][2][2][1] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[1][3][1][3] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[1][3][3][1] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[1][4][1][4] = dkappa;
  Curvature_Higgs_CT_L4[1][4][4][1] = dkappa;
  Curvature_Higgs_CT_L4[1][5][1][5] = dkappa;
  Curvature_Higgs_CT_L4[1][5][5][1] = dkappa;
  Curvature_Higgs_CT_L4[2][0][0][2] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[2][0][2][0] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[2][1][1][2] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[2][1][2][1] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[2][2][0][0] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[2][2][1][1] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[2][2][2][2] = 6 * dlambdaH;
  Curvature_Higgs_CT_L4[2][2][3][3] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[2][2][4][4] = dkappa;
  Curvature_Higgs_CT_L4[2][2][5][5] = dkappa;
  Curvature_Higgs_CT_L4[2][3][2][3] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[2][3][3][2] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[2][4][2][4] = dkappa;
  Curvature_Higgs_CT_L4[2][4][4][2] = dkappa;
  Curvature_Higgs_CT_L4[2][5][2][5] = dkappa;
  Curvature_Higgs_CT_L4[2][5][5][2] = dkappa;
  Curvature_Higgs_CT_L4[3][0][0][3] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[3][0][3][0] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[3][1][1][3] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[3][1][3][1] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[3][2][2][3] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[3][2][3][2] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[3][3][0][0] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[3][3][1][1] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[3][3][2][2] = 2 * dlambdaH;
  Curvature_Higgs_CT_L4[3][3][3][3] = 6 * dlambdaH;
  Curvature_Higgs_CT_L4[3][3][4][4] = dkappa;
  Curvature_Higgs_CT_L4[3][3][5][5] = dkappa;
  Curvature_Higgs_CT_L4[3][4][3][4] = dkappa;
  Curvature_Higgs_CT_L4[3][4][4][3] = dkappa;
  Curvature_Higgs_CT_L4[3][5][3][5] = dkappa;
  Curvature_Higgs_CT_L4[3][5][5][3] = dkappa;
  Curvature_Higgs_CT_L4[4][0][0][4] = dkappa;
  Curvature_Higgs_CT_L4[4][0][4][0] = dkappa;
  Curvature_Higgs_CT_L4[4][1][1][4] = dkappa;
  Curvature_Higgs_CT_L4[4][1][4][1] = dkappa;
  Curvature_Higgs_CT_L4[4][2][2][4] = dkappa;
  Curvature_Higgs_CT_L4[4][2][4][2] = dkappa;
  Curvature_Higgs_CT_L4[4][3][3][4] = dkappa;
  Curvature_Higgs_CT_L4[4][3][4][3] = dkappa;
  Curvature_Higgs_CT_L4[4][4][0][0] = dkappa;
  Curvature_Higgs_CT_L4[4][4][1][1] = dkappa;
  Curvature_Higgs_CT_L4[4][4][2][2] = dkappa;
  Curvature_Higgs_CT_L4[4][4][3][3] = dkappa;
  Curvature_Higgs_CT_L4[4][4][4][4] = 6 * dlambdaS;
  Curvature_Higgs_CT_L4[4][4][5][5] = 2 * dlambdaS;
  Curvature_Higgs_CT_L4[4][5][4][5] = 2 * dlambdaS;
  Curvature_Higgs_CT_L4[4][5][5][4] = 2 * dlambdaS;
  Curvature_Higgs_CT_L4[5][0][0][5] = dkappa;
  Curvature_Higgs_CT_L4[5][0][5][0] = dkappa;
  Curvature_Higgs_CT_L4[5][1][1][5] = dkappa;
  Curvature_Higgs_CT_L4[5][1][5][1] = dkappa;
  Curvature_Higgs_CT_L4[5][2][2][5] = dkappa;
  Curvature_Higgs_CT_L4[5][2][5][2] = dkappa;
  Curvature_Higgs_CT_L4[5][3][3][5] = dkappa;
  Curvature_Higgs_CT_L4[5][3][5][3] = dkappa;
  Curvature_Higgs_CT_L4[5][4][4][5] = 2 * dlambdaS;
  Curvature_Higgs_CT_L4[5][4][5][4] = 2 * dlambdaS;
  Curvature_Higgs_CT_L4[5][5][0][0] = dkappa;
  Curvature_Higgs_CT_L4[5][5][1][1] = dkappa;
  Curvature_Higgs_CT_L4[5][5][2][2] = dkappa;
  Curvature_Higgs_CT_L4[5][5][3][3] = dkappa;
  Curvature_Higgs_CT_L4[5][5][4][4] = 2 * dlambdaS;
  Curvature_Higgs_CT_L4[5][5][5][5] = 6 * dlambdaS;
}

/**
 * console output of all Parameters
 */
void Class_VDM::write() const
{

  std::stringstream ss;
  ss << "Model = " << Model << std::endl;

  ss << "The parameters are : " << std::endl;
  ss << "\tmuHsq = " << muHsq << std::endl
     << "\tlambdaH = " << lambdaH << std::endl
     << "\tkappa = " << kappa << std::endl
     << "\tmuSsq = " << muSsq << std::endl
     << "\tlambdaS = " << lambdaS << std::endl
     << "\tv = " << v << "\n"
     << "\tvs = " << vs << "\n";

  ss << "The counterterm parameters are : " << std::endl;
  ss << "\tdmuHsq= " << dmuHsq << std::endl
     << "\tdlambdaH = " << dlambdaH << std::endl
     << "\tdkappa = " << dkappa << std::endl
     << "\tdmuSsq = " << dmuSsq << "\n"
     << "\tdlambdaS = " << dlambdaS << "\n"
     << "\tdT1 = " << dT1 << "\n"
     << "\tdT2 = " << dT2 << "\n"
     << "\tdT3 = " << dT3 << "\n"
     << "\tdT4 = " << dT4 << "\n"
     << "\tdT5 = " << dT5 << "\n"
     << "\tdT6 = " << dT6 << std::endl;

  ss << "The scale is given by mu = " << scale << " GeV " << std::endl;
  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_VDM::calc_CT() const
{

  std::vector<double> parCT;
  parCT.resize(nParCT);

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsdone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  // Here you have to use your formulae for the counterterm scheme
  parCT[0] = (-HesseWeinberg(2, 2) * v - HesseWeinberg(2, 4) * vs +
              3 * HesseWeinberg(3, 3) * v) /
             v / 2;
  parCT[1] = (-HesseWeinberg(2, 2) + HesseWeinberg(3, 3)) * pow(v, -2) / 2;
  parCT[2] = -HesseWeinberg(2, 4) / v / vs;
  parCT[3] = (-HesseWeinberg(2, 4) * v -
              vs * (HesseWeinberg(4, 4) - 3 * HesseWeinberg(5, 5))) /
             vs / 2;
  parCT[4]  = (-HesseWeinberg(4, 4) + HesseWeinberg(5, 5)) * pow(vs, -2) / 2;
  parCT[5]  = -NablaWeinberg(0);
  parCT[6]  = -NablaWeinberg(1);
  parCT[7]  = HesseWeinberg(3, 3) * v - NablaWeinberg(2);
  parCT[8]  = -NablaWeinberg(3);
  parCT[9]  = HesseWeinberg(5, 5) * vs - NablaWeinberg(4);
  parCT[10] = -NablaWeinberg(5);

  return parCT;
}

void Class_VDM::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  std::vector<double> HiggsOrder(NHiggs);
  // Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] =
  // 5 you always want your 6th lightest particle to be the first particle in
  // the vector (which has the index 5 because they are sorted by mass)

  // example for keeping the mass order
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsOrder[i] = i;
  }

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        GaugeBasis[i][j][k] =
            TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
      }
    }
  }

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  MatrixXd HiggsRotSort(NHiggs, NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
  }

  TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
  TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
  TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        TripleHiggsCorrectionsCWPhysical[i][j][k]   = 0;
        TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
        TripleHiggsCorrectionsCTPhysical[i][j][k]   = 0;
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          for (std::size_t m = 0; m < NHiggs; m++)
          {
            for (std::size_t n = 0; n < NHiggs; n++)
            {
              double RotFac =
                  HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
              TripleHiggsCorrectionsCWPhysical[i][j][k] +=
                  RotFac * GaugeBasis[l][m][n];
              TripleHiggsCorrectionsTreePhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3[l][m][n];
              TripleHiggsCorrectionsCTPhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3_CT[l][m][n];
            }
          }
        }
      }
    }
  }
}

void Class_VDM::SetCurvatureArrays()
{
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

  initVectors();

  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  Curvature_Higgs_L2[0][0] = -muHsq;
  Curvature_Higgs_L2[1][1] = -muHsq;
  Curvature_Higgs_L2[2][2] = -muHsq;
  Curvature_Higgs_L2[3][3] = -muHsq;
  Curvature_Higgs_L2[4][4] = -muSsq;
  Curvature_Higgs_L2[5][5] = -muSsq;

  Curvature_Higgs_L4[0][0][0][0] = 0.3e1 * 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][0][1][1] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][0][2][2] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][0][3][3] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][0][4][4] = kappa;
  Curvature_Higgs_L4[0][0][5][5] = kappa;
  Curvature_Higgs_L4[0][1][0][1] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][1][1][0] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][2][0][2] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][2][2][0] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][3][0][3] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][3][3][0] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[0][4][0][4] = kappa;
  Curvature_Higgs_L4[0][4][4][0] = kappa;
  Curvature_Higgs_L4[0][5][0][5] = kappa;
  Curvature_Higgs_L4[0][5][5][0] = kappa;
  Curvature_Higgs_L4[1][0][0][1] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][0][1][0] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][1][0][0] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][1][1][1] = 0.3e1 * 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][1][2][2] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][1][3][3] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][1][4][4] = kappa;
  Curvature_Higgs_L4[1][1][5][5] = kappa;
  Curvature_Higgs_L4[1][2][1][2] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][2][2][1] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][3][1][3] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][3][3][1] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[1][4][1][4] = kappa;
  Curvature_Higgs_L4[1][4][4][1] = kappa;
  Curvature_Higgs_L4[1][5][1][5] = kappa;
  Curvature_Higgs_L4[1][5][5][1] = kappa;
  Curvature_Higgs_L4[2][0][0][2] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][0][2][0] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][1][1][2] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][1][2][1] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][2][0][0] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][2][1][1] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][2][2][2] = 0.3e1 * 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][2][3][3] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][2][4][4] = kappa;
  Curvature_Higgs_L4[2][2][5][5] = kappa;
  Curvature_Higgs_L4[2][3][2][3] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][3][3][2] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[2][4][2][4] = kappa;
  Curvature_Higgs_L4[2][4][4][2] = kappa;
  Curvature_Higgs_L4[2][5][2][5] = kappa;
  Curvature_Higgs_L4[2][5][5][2] = kappa;
  Curvature_Higgs_L4[3][0][0][3] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][0][3][0] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][1][1][3] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][1][3][1] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][2][2][3] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][2][3][2] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][3][0][0] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][3][1][1] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][3][2][2] = 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][3][3][3] = 0.3e1 * 0.2e1 * lambdaH;
  Curvature_Higgs_L4[3][3][4][4] = kappa;
  Curvature_Higgs_L4[3][3][5][5] = kappa;
  Curvature_Higgs_L4[3][4][3][4] = kappa;
  Curvature_Higgs_L4[3][4][4][3] = kappa;
  Curvature_Higgs_L4[3][5][3][5] = kappa;
  Curvature_Higgs_L4[3][5][5][3] = kappa;
  Curvature_Higgs_L4[4][0][0][4] = kappa;
  Curvature_Higgs_L4[4][0][4][0] = kappa;
  Curvature_Higgs_L4[4][1][1][4] = kappa;
  Curvature_Higgs_L4[4][1][4][1] = kappa;
  Curvature_Higgs_L4[4][2][2][4] = kappa;
  Curvature_Higgs_L4[4][2][4][2] = kappa;
  Curvature_Higgs_L4[4][3][3][4] = kappa;
  Curvature_Higgs_L4[4][3][4][3] = kappa;
  Curvature_Higgs_L4[4][4][0][0] = kappa;
  Curvature_Higgs_L4[4][4][1][1] = kappa;
  Curvature_Higgs_L4[4][4][2][2] = kappa;
  Curvature_Higgs_L4[4][4][3][3] = kappa;
  Curvature_Higgs_L4[4][4][4][4] = 0.3e1 * 0.2e1 * lambdaS;
  Curvature_Higgs_L4[4][4][5][5] = 0.2e1 * lambdaS;
  Curvature_Higgs_L4[4][5][4][5] = 0.2e1 * lambdaS;
  Curvature_Higgs_L4[4][5][5][4] = 0.2e1 * lambdaS;
  Curvature_Higgs_L4[5][0][0][5] = kappa;
  Curvature_Higgs_L4[5][0][5][0] = kappa;
  Curvature_Higgs_L4[5][1][1][5] = kappa;
  Curvature_Higgs_L4[5][1][5][1] = kappa;
  Curvature_Higgs_L4[5][2][2][5] = kappa;
  Curvature_Higgs_L4[5][2][5][2] = kappa;
  Curvature_Higgs_L4[5][3][3][5] = kappa;
  Curvature_Higgs_L4[5][3][5][3] = kappa;
  Curvature_Higgs_L4[5][4][4][5] = 0.2e1 * lambdaS;
  Curvature_Higgs_L4[5][4][5][4] = 0.2e1 * lambdaS;
  Curvature_Higgs_L4[5][5][0][0] = kappa;
  Curvature_Higgs_L4[5][5][1][1] = kappa;
  Curvature_Higgs_L4[5][5][2][2] = kappa;
  Curvature_Higgs_L4[5][5][3][3] = kappa;
  Curvature_Higgs_L4[5][5][4][4] = 0.2e1 * lambdaS;
  Curvature_Higgs_L4[5][5][5][5] = 0.3e1 * 0.2e1 * lambdaS;

  Curvature_Gauge_G2H2[0][0][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][0][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][0][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][2] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][1] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][0][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][2] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][1] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][3][0][0] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][1][1] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][2][2] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][3][3] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[4][4][4][4] = 2 * C_gX * C_gX;
  Curvature_Gauge_G2H2[4][4][5][5] = 2 * C_gX * C_gX;

  std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
  V11 = SMConstants.C_Vud;
  V12 = SMConstants.C_Vus;
  V13 = SMConstants.C_Vub;
  V21 = SMConstants.C_Vcd;
  V22 = SMConstants.C_Vcs;
  V23 = SMConstants.C_Vcb;
  V31 = SMConstants.C_Vtd;
  V32 = SMConstants.C_Vts;
  V33 = SMConstants.C_Vtb;

  Curvature_Lepton_F2H1[0][1][2] = 0.1e1 / v * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[0][1][3] = II / v * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][2] = 0.1e1 / v * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][3] = II / v * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / v * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][1] = II / v * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[2][3][2] = 0.1e1 / v * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[2][3][3] = II / v * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][2] = 0.1e1 / v * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][3] = II / v * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / v * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][1] = II / v * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[4][5][2] = 0.1e1 / v * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[4][5][3] = II / v * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][2] = 0.1e1 / v * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][3] = II / v * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / v * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][1] = II / v * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / v * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[6][1][1] = II / v * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / v * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[7][3][1] = II / v * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / v * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[8][5][1] = II / v * SMConstants.C_MassTau;

  Curvature_Quark_F2H1[0][6][2] = 0.1e1 / v * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][6][3] = -II / v * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][9][0] = -0.1e1 / v * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][9][1] = II / v * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][10][0] =
      -0.1e1 / v * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][10][1] = II / v * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][11][0] =
      -0.1e1 / v * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[0][11][1] = II / v * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[1][7][2]  = 0.1e1 / v * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[1][7][3]  = -II / v * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[1][9][0] =
      -0.1e1 / v * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][9][1] = II / v * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][10][0] =
      -0.1e1 / v * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][10][1] = II / v * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][11][0] =
      -0.1e1 / v * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[1][11][1] = II / v * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[2][8][2]  = 0.1e1 / v * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[2][8][3]  = -II / v * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[2][9][0] =
      -0.1e1 / v * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][9][1] = II / v * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][10][0] =
      -0.1e1 / v * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][10][1] = II / v * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][11][0] =
      -0.1e1 / v * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[2][11][1] = II / v * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[3][6][0]  = V11 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][6][1]  = II * V11 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][7][0]  = V21 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][7][1]  = II * V21 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][8][0]  = V31 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][8][1]  = II * V31 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][9][2]  = 0.1e1 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][9][3]  = II / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[4][6][0]  = V12 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][6][1]  = II * V12 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][7][0]  = V22 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][7][1]  = II * V22 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][8][0]  = V32 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][8][1]  = II * V32 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][10][2] = 0.1e1 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][10][3] = II / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[5][6][0]  = V13 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][6][1]  = II * V13 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][7][0]  = V23 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][7][1]  = II * V23 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][8][0]  = V33 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][8][1]  = II * V33 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][11][2] = 0.1e1 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][11][3] = II / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[6][0][2]  = 0.1e1 / v * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[6][0][3]  = -II / v * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[6][3][0]  = V11 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[6][3][1]  = II * V11 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[6][4][0]  = V12 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[6][4][1]  = II * V12 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[6][5][0]  = V13 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[6][5][1]  = II * V13 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[7][1][2]  = 0.1e1 / v * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[7][1][3]  = -II / v * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[7][3][0]  = V21 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[7][3][1]  = II * V21 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[7][4][0]  = V22 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][4][1]  = II * V22 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][5][0]  = V23 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[7][5][1]  = II * V23 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[8][2][2]  = 0.1e1 / v * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][2][3]  = -II / v * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][3][0]  = V31 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[8][3][1]  = II * V31 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[8][4][0]  = V32 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[8][4][1]  = II * V32 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[8][5][0]  = V33 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[8][5][1]  = II * V33 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[9][0][0] = -0.1e1 / v * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][0][1] = II / v * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][1][0] =
      -0.1e1 / v * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][1][1] = II / v * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][2][0] =
      -0.1e1 / v * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][2][1] = II / v * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][3][2] = 0.1e1 / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[9][3][3] = II / v * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[10][0][0] =
      -0.1e1 / v * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][0][1] = II / v * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][1][0] =
      -0.1e1 / v * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][1][1] = II / v * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][2][0] =
      -0.1e1 / v * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][2][1] = II / v * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][4][2] = 0.1e1 / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[10][4][3] = II / v * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[11][0][0] =
      -0.1e1 / v * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][0][1] = II / v * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][1][0] =
      -0.1e1 / v * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][1][1] = II / v * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][2][0] =
      -0.1e1 / v * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][2][1] = II / v * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][5][2] = 0.1e1 / v * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][5][3] = II / v * SMConstants.C_MassBottom;

  SetCurvatureDone = true;
}

bool Class_VDM::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_VDM::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */

  return false;
}
double Class_VDM::VTreeSimplified(const std::vector<double> &vev) const
{
  (void)vev;
  double res = 0;

  // double vIn = v[0];
  // res = 0.5 * ms * std::pow(vIn, 2) + 1.0 / 24.0 * lambda * std::pow(vIn, 4);

  return res;
}

double Class_VDM::VCounterSimplified(const std::vector<double> &vev) const
{
  (void)vev;
  double res = 0;

  // double vIn = v[0];
  // res = 0.5 * dms * std::pow(vIn, 2) + 1.0 / 24.0 * dlambda * std::pow(vIn,
  // 4) + dT * vIn;

  return res;
}

void Class_VDM::Debugging(const std::vector<double> &input,
                          std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
