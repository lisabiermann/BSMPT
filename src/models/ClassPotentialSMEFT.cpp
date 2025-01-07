// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialSMEFT.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

/**
 * @file
 * Implementation of the SMEFT
 */

namespace BSMPT
{
namespace Models
{

Class_Potential_SMEFT::Class_Potential_SMEFT()
{
  Model         = ModelID::ModelIDs::SMEFT;
  NNeutralHiggs = 2; // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 2; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 2;  // number of parameters in the tree-Level Lagrangian
  nParCT = 24; // number of parameters in the counterterm potential

  nVEV = 1; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  VevOrder[0] = 2;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_SMEFT::~Class_Potential_SMEFT()
{
}

std::vector<std::string> Class_Potential_SMEFT::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dmuSq");
  labels.push_back("dlambda");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  labels.push_back("dZtL");
  labels.push_back("dZbL");
  labels.push_back("dZtR");
  labels.push_back("dZtLv2");
  labels.push_back("dZbLv2");
  labels.push_back("dZtRv2");
  labels.push_back("dZtLyv2");
  labels.push_back("dZbLyv2");
  labels.push_back("dZtRyv2");
  labels.push_back("dZtLv3");
  labels.push_back("dZbLv3");
  labels.push_back("dZtRv3");
  labels.push_back("dZtLy1v3");
  labels.push_back("dZbLy1v3");
  labels.push_back("dZtRy1v3");
  labels.push_back("dZtLy2v3");
  labels.push_back("dZbLy2v3");
  labels.push_back("dZtRy2v3");
  return labels;
}

/**
 * chronological order of the EFT parameters
 */
std::vector<std::string> Class_Potential_SMEFT::addLegendEFT() const
{
  std::vector<std::string> labels;
  labels.push_back("Ouphi");
  return labels;
}

/**
 * numerical values of the EFT parameters
 */
std::vector<double> Class_Potential_SMEFT::getParamsEFT() const
{
  std::vector<double> valsEFT;
  valsEFT.push_back(Ouphi);
  return valsEFT;
}

std::vector<std::string> Class_Potential_SMEFT::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");
  labels.push_back("v_c");
  labels.push_back("omega_c/T_c");
  labels.push_back("omega_c");
  return labels;
}

std::vector<std::string> Class_Potential_SMEFT::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);

  particles[0] = "G+";
  particles[1] = "G-";
  particles[2] = "G0";
  particles[3] = "H";

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

std::vector<std::string> Class_Potential_SMEFT::addLegendVEV() const
{
  std::vector<std::string> labels;
  labels.push_back("omega");
  return labels;
}

void Class_Potential_SMEFT::ReadAndSet(const std::string &linestr,
                                       std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 3; k++)
  {
    ss >> tmp;

    if (k == 1)
      par[0] = tmp; // muSq
    else if (k == 2)
      par[1] = tmp; // lambda
    else if (k == 3)
      par[2] = tmp; // Ouphi
  }

  set_gen(par);
  return;
}

void Class_Potential_SMEFT::set_gen(const std::vector<double> &par)
{
  v0 = C_vev0;

  // uncomment if you want to load muSq from data
  // muSq = par[0];

  // muSq calculated from SMEFT Higgs Mass directly
  (void)par;
  muSq = -std::pow(C_MassSMHiggs, 0.2e1) / 2;

  lambda = -muSq * std::pow(v0, -0.2e1);

  Ouphi = par[2];

  scale = v0;

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);

  vevTreeMin[0] = v0;

  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

void Class_Potential_SMEFT::set_CT_Pot_Par(const std::vector<double> &par)
{
  dmuSq    = par[0];
  dlambda  = par[1];
  dT1      = par[2];
  dT2      = par[3];
  dT3      = par[4];
  dT4      = par[5];
  dZtL     = par[6];
  dZbL     = par[7];
  dZtR     = par[8];
  dZtLv2   = par[9];
  dZbLv2   = par[10];
  dZtRv2   = par[11];
  dZtLyv2  = par[12];
  dZbLyv2  = par[13];
  dZtRyv2  = par[14];
  dZtLv3   = par[15];
  dZbLv3   = par[16];
  dZtRv3   = par[17];
  dZtLy1v3 = par[18];
  dZbLy1v3 = par[19];
  dZtRy1v3 = par[20];
  dZtLy2v3 = par[21];
  dZbLy2v3 = par[22];
  dZtRy2v3 = par[23];

  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;

  Curvature_Higgs_CT_L2[0][0] = dmuSq;
  Curvature_Higgs_CT_L2[1][1] = dmuSq;
  Curvature_Higgs_CT_L2[2][2] = dmuSq;
  Curvature_Higgs_CT_L2[3][3] = dmuSq;

  Curvature_Higgs_CT_L3[0][0][0] = -0.3e1 * std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[0][0][2] = std::sqrt(0.2e1) * (dZtL + dZtR);
  Curvature_Higgs_CT_L3[0][1][1] = -std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[0][2][0] = std::sqrt(0.2e1) * (dZtL + dZtR);
  Curvature_Higgs_CT_L3[0][2][2] = -std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[0][3][3] = -std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[1][0][1] = -std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[1][1][0] = -std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[1][1][2] = std::sqrt(0.2e1) * (dZtL + dZtR);
  Curvature_Higgs_CT_L3[1][2][1] = std::sqrt(0.2e1) * (dZtL + dZtR);
  Curvature_Higgs_CT_L3[2][0][0] = std::sqrt(0.2e1) * (dZtL + dZtR);
  Curvature_Higgs_CT_L3[2][0][2] = -std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[2][1][1] = std::sqrt(0.2e1) * (dZtL + dZtR);
  Curvature_Higgs_CT_L3[2][2][0] = -std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[2][2][2] = 0.3e1 * std::sqrt(0.2e1) * (dZtL + dZtR);
  Curvature_Higgs_CT_L3[2][3][3] = std::sqrt(0.2e1) * (dZtL + dZtR);
  Curvature_Higgs_CT_L3[3][0][3] = -std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[3][2][3] = std::sqrt(0.2e1) * (dZtL + dZtR);
  Curvature_Higgs_CT_L3[3][3][0] = -std::sqrt(0.2e1) * (dZbL + dZtR);
  Curvature_Higgs_CT_L3[3][3][2] = std::sqrt(0.2e1) * (dZtL + dZtR);

  sym3Dim(Curvature_Higgs_CT_L3, NHiggs, NHiggs, NHiggs);

  Curvature_Higgs_CT_L4[0][0][0][0] =
      6 * dlambda - 24 * dZbLv2 - 24 * dZtLyv2 - 24 * dZtRv2 - 24 * dZtRyv2;
  Curvature_Higgs_CT_L4[0][0][0][2] =
      -6 * dZbLv2 - 6 * dZbLyv2 + 6 * dZtLv2 + 6 * dZtLyv2;
  Curvature_Higgs_CT_L4[0][0][1][1] =
      2 * dlambda - 4 * dZbLv2 - 4 * dZtLyv2 - 4 * dZtRv2 - 4 * dZtRyv2;
  Curvature_Higgs_CT_L4[0][0][2][2] =
      2 * dlambda + 4 * dZbLyv2 + 4 * dZtLv2 - 4 * dZbLv2 - 4 * dZtLyv2;
  Curvature_Higgs_CT_L4[0][0][3][3] =
      2 * dlambda - 4 * dZbLv2 - 4 * dZtLyv2 - 4 * dZtRv2 - 4 * dZtRyv2;
  Curvature_Higgs_CT_L4[0][1][1][2] =
      -2 * dZbLv2 + 2 * dZtLv2 - 2 * dZbLyv2 + 2 * dZtLyv2;
  Curvature_Higgs_CT_L4[0][2][2][2] =
      -6 * dZbLv2 - 6 * dZbLyv2 + 6 * dZtLv2 + 6 * dZtLyv2;
  Curvature_Higgs_CT_L4[0][2][3][3] =
      -2 * dZbLv2 + 2 * dZtLv2 - 2 * dZbLyv2 + 2 * dZtLyv2;
  Curvature_Higgs_CT_L4[1][1][1][1] = 6 * dlambda;
  Curvature_Higgs_CT_L4[1][1][2][2] =
      2 * dlambda + 4 * dZbLyv2 + 4 * dZtLv2 + 4 * dZtRv2 + 4 * dZtRyv2;
  Curvature_Higgs_CT_L4[1][1][3][3] = 2 * dlambda;
  Curvature_Higgs_CT_L4[2][2][2][2] =
      6 * dlambda + 24 * dZbLyv2 + 24 * dZtLv2 + 24 * dZtRv2 + 24 * dZtRyv2;
  Curvature_Higgs_CT_L4[2][2][3][3] =
      2 * dlambda + 4 * dZbLyv2 + 4 * dZtLv2 + 4 * dZtRv2 + 4 * dZtRyv2;
  Curvature_Higgs_CT_L4[3][3][3][3] = 6 * dlambda;

  sym4Dim(Curvature_Higgs_CT_L4, NHiggs, NHiggs, NHiggs, NHiggs);

  Curvature_Higgs_CT_L5[0][0][0][0][0] =
      -0.120e3 * std::sqrt(0.2e1) *
      (dZbLv3 + dZtLy1v3 + dZtLy2v3 + dZtRv3 + dZtRy1v3 + dZtRy2v3);
  Curvature_Higgs_CT_L5[0][0][0][0][2] =
      -0.24e2 * std::sqrt(0.2e1) *
      (double)(2 * dZbLv3 + dZtRv3 + dZtRy1v3 + dZtRy2v3 + dZbLy1v3 + dZbLy2v3 -
               dZtLv3);
  Curvature_Higgs_CT_L5[0][0][0][1][1] =
      -0.12e2 * std::sqrt(0.2e1) *
      (dZbLv3 + dZtLy1v3 + dZtLy2v3 + dZtRv3 + dZtRy1v3 + dZtRy2v3);
  Curvature_Higgs_CT_L5[0][0][0][2][2] =
      -0.24e2 * std::sqrt(0.2e1) * (dZbLv3 - dZtLv3);
  Curvature_Higgs_CT_L5[0][0][0][3][3] =
      -0.12e2 * std::sqrt(0.2e1) *
      (dZbLv3 + dZtLy1v3 + dZtLy2v3 + dZtRv3 + dZtRy1v3 + dZtRy2v3);
  Curvature_Higgs_CT_L5[0][0][1][1][2] =
      std::sqrt(0.2e1) *
      (double)(-8 * dZbLv3 - 4 * dZbLy1v3 - 4 * dZbLy2v3 + 4 * dZtLv3 -
               4 * dZtRv3 - 4 * dZtRy1v3 - 4 * dZtRy2v3);
  Curvature_Higgs_CT_L5[0][0][2][2][2] =
      -0.24e2 * std::sqrt(0.2e1) * (dZbLv3 - dZtLv3);
  Curvature_Higgs_CT_L5[0][0][2][3][3] =
      std::sqrt(0.2e1) *
      (double)(-8 * dZbLv3 - 4 * dZbLy1v3 - 4 * dZbLy2v3 + 4 * dZtLv3 -
               4 * dZtRv3 - 4 * dZtRy1v3 - 4 * dZtRy2v3);
  Curvature_Higgs_CT_L5[0][1][1][2][2] =
      -0.4e1 * std::sqrt(0.2e1) *
      (double)(-2 * dZtLv3 - dZtRv3 - dZtRy1v3 - dZtRy2v3 + dZbLv3 - dZtLy1v3 -
               dZtLy2v3);
  Curvature_Higgs_CT_L5[0][2][2][2][2] =
      -0.24e2 * std::sqrt(0.2e1) *
      (double)(-2 * dZtLv3 - dZtRv3 - dZtRy1v3 - dZtRy2v3 + dZbLv3 - dZtLy1v3 -
               dZtLy2v3);
  Curvature_Higgs_CT_L5[0][2][2][3][3] =
      -0.4e1 * std::sqrt(0.2e1) *
      (double)(-2 * dZtLv3 - dZtRv3 - dZtRy1v3 - dZtRy2v3 + dZbLv3 - dZtLy1v3 -
               dZtLy2v3);
  Curvature_Higgs_CT_L5[1][1][2][2][2] =
      0.12e2 * std::sqrt(0.2e1) *
      (dZbLy1v3 + dZbLy2v3 + dZtRy1v3 + dZtRy2v3 + dZtLv3 + dZtRv3);
  Curvature_Higgs_CT_L5[2][2][2][2][2] =
      0.120e3 * std::sqrt(0.2e1) *
      (dZbLy1v3 + dZbLy2v3 + dZtRy1v3 + dZtRy2v3 + dZtLv3 + dZtRv3);
  Curvature_Higgs_CT_L5[2][2][2][3][3] =
      0.12e2 * std::sqrt(0.2e1) *
      (dZbLy1v3 + dZbLy2v3 + dZtRy1v3 + dZtRy2v3 + dZtLv3 + dZtRv3);

  sym5Dim(Curvature_Higgs_CT_L5, NHiggs, NHiggs, NHiggs, NHiggs, NHiggs);

  return;
}

void Class_Potential_SMEFT::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "The parameters are : "
     << "\n"
     << "\tmuSq = " << muSq << "\n"
     << "\tlambda = " << lambda << "\n"
     << "\tv0 = " << v0 << "\n"
     << "\tOuphi = " << Ouphi << "\n";

  ss << "The counterterm parameters are : "
     << "\n";
  ss << "\tdmuSq = " << dmuSq << "\n"
     << "\tdlambda = " << dlambda << "\n"
     << "\tdT1 = " << dT1 << "\n"
     << "\tdT2 = " << dT2 << "\n"
     << "\tdT3 = " << dT3 << "\n"
     << "\tdT4 = " << dT4 << "\n"
     << "\tdZtL = " << dZtL << "\n"
     << "\tdZbL = " << dZbL << "\n"
     << "\tdZtR = " << dZtR << "\n"
     << "\tdZtLv2 = " << dZtLv2 << "\n"
     << "\tdZbLv2 = " << dZbLv2 << "\n"
     << "\tdZtRv2 = " << dZtRv2 << "\n"
     << "\tdZtLyv2 = " << dZtLyv2 << "\n"
     << "\tdZbLyv2 = " << dZbLyv2 << "\n"
     << "\tdZtRyv2 = " << dZtRyv2 << "\n"
     << "\tdZtLv3 = " << dZtLv3 << "\n"
     << "\tdZbLv3 = " << dZbLv3 << "\n"
     << "\tdZtRv3 = " << dZtRv3 << "\n"
     << "\tdZtLy1v3 = " << dZtLy1v3 << "\n"
     << "\tdZbLy1v3 = " << dZbLy1v3 << "\n"
     << "\tdZtRy1v3 = " << dZtRy1v3 << "\n"
     << "\tdZtLy2v3 = " << dZtLy2v3 << "\n"
     << "\tdZbLy2v3 = " << dZbLy2v3 << "\n"
     << "\tdZtRy2v3 = " << dZtRy2v3 << "\n";

  ss << "The scale is given by mu = " << scale << " GeV "
     << "\n";

  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  ss << "The mass spectrum is given by :\n";
  ss << "m_{G^+}^2 = " << HiggsMasses[0] << " GeV^2 \n"
     << "m_{G^-}^2 = " << HiggsMasses[1] << " GeV^2 \n"
     << "m_{G^0}^2 = " << HiggsMasses[2] << " GeV^2 \n"
     << "m_{H_SMEFT} = " << std::sqrt(HiggsMasses[3]) << " GeV \n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

std::vector<double> Class_Potential_SMEFT::calc_CT() const
{
  std::vector<double> parCT;

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

  parCT.push_back((-2 * HesseWeinberg(0, 2) * v0 + 6 * NablaWeinberg(0)) / v0 +
                  0.3e1 / 0.2e1 * HesseWeinberg(0, 0) - HesseWeinberg(0, 2) +
                  HesseWeinberg(2, 2) / 2 - 3 * HesseWeinberg(3, 3)); // dmuSq
  parCT.push_back((2 * HesseWeinberg(0, 2) * v0 - HesseWeinberg(0, 0) * v0 -
                   HesseWeinberg(2, 2) * v0 - 2 * HesseWeinberg(3, 3) * v0 -
                   4 * NablaWeinberg(0) + 4 * NablaWeinberg(2)) *
                  std::pow(v0, -3) / 2); // dlambda
  parCT.push_back(0);                    // dT1
  parCT.push_back(-NablaWeinberg(1));    // dT2
  parCT.push_back(0);                    // dT3
  parCT.push_back(-NablaWeinberg(3));    // dT4
  parCT.push_back(-std::sqrt(2) *
                  (HesseWeinberg(0, 0) * v0 - 3 * HesseWeinberg(0, 2) * v0 -
                   2 * HesseWeinberg(3, 3) * v0 + 7 * NablaWeinberg(0) +
                   NablaWeinberg(2)) *
                  std::pow(v0, -2)); // dZtL
  parCT.push_back(0);                // dZbL
  parCT.push_back(-std::sqrt(2) *
                  (HesseWeinberg(0, 2) * v0 - 3 * NablaWeinberg(0)) *
                  std::pow(v0, -2)); // dZtR
  parCT.push_back((-HesseWeinberg(0, 2) * v0 + 2 * NablaWeinberg(0)) *
                  std::pow(v0, -3)); // dZtLv2
  parCT.push_back(0);                // dZbLv2
  parCT.push_back(std::pow(v0, -2) *
                  (HesseWeinberg(0, 0) - HesseWeinberg(3, 3)) / 2); // dZtRv2
  parCT.push_back(0);                                               // dZtLyv2
  parCT.push_back(0);                                               // dZbLyv2
  parCT.push_back(0);                                               // dZtRyv2
  parCT.push_back(0);                                               // dZtLv3
  parCT.push_back(0);                                               // dZbLv3
  parCT.push_back(0);                                               // dZtRv3
  parCT.push_back(0);                                               // dZtLy1v3
  parCT.push_back(0);                                               // dZbLy1v3
  parCT.push_back(0);                                               // dZtRy1v3
  parCT.push_back(0);                                               // dZtLy2v3
  parCT.push_back(0);                                               // dZbLy2v3
  parCT.push_back(0);                                               // dZtRy2v3

  return parCT;
}

void Class_Potential_SMEFT::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  std::vector<double> HiggsOrder(NHiggs);

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

void Class_Potential_SMEFT::SetCurvatureArrays()
{
  initVectors();

  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  Curvature_Higgs_L2[0][0] = muSq;
  Curvature_Higgs_L2[1][1] = muSq;
  Curvature_Higgs_L2[2][2] = muSq;
  Curvature_Higgs_L2[3][3] = muSq;

  Curvature_Higgs_L4[0][0][0][0] = 6 * lambda;
  Curvature_Higgs_L4[0][0][1][1] = 2 * lambda;
  Curvature_Higgs_L4[0][0][2][2] = 2 * lambda;
  Curvature_Higgs_L4[0][0][3][3] = 2 * lambda;
  Curvature_Higgs_L4[1][1][1][1] = 6 * lambda;
  Curvature_Higgs_L4[1][1][2][2] = 2 * lambda;
  Curvature_Higgs_L4[1][1][3][3] = 2 * lambda;
  Curvature_Higgs_L4[2][2][2][2] = 6 * lambda;
  Curvature_Higgs_L4[2][2][3][3] = 2 * lambda;
  Curvature_Higgs_L4[3][3][3][3] = 6 * lambda;

  sym4Dim(Curvature_Higgs_L4, NHiggs, NHiggs, NHiggs, NHiggs);

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

  MatrixXcd YIJQc1(NQuarks, NQuarks), YIJQc2(NQuarks, NQuarks),
      YIJQc2OI(NQuarks, NQuarks), YIJQg0(NQuarks, NQuarks),
      YIJQg0OI(NQuarks, NQuarks), YIJQh1(NQuarks, NQuarks),
      YIJQh2(NQuarks, NQuarks), YIJQh3(NQuarks, NQuarks);

  std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
  V11 = C_Vud;
  V12 = C_Vus;
  V13 = C_Vub;
  V21 = C_Vcd;
  V22 = C_Vcs;
  V23 = C_Vcb;
  V31 = C_Vtd;
  V32 = C_Vts;
  V33 = C_Vtb;

  std::complex<double> VC11, VC12, VC13, VC21, VC22, VC23, VC31, VC32, VC33;
  VC11 = std::conj(C_Vud);
  VC12 = std::conj(C_Vus);
  VC13 = std::conj(C_Vub);
  VC21 = std::conj(C_Vcd);
  VC22 = std::conj(C_Vcs);
  VC23 = std::conj(C_Vcb);
  VC31 = std::conj(C_Vtd);
  VC32 = std::conj(C_Vts);
  VC33 = std::conj(C_Vtb);

  Curvature_Lepton_F2H1[0][1][2] = 0.1e1 / v0 * C_MassElectron;
  Curvature_Lepton_F2H1[0][1][3] = II / v0 * C_MassElectron;
  Curvature_Lepton_F2H1[1][0][2] = 0.1e1 / v0 * C_MassElectron;
  Curvature_Lepton_F2H1[1][0][3] = II / v0 * C_MassElectron;
  Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / v0 * C_MassElectron;
  Curvature_Lepton_F2H1[1][6][1] = II / v0 * C_MassElectron;
  Curvature_Lepton_F2H1[2][3][2] = 0.1e1 / v0 * C_MassMu;
  Curvature_Lepton_F2H1[2][3][3] = II / v0 * C_MassMu;
  Curvature_Lepton_F2H1[3][2][2] = 0.1e1 / v0 * C_MassMu;
  Curvature_Lepton_F2H1[3][2][3] = II / v0 * C_MassMu;
  Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / v0 * C_MassMu;
  Curvature_Lepton_F2H1[3][7][1] = II / v0 * C_MassMu;
  Curvature_Lepton_F2H1[4][5][2] = 0.1e1 / v0 * C_MassTau;
  Curvature_Lepton_F2H1[4][5][3] = II / v0 * C_MassTau;
  Curvature_Lepton_F2H1[5][4][2] = 0.1e1 / v0 * C_MassTau;
  Curvature_Lepton_F2H1[5][4][3] = II / v0 * C_MassTau;
  Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / v0 * C_MassTau;
  Curvature_Lepton_F2H1[5][8][1] = II / v0 * C_MassTau;
  Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / v0 * C_MassElectron;
  Curvature_Lepton_F2H1[6][1][1] = II / v0 * C_MassElectron;
  Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / v0 * C_MassMu;
  Curvature_Lepton_F2H1[7][3][1] = II / v0 * C_MassMu;
  Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / v0 * C_MassTau;
  Curvature_Lepton_F2H1[8][5][1] = II / v0 * C_MassTau;

  Curvature_Quark_F2H1[0][6][2]  = 0.1e1 / v0 * C_MassUp;
  Curvature_Quark_F2H1[0][6][3]  = -II / v0 * C_MassUp;
  Curvature_Quark_F2H1[0][9][0]  = -0.1e1 / v0 * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][9][1]  = II / v0 * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][10][0] = -0.1e1 / v0 * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][10][1] = II / v0 * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][11][0] = -0.1e1 / v0 * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[0][11][1] = II / v0 * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[1][7][2]  = 0.1e1 / v0 * C_MassCharm;
  Curvature_Quark_F2H1[1][7][3]  = -II / v0 * C_MassCharm;
  Curvature_Quark_F2H1[1][9][0]  = -0.1e1 / v0 * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][9][1]  = II / v0 * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][10][0] = -0.1e1 / v0 * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][10][1] = II / v0 * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][11][0] = -0.1e1 / v0 * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[1][11][1] = II / v0 * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[2][8][2]  = 0.1e1 / v0 * C_MassTop;
  Curvature_Quark_F2H1[2][8][3]  = -II / v0 * C_MassTop;
  Curvature_Quark_F2H1[2][9][0]  = -0.1e1 / v0 * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][9][1]  = II / v0 * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][10][0] = -0.1e1 / v0 * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][10][1] = II / v0 * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][11][0] = -0.1e1 / v0 * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[2][11][1] = II / v0 * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[3][6][0]  = 0.1e1 / v0 * C_MassDown * V11;
  Curvature_Quark_F2H1[3][6][1]  = II / v0 * C_MassDown * V11;
  Curvature_Quark_F2H1[3][7][0]  = V21 / v0 * C_MassDown;
  Curvature_Quark_F2H1[3][7][1]  = II * V21 / v0 * C_MassDown;
  Curvature_Quark_F2H1[3][8][0]  = 0.1e1 / v0 * C_MassDown * V31;
  Curvature_Quark_F2H1[3][8][1]  = II / v0 * C_MassDown * V31;
  Curvature_Quark_F2H1[3][9][2]  = 0.1e1 / v0 * C_MassDown;
  Curvature_Quark_F2H1[3][9][3]  = II / v0 * C_MassDown;
  Curvature_Quark_F2H1[4][6][0]  = 0.1e1 / v0 * C_MassStrange * V12;
  Curvature_Quark_F2H1[4][6][1]  = II / v0 * C_MassStrange * V12;
  Curvature_Quark_F2H1[4][7][0]  = V22 / v0 * C_MassStrange;
  Curvature_Quark_F2H1[4][7][1]  = II * V22 / v0 * C_MassStrange;
  Curvature_Quark_F2H1[4][8][0]  = 0.1e1 / v0 * C_MassStrange * V32;
  Curvature_Quark_F2H1[4][8][1]  = II / v0 * C_MassStrange * V32;
  Curvature_Quark_F2H1[4][10][2] = 0.1e1 / v0 * C_MassStrange;
  Curvature_Quark_F2H1[4][10][3] = II / v0 * C_MassStrange;
  Curvature_Quark_F2H1[5][6][0]  = V13 / v0 * C_MassBottom;
  Curvature_Quark_F2H1[5][6][1]  = II / v0 * C_MassBottom * V13;
  Curvature_Quark_F2H1[5][7][0]  = V23 / v0 * C_MassBottom;
  Curvature_Quark_F2H1[5][7][1]  = II / v0 * C_MassBottom * V23;
  Curvature_Quark_F2H1[5][8][0]  = V33 / v0 * C_MassBottom;
  Curvature_Quark_F2H1[5][8][1]  = II / v0 * C_MassBottom * V33;
  Curvature_Quark_F2H1[5][11][2] = 0.1e1 / v0 * C_MassBottom;
  Curvature_Quark_F2H1[5][11][3] = II / v0 * C_MassBottom;
  Curvature_Quark_F2H1[6][0][2]  = 0.1e1 / v0 * C_MassUp;
  Curvature_Quark_F2H1[6][0][3]  = -II / v0 * C_MassUp;
  Curvature_Quark_F2H1[6][3][0]  = 0.1e1 / v0 * C_MassDown * V11;
  Curvature_Quark_F2H1[6][3][1]  = II / v0 * C_MassDown * V11;
  Curvature_Quark_F2H1[6][4][0]  = 0.1e1 / v0 * C_MassStrange * V12;
  Curvature_Quark_F2H1[6][4][1]  = II / v0 * C_MassStrange * V12;
  Curvature_Quark_F2H1[6][5][0]  = V13 / v0 * C_MassBottom;
  Curvature_Quark_F2H1[6][5][1]  = II / v0 * C_MassBottom * V13;
  Curvature_Quark_F2H1[7][1][2]  = 0.1e1 / v0 * C_MassCharm;
  Curvature_Quark_F2H1[7][1][3]  = -II / v0 * C_MassCharm;
  Curvature_Quark_F2H1[7][3][0]  = V21 / v0 * C_MassDown;
  Curvature_Quark_F2H1[7][3][1]  = II * V21 / v0 * C_MassDown;
  Curvature_Quark_F2H1[7][4][0]  = V22 / v0 * C_MassStrange;
  Curvature_Quark_F2H1[7][4][1]  = II * V22 / v0 * C_MassStrange;
  Curvature_Quark_F2H1[7][5][0]  = V23 / v0 * C_MassBottom;
  Curvature_Quark_F2H1[7][5][1]  = II / v0 * C_MassBottom * V23;
  Curvature_Quark_F2H1[8][2][2]  = 0.1e1 / v0 * C_MassTop;
  Curvature_Quark_F2H1[8][2][3]  = -II / v0 * C_MassTop;
  Curvature_Quark_F2H1[8][3][0]  = 0.1e1 / v0 * C_MassDown * V31;
  Curvature_Quark_F2H1[8][3][1]  = II / v0 * C_MassDown * V31;
  Curvature_Quark_F2H1[8][4][0]  = 0.1e1 / v0 * C_MassStrange * V32;
  Curvature_Quark_F2H1[8][4][1]  = II / v0 * C_MassStrange * V32;
  Curvature_Quark_F2H1[8][5][0]  = V33 / v0 * C_MassBottom;
  Curvature_Quark_F2H1[8][5][1]  = II / v0 * C_MassBottom * V33;
  Curvature_Quark_F2H1[9][0][0]  = -0.1e1 / v0 * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][0][1]  = II / v0 * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][1][0]  = -0.1e1 / v0 * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][1][1]  = II / v0 * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][2][0]  = -0.1e1 / v0 * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][2][1]  = II / v0 * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][3][2]  = 0.1e1 / v0 * C_MassDown;
  Curvature_Quark_F2H1[9][3][3]  = II / v0 * C_MassDown;
  Curvature_Quark_F2H1[10][0][0] = -0.1e1 / v0 * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][0][1] = II / v0 * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][1][0] = -0.1e1 / v0 * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][1][1] = II / v0 * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][2][0] = -0.1e1 / v0 * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][2][1] = II / v0 * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][4][2] = 0.1e1 / v0 * C_MassStrange;
  Curvature_Quark_F2H1[10][4][3] = II / v0 * C_MassStrange;
  Curvature_Quark_F2H1[11][0][0] = -0.1e1 / v0 * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][0][1] = II / v0 * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][1][0] = -0.1e1 / v0 * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][1][1] = II / v0 * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][2][0] = -0.1e1 / v0 * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][2][1] = II / v0 * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][5][2] = 0.1e1 / v0 * C_MassBottom;
  Curvature_Quark_F2H1[11][5][3] = II / v0 * C_MassBottom;

  Curvature_Quark_F2H3[0][8][0][0][2] =
      Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][8][0][0][3] =
      -II * Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][8][1][1][2] =
      Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][8][1][1][3] =
      -II * Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][8][2][2][2] =
      0.3e1 / 0.2e1 * Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Quark_F2H3[0][8][2][2][3] =
      -II * Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][8][2][3][3] =
      Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][8][3][3][3] = -0.3e1 / 0.2e1 * II * Ouphi *
                                        std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1);
  Curvature_Quark_F2H3[0][11][0][0][0] =
      -0.3e1 / 0.2e1 * Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Quark_F2H3[0][11][0][0][1] =
      II * Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][11][0][1][1] =
      -Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][11][0][2][2] =
      -Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][11][0][3][3] =
      -Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][11][1][1][1] = 0.3e1 / 0.2e1 * II * Ouphi *
                                         std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1);
  Curvature_Quark_F2H3[0][11][1][2][2] =
      II * Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[0][11][1][3][3] =
      II * Ouphi * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;

  sym5Dim(Curvature_Quark_F2H3, NQuarks, NQuarks, NHiggs, NHiggs, NHiggs);

  CorrectQuarkTensorsDim6();

  SetCurvatureDone = true;
}

void Class_Potential_SMEFT::PerformVCTShift()
{
  std::vector<double> WeinbergNabla, WeinbergHesse;

  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaCW(NHiggs), NablaCT(NHiggs);
  MatrixXd HesseCW(NHiggs, NHiggs), HesseCT(NHiggs, NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)

  {
    NablaCW[i] = WeinbergNabla[i];
    NablaCT(i) = Curvature_Higgs_CT_L1[i];

    for (std::size_t j = 0; j < NHiggs; j++)

    {
      HesseCW(i, j) = WeinbergHesse.at(j * NHiggs + i);
      NablaCT(i) += Curvature_Higgs_CT_L2[i][j] * vevTree[j];
      HesseCT(i, j) = Curvature_Higgs_CT_L2[i][j];

      for (size_t k = 0; k < NHiggs; k++)

      {
        NablaCT(i) +=
            0.5 * Curvature_Higgs_CT_L3[i][j][k] * vevTree[j] * vevTree[k];
        HesseCT(i, j) += Curvature_Higgs_CT_L3[i][j][k] * vevTree[k];

        for (std::size_t l = 0; l < NHiggs; l++)
        {
          NablaCT(i) += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][l] *
                        vevTree[j] * vevTree[k] * vevTree[l];
          HesseCT(i, j) +=
              0.5 * Curvature_Higgs_CT_L4[i][j][k][l] * vevTree[k] * vevTree[l];
        }
      }
    }
  }

  // correct HesseCT
  double delta;
  for (size_t i = 0; i < NHiggs; i++)
  {
    for (size_t j = 0; j < NHiggs; j++)
    {
      delta = HesseCW(i, j) + HesseCT(i, j);
      if (abs(delta) > 1e-6)
      {
        std::cout << "Correct HesseCT in direction ( " << i << ", " << j
                  << " ) with delta = " << delta << std::endl;
        Curvature_Higgs_CT_L2[i][j] -= delta;
      }
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)

  {
    NablaCW[i] = WeinbergNabla[i];
    NablaCT(i) = Curvature_Higgs_CT_L1[i];

    for (std::size_t j = 0; j < NHiggs; j++)

    {
      HesseCW(i, j) = WeinbergHesse.at(j * NHiggs + i);
      NablaCT(i) += Curvature_Higgs_CT_L2[i][j] * vevTree[j];
      HesseCT(i, j) = Curvature_Higgs_CT_L2[i][j];

      for (size_t k = 0; k < NHiggs; k++)

      {
        NablaCT(i) +=
            0.5 * Curvature_Higgs_CT_L3[i][j][k] * vevTree[j] * vevTree[k];
        HesseCT(i, j) += Curvature_Higgs_CT_L3[i][j][k] * vevTree[k];

        for (std::size_t l = 0; l < NHiggs; l++)
        {
          NablaCT(i) += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][l] *
                        vevTree[j] * vevTree[k] * vevTree[l];
          HesseCT(i, j) +=
              0.5 * Curvature_Higgs_CT_L4[i][j][k][l] * vevTree[k] * vevTree[l];
        }
      }
    }
  }

  // correct NablaCT taking into account redefined L2
  for (size_t i = 0; i < NHiggs; i++)
  {
    delta = NablaCW(i) + NablaCT(i);
    if (abs(delta) > 1e-6)
    {
      std::cout << "Correct NablaCT in direction ( " << i
                << " ) with delta = " << delta << std::endl;
      std::cout << "tadpole " << i << " before correction "
                << Curvature_Higgs_CT_L1[i] << std::endl;
      Curvature_Higgs_CT_L1[i] -= delta;
      std::cout << "tadpole " << i << " after correction "
                << Curvature_Higgs_CT_L1[i] << std::endl;
    }
  }

  // check if corrected tadpoles are zero
  for (size_t i = 0; i < NHiggs; i++)
  {
    if (abs(Curvature_Higgs_CT_L1[i]) > 1e-3)
    {
      std::cout << "Non-zero tadpoles detected with: "
                << Curvature_Higgs_CT_L1[i] << std::endl;
      throw std::runtime_error(
          "Non-zero tadpole CT not allowed in Landau gauge!");
    }
  }
}

double Class_Potential_SMEFT::SymFac_Higgs_OneLoop(
    const int &i,
    const int &j,
    const std::vector<double> &point) const
{
  if (i == 2 and j == 2) // zetazeta
  {
    return (-Ouphi * C_MassTop) /
           (std::sqrt(2) * std::pow(LambdaEFT, 2) * C_vev0) *
           (point[2] * point[2] - C_vev0 * C_vev0); // zetazeta}
  }
  else
  {
    return 0;
  }
}

double Class_Potential_SMEFT::SymFac_Higgs_TwoLoop(const int &i,
                                                   const int &j) const
{
  if (i == 0 and j == 0)
  {
    return (Ouphi * C_MassTop) /
           (4. * std::sqrt(2) * std::pow(LambdaEFT, 2) * v0); // rhorho
  }
  else if (i == 1 and j == 1)
  {
    return (Ouphi * C_MassTop) /
           (4. * std::sqrt(2) * std::pow(LambdaEFT, 2) * v0); // etaeta
  }
  else if (i == 2 and j == 2)
  {
    return (Ouphi * C_MassTop) /
           (4. * std::sqrt(2) * std::pow(LambdaEFT, 2) * v0); // zetazeta
  }
  else if (i == 3 and j == 3)
  {
    return (Ouphi * C_MassTop) /
           (4. * std::sqrt(2) * std::pow(LambdaEFT, 2) * v0); // psipsi
  }
  else
  {
    return 0;
  }
}

bool Class_Potential_SMEFT::CalculateDebyeSimplified()
{
  return false;
}

bool Class_Potential_SMEFT::CalculateDebyeGaugeSimplified()
{
  return false;
}
double
Class_Potential_SMEFT::VTreeSimplified(const std::vector<double> &v) const
{
  double res = 0;
  (void)v;
  // not implemented
  return res;
}

double
Class_Potential_SMEFT::VCounterSimplified(const std::vector<double> &v) const
{
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  (void)v;
  // not implemented
  return res;
}

void Class_Potential_SMEFT::Debugging(const std::vector<double> &input,
                                      std::vector<double> &output) const
{
  (void)input;

  std::cout
      << "\n------------------------------------------------------------"
      << "\nDebugging output - First and second derivative of VTree, VCT, "
         "VCW:\n"
      << std::endl;

  std::vector<double> WeinbergNabla, WeinbergHesse;

  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs), NablaVCT(NHiggs), NablaTree(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HesseVCT(NHiggs, NHiggs),
      HesseTree(NHiggs, NHiggs);
  VectorXd DeltaNabla(NHiggs);
  MatrixXd DeltaHesse(NHiggs, NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)

  {
    NablaWeinberg[i] = WeinbergNabla[i];
    NablaTree(i)     = Curvature_Higgs_L1[i];
    NablaVCT(i)      = Curvature_Higgs_CT_L1[i];

    for (std::size_t j = 0; j < NHiggs; j++)

    {
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
      NablaTree(i) += Curvature_Higgs_L2[i][j] * vevTree[j];
      HesseTree(i, j) = Curvature_Higgs_L2[i][j];
      NablaVCT(i) += Curvature_Higgs_CT_L2[i][j] * vevTree[j];
      HesseVCT(i, j) = Curvature_Higgs_CT_L2[i][j];

      for (size_t k = 0; k < NHiggs; k++)

      {
        NablaTree(i) +=
            0.5 * Curvature_Higgs_L3[i][j][k] * vevTree[j] * vevTree[k];
        HesseTree(i, j) += Curvature_Higgs_L3[i][j][k] * vevTree[k];
        NablaVCT(i) +=
            0.5 * Curvature_Higgs_CT_L3[i][j][k] * vevTree[j] * vevTree[k];
        HesseVCT(i, j) += Curvature_Higgs_CT_L3[i][j][k] * vevTree[k];

        for (std::size_t l = 0; l < NHiggs; l++)
        {
          NablaTree(i) += 1.0 / 6.0 * Curvature_Higgs_L4[i][j][k][l] *
                          vevTree[j] * vevTree[k] * vevTree[l];
          HesseTree(i, j) +=
              0.5 * Curvature_Higgs_L4[i][j][k][l] * vevTree[k] * vevTree[l];
          NablaVCT(i) += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][l] *
                         vevTree[j] * vevTree[k] * vevTree[l];
          HesseVCT(i, j) +=
              0.5 * Curvature_Higgs_CT_L4[i][j][k][l] * vevTree[k] * vevTree[l];
        }
      }
    }
  }

  double thres_coupling = 1e-12;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    double nablaCW   = NablaWeinberg[i];
    double nablaCT   = NablaVCT[i];
    double nablaTree = NablaTree[i];

    if (std::abs(nablaCW) < thres_coupling and
        std::abs(nablaCT) < thres_coupling and
        std::abs(nablaTree) < thres_coupling)
    {
      DeltaNabla[i] = 0;
    }
    else if (std::abs(nablaTree) < thres_coupling and
             std::abs(nablaCT) < thres_coupling and
             std::abs(nablaCW) > thres_coupling)
    {
      DeltaNabla[i] = 1;
    }
    else if (std::abs(nablaTree) < thres_coupling and
             (std::abs(nablaCW + nablaCT) < thres_coupling))
    {
      DeltaNabla[i] = 0;
    }
    else if (std::abs(nablaTree) < thres_coupling and
             (std::abs(nablaCW + nablaCT) > thres_coupling))
    {
      DeltaNabla[i] = nablaCW + nablaCT;
    }
    else
    {
      DeltaNabla[i] = (nablaCW + nablaCT) / std::abs(nablaCW);

      if (std::abs(DeltaNabla[i]) > 0)
      {
        std::cout << std::cout.precision(10) << "field comp.: (" << i
                  << "): rel. diff nabla = " << DeltaNabla[i]
                  << " with NTree = " << NablaTree[i]
                  << " - NCW = " << NablaWeinberg[i]
                  << " - NCT = " << NablaVCT[i] << std::endl;
      }
      else
      {
        DeltaNabla[i] = 0;
      }
    }
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      double hesseCW   = HesseWeinberg(i, j);
      double hesseCT   = HesseVCT(i, j);
      double hesseTree = HesseTree(i, j);

      if (std::abs(hesseCW) < thres_coupling and
          std::abs(hesseCT) < thres_coupling and
          std::abs(hesseTree) < thres_coupling)
      {
        DeltaHesse(i, j) = 0;
      }
      else if (std::abs(hesseTree) < thres_coupling and
               std::abs(hesseCT) < thres_coupling and
               std::abs(hesseCW) > thres_coupling)
      {
        DeltaHesse(i, j) = 1;
      }
      else if (std::abs(hesseTree) < thres_coupling and
               (std::abs(hesseCW + hesseCT) < thres_coupling))
      {
        DeltaHesse(i, j) = 0;
      }
      else if (std::abs(hesseTree) < thres_coupling and
               (std::abs(hesseCW + hesseCT) > thres_coupling))
      {
        DeltaHesse(i, j) = hesseCW + hesseCT;
      }
      else
      {
        DeltaHesse(i, j) = (hesseCW + hesseCT) / std::abs(hesseCW);
        if (std::abs(DeltaHesse(i, j)) > 0)
        {
          std::cout << std::cout.precision(10) << "field comp.: (" << i << " , "
                    << j << "): rel. diff hesse = " << DeltaHesse(i, j)
                    << " with HTree = " << HesseTree(i, j)
                    << " - HCW = " << HesseWeinberg(i, j)
                    << " - HCT = " << HesseVCT(i, j) << std::endl;
        }
        else
        {
          DeltaHesse(i, j) = 0;
        }
      }
    }
  }

  double deltanablamax{0.0}, deltahessemax{0.0};
  double absnablamax{0.0}, abshessemax{0.0};

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    if (std::abs(DeltaNabla[i]) > std::abs(deltanablamax))
    {
      deltanablamax = DeltaNabla[i];
      absnablamax   = std::abs(NablaWeinberg[i] + NablaVCT[i]);
    }
    for (std::size_t j = i; j < NHiggs; j++)
    {
      if (std::abs(DeltaHesse(i, j)) > std::abs(deltahessemax))
      {
        deltahessemax = DeltaHesse(i, j);
        abshessemax   = std::abs(HesseWeinberg(i, j) + HesseVCT(i, j));
      }
    }
  }

  double NablaCWCTnotequal = deltanablamax;
  double HesseCWCTnotequal = deltahessemax;
  double NCWnotZero        = absnablamax;
  double HCWnotZero        = abshessemax;

  std::cout << "Compare couplings between LO and NLO:\n"
            << "max rel diff in first deriv   = " << NablaCWCTnotequal << "\n"
            << "corresponds to abs diff of    = " << NCWnotZero << "\n"
            << "max rel diff in second deriv  = " << HesseCWCTnotequal << "\n"
            << "corresponds to abs diff of    = " << HCWnotZero << "\n"
            << "------------------------------------------------------------"
            << std::endl;

  if (std::abs(NablaCWCTnotequal) + std::abs(HesseCWCTnotequal) > 1e-5)
  {
    std::cout << "First and second derivative of (VCW + VCT) are not zero, "
                 "check implementation of counterterm potential."
              << std::endl;
  }
  else
  {
    std::cout << "First and second derivative of (VCW + VCT) are zero, "
                 "masses and angles are fixed."
              << std::endl;
  }

  (void)output;
}
} // namespace Models
} // namespace BSMPT
