// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialR2HDMEFTPHI6_PSI2PHI3.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>

namespace BSMPT
{
namespace Models
{
Class_Potential_R2HDMEFTPHI6_PSI2PHI3::Class_Potential_R2HDMEFTPHI6_PSI2PHI3()
{
  // TODO Auto-generated constructor stub
  Model         = ModelID::ModelIDs::R2HDMEFTPHI6_PSI2PHI3;
  NNeutralHiggs = 4;
  NChargedHiggs = 4;

  NHiggs  = NNeutralHiggs + NChargedHiggs;
  NGauge  = 4;
  NLepton = 9;
  NQuarks = 12;

  nPar   = 10;
  nParCT = 31;

  nVEV = 4;

  VevOrder.resize(nVEV);
  VevOrder[0] = 2;
  VevOrder[1] = 4;
  VevOrder[2] = 6;
  VevOrder[3] = 7;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_R2HDMEFTPHI6_PSI2PHI3::~Class_Potential_R2HDMEFTPHI6_PSI2PHI3()
{
  // TODO Auto-generated destructor stub
}

/**
 * chronological order of the counterterms
 */
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6_PSI2PHI3::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dm11Sq");
  labels.push_back("dm22Sq");
  labels.push_back("dm12Sq");
  labels.push_back("dL1");
  labels.push_back("dL2");
  labels.push_back("dL3");
  labels.push_back("dL4");
  labels.push_back("dL5");
  labels.push_back("dL6");
  labels.push_back("dL7");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  labels.push_back("dT5");
  labels.push_back("dT6");
  labels.push_back("dT7");
  labels.push_back("dT8");
  labels.push_back("dOp6_111111");
  labels.push_back("dOp6_111122");
  labels.push_back("dOp6_122111");
  labels.push_back("dOp6_121211");
  labels.push_back("dOp6_111112");
  labels.push_back("dOp6_121221");
  labels.push_back("dOp6_112212");
  labels.push_back("dOp6_222222");
  labels.push_back("dOp6_112222");
  labels.push_back("dOp6_122122");
  labels.push_back("dOp6_121222");
  labels.push_back("dOp6_122222");
  labels.push_back("dOp6_121212");

  return labels;
}

/**
 * chronological order of the VEVs and the critical temperature
 */
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6_PSI2PHI3::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");
  labels.push_back("omega_c");
  labels.push_back("omega_c/T_c");

  labels.push_back("omega_CB(T_c)");
  labels.push_back("omega_1(T_c)");
  labels.push_back("omega_2(T_c)");
  labels.push_back("omega_CP(T_c)");
  return labels;
}

/**
 * chronological order of the VEVs
 */
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6_PSI2PHI3::addLegendVEV() const
{
  std::vector<std::string> labels;
  labels.push_back("omega_CB");
  labels.push_back("omega_1");
  labels.push_back("omega_2");
  labels.push_back("omega_CP");
  return labels;
}

/**
 * chronological order of the EFT parameters
 */
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6_PSI2PHI3::addLegendEFT() const
{
  std::vector<std::string> labels;
  labels.push_back("Op6_111111");
  labels.push_back("Op6_111122");
  labels.push_back("Op6_122111");
  labels.push_back("Op6_121211");
  labels.push_back("Op6_111112");
  labels.push_back("Op6_121221");
  labels.push_back("Op6_112212");
  labels.push_back("Op6_222222");
  labels.push_back("Op6_112222");
  labels.push_back("Op6_122122");
  labels.push_back("Op6_121222");
  labels.push_back("Op6_122222");
  labels.push_back("Op6_121212");
  labels.push_back("OQu_1b11b");
  labels.push_back("OQu_1b12b");
  labels.push_back("OQu_1b21b");
  labels.push_back("OQu_1b22b");
  labels.push_back("OQu_2b11b");
  labels.push_back("OQu_2b12b");
  labels.push_back("OQu_2b21b");
  labels.push_back("OQu_2b22b");
  labels.push_back("OQd_1b11b");
  labels.push_back("OQd_1b12b");
  labels.push_back("OQd_1b21b");
  labels.push_back("OQd_1b22b");
  labels.push_back("OQd_2b11b");
  labels.push_back("OQd_2b12b");
  labels.push_back("OQd_2b21b");
  labels.push_back("OQd_2b22b");
  labels.push_back("OL_1b11b");
  labels.push_back("OL_1b12b");
  labels.push_back("OL_1b21b");
  labels.push_back("OL_1b22b");
  labels.push_back("OL_2b11b");
  labels.push_back("OL_2b12b");
  labels.push_back("OL_2b21b");
  labels.push_back("OL_2b22b");
  return labels;
}

/**
 * numerical values of the EFT parameters
 */
std::vector<double> Class_Potential_R2HDMEFTPHI6_PSI2PHI3::getParamsEFT() const
{
  std::vector<double> valsEFT;
  valsEFT.push_back(Op6_111111);
  valsEFT.push_back(Op6_111122);
  valsEFT.push_back(Op6_122111);
  valsEFT.push_back(Op6_121211);
  valsEFT.push_back(Op6_111112);
  valsEFT.push_back(Op6_121221);
  valsEFT.push_back(Op6_112212);
  valsEFT.push_back(Op6_222222);
  valsEFT.push_back(Op6_112222);
  valsEFT.push_back(Op6_122122);
  valsEFT.push_back(Op6_121222);
  valsEFT.push_back(Op6_122222);
  valsEFT.push_back(Op6_121212);
  valsEFT.push_back(OQu_1b11b);
  valsEFT.push_back(OQu_1b12b);
  valsEFT.push_back(OQu_1b21b);
  valsEFT.push_back(OQu_1b22b);
  valsEFT.push_back(OQu_2b11b);
  valsEFT.push_back(OQu_2b12b);
  valsEFT.push_back(OQu_2b21b);
  valsEFT.push_back(OQu_2b22b);
  valsEFT.push_back(OQd_1b11b);
  valsEFT.push_back(OQd_1b12b);
  valsEFT.push_back(OQd_1b21b);
  valsEFT.push_back(OQd_1b22b);
  valsEFT.push_back(OQd_2b11b);
  valsEFT.push_back(OQd_2b12b);
  valsEFT.push_back(OQd_2b21b);
  valsEFT.push_back(OQd_2b22b);
  valsEFT.push_back(OL_1b11b);
  valsEFT.push_back(OL_1b12b);
  valsEFT.push_back(OL_1b21b);
  valsEFT.push_back(OL_1b22b);
  valsEFT.push_back(OL_2b11b);
  valsEFT.push_back(OL_2b12b);
  valsEFT.push_back(OL_2b21b);
  valsEFT.push_back(OL_2b22b);
  return valsEFT;
}

/**
 * chronological order of the triple higgs couplings
 */
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6_PSI2PHI3::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;

  particles.push_back("G^+");
  particles.push_back("G^-");
  particles.push_back("H^+");
  particles.push_back("H^-");
  particles.push_back("G^0");
  particles.push_back("A");
  particles.push_back("h");
  particles.push_back("H");

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

void Class_Potential_R2HDMEFTPHI6_PSI2PHI3::ReadAndSet(
    const std::string &linestr,
    std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;
  double lType = 0, lL1 = 0, lL2 = 0, lL3 = 0, lL4 = 0, lL5 = 0, lm12Sq = 0,
         lTanBeta = 0;
  //   double lm11Sq = 0, lm22Sq = 0;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 47; k++)
  {
    ss >> tmp;
    if (k == 1)
      lType = tmp;
    else if (k == 2)
      lL1 = tmp;
    else if (k == 3)
      lL2 = tmp;
    else if (k == 4)
      lL3 = tmp;
    else if (k == 5)
      lL4 = tmp;
    else if (k == 6)
      lL5 = tmp;
    // else if (k == 7)
    //   lm11Sq = tmp;
    // else if (k == 8)
    //   lm22Sq = tmp;
    else if (k == 9)
      lm12Sq = tmp;
    else if (k == 10)
      lTanBeta = tmp;
    else if (k == 11)
      Op6_111111 = tmp;
    else if (k == 12)
      Op6_111122 = tmp;
    else if (k == 13)
      Op6_122111 = tmp;
    else if (k == 14)
      Op6_121211 = tmp;
    else if (k == 15)
      Op6_111112 = tmp;
    else if (k == 16)
      Op6_121221 = tmp;
    else if (k == 17)
      Op6_112212 = tmp;
    else if (k == 18)
      Op6_222222 = tmp;
    else if (k == 19)
      Op6_112222 = tmp;
    else if (k == 20)
      Op6_122122 = tmp;
    else if (k == 21)
      Op6_121222 = tmp;
    else if (k == 22)
      Op6_122222 = tmp;
    else if (k == 23)
      Op6_121212 = tmp;
    else if (k == 24)
      OQu_1b11b = tmp; // Z2-violating
    else if (k == 25)
      OQu_1b12b = tmp;
    else if (k == 26)
      OQu_1b21b = tmp;
    else if (k == 27)
      OQu_1b22b = tmp; // Z2-violating
    else if (k == 28)
      OQu_2b11b = tmp;
    else if (k == 29)
      OQu_2b12b = tmp; // Z2-violating
    else if (k == 30)
      OQu_2b21b = tmp; // Z2-violating
    else if (k == 31)
      OQu_2b22b = tmp;
    else if (k == 32)
      OQd_1b11b = tmp; // Z2-violating
    else if (k == 33)
      OQd_1b12b = tmp;
    else if (k == 34)
      OQd_1b21b = tmp;
    else if (k == 35)
      OQd_1b22b = tmp; // Z2-violating
    else if (k == 36)
      OQd_2b11b = tmp;
    else if (k == 37)
      OQd_2b12b = tmp; // Z2-violating
    else if (k == 38)
      OQd_2b21b = tmp; // Z2-violating
    else if (k == 39)
      OQd_2b22b = tmp;
    else if (k == 40)
      OL_1b11b = tmp; // Z2-violating
    else if (k == 41)
      OL_1b12b = tmp;
    else if (k == 42)
      OL_1b21b = tmp;
    else if (k == 43)
      OL_1b22b = tmp; // Z2-violating
    else if (k == 44)
      OL_2b11b = tmp;
    else if (k == 45)
      OL_2b12b = tmp; // Z2-violating
    else if (k == 46)
      OL_2b21b = tmp; // Z2-violating
    else if (k == 47)
      OL_2b22b = tmp;
  }

  par[0] = lL1;
  par[1] = lL2;
  par[2] = lL3;
  par[3] = lL4;
  par[4] = lL5;
  par[5] = lm12Sq;
  par[6] = lTanBeta;
  par[7] = lType;

  set_gen(par);
  return;
}

/**
 * Set Class Object with an CP-Conserving Point
 */
void Class_Potential_R2HDMEFTPHI6_PSI2PHI3::set_gen(
    const std::vector<double> &par)
{
  // double *p = (double *)par;
  scale            = C_vev0;
  L1               = par[0];
  L2               = par[1];
  L3               = par[2];
  L4               = par[3];
  L5               = par[4];
  m12Sq            = par[5];
  TanBeta          = par[6];
  beta             = std::atan(TanBeta);
  Type             = static_cast<int>(par[7]);
  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = std::sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = std::sqrt(C_SinBetaSquared);

  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;

  m11Sq =
      -(-(double)(4 * m12Sq * v2 * LambdaEFT * LambdaEFT) +
        0.4e1 * L1 * std::pow(v1, 0.3e1) * (double)(LambdaEFT * LambdaEFT) +
        0.2e1 * L3 * v1 * (double)(v2 * v2) * (double)(LambdaEFT * LambdaEFT) +
        0.2e1 * L4 * (double)(v2 * v2) * v1 * (double)(LambdaEFT * LambdaEFT) +
        0.2e1 * L5 * (double)(v2 * v2) * v1 * (double)(LambdaEFT * LambdaEFT) +
        0.6e1 * L6 * v1 * v1 * (double)v2 * (double)(LambdaEFT * LambdaEFT) +
        (double)(2 * L7 * std::pow((double)v2, (double)3) * LambdaEFT *
                 LambdaEFT) -
        0.3e1 * Op6_111111 * std::pow(v1, 0.5e1) -
        0.5e1 * Op6_111112 * std::pow(v1, 0.4e1) * (double)v2 -
        0.2e1 * Op6_111122 * std::pow(v1, 0.3e1) * (double)(v2 * v2) -
        0.4e1 * Op6_121211 * std::pow(v1, 0.3e1) * (double)(v2 * v2) -
        0.2e1 * Op6_122111 * (double)(v2 * v2) * std::pow(v1, 0.3e1) -
        0.3e1 * Op6_112212 * v1 * v1 * std::pow((double)v2, (double)3) -
        0.3e1 * Op6_121212 * std::pow((double)v2, (double)3) * v1 * v1 -
        0.3e1 * Op6_121221 * std::pow((double)v2, (double)3) * v1 * v1 -
        Op6_112222 * v1 * std::pow((double)v2, (double)4) -
        0.2e1 * Op6_121222 * std::pow((double)v2, (double)4) * v1 -
        Op6_122122 * std::pow((double)v2, (double)4) * v1 -
        (double)(Op6_122222 * std::pow((double)v2, (double)5))) /
      v1 * std::pow((double)LambdaEFT, (double)(-2)) / 0.4e1;

  m22Sq =
      -(-(double)(4 * m12Sq * v1 * LambdaEFT * LambdaEFT) +
        0.4e1 * L2 * std::pow(v2, 0.3e1) * (double)(LambdaEFT * LambdaEFT) +
        0.2e1 * L3 * (double)(v1 * v1) * v2 * (double)(LambdaEFT * LambdaEFT) +
        0.2e1 * L4 * (double)(v1 * v1) * v2 * (double)(LambdaEFT * LambdaEFT) +
        0.2e1 * L5 * v2 * (double)(v1 * v1) * (double)(LambdaEFT * LambdaEFT) +
        0.6e1 * L7 * v2 * v2 * (double)v1 * (double)(LambdaEFT * LambdaEFT) +
        (double)(2 * L6 * std::pow((double)v1, (double)3) * LambdaEFT *
                 LambdaEFT) -
        (double)(Op6_111112 * std::pow((double)v1, (double)5)) -
        Op6_111122 * std::pow((double)v1, (double)4) * v2 -
        0.2e1 * Op6_121211 * std::pow((double)v1, (double)4) * v2 -
        Op6_122111 * std::pow((double)v1, (double)4) * v2 -
        0.3e1 * Op6_112212 * std::pow((double)v1, (double)3) * v2 * v2 -
        0.3e1 * Op6_121212 * v2 * v2 * std::pow((double)v1, (double)3) -
        0.3e1 * Op6_121221 * v2 * v2 * std::pow((double)v1, (double)3) -
        0.2e1 * Op6_112222 * (double)(v1 * v1) * std::pow(v2, 0.3e1) -
        0.4e1 * Op6_121222 * std::pow(v2, 0.3e1) * (double)(v1 * v1) -
        0.2e1 * Op6_122122 * (double)(v1 * v1) * std::pow(v2, 0.3e1) -
        0.5e1 * Op6_122222 * (double)v1 * std::pow(v2, 0.4e1) -
        0.3e1 * Op6_222222 * std::pow(v2, 0.5e1)) /
      v2 * std::pow((double)LambdaEFT, (double)(-2)) / 0.4e1;

  double cb = 0;

  if (Type == 1 or Type == 3) // Type I 2HDM oder Lepton Specific
  {
    cb = std::sqrt(2) * C_MassBottom / (C_vev0 * C_SinBeta);
  }
  if (Type == 2 or Type == 4) // Type II 2HDM oder Flipped
  {
    cb = std::sqrt(2) * C_MassBottom / (C_vev0 * C_CosBeta);
  }
  CTempC1 = 1.0 / 48 *
            (12 * L1 + 8 * L3 + 4 * L4 + 3 * (3 * C_g * C_g + C_gs * C_gs));
  double ct = std::sqrt(2) * C_MassTop / (C_vev0 * C_SinBeta);
  CTempC2   = 1.0 / 48 *
            (12 * L2 + 8 * L3 + 4 * L4 + 3 * (3 * C_g * C_g + C_gs * C_gs) +
             12 * ct * ct);

  if (Type == 1 or Type == 3)
  {
    CTempC2 += 12.0 / 48.0 * cb * cb;
  }
  else
  {
    CTempC1 += 12.0 / 48.0 * cb * cb;
  }

  vevTreeMin.resize(nVEV);
  vevTreeMin[0] = 0;
  vevTreeMin[1] = C_vev0 * C_CosBeta;
  vevTreeMin[2] = C_vev0 * C_SinBeta;
  vevTreeMin[3] = 0;
  vevTree.resize(NHiggs);
  vevTree = MinimizeOrderVEV(vevTreeMin);
}

void Class_Potential_R2HDMEFTPHI6_PSI2PHI3::set_CT_Pot_Par(
    const std::vector<double> &p)
{
  //	double *p = (double *)par;
  dm11Sq      = p[0];
  dm22Sq      = p[1];
  dm12Sq      = p[2];
  dL1         = p[3];
  dL2         = p[4];
  dL3         = p[5];
  dL4         = p[6];
  dL5         = p[7];
  dL6         = p[8];
  dL7         = p[9];
  dT1         = p[10];
  dT2         = p[11];
  dT3         = p[12];
  dT4         = p[13];
  dT5         = p[14];
  dT6         = p[15];
  dT7         = p[16];
  dT8         = p[17];
  dOp6_111111 = p[18];
  dOp6_111122 = p[19];
  dOp6_122111 = p[20];
  dOp6_121211 = p[21];
  dOp6_111112 = p[22];
  dOp6_121221 = p[23];
  dOp6_112212 = p[24];
  dOp6_222222 = p[25];
  dOp6_112222 = p[26];
  dOp6_122122 = p[27];
  dOp6_121222 = p[28];
  dOp6_122222 = p[29];
  dOp6_121212 = p[30];

  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;
  Curvature_Higgs_CT_L1[4] = dT5;
  Curvature_Higgs_CT_L1[5] = dT6;
  Curvature_Higgs_CT_L1[6] = dT7;
  Curvature_Higgs_CT_L1[7] = dT8;

  Curvature_Higgs_CT_L2[0][0] = dm11Sq;
  Curvature_Higgs_CT_L2[0][2] = -dm12Sq;
  Curvature_Higgs_CT_L2[1][1] = dm11Sq;
  Curvature_Higgs_CT_L2[1][3] = -dm12Sq;
  Curvature_Higgs_CT_L2[2][2] = dm22Sq;
  Curvature_Higgs_CT_L2[3][3] = dm22Sq;
  Curvature_Higgs_CT_L2[4][4] = dm11Sq;
  Curvature_Higgs_CT_L2[4][6] = -dm12Sq;
  Curvature_Higgs_CT_L2[5][5] = dm11Sq;
  Curvature_Higgs_CT_L2[5][7] = -dm12Sq;
  Curvature_Higgs_CT_L2[6][6] = dm22Sq;
  Curvature_Higgs_CT_L2[7][7] = dm22Sq;

  sym2Dim(Curvature_Higgs_CT_L2, NHiggs, NHiggs);

  Curvature_Higgs_CT_L4[0][0][0][0] = 6 * dL1;
  Curvature_Higgs_CT_L4[0][0][0][2] = 3 * dL6;
  Curvature_Higgs_CT_L4[0][0][1][1] = 2 * dL1;
  Curvature_Higgs_CT_L4[0][0][1][3] = dL6;
  Curvature_Higgs_CT_L4[0][0][2][2] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[0][0][3][3] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[0][0][4][4] = 2 * dL1;
  Curvature_Higgs_CT_L4[0][0][4][6] = dL6;
  Curvature_Higgs_CT_L4[0][0][5][5] = 2 * dL1;
  Curvature_Higgs_CT_L4[0][0][5][7] = dL6;
  Curvature_Higgs_CT_L4[0][0][6][6] = dL3;
  Curvature_Higgs_CT_L4[0][0][7][7] = dL3;
  Curvature_Higgs_CT_L4[0][1][1][2] = dL6;
  Curvature_Higgs_CT_L4[0][1][2][3] = dL5;
  Curvature_Higgs_CT_L4[0][2][2][2] = 3 * dL7;
  Curvature_Higgs_CT_L4[0][2][3][3] = dL7;
  Curvature_Higgs_CT_L4[0][2][4][4] = dL6;
  Curvature_Higgs_CT_L4[0][2][4][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][5][5] = dL6;
  Curvature_Higgs_CT_L4[0][2][5][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][6][6] = dL7;
  Curvature_Higgs_CT_L4[0][2][7][7] = dL7;
  Curvature_Higgs_CT_L4[0][3][4][7] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][5][6] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][1][1] = 6 * dL1;
  Curvature_Higgs_CT_L4[1][1][1][3] = 3 * dL6;
  Curvature_Higgs_CT_L4[1][1][2][2] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[1][1][3][3] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[1][1][4][4] = 2 * dL1;
  Curvature_Higgs_CT_L4[1][1][4][6] = dL6;
  Curvature_Higgs_CT_L4[1][1][5][5] = 2 * dL1;
  Curvature_Higgs_CT_L4[1][1][5][7] = dL6;
  Curvature_Higgs_CT_L4[1][1][6][6] = dL3;
  Curvature_Higgs_CT_L4[1][1][7][7] = dL3;
  Curvature_Higgs_CT_L4[1][2][2][3] = dL7;
  Curvature_Higgs_CT_L4[1][2][4][7] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][5][6] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][3][3] = 3 * dL7;
  Curvature_Higgs_CT_L4[1][3][4][4] = dL6;
  Curvature_Higgs_CT_L4[1][3][4][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][5][5] = dL6;
  Curvature_Higgs_CT_L4[1][3][5][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][6][6] = dL7;
  Curvature_Higgs_CT_L4[1][3][7][7] = dL7;
  Curvature_Higgs_CT_L4[2][2][2][2] = 6 * dL2;
  Curvature_Higgs_CT_L4[2][2][3][3] = 2 * dL2;
  Curvature_Higgs_CT_L4[2][2][4][4] = dL3;
  Curvature_Higgs_CT_L4[2][2][4][6] = dL7;
  Curvature_Higgs_CT_L4[2][2][5][5] = dL3;
  Curvature_Higgs_CT_L4[2][2][5][7] = dL7;
  Curvature_Higgs_CT_L4[2][2][6][6] = 2 * dL2;
  Curvature_Higgs_CT_L4[2][2][7][7] = 2 * dL2;
  Curvature_Higgs_CT_L4[3][3][3][3] = 6 * dL2;
  Curvature_Higgs_CT_L4[3][3][4][4] = dL3;
  Curvature_Higgs_CT_L4[3][3][4][6] = dL7;
  Curvature_Higgs_CT_L4[3][3][5][5] = dL3;
  Curvature_Higgs_CT_L4[3][3][5][7] = dL7;
  Curvature_Higgs_CT_L4[3][3][6][6] = 2 * dL2;
  Curvature_Higgs_CT_L4[3][3][7][7] = 2 * dL2;
  Curvature_Higgs_CT_L4[4][4][4][4] = 6 * dL1;
  Curvature_Higgs_CT_L4[4][4][4][6] = 3 * dL6;
  Curvature_Higgs_CT_L4[4][4][5][5] = 2 * dL1;
  Curvature_Higgs_CT_L4[4][4][5][7] = dL6;
  Curvature_Higgs_CT_L4[4][4][6][6] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[4][4][7][7] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[4][5][5][6] = dL6;
  Curvature_Higgs_CT_L4[4][5][6][7] = dL5;
  Curvature_Higgs_CT_L4[4][6][6][6] = 3 * dL7;
  Curvature_Higgs_CT_L4[4][6][7][7] = dL7;
  Curvature_Higgs_CT_L4[5][5][5][5] = 6 * dL1;
  Curvature_Higgs_CT_L4[5][5][5][7] = 3 * dL6;
  Curvature_Higgs_CT_L4[5][5][6][6] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[5][5][7][7] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[5][6][6][7] = dL7;
  Curvature_Higgs_CT_L4[5][7][7][7] = 3 * dL7;
  Curvature_Higgs_CT_L4[6][6][6][6] = 6 * dL2;
  Curvature_Higgs_CT_L4[6][6][7][7] = 2 * dL2;
  Curvature_Higgs_CT_L4[7][7][7][7] = 6 * dL2;

  sym4Dim(Curvature_Higgs_CT_L4, NHiggs, NHiggs, NHiggs, NHiggs);

  Curvature_Higgs_CT_L6[0][0][0][0][0][0] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][0][0][0][2] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][0][0][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][0][0][1][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][0][0][2][2] =
      (double)(-6 * dOp6_111122 - 12 * dOp6_121211 - 6 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][0][0][3][3] =
      (double)(-6 * dOp6_122111 + 12 * dOp6_121211 - 6 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][0][0][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][0][0][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][0][0][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][0][0][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][0][0][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][0][0][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][0][1][1][2] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][0][1][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][0][0][2][2][2] =
      (double)(-9 * dOp6_112212 - 9 * dOp6_121212 - 9 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][0][2][3][3] =
      (double)(9 * dOp6_121212 - 3 * dOp6_121221 - 3 * dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][0][2][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][0][2][4][6] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][0][2][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][0][2][5][7] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][0][2][6][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][0][0][2][7][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][0][0][3][4][7] =
      (double)(-3 * dOp6_122111 + 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][0][3][5][6] =
      (double)(3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][1][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][1][1][1][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][1][1][2][2] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][1][1][3][3] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][1][1][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][1][1][4][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][1][1][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][1][1][5][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][1][1][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][1][1][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][1][2][2][3] =
      (double)(-9 * dOp6_121212 - dOp6_121221 - dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][1][2][4][7] =
      (double)(dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][2][5][6] =
      (double)(-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][3][3][3] =
      (double)(9 * dOp6_121212 - 3 * dOp6_121221 - 3 * dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][1][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][1][3][4][6] =
      (double)(-2 * dOp6_121211 - dOp6_122111) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][1][3][5][7] =
      (double)(-2 * dOp6_121211 - dOp6_122111) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][3][6][6] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][0][1][3][7][7] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][0][2][2][2][2] =
      (double)(-6 * dOp6_112222 - 12 * dOp6_121222 - 6 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][3][3] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][4][4] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_121211 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][4][6] =
      (double)(-dOp6_112212 - 3 * dOp6_121212 - 3 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][5][5] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_121211 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][5][7] =
      (double)(-dOp6_112212 - 3 * dOp6_121212 - 3 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][6][6] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][7][7] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][3][4][7] =
      (double)(3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][3][5][6] =
      (double)(-3 * dOp6_121212 + dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][3][3] =
      (double)(-6 * dOp6_112222 + 12 * dOp6_121222 - 6 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][4][4] =
      (double)(2 * dOp6_121211 - 2 * dOp6_111122 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][4][6] =
      (double)(-dOp6_112212 - dOp6_121221 + 3 * dOp6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][5][5] =
      (double)(2 * dOp6_121211 - 2 * dOp6_111122 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][5][7] =
      (double)(-dOp6_112212 - dOp6_121221 + 3 * dOp6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][6][6] =
      (double)(-2 * dOp6_112222 + 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][7][7] =
      (double)(-2 * dOp6_112222 + 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][4][4][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][4][4][5][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][4][4][6][6] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_121211 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][4][4][7][7] =
      (double)(2 * dOp6_121211 - 2 * dOp6_111122 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][4][5][5][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][0][4][6][6][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][0][4][6][7][7] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][0][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][5][5][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][0][5][5][6][6] =
      (double)(2 * dOp6_121211 - 2 * dOp6_111122 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][5][5][7][7] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_121211 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][5][6][6][7] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][0][5][7][7][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][0][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[0][0][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[0][0][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[0][1][1][1][1][2] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][1][1][1][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][1][1][2][2][2] =
      (double)(9 * dOp6_121212 - 3 * dOp6_121221 - 3 * dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][1][1][2][3][3] =
      (double)(-9 * dOp6_121212 - dOp6_121221 - dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][1][1][2][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][1][1][2][4][6] =
      (double)(-2 * dOp6_121211 - dOp6_122111) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][1][1][2][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][1][1][2][5][7] =
      (double)(-2 * dOp6_121211 - dOp6_122111) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][1][1][2][6][6] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][1][1][2][7][7] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][1][1][3][4][7] =
      (double)(-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][1][1][3][5][6] =
      (double)(dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][1][2][2][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[0][1][2][2][4][7] =
      (double)(-3 * dOp6_121212 + dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][1][2][2][5][6] =
      (double)(3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][1][2][3][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[0][1][2][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][1][2][3][4][6] =
      (double)(-3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][1][2][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][1][2][3][5][7] =
      (double)(-3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][1][2][3][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[0][1][2][3][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[0][1][3][3][4][7] =
      (double)(3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][1][3][3][5][6] =
      (double)(-3 * dOp6_121212 + dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][2][2][2][2][2] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][2][2][2][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][2][2][2][4][4] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][2][2][2][4][6] =
      (double)(-6 * dOp6_121222 - 3 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][2][2][5][5] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][2][2][2][5][7] =
      (double)(-6 * dOp6_121222 - 3 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][2][2][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][2][2][2][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][2][2][3][4][7] =
      (double)(-dOp6_122122 + 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][2][3][5][6] =
      (double)(dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][3][3][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][2][3][3][4][4] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][2][3][3][4][6] =
      (double)(-2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][3][3][5][5] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[0][2][3][3][5][7] =
      (double)(-2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][3][3][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][2][3][3][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][2][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][2][4][4][4][6] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][2][4][4][5][7] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][4][6][6] =
      (double)(-dOp6_112212 - 3 * dOp6_121212 - 3 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][2][4][4][7][7] =
      (double)(-dOp6_112212 - dOp6_121221 + 3 * dOp6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][2][4][5][5][6] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][5][6][7] =
      (double)(-3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][2][4][6][6][6] =
      (double)(-6 * dOp6_121222 - 3 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][6][7][7] =
      (double)(-2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[0][2][5][5][5][7] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][5][5][6][6] =
      (double)(-dOp6_112212 - dOp6_121221 + 3 * dOp6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][2][5][5][7][7] =
      (double)(-dOp6_112212 - 3 * dOp6_121212 - 3 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][2][5][6][6][7] =
      (double)(-2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][5][7][7][7] =
      (double)(-6 * dOp6_121222 - 3 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][2][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][2][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[0][3][3][3][4][7] =
      (double)(-3 * dOp6_122122 + 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][3][3][5][6] =
      (double)(3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][4][4][7] =
      (double)(-3 * dOp6_122111 + 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][4][5][6] =
      (double)(dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][4][6][7] =
      (double)(3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][3][4][5][5][7] =
      (double)(-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][5][6][6] =
      (double)(-3 * dOp6_121212 + dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][3][4][5][7][7] =
      (double)(3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][3][4][6][6][7] =
      (double)(-dOp6_122122 + 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][7][7][7] =
      (double)(-3 * dOp6_122122 + 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][5][5][5][6] =
      (double)(3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][5][5][6][7] =
      (double)(-3 * dOp6_121212 + dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][3][5][6][6][6] =
      (double)(3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][5][6][7][7] =
      (double)(dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][1][1][1] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][1][1][1][3] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][1][1][1][2][2] =
      (double)(12 * dOp6_121211 - 6 * dOp6_111122 - 6 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][1][1][3][3] =
      (double)(-6 * dOp6_111122 - 12 * dOp6_121211 - 6 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][1][1][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][1][1][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][1][1][1][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][1][1][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][1][1][1][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[1][1][1][1][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[1][1][1][2][2][3] =
      (double)(-3 * dOp6_121221 + 9 * dOp6_121212 - 3 * dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][1][2][4][7] =
      (double)(3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][2][5][6] =
      (double)(-3 * dOp6_122111 + 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][3][3][3] =
      (double)(-9 * dOp6_112212 - 9 * dOp6_121212 - 9 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][1][3][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][1][1][3][4][6] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][3][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][1][1][3][5][7] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][3][6][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][1][1][3][7][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][1][2][2][2][2] =
      (double)(-6 * dOp6_112222 + 12 * dOp6_121222 - 6 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][3][3] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][4][4] =
      (double)(2 * dOp6_121211 - 2 * dOp6_111122 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][4][6] =
      (double)(-dOp6_112212 - dOp6_121221 + 3 * dOp6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][5][5] =
      (double)(2 * dOp6_121211 - 2 * dOp6_111122 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][5][7] =
      (double)(-dOp6_112212 - dOp6_121221 + 3 * dOp6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][6][6] =
      (double)(-2 * dOp6_112222 + 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][7][7] =
      (double)(-2 * dOp6_112222 + 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][3][4][7] =
      (double)(-3 * dOp6_121212 + dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][3][5][6] =
      (double)(3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][3][3] =
      (double)(-6 * dOp6_112222 - 12 * dOp6_121222 - 6 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][4][4] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_121211 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][4][6] =
      (double)(-dOp6_112212 - 3 * dOp6_121212 - 3 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][5][5] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_121211 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][5][7] =
      (double)(-dOp6_112212 - 3 * dOp6_121212 - 3 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][6][6] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][7][7] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][4][4][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][1][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][4][4][5][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][1][4][4][6][6] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_121211 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][4][4][7][7] =
      (double)(2 * dOp6_121211 - 2 * dOp6_111122 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][4][5][5][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][1][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[1][1][4][6][6][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][1][4][6][7][7] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][1][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][5][5][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][1][5][5][6][6] =
      (double)(2 * dOp6_121211 - 2 * dOp6_111122 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][5][5][7][7] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_121211 - dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][5][6][6][7] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][1][5][7][7][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][1][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[1][1][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[1][1][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[1][2][2][2][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[1][2][2][2][4][7] =
      (double)(3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][2][5][6] =
      (double)(-3 * dOp6_122122 + 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][3][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[1][2][2][3][4][4] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][2][2][3][4][6] =
      (double)(-2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][3][5][5] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][2][2][3][5][7] =
      (double)(-2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][3][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[1][2][2][3][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[1][2][3][3][4][7] =
      (double)(dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][3][3][5][6] =
      (double)(-dOp6_122122 + 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][4][4][7] =
      (double)(3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][4][5][6] =
      (double)(-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][4][6][7] =
      (double)(-3 * dOp6_121212 + dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][2][4][5][5][7] =
      (double)(dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][5][6][6] =
      (double)(3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][2][4][5][7][7] =
      (double)(-3 * dOp6_121212 + dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][2][4][6][6][7] =
      (double)(dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][7][7][7] =
      (double)(3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][5][5][5][6] =
      (double)(-3 * dOp6_122111 + 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][5][5][6][7] =
      (double)(3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][2][5][6][6][6] =
      (double)(-3 * dOp6_122122 + 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][5][6][7][7] =
      (double)(-dOp6_122122 + 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][3][3][3][3] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[1][3][3][3][4][4] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][3][3][3][4][6] =
      (double)(-6 * dOp6_121222 - 3 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][3][3][5][5] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[1][3][3][3][5][7] =
      (double)(-6 * dOp6_121222 - 3 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][3][3][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[1][3][3][3][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[1][3][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][3][4][4][4][6] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][3][4][4][5][7] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][4][6][6] =
      (double)(-dOp6_112212 - 3 * dOp6_121212 - 3 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][3][4][4][7][7] =
      (double)(-dOp6_112212 - dOp6_121221 + 3 * dOp6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][3][4][5][5][6] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][5][6][7] =
      (double)(-3 * dOp6_121212 - dOp6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][3][4][6][6][6] =
      (double)(-6 * dOp6_121222 - 3 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][6][7][7] =
      (double)(-2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[1][3][5][5][5][7] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][5][5][6][6] =
      (double)(-dOp6_112212 - dOp6_121221 + 3 * dOp6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][3][5][5][7][7] =
      (double)(-dOp6_112212 - 3 * dOp6_121212 - 3 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][3][5][6][6][7] =
      (double)(-2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][5][7][7][7] =
      (double)(-6 * dOp6_121222 - 3 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[1][3][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[1][3][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[2][2][2][2][2][2] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][2][2][3][3] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][2][2][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[2][2][2][2][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[2][2][2][2][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[2][2][2][2][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[2][2][2][2][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][2][2][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][3][3][3][3] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][3][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[2][2][3][3][4][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[2][2][3][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[2][2][3][3][5][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[2][2][3][3][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][3][3][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[2][2][4][4][4][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[2][2][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[2][2][4][4][5][7] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[2][2][4][4][6][6] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][4][4][7][7] =
      (double)(2 * dOp6_121222 - 2 * dOp6_112222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][4][5][5][6] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[2][2][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[2][2][4][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[2][2][4][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[2][2][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[2][2][5][5][5][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[2][2][5][5][6][6] =
      (double)(2 * dOp6_121222 - 2 * dOp6_112222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][5][5][7][7] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][5][6][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[2][2][5][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[2][2][6][6][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][6][6][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][7][7][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][3][3][3][3] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][3][3][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[3][3][3][3][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[3][3][3][3][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[3][3][3][3][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[3][3][3][3][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][3][3][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[3][3][4][4][4][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[3][3][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[3][3][4][4][5][7] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[3][3][4][4][6][6] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][4][4][7][7] =
      (double)(2 * dOp6_121222 - 2 * dOp6_112222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][4][5][5][6] =
      -std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[3][3][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[3][3][4][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[3][3][4][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[3][3][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[3][3][5][5][5][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112212;
  Curvature_Higgs_CT_L6[3][3][5][5][6][6] =
      (double)(2 * dOp6_121222 - 2 * dOp6_112222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][5][5][7][7] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_121222 - dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][5][6][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[3][3][5][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[3][3][6][6][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][6][6][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][7][7][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[4][4][4][4][4][4] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[4][4][4][4][4][6] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[4][4][4][4][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[4][4][4][4][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[4][4][4][4][6][6] =
      (double)(-6 * dOp6_111122 - 12 * dOp6_121211 - 6 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][4][4][7][7] =
      (double)(12 * dOp6_121211 - 6 * dOp6_122111 - 6 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][4][5][5][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[4][4][4][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[4][4][4][6][6][6] =
      (double)(-9 * dOp6_112212 - 9 * dOp6_121212 - 9 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][4][6][7][7] =
      (double)(-3 * dOp6_121221 + 9 * dOp6_121212 - 3 * dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[4][4][5][5][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[4][4][5][5][6][6] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][5][5][7][7] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][5][6][6][7] =
      (double)(-9 * dOp6_121212 - dOp6_121221 - dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][5][7][7][7] =
      (double)(-3 * dOp6_121221 + 9 * dOp6_121212 - 3 * dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][6][6][6][6] =
      (double)(-6 * dOp6_112222 - 12 * dOp6_121222 - 6 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][6][6][7][7] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][7][7][7][7] =
      (double)(12 * dOp6_121222 - 6 * dOp6_112222 - 6 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][5][5][5][5][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[4][5][5][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[4][5][5][6][6][6] =
      (double)(-3 * dOp6_121221 + 9 * dOp6_121212 - 3 * dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][5][5][6][7][7] =
      (double)(-9 * dOp6_121212 - dOp6_121221 - dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][5][6][6][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[4][5][6][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[4][6][6][6][6][6] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[4][6][6][6][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[4][6][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[5][5][5][5][5][5] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[5][5][5][5][5][7] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111112;
  Curvature_Higgs_CT_L6[5][5][5][5][6][6] =
      (double)(12 * dOp6_121211 - 6 * dOp6_122111 - 6 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][5][5][7][7] =
      (double)(-6 * dOp6_111122 - 12 * dOp6_121211 - 6 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][5][6][6][7] =
      (double)(-3 * dOp6_121221 + 9 * dOp6_121212 - 3 * dOp6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][5][7][7][7] =
      (double)(-9 * dOp6_112212 - 9 * dOp6_121212 - 9 * dOp6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][6][6][6][6] =
      (double)(12 * dOp6_121222 - 6 * dOp6_112222 - 6 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][6][6][7][7] =
      (double)(-2 * dOp6_112222 - 2 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][7][7][7][7] =
      (double)(-6 * dOp6_112222 - 12 * dOp6_121222 - 6 * dOp6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][6][6][6][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[5][6][6][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[5][7][7][7][7][7] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_122222;
  Curvature_Higgs_CT_L6[6][6][6][6][6][6] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[6][6][6][6][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[6][6][7][7][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[7][7][7][7][7][7] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;

  sym6Dim(
      Curvature_Higgs_CT_L6, NHiggs, NHiggs, NHiggs, NHiggs, NHiggs, NHiggs);

  return;
}

/**
 * Console-Output of all Parameters
 */
void Class_Potential_R2HDMEFTPHI6_PSI2PHI3::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "scale = " << scale << "\n";

  ss << "The parameters are : \n";
  ss << "Model = " << Model << "\n";
  ss << "v1 = " << C_vev0 * C_CosBeta << "\n";
  ss << "v2 = " << C_vev0 * C_SinBeta << "\n";
  ss << "Type = " << Type << "\n";

  ss << "beta = " << beta << "\n";
  ss << "tan(beta) = " << TanBeta << "\n";
  ss << "L1 = " << L1 << "\n";
  ss << "L2 = " << L2 << "\n";
  ss << "L3 = " << L3 << "\n";
  ss << "L4 = " << L4 << "\n";
  ss << "L5 = " << L5 << "\n";
  ss << "L6 = " << L6 << "\n";
  ss << "L7 = " << L7 << "\n";
  ss << "m12Sq = " << m12Sq << "\n";
  ss << "m11Sq = " << m11Sq << "\n";
  ss << "m22Sq = " << m22Sq << "\n";

  ss << "The counterterms are :\n";

  ss << "dL1 := " << dL1 << ";\n";
  ss << "dL2 := " << dL2 << ";\n";
  ss << "dL3 := " << dL3 << ";\n";
  ss << "dL4 := " << dL4 << ";\n";
  ss << "dL5 := " << dL5 << ";\n";
  ss << "dL6 := " << dL6 << ";\n";
  ss << "dL7 := " << dL7 << ";\n";
  ss << "dm11Sq := " << dm11Sq << ";\n";
  ss << "dm22Sq := " << dm22Sq << ";\n";
  ss << "dm12Sq := " << dm12Sq << ";\n";

  ss << "dT1 := " << dT1 << ";\n";
  ss << "dT2 := " << dT2 << ";\n";
  ss << "dT3 := " << dT3 << ";\n";
  ss << "dT4 := " << dT4 << ";\n";
  ss << "dT5 := " << dT5 << ";\n";
  ss << "dT6 := " << dT6 << ";\n";
  ss << "dT7 := " << dT7 << ";\n";
  ss << "dT8 := " << dT8 << ";\n";

  ss << "dOp6_111111 := " << dOp6_111111 << ";\n";
  ss << "dOp6_111122 := " << dOp6_111122 << ";\n";
  ss << "dOp6_122111 := " << dOp6_122111 << ";\n";
  ss << "dOp6_121211 := " << dOp6_121211 << ";\n";
  ss << "dOp6_111112 := " << dOp6_111112 << ";\n";
  ss << "dOp6_121221 := " << dOp6_121221 << ";\n";
  ss << "dOp6_112212 := " << dOp6_112212 << ";\n";
  ss << "dOp6_222222 := " << dOp6_222222 << ";\n";
  ss << "dOp6_112222 := " << dOp6_112222 << ";\n";
  ss << "dOp6_122122 := " << dOp6_122122 << ";\n";
  ss << "dOp6_121222 := " << dOp6_121222 << ";\n";
  ss << "dOp6_122222 := " << dOp6_122222 << ";\n";
  ss << "dOp6_121212 := " << dOp6_121212 << ";\n";

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  int posMHCS1 = 0;
  int posN[2];
  int countposN = 0;
  int posA      = 0;
  int posG1 = 0, posG0 = 0;
  double testsum             = 0;
  const double ZeroThreshold = 1e-5;
  for (std::size_t i = 0; i < 3; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
    if (testsum > ZeroThreshold) posG1 = i;
    //		testsum = std::abs(HiggsRot(i,1)) + std::abs(HiggsRot(i,3));
    //		if(testsum > ZeroThreshold) posG2 = i;
    testsum = std::abs(HiggsRot(i, 5)) + std::abs(HiggsRot(i, 7));
    if (testsum > ZeroThreshold) posG0 = i;
  }
  for (std::size_t i = 3; i < NHiggs; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
    if (testsum > ZeroThreshold) posMHCS1 = i;
    //		testsum = std::abs(HiggsRot(i,1)) + std::abs(HiggsRot(i,3));
    //		if(testsum > ZeroThreshold) posMHCS2 = i;
    testsum = std::abs(HiggsRot(i, 4)) + std::abs(HiggsRot(i, 6));
    if (testsum > ZeroThreshold)
    {
      posN[countposN] = i;
      countposN++;
    }
    testsum = std::abs(HiggsRot(i, 5)) + std::abs(HiggsRot(i, 7));
    if (testsum > ZeroThreshold) posA = i;
  }

  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  ss << "The mass spectrum is given by :\n";
  ss << "m_{G^+}^2 = " << HiggsMasses[posG1] << " GeV^2 \n";
  ss << "m_{G^0}^2 = " << HiggsMasses[posG0] << " GeV^2 \n";
  ss << "m_{H^+} = " << std::sqrt(HiggsMasses[posMHCS1]) << " GeV \n"
     << "m_h = " << std::sqrt(HiggsMasses[posN[0]]) << " GeV \n"
     << "m_H = " << std::sqrt(HiggsMasses[posN[1]]) << " GeV \n"
     << "m_A = " << std::sqrt(HiggsMasses[posA]) << " GeV \n";

  ss << "The neutral mixing Matrix is given by :\n";
  ss << "h = " << HiggsRot(posN[0], 4) << " zeta_1 ";
  bool IsNegative = HiggsRot(posN[0], 6) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(HiggsRot(posN[0], 6)) << " zeta_2\n"
     << "H = " << HiggsRot(posN[1], 4) << " zeta_1 ";
  IsNegative = HiggsRot(posN[1], 6) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(HiggsRot(posN[1], 6)) << " zeta_2\n";

  ss << "The dim-6 operator are set to:"
     << "\n"
     << "purely-scalar:\n";
  ss << "Op6_111111 = " << Op6_111111 << "\n";
  ss << "Op6_111122 = " << Op6_111122 << "\n";
  ss << "Op6_122111 = " << Op6_122111 << "\n";
  ss << "Op6_121211 = " << Op6_121211 << "\n";
  ss << "Op6_111112 = " << Op6_111112 << "\n";
  ss << "Op6_121221 = " << Op6_121221 << "\n";
  ss << "Op6_112212 = " << Op6_112212 << "\n";
  ss << "Op6_222222 = " << Op6_222222 << "\n";
  ss << "Op6_112222 = " << Op6_112222 << "\n";
  ss << "Op6_122122 = " << Op6_122122 << "\n";
  ss << "Op6_121222 = " << Op6_121222 << "\n";
  ss << "Op6_122222 = " << Op6_122222 << "\n";
  ss << "Op6_121212 = " << Op6_121212 << "\n";
  ss << "up-type quarks:\n";
  ss << "OQu_1b11b = " << OQu_1b11b << "\n";
  ss << "OQu_1b12b = " << OQu_1b12b << "\n";
  ss << "OQu_1b21b = " << OQu_1b21b << "\n";
  ss << "OQu_1b22b = " << OQu_1b22b << "\n";
  ss << "OQu_2b11b = " << OQu_2b11b << "\n";
  ss << "OQu_2b12b = " << OQu_2b12b << "\n";
  ss << "OQu_2b21b = " << OQu_2b21b << "\n";
  ss << "OQu_2b22b = " << OQu_2b22b << "\n";
  ss << "down-type quarks:\n";
  ss << "OQd_1b11b = " << OQd_1b11b << "\n";
  ss << "OQd_1b12b = " << OQd_1b12b << "\n";
  ss << "OQd_1b21b = " << OQd_1b21b << "\n";
  ss << "OQd_1b22b = " << OQd_1b22b << "\n";
  ss << "OQd_2b11b = " << OQd_2b11b << "\n";
  ss << "OQd_2b12b = " << OQd_2b12b << "\n";
  ss << "OQd_2b21b = " << OQd_2b21b << "\n";
  ss << "OQd_2b21b = " << OQd_2b21b << "\n";
  ss << "leptons:\n";
  ss << "OL_1b11b = " << OL_1b11b << "\n";
  ss << "OL_1b12b = " << OL_1b12b << "\n";
  ss << "OL_1b21b = " << OL_1b21b << "\n";
  ss << "OL_1b22b = " << OL_1b22b << "\n";
  ss << "OL_2b11b = " << OL_2b11b << "\n";
  ss << "OL_2b12b = " << OL_2b12b << "\n";
  ss << "OL_2b21b = " << OL_2b21b << "\n";
  ss << "OL_2b22b = " << OL_2b22b << "\n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms for the R2HDM-EFT
 */
std::vector<double> Class_Potential_R2HDMEFTPHI6_PSI2PHI3::calc_CT() const
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

  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;

  VectorXd NablaWeinberg(8);
  MatrixXd HesseWeinberg(8, 8), HiggsRot(8, 8);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg(i) = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HesseWeinberg(i, j) = WeinbergHesse.at((j)*NHiggs + i);
      if (std::abs(HesseWeinberg(i, j)) <= 1e-3) HesseWeinberg(i, j) = 0;
    }
  }

  // dm11Sq
  parCT.push_back(
      (-4 * v2 * v2 * HesseWeinberg(3, 3) + HesseWeinberg(4, 4) * v1 * v1 -
       3 * HesseWeinberg(5, 5) * v1 * v1 + 4 * HesseWeinberg(7, 7) * v2 * v2 +
       HesseWeinberg(4, 6) * v2 * v1 - HesseWeinberg(5, 7) * v2 * v1) *
      std::pow(v1, -2) / 2);
  // dm22Sq
  parCT.push_back((-4 * HesseWeinberg(3, 3) * v2 + HesseWeinberg(4, 6) * v1 -
                   HesseWeinberg(5, 7) * v1 +
                   v2 * (HesseWeinberg(6, 6) + HesseWeinberg(7, 7))) /
                  v2 / 2);
  // dm12Sq
  parCT.push_back((-2 * HesseWeinberg(3, 3) * v2 + HesseWeinberg(5, 7) * v1 +
                   2 * HesseWeinberg(7, 7) * v2) /
                  v1);
  // dL1
  parCT.push_back(
      (2 * v2 * v2 * HesseWeinberg(3, 3) - HesseWeinberg(4, 4) * v1 * v1 +
       HesseWeinberg(5, 5) * v1 * v1 - 2 * HesseWeinberg(7, 7) * v2 * v2) *
      std::pow(v1, -4) / 2);
  // dL2
  parCT.push_back(
      (2 * HesseWeinberg(3, 3) - HesseWeinberg(6, 6) - HesseWeinberg(7, 7)) *
      std::pow(v2, -2) / 2);
  // dL3
  parCT.push_back((-HesseWeinberg(4, 6) + HesseWeinberg(5, 7)) / v2 / v1);
  // dL4
  parCT.push_back(0);
  // dL5
  parCT.push_back((-2 * HesseWeinberg(3, 3) + 2 * HesseWeinberg(7, 7)) *
                  std::pow(v1, -2));
  // dL6
  parCT.push_back(0);
  // dL7
  parCT.push_back(0);
  // dT1
  parCT.push_back(-NablaWeinberg(0));
  // dT2
  parCT.push_back(-NablaWeinberg(1));
  // dT3
  parCT.push_back(-NablaWeinberg(2));
  // dT4
  parCT.push_back(-NablaWeinberg(3));
  // dT5
  parCT.push_back(HesseWeinberg(5, 5) * v1 + HesseWeinberg(5, 7) * v2 -
                  NablaWeinberg(4));
  // dT6
  parCT.push_back(-NablaWeinberg(5));
  // dT7
  parCT.push_back(HesseWeinberg(5, 7) * v1 + HesseWeinberg(7, 7) * v2 -
                  NablaWeinberg(6));
  // dT8
  parCT.push_back(-NablaWeinberg(7));
  // dOp6_111111
  parCT.push_back(0);
  // dOp6_111122
  parCT.push_back(0);
  // dOp6_122111
  parCT.push_back(0);
  // dOp6_121211
  parCT.push_back(0);
  // dOp6_111112
  parCT.push_back(0);
  // dOp6_121221
  parCT.push_back(0);
  // dOp6_112212
  parCT.push_back(0);
  // dOp6_222222
  parCT.push_back(0);
  // dOp6_112222
  parCT.push_back(0);
  // dOp6_122122
  parCT.push_back(0);
  // dOp6_121222
  parCT.push_back(0);
  // dOp6_122222
  parCT.push_back(0);
  // dOp6_121212
  parCT.push_back(0);

  return parCT;
}

/**
 * Calculates the corrections to the Triple higgs couplings in the mass basis.
 */
void Class_Potential_R2HDMEFTPHI6_PSI2PHI3::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

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
  int posMHCS1 = 0, posMHCS2 = 0;
  int posN[2]   = {-1, -1};
  int countposN = 0;
  int posG1 = 0, posG2 = 0, posG0 = 0;
  int posA = 0, posh = 0, posH = 0;
  double testsum             = 0;
  const double ZeroThreshold = 1e-5;
  for (std::size_t i = 0; i < 3; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
    if (testsum > ZeroThreshold) posG1 = i;
    testsum = std::abs(HiggsRot(i, 1)) + std::abs(HiggsRot(i, 3));
    if (testsum > ZeroThreshold) posG2 = i;
    testsum = std::abs(HiggsRot(i, 5)) + std::abs(HiggsRot(i, 7));
    if (testsum > ZeroThreshold) posG0 = i;
  }
  for (std::size_t i = 3; i < NHiggs; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
    if (testsum > ZeroThreshold) posMHCS1 = i;
    testsum = std::abs(HiggsRot(i, 1)) + std::abs(HiggsRot(i, 3));
    if (testsum > ZeroThreshold) posMHCS2 = i;
    testsum = std::abs(HiggsRot(i, 5)) + std::abs(HiggsRot(i, 7));
    if (testsum > ZeroThreshold) posA = i;
    testsum = 0;
    testsum = std::abs(HiggsRot(i, 4)) + std::abs(HiggsRot(i, 6));
    if (testsum > ZeroThreshold)
    {
      posN[countposN] = i;
      countposN++;
    }
  }

  posh = posN[0];
  posH = posN[1];

  std::vector<double> HiggsOrder(NHiggs);
  HiggsOrder[0] = posG1;
  HiggsOrder[1] = posG2;
  HiggsOrder[2] = posMHCS1;
  HiggsOrder[3] = posMHCS2;
  HiggsOrder[4] = posG0;
  HiggsOrder[5] = posA;
  HiggsOrder[6] = posh;
  HiggsOrder[7] = posH;

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
              //  			  double RotFac =
              //  (HiggsRot(i,l)*HiggsRot(j,m)*HiggsRot(k,n)).real();
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

void Class_Potential_R2HDMEFTPHI6_PSI2PHI3::SetCurvatureArrays()
{
  initVectors();

  HiggsVev[4] = C_vev0 * C_CosBeta;
  HiggsVev[6] = C_vev0 * C_SinBeta;

  Curvature_Higgs_L2[0][0] = m11Sq;
  Curvature_Higgs_L2[0][2] = -m12Sq;
  Curvature_Higgs_L2[1][1] = m11Sq;
  Curvature_Higgs_L2[1][3] = -m12Sq;
  Curvature_Higgs_L2[2][2] = m22Sq;
  Curvature_Higgs_L2[3][3] = m22Sq;
  Curvature_Higgs_L2[4][4] = m11Sq;
  Curvature_Higgs_L2[4][6] = -m12Sq;
  Curvature_Higgs_L2[5][5] = m11Sq;
  Curvature_Higgs_L2[5][7] = -m12Sq;
  Curvature_Higgs_L2[6][6] = m22Sq;
  Curvature_Higgs_L2[7][7] = m22Sq;

  sym2Dim(Curvature_Higgs_L2, NHiggs, NHiggs);

  Curvature_Higgs_L4[0][0][0][0] = 6 * L1;
  Curvature_Higgs_L4[0][0][0][2] = 3 * L6;
  Curvature_Higgs_L4[0][0][1][1] = 2 * L1;
  Curvature_Higgs_L4[0][0][1][3] = L6;
  Curvature_Higgs_L4[0][0][2][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][0][4][4] = 2 * L1;
  Curvature_Higgs_L4[0][0][4][6] = L6;
  Curvature_Higgs_L4[0][0][5][5] = 2 * L1;
  Curvature_Higgs_L4[0][0][5][7] = L6;
  Curvature_Higgs_L4[0][0][6][6] = L3;
  Curvature_Higgs_L4[0][0][7][7] = L3;
  Curvature_Higgs_L4[0][1][1][2] = L6;
  Curvature_Higgs_L4[0][1][2][3] = L5;
  Curvature_Higgs_L4[0][2][2][2] = 3 * L7;
  Curvature_Higgs_L4[0][2][3][3] = L7;
  Curvature_Higgs_L4[0][2][4][4] = L6;
  Curvature_Higgs_L4[0][2][4][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][2][5][5] = L6;
  Curvature_Higgs_L4[0][2][5][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][2][6][6] = L7;
  Curvature_Higgs_L4[0][2][7][7] = L7;
  Curvature_Higgs_L4[0][3][4][7] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][6] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][1][1][1] = 6 * L1;
  Curvature_Higgs_L4[1][1][1][3] = 3 * L6;
  Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][1][3][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][1][4][4] = 2 * L1;
  Curvature_Higgs_L4[1][1][4][6] = L6;
  Curvature_Higgs_L4[1][1][5][5] = 2 * L1;
  Curvature_Higgs_L4[1][1][5][7] = L6;
  Curvature_Higgs_L4[1][1][6][6] = L3;
  Curvature_Higgs_L4[1][1][7][7] = L3;
  Curvature_Higgs_L4[1][2][2][3] = L7;
  Curvature_Higgs_L4[1][2][4][7] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][6] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[1][3][3][3] = 3 * L7;
  Curvature_Higgs_L4[1][3][4][4] = L6;
  Curvature_Higgs_L4[1][3][4][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][3][5][5] = L6;
  Curvature_Higgs_L4[1][3][5][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][3][6][6] = L7;
  Curvature_Higgs_L4[1][3][7][7] = L7;
  Curvature_Higgs_L4[2][2][2][2] = 6 * L2;
  Curvature_Higgs_L4[2][2][3][3] = 2 * L2;
  Curvature_Higgs_L4[2][2][4][4] = L3;
  Curvature_Higgs_L4[2][2][4][6] = L7;
  Curvature_Higgs_L4[2][2][5][5] = L3;
  Curvature_Higgs_L4[2][2][5][7] = L7;
  Curvature_Higgs_L4[2][2][6][6] = 2 * L2;
  Curvature_Higgs_L4[2][2][7][7] = 2 * L2;
  Curvature_Higgs_L4[3][3][3][3] = 6 * L2;
  Curvature_Higgs_L4[3][3][4][4] = L3;
  Curvature_Higgs_L4[3][3][4][6] = L7;
  Curvature_Higgs_L4[3][3][5][5] = L3;
  Curvature_Higgs_L4[3][3][5][7] = L7;
  Curvature_Higgs_L4[3][3][6][6] = 2 * L2;
  Curvature_Higgs_L4[3][3][7][7] = 2 * L2;
  Curvature_Higgs_L4[4][4][4][4] = 6 * L1;
  Curvature_Higgs_L4[4][4][4][6] = 3 * L6;
  Curvature_Higgs_L4[4][4][5][5] = 2 * L1;
  Curvature_Higgs_L4[4][4][5][7] = L6;
  Curvature_Higgs_L4[4][4][6][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][5][5][6] = L6;
  Curvature_Higgs_L4[4][5][6][7] = L5;
  Curvature_Higgs_L4[4][6][6][6] = 3 * L7;
  Curvature_Higgs_L4[4][6][7][7] = L7;
  Curvature_Higgs_L4[5][5][5][5] = 6 * L1;
  Curvature_Higgs_L4[5][5][5][7] = 3 * L6;
  Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][5][7][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[5][6][6][7] = L7;
  Curvature_Higgs_L4[5][7][7][7] = 3 * L7;
  Curvature_Higgs_L4[6][6][6][6] = 6 * L2;
  Curvature_Higgs_L4[6][6][7][7] = 2 * L2;
  Curvature_Higgs_L4[7][7][7][7] = 6 * L2;

  sym4Dim(Curvature_Higgs_L4, NHiggs, NHiggs, NHiggs, NHiggs);

  Curvature_Higgs_L6[0][0][0][0][0][0] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][0][0][0][2] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][0][0][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][0][0][1][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][0][0][2][2] =
      (double)(-6 * Op6_122111 - 12 * Op6_121211 - 6 * Op6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][0][0][3][3] =
      (double)(-6 * Op6_122111 + 12 * Op6_121211 - 6 * Op6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][0][0][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][0][0][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][0][0][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][0][0][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][0][0][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][0][0][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][0][1][1][2] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][0][1][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][0][0][2][2][2] =
      (double)(-9 * Op6_112212 - 9 * Op6_121212 - 9 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][0][2][3][3] =
      (double)(-3 * Op6_112212 + 9 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][0][2][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][0][2][4][6] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][0][2][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][0][2][5][7] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][0][2][6][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][0][0][2][7][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][0][0][3][4][7] =
      (double)(-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][0][3][5][6] =
      (double)(3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][1][1][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][1][1][1][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][1][1][2][2] =
      (double)(-2 * Op6_122111 - 2 * Op6_111122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][1][1][3][3] =
      (double)(-2 * Op6_122111 - 2 * Op6_111122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][1][1][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][1][1][4][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][1][1][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][1][1][5][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][1][1][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][1][1][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][1][2][2][3] =
      (double)(-Op6_112212 - Op6_121221 - 9 * Op6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][1][2][4][7] = (double)(Op6_122111 - 2 * Op6_121211) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][1][2][5][6] =
      (double)(-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][1][3][3][3] =
      (double)(-3 * Op6_112212 + 9 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][1][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][1][3][4][6] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][1][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][1][3][5][7] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][1][3][6][6] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][0][1][3][7][7] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][0][2][2][2][2] =
      (double)(-6 * Op6_122122 - 12 * Op6_121222 - 6 * Op6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][3][3] =
      (double)(-2 * Op6_122122 - 2 * Op6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][4][4] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][4][6] =
      (double)(-Op6_112212 - 3 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][5][5] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][5][7] =
      (double)(-Op6_112212 - 3 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][6][6] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][3][4][7] =
      (double)(-Op6_121221 + 3 * Op6_121212) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][3][5][6] =
      std::pow(LambdaEFT, -0.2e1) * (double)(Op6_121221 - 3 * Op6_121212);
  Curvature_Higgs_L6[0][0][3][3][3][3] =
      (double)(12 * Op6_121222 - 6 * Op6_122122 - 6 * Op6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][4][4] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][4][6] =
      (double)(-Op6_112212 + 3 * Op6_121212 - Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][5][5] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][5][7] =
      (double)(-Op6_112212 + 3 * Op6_121212 - Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][6][6] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][7][7] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][4][4][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][4][4][5][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][4][4][6][6] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][4][4][7][7] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][4][5][5][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][0][4][6][6][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][0][4][6][7][7] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][0][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][5][5][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][0][5][5][6][6] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][5][5][7][7] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][5][6][6][7] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][0][5][7][7][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][0][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[0][0][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[0][0][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[0][1][1][1][1][2] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][1][1][1][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][1][1][2][2][2] =
      (double)(-3 * Op6_112212 + 9 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][1][1][2][3][3] =
      (double)(-Op6_112212 - Op6_121221 - 9 * Op6_121212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][1][1][2][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][1][1][2][4][6] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][1][1][2][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][1][1][2][5][7] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][1][1][2][6][6] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][1][1][2][7][7] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][1][1][3][4][7] =
      (double)(-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][1][1][3][5][6] = (double)(Op6_122111 - 2 * Op6_121211) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][1][2][2][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[0][1][2][2][4][7] =
      std::pow(LambdaEFT, -0.2e1) * (double)(Op6_121221 - 3 * Op6_121212);
  Curvature_Higgs_L6[0][1][2][2][5][6] =
      (double)(-Op6_121221 + 3 * Op6_121212) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][1][2][3][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[0][1][2][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][1][2][3][4][6] =
      (double)(-3 * Op6_121212 - Op6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][1][2][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][1][2][3][5][7] =
      (double)(-3 * Op6_121212 - Op6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][1][2][3][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[0][1][2][3][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[0][1][3][3][4][7] =
      (double)(-Op6_121221 + 3 * Op6_121212) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][1][3][3][5][6] =
      std::pow(LambdaEFT, -0.2e1) * (double)(Op6_121221 - 3 * Op6_121212);
  Curvature_Higgs_L6[0][2][2][2][2][2] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][2][2][2][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][2][2][2][4][4] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][2][2][2][4][6] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][2][2][5][5] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][2][2][2][5][7] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][2][2][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][2][2][2][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][2][2][3][4][7] =
      (double)(-Op6_122122 + 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][2][3][5][6] = (double)(Op6_122122 - 2 * Op6_121222) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][3][3][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][2][3][3][4][4] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][2][3][3][4][6] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][3][3][5][5] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[0][2][3][3][5][7] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][3][3][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][2][3][3][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][2][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][2][4][4][4][6] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][2][4][4][5][7] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][4][4][6][6] =
      (double)(-Op6_112212 - 3 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][2][4][4][7][7] =
      (double)(-Op6_112212 + 3 * Op6_121212 - Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][2][4][5][5][6] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][4][5][6][7] =
      (double)(-3 * Op6_121212 - Op6_121221) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][2][4][6][6][6] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][4][6][7][7] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[0][2][5][5][5][7] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][5][5][6][6] =
      (double)(-Op6_112212 + 3 * Op6_121212 - Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][2][5][5][7][7] =
      (double)(-Op6_112212 - 3 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][2][5][6][6][7] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][5][7][7][7] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][2][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][2][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[0][3][3][3][4][7] =
      (double)(-3 * Op6_122122 + 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][3][3][5][6] =
      (double)(3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][4][4][4][7] =
      (double)(-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][4][4][5][6] = (double)(Op6_122111 - 2 * Op6_121211) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][4][4][6][7] =
      (double)(-Op6_121221 + 3 * Op6_121212) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][3][4][5][5][7] =
      (double)(-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][4][5][6][6] =
      std::pow(LambdaEFT, -0.2e1) * (double)(Op6_121221 - 3 * Op6_121212);
  Curvature_Higgs_L6[0][3][4][5][7][7] =
      (double)(-Op6_121221 + 3 * Op6_121212) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][3][4][6][6][7] =
      (double)(-Op6_122122 + 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][4][7][7][7] =
      (double)(-3 * Op6_122122 + 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][5][5][5][6] =
      (double)(3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][5][5][6][7] =
      std::pow(LambdaEFT, -0.2e1) * (double)(Op6_121221 - 3 * Op6_121212);
  Curvature_Higgs_L6[0][3][5][6][6][6] =
      (double)(3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][5][6][7][7] = (double)(Op6_122122 - 2 * Op6_121222) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][1][1][1][1][1] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][1][1][1][3] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][1][1][1][2][2] =
      (double)(-6 * Op6_122111 + 12 * Op6_121211 - 6 * Op6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][1][1][3][3] =
      (double)(-6 * Op6_122111 - 12 * Op6_121211 - 6 * Op6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][1][1][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][1][1][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][1][1][1][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][1][1][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][1][1][1][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[1][1][1][1][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[1][1][1][2][2][3] =
      (double)(-3 * Op6_112212 + 9 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][1][2][4][7] =
      (double)(3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][1][1][2][5][6] =
      (double)(-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][1][1][3][3][3] =
      (double)(-9 * Op6_121221 - 9 * Op6_121212 - 9 * Op6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][1][3][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][1][1][3][4][6] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][1][1][3][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][1][1][3][5][7] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][1][1][3][6][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][1][1][3][7][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][1][2][2][2][2] =
      (double)(12 * Op6_121222 - 6 * Op6_122122 - 6 * Op6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][3][3] =
      (double)(-2 * Op6_122122 - 2 * Op6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][4][4] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][4][6] =
      (double)(-Op6_112212 + 3 * Op6_121212 - Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][5][5] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][5][7] =
      (double)(-Op6_112212 + 3 * Op6_121212 - Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][6][6] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][7][7] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][3][4][7] =
      std::pow(LambdaEFT, -0.2e1) * (double)(Op6_121221 - 3 * Op6_121212);
  Curvature_Higgs_L6[1][1][2][3][5][6] =
      (double)(-Op6_121221 + 3 * Op6_121212) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][3][3] =
      (double)(-6 * Op6_122122 - 12 * Op6_121222 - 6 * Op6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][4][4] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][4][6] =
      (double)(-Op6_112212 - 3 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][5][5] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][5][7] =
      (double)(-Op6_112212 - 3 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][6][6] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][4][4][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][1][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][4][4][5][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][1][4][4][6][6] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][4][4][7][7] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][4][5][5][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][1][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[1][1][4][6][6][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][1][4][6][7][7] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][1][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][5][5][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][1][5][5][6][6] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][5][5][7][7] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][5][6][6][7] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][1][5][7][7][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][1][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[1][1][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[1][1][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[1][2][2][2][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[1][2][2][2][4][7] =
      (double)(3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][2][2][5][6] =
      (double)(-3 * Op6_122122 + 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][2][3][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[1][2][2][3][4][4] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][2][2][3][4][6] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][2][3][5][5] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][2][2][3][5][7] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][2][3][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[1][2][2][3][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[1][2][3][3][4][7] = (double)(Op6_122122 - 2 * Op6_121222) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][3][3][5][6] =
      (double)(-Op6_122122 + 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][4][4][4][7] =
      (double)(-6 * Op6_121211 + 3 * Op6_122111) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][4][4][5][6] = (double)(2 * Op6_121211 - Op6_122111) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][4][4][6][7] =
      std::pow(LambdaEFT, -0.2e1) * (double)(Op6_121221 - 3 * Op6_121212);
  Curvature_Higgs_L6[1][2][4][5][5][7] =
      (double)(-2 * Op6_121211 + Op6_122111) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][4][5][6][6] =
      (double)(-Op6_121221 + 3 * Op6_121212) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][2][4][5][7][7] =
      std::pow(LambdaEFT, -0.2e1) * (double)(Op6_121221 - 3 * Op6_121212);
  Curvature_Higgs_L6[1][2][4][6][6][7] = (double)(Op6_122122 - 2 * Op6_121222) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][4][7][7][7] =
      (double)(3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][5][5][5][6] =
      (double)(6 * Op6_121211 - 3 * Op6_122111) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][5][5][6][7] =
      (double)(-Op6_121221 + 3 * Op6_121212) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][2][5][6][6][6] =
      (double)(-3 * Op6_122122 + 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][5][6][7][7] =
      (double)(-Op6_122122 + 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][3][3][3][3] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[1][3][3][3][4][4] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][3][3][3][4][6] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][3][3][5][5] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[1][3][3][3][5][7] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][3][3][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[1][3][3][3][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[1][3][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][3][4][4][4][6] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][3][4][4][5][7] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][4][4][6][6] =
      (double)(-Op6_112212 - 3 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][3][4][4][7][7] =
      (double)(-Op6_112212 + 3 * Op6_121212 - Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][3][4][5][5][6] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][4][5][6][7] =
      (double)(-Op6_121221 - 3 * Op6_121212) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][3][4][6][6][6] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][4][6][7][7] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[1][3][5][5][5][7] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][5][5][6][6] =
      (double)(-Op6_112212 + 3 * Op6_121212 - Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][3][5][5][7][7] =
      (double)(-Op6_112212 - 3 * Op6_121212 - 3 * Op6_121221) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][3][5][6][6][7] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][5][7][7][7] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[1][3][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[1][3][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[2][2][2][2][2][2] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][2][2][3][3] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][2][2][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[2][2][2][2][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[2][2][2][2][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[2][2][2][2][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[2][2][2][2][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][2][2][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][3][3][3][3] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][3][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[2][2][3][3][4][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[2][2][3][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[2][2][3][3][5][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[2][2][3][3][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][3][3][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[2][2][4][4][4][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[2][2][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[2][2][4][4][5][7] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[2][2][4][4][6][6] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][4][4][7][7] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][4][5][5][6] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[2][2][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[2][2][4][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[2][2][4][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[2][2][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[2][2][5][5][5][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[2][2][5][5][6][6] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][5][5][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][5][6][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[2][2][5][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[2][2][6][6][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][6][6][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][7][7][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][3][3][3][3] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][3][3][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[3][3][3][3][4][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[3][3][3][3][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[3][3][3][3][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[3][3][3][3][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][3][3][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[3][3][4][4][4][6] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[3][3][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[3][3][4][4][5][7] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[3][3][4][4][6][6] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][4][4][7][7] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][4][5][5][6] =
      -std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[3][3][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[3][3][4][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[3][3][4][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[3][3][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[3][3][5][5][5][7] =
      -0.3e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112212;
  Curvature_Higgs_L6[3][3][5][5][6][6] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][5][5][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][5][6][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[3][3][5][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[3][3][6][6][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][6][6][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][7][7][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[4][4][4][4][4][4] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[4][4][4][4][4][6] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[4][4][4][4][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[4][4][4][4][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[4][4][4][4][6][6] =
      (double)(-12 * Op6_121211 - 6 * Op6_122111 - 6 * Op6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][4][4][7][7] =
      (double)(12 * Op6_121211 - 6 * Op6_122111 - 6 * Op6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][4][5][5][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[4][4][4][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[4][4][4][6][6][6] =
      (double)(-9 * Op6_121221 - 9 * Op6_121212 - 9 * Op6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][4][6][7][7] =
      (double)(-3 * Op6_121221 + 9 * Op6_121212 - 3 * Op6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[4][4][5][5][5][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[4][4][5][5][6][6] =
      (double)(-2 * Op6_122111 - 2 * Op6_111122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][5][5][7][7] =
      (double)(-2 * Op6_122111 - 2 * Op6_111122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][5][6][6][7] =
      (double)(-Op6_121221 - 9 * Op6_121212 - Op6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][5][7][7][7] =
      (double)(-3 * Op6_121221 + 9 * Op6_121212 - 3 * Op6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][6][6][6][6] =
      (double)(-6 * Op6_122122 - 6 * Op6_112222 - 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][6][6][7][7] =
      (double)(-2 * Op6_122122 - 2 * Op6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][7][7][7][7] =
      (double)(12 * Op6_121222 - 6 * Op6_122122 - 6 * Op6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][5][5][5][5][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[4][5][5][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[4][5][5][6][6][6] =
      (double)(-3 * Op6_121221 + 9 * Op6_121212 - 3 * Op6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][5][5][6][7][7] =
      (double)(-Op6_121221 - 9 * Op6_121212 - Op6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][5][6][6][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[4][5][6][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[4][6][6][6][6][6] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[4][6][6][6][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[4][6][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[5][5][5][5][5][5] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[5][5][5][5][5][7] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111112;
  Curvature_Higgs_L6[5][5][5][5][6][6] =
      (double)(12 * Op6_121211 - 6 * Op6_122111 - 6 * Op6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][5][5][7][7] =
      (double)(-12 * Op6_121211 - 6 * Op6_122111 - 6 * Op6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][5][6][6][7] =
      (double)(-3 * Op6_121221 + 9 * Op6_121212 - 3 * Op6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][5][7][7][7] =
      (double)(-9 * Op6_121221 - 9 * Op6_121212 - 9 * Op6_112212) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][6][6][6][6] =
      (double)(12 * Op6_121222 - 6 * Op6_122122 - 6 * Op6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][6][6][7][7] =
      (double)(-2 * Op6_122122 - 2 * Op6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][7][7][7][7] =
      (double)(-6 * Op6_122122 - 6 * Op6_112222 - 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][6][6][6][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[5][6][6][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[5][7][7][7][7][7] =
      -0.30e2 * std::pow(LambdaEFT, -0.2e1) * Op6_122222;
  Curvature_Higgs_L6[6][6][6][6][6][6] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[6][6][6][6][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[6][6][7][7][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[7][7][7][7][7][7] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;

  sym6Dim(Curvature_Higgs_L6, NHiggs, NHiggs, NHiggs, NHiggs, NHiggs, NHiggs);

  Curvature_Gauge_G2H2[0][0][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][4][4] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][5][5] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][6][6] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][7][7] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][0][4] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][5] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][6] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][7] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][4][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][5][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][6][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][7][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][4][4] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][5][5] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][6][6] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][7][7] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][0][5] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][4] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][7] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][6] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][4][1] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][5][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][6][3] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][7][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][4][4] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][5][5] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][6][6] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][7][7] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][0][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][4][4] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][5][5] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][6][6] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][7][7] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][4] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][5] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][6] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][7] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][4][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][5][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][6][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][7][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][5] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][4] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][7] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][6] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][4][1] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][5][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][6][3] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][7][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][4][4] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][5][5] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][6][6] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][7][7] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][3][0][0] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][1][1] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][2][2] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][3][3] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][4][4] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][5][5] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][6][6] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][7][7] = C_gs * C_gs / 0.2e1;

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

  MatrixXcd YIJR2(NQuarks, NQuarks), YIJE2(NQuarks, NQuarks),
      YIJS2(NQuarks, NQuarks), YIJP2(NQuarks, NQuarks), YIJRD(NQuarks, NQuarks),
      YIJED(NQuarks, NQuarks), YIJSD(NQuarks, NQuarks), YIJPD(NQuarks, NQuarks);
  MatrixXcd YIJRL(NLepton, NLepton), YIJEL(NLepton, NLepton),
      YIJSL(NLepton, NLepton), YIJPL(NLepton, NLepton);
  YIJR2 = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJE2 = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJS2 = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJP2 = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJRD = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJED = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJSD = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJPD = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJRL = MatrixXcd::Zero(NLepton, NLepton);
  YIJEL = MatrixXcd::Zero(NLepton, NLepton);
  YIJSL = MatrixXcd::Zero(NLepton, NLepton);
  YIJPL = MatrixXcd::Zero(NLepton, NLepton);

  std::complex<double> II(0, 1);

  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;
  double vL = v2;
  double vD = v2;
  if (Type == 2)
  {
    vL = v1;
    vD = v1;
  }
  else if (Type == 3)
    vL = v1;
  else if (Type == 4)
    vD = v1;

  YIJR2(0, 9)  = -std::conj(V11) * C_MassUp / v2;
  YIJR2(0, 10) = -std::conj(V12) * C_MassUp / v2;
  YIJR2(0, 11) = -std::conj(V13) * C_MassUp / v2;

  YIJR2(1, 9)  = -std::conj(V21) * C_MassCharm / v2;
  YIJR2(1, 10) = -std::conj(V22) * C_MassCharm / v2;
  YIJR2(1, 11) = -std::conj(V23) * C_MassCharm / v2;

  YIJR2(2, 9)  = -std::conj(V31) * C_MassTop / v2;
  YIJR2(2, 10) = -std::conj(V32) * C_MassTop / v2;
  YIJR2(2, 11) = -std::conj(V33) * C_MassTop / v2;

  YIJS2(0, 6) = C_MassUp / v2;
  YIJS2(1, 7) = C_MassCharm / v2;
  YIJS2(2, 8) = C_MassTop / v2;

  YIJSD(3, 9)  = C_MassDown / vD;
  YIJSD(4, 10) = C_MassStrange / vD;
  YIJSD(5, 11) = C_MassBottom / vD;

  YIJRD(3, 6) = V11 * C_MassDown / vD;
  YIJRD(3, 7) = V21 * C_MassDown / vD;
  YIJRD(3, 8) = V31 * C_MassDown / vD;
  YIJRD(4, 6) = V12 * C_MassStrange / vD;
  YIJRD(4, 7) = V22 * C_MassStrange / vD;
  YIJRD(4, 8) = V32 * C_MassStrange / vD;
  YIJRD(5, 6) = V13 * C_MassBottom / vD;
  YIJRD(5, 7) = V23 * C_MassBottom / vD;
  YIJRD(5, 8) = V33 * C_MassBottom / vD;

  YIJRL(1, 6) = C_MassElectron / vL;
  YIJRL(3, 7) = C_MassMu / vL;
  YIJRL(5, 8) = C_MassTau / vL;

  YIJSL(0, 1) = C_MassElectron / vL;
  YIJSL(2, 3) = C_MassMu / vL;
  YIJSL(4, 5) = C_MassTau / vL;

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    for (std::size_t j = 0; j < i; j++)
    {
      YIJR2(i, j) = YIJR2(j, i);
      YIJS2(i, j) = YIJS2(j, i);
      YIJRD(i, j) = YIJRD(j, i);
      YIJSD(i, j) = YIJSD(j, i);
    }
  }
  for (std::size_t i = 0; i < NLepton; i++)
  {
    for (std::size_t j = 0; j < i; j++)
    {
      YIJRL(i, j) = YIJRL(j, i);
      YIJSL(i, j) = YIJSL(j, i);
    }
  }

  YIJP2 = std::complex<double>(-1, 0) * II * YIJS2;
  YIJE2 = std::complex<double>(-1, 0) * II * YIJR2;

  YIJPD = II * YIJSD;
  YIJED = II * YIJRD;

  YIJPL = II * YIJSL;
  YIJEL = II * YIJRL;

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    for (std::size_t j = 0; j < NQuarks; j++)
    {
      Curvature_Quark_F2H1[i][j][0] = 0;
      Curvature_Quark_F2H1[i][j][1] = 0;
      Curvature_Quark_F2H1[i][j][2] = YIJR2(i, j);
      Curvature_Quark_F2H1[i][j][3] = YIJE2(i, j);
      Curvature_Quark_F2H1[i][j][4] = 0;
      Curvature_Quark_F2H1[i][j][5] = 0;
      Curvature_Quark_F2H1[i][j][6] = YIJS2(i, j);
      Curvature_Quark_F2H1[i][j][7] = YIJP2(i, j);

      if (Type == 1 or Type == 3)
      {
        Curvature_Quark_F2H1[i][j][2] += YIJRD(i, j);
        Curvature_Quark_F2H1[i][j][3] += YIJED(i, j);
        Curvature_Quark_F2H1[i][j][6] += YIJSD(i, j);
        Curvature_Quark_F2H1[i][j][7] += YIJPD(i, j);
      }
      else
      {
        Curvature_Quark_F2H1[i][j][0] += YIJRD(i, j);
        Curvature_Quark_F2H1[i][j][1] += YIJED(i, j);
        Curvature_Quark_F2H1[i][j][4] += YIJSD(i, j);
        Curvature_Quark_F2H1[i][j][5] += YIJPD(i, j);
      }
    }
  }

  for (std::size_t i = 0; i < NLepton; i++)
  {
    for (std::size_t j = 0; j < NLepton; j++)
    {
      if (Type == 1 or Type == 4)
      {
        Curvature_Lepton_F2H1[i][j][0] = 0;
        Curvature_Lepton_F2H1[i][j][1] = 0;
        Curvature_Lepton_F2H1[i][j][2] = YIJRL(i, j);
        Curvature_Lepton_F2H1[i][j][3] = YIJEL(i, j);
        Curvature_Lepton_F2H1[i][j][4] = 0;
        Curvature_Lepton_F2H1[i][j][5] = 0;
        Curvature_Lepton_F2H1[i][j][6] = YIJSL(i, j);
        Curvature_Lepton_F2H1[i][j][7] = YIJPL(i, j);
      }
      else
      {
        Curvature_Lepton_F2H1[i][j][2] = 0;
        Curvature_Lepton_F2H1[i][j][3] = 0;
        Curvature_Lepton_F2H1[i][j][0] = YIJRL(i, j);
        Curvature_Lepton_F2H1[i][j][1] = YIJEL(i, j);
        Curvature_Lepton_F2H1[i][j][6] = 0;
        Curvature_Lepton_F2H1[i][j][7] = 0;
        Curvature_Lepton_F2H1[i][j][4] = YIJSL(i, j);
        Curvature_Lepton_F2H1[i][j][5] = YIJPL(i, j);
      }
    }
  }

  Curvature_Lepton_F2H3[4][5][0][0][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][0][0][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][0][0][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][0][0][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][0][2][4] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b + OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][2][5] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b + OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][2][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b + OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][2][7] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b + OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][3][4] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][3][5] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][3][6] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b - OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][3][7] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b - OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][1][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][1][1][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][1][1][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][1][1][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][1][2][4] = -II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][2][5] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][2][6] = -II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b - OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][2][7] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b - OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][3][4] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b + OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][3][5] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b + OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][3][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b + OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][3][7] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b + OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][2][2][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][2][2][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][2][2][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][2][2][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][3][3][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][3][3][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][3][3][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][3][3][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][4][4] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b;
  Curvature_Lepton_F2H3[4][5][4][4][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][4][6] = std::sqrt(0.2e1) *
                                         (OL_1b12b + OL_1b21b + OL_2b11b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][4][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][5][5] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][5][6] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][5][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][6][6] = std::sqrt(0.2e1) *
                                         (OL_1b22b + OL_2b12b + OL_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][6][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b12b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][7][7] = std::sqrt(0.2e1) *
                                         (OL_1b22b - OL_2b12b + OL_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][5][5][5] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) * OL_1b11b;
  Curvature_Lepton_F2H3[4][5][5][5][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][5][5][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OL_1b12b + OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][5][6][6] = std::sqrt(0.2e1) * II *
                                         (OL_1b22b - OL_2b12b + OL_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][5][6][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b12b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][5][7][7] = std::sqrt(0.2e1) * II *
                                         (OL_1b22b + OL_2b12b + OL_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[4][5][6][6][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][6][7][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][7][7][7] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[5][8][0][0][0] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b;
  Curvature_Lepton_F2H3[5][8][0][0][1] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][0][2] = std::sqrt(0.2e1) *
                                         (OL_1b12b + OL_1b21b + OL_2b11b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][0][3] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][1][1] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][1][2] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][1][3] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][2][2] = std::sqrt(0.2e1) *
                                         (OL_1b22b + OL_2b12b + OL_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][2][3] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b12b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][3][3] = std::sqrt(0.2e1) *
                                         (OL_1b22b - OL_2b12b + OL_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][4][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][4][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b + OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][0][4][7] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][0][5][5] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][5][6] = -II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][0][5][7] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b + OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][0][6][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][7][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][1][1] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) * OL_1b11b;
  Curvature_Lepton_F2H3[5][8][1][1][2] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][1][3] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OL_1b12b + OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][2][2] = std::sqrt(0.2e1) * II *
                                         (OL_1b22b - OL_2b12b + OL_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][2][3] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b12b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][3][3] = std::sqrt(0.2e1) * II *
                                         (OL_1b22b + OL_2b12b + OL_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][4][4] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][4][6] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b + OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][1][4][7] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][1][5][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][5][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][1][5][7] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_1b12b + OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][1][6][6] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][7][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_1b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][2][2] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[5][8][2][2][3] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][3][3] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][4][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][4][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b + OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][2][4][7] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b - OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][2][5][5] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][5][6] = -II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b - OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][2][5][7] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b + OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][2][6][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][7][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][3][3][3] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[5][8][3][4][4] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][3][4][6] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b + OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][3][4][7] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b - OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][3][5][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][3][5][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b - OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][3][5][7] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OL_2b12b + OL_2b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][3][6][6] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][3][7][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OL_2b22b / 0.2e1;

  sym5Dim(Curvature_Lepton_F2H3, NLepton, NLepton, NHiggs, NHiggs, NHiggs);

  CorrectLeptonTensorsDim6();

  Curvature_Quark_F2H3[2][8][0][0][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][0][0][5] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][0][0][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][0][0][7] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][0][2][4] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_1b12b + OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][2][5] = -II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_1b12b + OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][2][6] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_2b12b + OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][2][7] = -II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_2b12b + OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][3][4] = II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_1b12b - OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][3][5] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_1b12b - OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][3][6] = II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_2b12b - OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][3][7] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_2b12b - OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][1][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][1][1][5] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][1][1][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][1][1][7] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][1][2][4] = -II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_1b12b - OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][2][5] = -std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_1b12b - OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][2][6] = -II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_2b12b - OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][2][7] = -std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_2b12b - OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][3][4] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_1b12b + OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][3][5] = -II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_1b12b + OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][3][6] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_2b12b + OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][3][7] = -II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQu_2b12b + OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][8][2][2][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][2][2][5] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][2][2][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][2][2][7] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][3][3][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][3][3][5] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][3][3][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][3][3][7] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][4][4] = 0.3e1 / 0.2e1 * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) * OQu_1b11b;
  Curvature_Quark_F2H3[2][8][4][4][5] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][4][6] = std::sqrt(0.2e1) *
                                        (OQu_1b12b + OQu_1b21b + OQu_2b11b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][4][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OQu_1b12b - OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][5][5] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][5][6] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b12b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][5][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b12b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][6][6] = std::sqrt(0.2e1) *
                                        (OQu_1b22b + OQu_2b12b + OQu_2b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][6][7] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][7][7] = std::sqrt(0.2e1) *
                                        (OQu_1b22b + OQu_2b12b - OQu_2b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][5][5] = -0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) * OQu_1b11b;
  Curvature_Quark_F2H3[2][8][5][5][6] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OQu_1b12b - OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][5][7] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OQu_1b12b + OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][6][6] = -std::sqrt(0.2e1) * II *
                                        (OQu_1b22b + OQu_2b12b - OQu_2b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][6][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][7][7] = -std::sqrt(0.2e1) * II *
                                        (OQu_1b22b + OQu_2b12b + OQu_2b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][6][6] = 0.3e1 / 0.2e1 * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[2][8][6][6][7] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][7][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][7][7] = -0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[2][11][0][0][0] = -0.3e1 / 0.2e1 * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         OQu_1b11b;
  Curvature_Quark_F2H3[2][11][0][0][1] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][0][2] = -std::sqrt(0.2e1) *
                                         (OQu_1b12b + OQu_1b21b + OQu_2b11b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][0][3] =
      -II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OQu_1b12b - OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][1][1] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][1][2] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b12b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][1][3] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b12b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][2][2] = -std::sqrt(0.2e1) *
                                         (OQu_1b22b + OQu_2b12b + OQu_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][2][3] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][3][3] = -std::sqrt(0.2e1) *
                                         (OQu_1b22b + OQu_2b12b - OQu_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][4][4] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][4][6] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_1b12b + OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][4][7] = -II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_1b12b - OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][5][5] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][5][6] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_1b12b - OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][5][7] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_1b12b + OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][6][6] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][7][7] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][1][1] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         OQu_1b11b;
  Curvature_Quark_F2H3[2][11][1][1][2] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OQu_1b12b - OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][1][3] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OQu_1b12b + OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][2][2] = std::sqrt(0.2e1) * II *
                                         (OQu_1b22b + OQu_2b12b - OQu_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][2][3] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][3][3] = std::sqrt(0.2e1) * II *
                                         (OQu_1b22b + OQu_2b12b + OQu_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][4][4] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][4][6] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_1b12b + OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][4][7] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_1b12b - OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][5][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][5][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_1b12b - OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][5][7] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_1b12b + OQu_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][6][6] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][7][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_1b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][2][2] = -0.3e1 / 0.2e1 * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         OQu_2b22b;
  Curvature_Quark_F2H3[2][11][2][2][3] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][3][3] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][4][4] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][4][6] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_2b12b + OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][2][4][7] = -II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_2b12b - OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][2][5][5] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][5][6] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_2b12b - OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][2][5][7] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_2b12b + OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][2][6][6] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][7][7] =
      -std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][3][3] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         OQu_2b22b;
  Curvature_Quark_F2H3[2][11][3][4][4] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][4][6] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_2b12b + OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][3][4][7] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_2b12b - OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][3][5][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][5][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_2b12b - OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][3][5][7] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQu_2b12b + OQu_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[2][11][3][6][6] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][7][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][0][0] = 0.3e1 / 0.2e1 * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) * OQd_1b11b;
  Curvature_Quark_F2H3[5][8][0][0][1] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][0][2] = std::sqrt(0.2e1) *
                                        (OQd_1b12b + OQd_1b21b + OQd_2b11b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][0][3] = std::sqrt(0.2e1) * II *
                                        (OQd_1b12b - OQd_1b21b + OQd_2b11b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][1][1] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][1][2] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][1][3] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][2][2] = std::sqrt(0.2e1) *
                                        (OQd_1b22b + OQd_2b12b + OQd_2b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][2][3] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b12b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][3][3] = std::sqrt(0.2e1) *
                                        (OQd_1b22b - OQd_2b12b + OQd_2b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][4][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][4][6] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQd_1b12b + OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][4][7] = std::sqrt(0.2e1) * II *
                                        (OQd_1b12b - OQd_1b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][5][5] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][5][6] = -std::sqrt(0.2e1) * II *
                                        (OQd_1b12b - OQd_1b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][5][7] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQd_1b12b + OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][6][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][7][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][1][1] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) * OQd_1b11b;
  Curvature_Quark_F2H3[5][8][1][1][2] = std::sqrt(0.2e1) *
                                        (OQd_1b12b - OQd_1b21b + OQd_2b11b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][1][3] = std::sqrt(0.2e1) * II *
                                        (OQd_1b12b + OQd_1b21b + OQd_2b11b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][2][2] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OQd_1b22b - OQd_2b12b + OQd_2b21b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][2][3] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b12b / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][3][3] = std::sqrt(0.2e1) * II *
                                        (OQd_1b22b + OQd_2b12b + OQd_2b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][4][4] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][4][6] = std::sqrt(0.2e1) * II *
                                        (OQd_1b12b + OQd_1b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][4][7] = -std::sqrt(0.2e1) *
                                        (OQd_1b12b - OQd_1b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][5][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][5][6] = std::sqrt(0.2e1) *
                                        (OQd_1b12b - OQd_1b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][5][7] = std::sqrt(0.2e1) * II *
                                        (OQd_1b12b + OQd_1b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][6][6] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][7][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][2][2] = 0.3e1 / 0.2e1 * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[5][8][2][2][3] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][3][3] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][4][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][4][6] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQd_2b12b + OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][2][4][7] = II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQd_2b12b - OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][2][5][5] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][5][6] = -II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQd_2b12b - OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][2][5][7] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQd_2b12b + OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][2][6][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][7][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][3][3] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[5][8][3][4][4] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][4][6] = std::sqrt(0.2e1) * II *
                                        (OQd_2b12b + OQd_2b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][8][3][4][7] = -std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQd_2b12b - OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][3][5][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][5][6] = std::sqrt(0.2e1) *
                                        std::pow(LambdaEFT, -0.2e1) *
                                        (OQd_2b12b - OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][3][5][7] = std::sqrt(0.2e1) * II *
                                        (OQd_2b12b + OQd_2b21b) *
                                        std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][8][3][6][6] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][7][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][0][0][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][0][0][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][0][0][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][0][0][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][0][2][4] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQd_1b12b + OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][2][5] = std::sqrt(0.2e1) * II *
                                         (OQd_1b12b + OQd_1b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][2][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQd_2b12b + OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][2][7] = std::sqrt(0.2e1) * II *
                                         (OQd_2b12b + OQd_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][3][4] = std::sqrt(0.2e1) * II *
                                         (OQd_1b12b - OQd_1b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][3][5] = -std::sqrt(0.2e1) *
                                         (OQd_1b12b - OQd_1b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][3][6] = II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQd_2b12b - OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][3][7] = -std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQd_2b12b - OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][1][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][1][1][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][1][1][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][1][1][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][1][2][4] = -std::sqrt(0.2e1) * II *
                                         (OQd_1b12b - OQd_1b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][2][5] = std::sqrt(0.2e1) *
                                         (OQd_1b12b - OQd_1b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][2][6] = -II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQd_2b12b - OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][2][7] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQd_2b12b - OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][3][4] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQd_1b12b + OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][3][5] = std::sqrt(0.2e1) * II *
                                         (OQd_1b12b + OQd_1b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][3][6] = std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         (OQd_2b12b + OQd_2b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][3][7] = std::sqrt(0.2e1) * II *
                                         (OQd_2b12b + OQd_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.4e1;
  Curvature_Quark_F2H3[5][11][2][2][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][2][2][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][2][2][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][2][2][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][3][3][4] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][3][3][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][3][3][6] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][3][3][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][4][4] = 0.3e1 / 0.2e1 * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         OQd_1b11b;
  Curvature_Quark_F2H3[5][11][4][4][5] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][4][6] = std::sqrt(0.2e1) *
                                         (OQd_1b12b + OQd_1b21b + OQd_2b11b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][4][7] = std::sqrt(0.2e1) * II *
                                         (OQd_1b12b - OQd_1b21b + OQd_2b11b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][5][5] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][5][6] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][5][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][6][6] = std::sqrt(0.2e1) *
                                         (OQd_1b22b + OQd_2b12b + OQd_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][6][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b12b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][7][7] = std::sqrt(0.2e1) *
                                         (OQd_1b22b - OQd_2b12b + OQd_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][5][5] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         OQd_1b11b;
  Curvature_Quark_F2H3[5][11][5][5][6] = std::sqrt(0.2e1) *
                                         (OQd_1b12b - OQd_1b21b + OQd_2b11b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][5][7] = std::sqrt(0.2e1) * II *
                                         (OQd_1b12b + OQd_1b21b + OQd_2b11b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][6][6] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) *
      (OQd_1b22b - OQd_2b12b + OQd_2b21b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][6][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b12b / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][7][7] = std::sqrt(0.2e1) * II *
                                         (OQd_1b22b + OQd_2b12b + OQd_2b21b) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][6][6] = 0.3e1 / 0.2e1 * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         OQd_2b22b;
  Curvature_Quark_F2H3[5][11][6][6][7] =
      II * std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][7][7] =
      std::sqrt(0.2e1) * std::pow(LambdaEFT, -0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][7][7] = 0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) *
                                         std::pow(LambdaEFT, -0.2e1) *
                                         OQd_2b22b;

  sym5Dim(Curvature_Quark_F2H3, NQuarks, NQuarks, NHiggs, NHiggs, NHiggs);

  CorrectQuarkTensorsDim6();

  SetCurvatureDone = true;

  GetYukawaModifier(alpha);
}

double Class_Potential_R2HDMEFTPHI6_PSI2PHI3::SymFac_Higgs_OneLoop(
    const int &i,
    const int &j,
    const std::vector<double> &point) const
{
  double vev1 = point[4];
  double vev2 = point[6];

  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;

  double vU = v2;
  double vD = v2;
  if (Type == 2)
  {
    vD = v1;
  }
  else if (Type == 4)
    vD = v1;

  double yt = std::sqrt(2) * C_MassTop / vU;
  double yb = std::sqrt(2) * C_MassBottom / vD;

  // rho1, rho2
  if (i == 0 and j == 0) // rho1rho1
  {
    return (std::pow(vev1, 2) * yb * OQd_1b11b) /
               (2. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 3) * yb * OQd_2b22b) /
               (4. * std::pow(LambdaEFT, 2) * vev1) +
           (vev1 * vev2 * (3 * yb * OQd_2b11b + 2 * yt * OQu_1b11b)) /
               (4. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) *
            (yb * (OQd_2b12b + OQd_2b21b) + 2 * yt * OQu_2b11b)) /
               (4. * std::pow(LambdaEFT, 2));
  }
  else if ((i == 0 and j == 2) or (i == 2 and j == 0)) // rho1rho2
  {
    return (std::pow(vev1, 2) *
            (yb * (2 * OQd_1b12b + 2 * OQd_1b21b - OQd_2b11b) -
             yt * OQu_1b11b)) /
               (8. * std::pow(LambdaEFT, 2)) +
           (vev1 * vev2 *
            (yb * (OQd_2b12b + OQd_2b21b) + yt * (OQu_1b12b + OQu_1b21b))) /
               (8. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) *
            (-(yb * OQd_2b22b) -
             yt * (OQu_1b22b - 2 * (OQu_2b12b + OQu_2b21b)))) /
               (8. * std::pow(LambdaEFT, 2));
  }
  else if (i == 2 and j == 2) // rho2rho2
  {
    return (std::pow(vev1, 3) * yt * OQu_1b11b) /
               (4. * std::pow(LambdaEFT, 2) * vev2) +
           (std::pow(vev1, 2) *
            (2 * yb * OQd_1b22b + yt * (OQu_1b12b + OQu_1b21b))) /
               (4. * std::pow(LambdaEFT, 2)) +
           (vev1 * vev2 * (2 * yb * OQd_2b22b + 3 * yt * OQu_1b22b)) /
               (4. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) * yt * OQu_2b22b) / (2. * std::pow(LambdaEFT, 2));
  }
  // eta1, eta2
  else if (i == 1 and j == 1) // eta1eta1
  {
    return (std::pow(vev1, 2) * yb * OQd_1b11b) /
               (2. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 3) * yb * OQd_2b22b) /
               (4. * std::pow(LambdaEFT, 2) * vev1) +
           (vev1 * vev2 * (3 * yb * OQd_2b11b + 2 * yt * OQu_1b11b)) /
               (4. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) *
            (yb * (OQd_2b12b + OQd_2b21b) + 2 * yt * OQu_2b11b)) /
               (4. * std::pow(LambdaEFT, 2));
  }
  else if ((i == 1 and j == 3) or (j == 3 and i == 1)) // eta1eta2
  {
    return (std::pow(vev1, 2) *
            (yb * (2 * OQd_1b12b + 2 * OQd_1b21b - OQd_2b11b) -
             yt * OQu_1b11b)) /
               (8. * std::pow(LambdaEFT, 2)) +
           (vev1 * vev2 *
            (yb * (OQd_2b12b + OQd_2b21b) + yt * (OQu_1b12b + OQu_1b21b))) /
               (8. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) *
            (-(yb * OQd_2b22b) -
             yt * (OQu_1b22b - 2 * (OQu_2b12b + OQu_2b21b)))) /
               (8. * std::pow(LambdaEFT, 2));
  }
  else if (i == 3 and j == 3) // eta2eta2
  {
    return (std::pow(vev1, 3) * yt * OQu_1b11b) /
               (4. * std::pow(LambdaEFT, 2) * vev2) +
           (std::pow(vev1, 2) *
            (2 * yb * OQd_1b22b + yt * (OQu_1b12b + OQu_1b21b))) /
               (4. * std::pow(LambdaEFT, 2)) +
           (vev1 * vev2 * (2 * yb * OQd_2b22b + 3 * yt * OQu_1b22b)) /
               (4. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) * yt * OQu_2b22b) / (2. * std::pow(LambdaEFT, 2));
  }
  // zeta1, zeta2
  else if (i == 4 and j == 4) // zeta1zeta1
  {
    return (std::pow(vev1, 2) * yb * OQd_1b11b) / std::pow(LambdaEFT, 2) +
           (std::pow(vev2, 3) * yb * OQd_2b22b) /
               (4. * std::pow(LambdaEFT, 2) * vev1) +
           (vev1 * vev2 *
            (yb * (OQd_1b12b + OQd_1b21b + OQd_2b11b) + 6 * yt * OQu_1b11b)) /
               (4. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) * yt * (OQu_1b12b + OQu_1b21b + OQu_2b11b)) /
               (2. * std::pow(LambdaEFT, 2));
  }
  else if ((i == 4 and j == 6) or (i == 6 and j == 6)) // zeta1zeta2
  {
    return (3 * std::pow(vev1, 2) *
            (yb * (OQd_1b12b + OQd_1b21b + OQd_2b11b) - yt * OQu_1b11b)) /
               (8. * std::pow(LambdaEFT, 2)) +
           (vev1 * vev2 *
            (yb * (OQd_1b22b + OQd_2b12b + OQd_2b21b) +
             yt * (OQu_1b12b + OQu_1b21b + OQu_2b11b))) /
               (4. * std::pow(LambdaEFT, 2)) +
           (3 * std::pow(vev2, 2) *
            (-(yb * OQd_2b22b) + yt * (OQu_1b22b + OQu_2b12b + OQu_2b21b))) /
               (8. * std::pow(LambdaEFT, 2));
  }
  else if (i == 6 and j == 6) // zeta2zeta2
  {
    return (std::pow(vev1, 2) * yb * (OQd_1b22b + OQd_2b12b + OQd_2b21b)) /
               (2. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev1, 3) * yt * OQu_1b11b) /
               (4. * std::pow(LambdaEFT, 2) * vev2) +
           (vev1 * vev2 *
            (6 * yb * OQd_2b22b + yt * (OQu_1b22b + OQu_2b12b + OQu_2b21b))) /
               (4. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) * yt * OQu_2b22b) / std::pow(LambdaEFT, 2);
  }
  // psi1, psi2
  else if (i == 5 and j == 5) // psi1psi1
  {
    return (std::pow(vev1, 2) * yb * OQd_1b11b) /
               (2. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 3) * yb * OQd_2b22b) /
               (4. * std::pow(LambdaEFT, 2) * vev1) +
           (vev1 * vev2 *
            (3 * yb * (OQd_1b12b - OQd_1b21b + OQd_2b11b) +
             2 * yt * OQu_1b11b)) /
               (4. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) *
            (yb * OQd_2b12b + yt * (-OQu_1b12b + OQu_1b21b + OQu_2b11b))) /
               (2. * std::pow(LambdaEFT, 2));
  }
  else if ((i == 5 and j == 7) or (i == 7 and i == 5)) // psi1psi2
  {
    return -0.125 *
               (std::pow(vev1, 2) *
                (yb * (OQd_1b12b - 5 * OQd_1b21b + OQd_2b11b) +
                 yt * OQu_1b11b)) /
               std::pow(LambdaEFT, 2) +
           (vev1 * vev2 * (yb * OQd_2b12b + yt * OQu_1b12b)) /
               (4. * std::pow(LambdaEFT, 2)) -
           (std::pow(vev2, 2) *
            (yb * OQd_2b22b + yt * (OQu_1b22b + OQu_2b12b - 5 * OQu_2b21b))) /
               (8. * std::pow(LambdaEFT, 2));
  }
  else if (i == 7 and j == 7) // psi2psi2
  {
    return (std::pow(vev1, 3) * yt * OQu_1b11b) /
               (4. * std::pow(LambdaEFT, 2) * vev2) +
           (std::pow(vev1, 2) *
            (yb * (OQd_1b22b - OQd_2b12b + OQd_2b21b) + yt * OQu_1b12b)) /
               (2. * std::pow(LambdaEFT, 2)) +
           (vev1 * vev2 *
            (2 * yb * OQd_2b22b +
             3 * yt * (OQu_1b22b + OQu_2b12b - OQu_2b21b))) /
               (4. * std::pow(LambdaEFT, 2)) +
           (std::pow(vev2, 2) * yt * OQu_2b22b) / (2. * std::pow(LambdaEFT, 2));
  }
  else
  {
    return 0;
  }
}

double
Class_Potential_R2HDMEFTPHI6_PSI2PHI3::SymFac_Higgs_TwoLoop(const int &i,
                                                            const int &j) const
{
  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;

  double vU = v2;
  double vD = v2;
  if (Type == 2)
  {
    vD = v1;
  }
  else if (Type == 4)
    vD = v1;

  double yt = std::sqrt(2) * C_MassTop / vU;
  double yb = std::sqrt(2) * C_MassBottom / vD;

  // rho1, rho2
  if (i == 0 and j == 0) // rho1rho1
  {
    return (yb * OQd_1b11b) / (8. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_1b22b) / (12. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b21b) / (24. * std::pow(LambdaEFT, 2));
  }
  else if ((i == 0 and j == 2) or (i == 2 and j == 0)) // rho1rho2
  {
    return (yb * OQd_1b12b) / (48. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b11b) / (24. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b22b) / (16. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_1b11b) / (16. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_1b22b) / (24. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b12b) / (48. * std::pow(LambdaEFT, 2));
  }
  else if (i == 2 and j == 2) // rho2rho2
  {
    return (yt * OQu_1b21b) / (24. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b11b) / (12. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b22b) / (8. * std::pow(LambdaEFT, 2));
  }
  // eta1, eta2
  else if (i == 1 and j == 1) // eta1eta1
  {
    return (yb * OQd_1b11b) / (8. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_1b22b) / (12. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b21b) / (24. * std::pow(LambdaEFT, 2));
  }
  else if ((i == 1 and j == 3) or (j == 3 and i == 1)) // eta1eta2
  {
    return (yb * OQd_1b12b) / (48. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b11b) / (24. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b22b) / (16. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_1b11b) / (16. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_1b22b) / (24. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b12b) / (48. * std::pow(LambdaEFT, 2));
  }
  else if (i == 3 and j == 3) // eta2eta2
  {
    return (yt * OQu_1b21b) / (24. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b11b) / (12. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b22b) / (8. * std::pow(LambdaEFT, 2));
  }
  // zeta1, zeta2
  else if (i == 4 and j == 4) // zeta1zeta1
  {
    return (yb * OQd_1b11b) / (8. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_1b22b) / (12. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b21b) / (24. * std::pow(LambdaEFT, 2));
  }
  else if ((i == 4 and j == 6) or (i == 6 and j == 6)) // zeta1zeta2
  {
    return (yb * OQd_1b12b) / (48. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b11b) / (24. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b22b) / (16. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_1b11b) / (16. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_1b22b) / (24. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b12b) / (48. * std::pow(LambdaEFT, 2));
  }
  else if (i == 6 and j == 6) // zeta2zeta2
  {
    return (yt * OQu_1b21b) / (24. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b11b) / (12. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b22b) / (8. * std::pow(LambdaEFT, 2));
  }
  // psi1, psi2
  else if (i == 5 and j == 5) // psi1psi1
  {
    return (yb * OQd_1b11b) / (8. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_1b22b) / (12. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b21b) / (24. * std::pow(LambdaEFT, 2));
  }
  else if ((i == 5 and j == 7) or (i == 7 and i == 5)) // psi1psi2
  {
    return (yb * OQd_1b12b) / (48. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b11b) / (24. * std::pow(LambdaEFT, 2)) +
           (yb * OQd_2b22b) / (16. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_1b11b) / (16. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_1b22b) / (24. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b12b) / (48. * std::pow(LambdaEFT, 2));
  }
  else if (i == 7 and j == 7) // psi2psi2
  {
    return (yt * OQu_1b21b) / (24. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b11b) / (12. * std::pow(LambdaEFT, 2)) +
           (yt * OQu_2b22b) / (8. * std::pow(LambdaEFT, 2));
  }
  else
  {
    return 0;
  }
}

bool Class_Potential_R2HDMEFTPHI6_PSI2PHI3::CalculateDebyeSimplified()
{
  // not implemented
  return false;
}

bool Class_Potential_R2HDMEFTPHI6_PSI2PHI3::CalculateDebyeGaugeSimplified()
{
  // not implemented
  return false;
}

double Class_Potential_R2HDMEFTPHI6_PSI2PHI3::VTreeSimplified(
    const std::vector<double> &v) const
{
  (void)v;
  double res = 0;
  // not implemented
  return res;
}

void Class_Potential_R2HDMEFTPHI6_PSI2PHI3::PerformVCTShift()
{
  // not implemented
}

double Class_Potential_R2HDMEFTPHI6_PSI2PHI3::VCounterSimplified(
    const std::vector<double> &v) const
{
  (void)v;
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  // not implemented
  return res;
}

void Class_Potential_R2HDMEFTPHI6_PSI2PHI3::Debugging(
    const std::vector<double> &input,
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

  std::cout << "identities check:\n" << std::endl;

  double v1   = C_vev0 * C_CosBeta;
  double v2   = C_vev0 * C_SinBeta;
  double v1sq = std::pow(v1, 2);
  double v2sq = std::pow(v2, 2);

  std::cout << "00: " << HesseWeinberg(0, 0) << " :: "
            << 1. / v1sq *
                   (v2sq * HesseWeinberg(3, 3) + v1sq * HesseWeinberg(5, 5) -
                    v2sq * HesseWeinberg(7, 7))
            << std::endl;
  std::cout << "02: " << HesseWeinberg(0, 2) << " :: "
            << -1. / v1 *
                   (v2 * HesseWeinberg(3, 3) - v1 * HesseWeinberg(5, 7) -
                    v2 * HesseWeinberg(7, 7))
            << std::endl;
  std::cout << "11: " << HesseWeinberg(1, 1) << " :: " << HesseWeinberg(0, 0)
            << std::endl;
  std::cout << "13: " << HesseWeinberg(1, 3) << " :: " << HesseWeinberg(0, 2)
            << std::endl;
  std::cout << "22: " << HesseWeinberg(2, 2) << " :: " << HesseWeinberg(3, 3)
            << std::endl;

  (void)output;
}

} // namespace Models
} // namespace BSMPT
