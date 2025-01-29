// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialR2HDMEFTPHI6.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>

namespace BSMPT
{
namespace Models
{
Class_Potential_R2HDMEFTPHI6::Class_Potential_R2HDMEFTPHI6(
    const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  // TODO Auto-generated constructor stub
  Model         = ModelID::ModelIDs::R2HDMEFTPHI6;
  NNeutralHiggs = 4;
  NChargedHiggs = 4;

  NHiggs  = NNeutralHiggs + NChargedHiggs;
  NGauge  = 4;
  NLepton = 9;
  NQuarks = 12;

  nPar   = 8;
  nParCT = 20;

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

Class_Potential_R2HDMEFTPHI6::~Class_Potential_R2HDMEFTPHI6()
{
  // TODO Auto-generated destructor stub
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given inputfile
 */
std::vector<std::string> Class_Potential_R2HDMEFTPHI6::addLegendCT() const
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
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dTCP");
  labels.push_back("dTCB");
  labels.push_back("dOp6_111111");
  labels.push_back("dOp6_111122");
  labels.push_back("dOp6_122111");
  labels.push_back("dOp6_121211");
  labels.push_back("dOp6_222222");
  labels.push_back("dOp6_112222");
  labels.push_back("dOp6_122122");
  labels.push_back("dOp6_121222");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * inputfile
 */
std::vector<std::string> Class_Potential_R2HDMEFTPHI6::addLegendTemp() const
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
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given inputfile
 */
std::vector<std::string> Class_Potential_R2HDMEFTPHI6::addLegendVEV() const
{
  std::vector<std::string> labels;
  labels.push_back("omega_CB");
  labels.push_back("omega_1");
  labels.push_back("omega_2");
  labels.push_back("omega_CP");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the EFT
 * parameters.
 */
std::vector<std::string> Class_Potential_R2HDMEFTPHI6::addLegendEFT() const
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
  labels.push_back("L1_shifted");
  labels.push_back("L2_shifted");
  labels.push_back("L4_shifted");
  labels.push_back("L5_shifted");
  labels.push_back("m12Sq_shifted");
  return labels;
}

/**
 * returns the numerical values of the EFT parameters. This has to be
 * specified in the model file.
 */
std::vector<double> Class_Potential_R2HDMEFTPHI6::getParamsEFT() const
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
  valsEFT.push_back(L1);
  valsEFT.push_back(L2);
  valsEFT.push_back(L4);
  valsEFT.push_back(L5);
  valsEFT.push_back(m12Sq);
  return valsEFT;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * higgs couplings. Use this to complement the legend of the given inputfile
 *
 */
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6::addLegendTripleCouplings() const
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

void Class_Potential_R2HDMEFTPHI6::ReadAndSet(const std::string &linestr,
                                              std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;
  double lType = 0, lL1 = 0, lL2 = 0, lL3 = 0, lL4 = 0, lL5 = 0, lm12Sq = 0,
         lTanBeta = 0;
  // double lm11Sq = 0, lm22Sq =0;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 23; k++)
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
    //  lm11Sq = tmp;
    // else if (k == 8)
    //  lm22Sq = tmp;
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
    // else if (k == 15)
    //   Op6_111112 = tmp;
    // else if (k == 16)
    //   Op6_121221 = tmp;
    // else if (k == 17)
    //   Op6_112212 = tmp;
    else if (k == 15)
      Op6_222222 = tmp;
    else if (k == 16)
      Op6_112222 = tmp;
    else if (k == 17)
      Op6_122122 = tmp;
    else if (k == 18)
      Op6_121222 = tmp;
    // else if (k == 22)
    //   Op6_122222 = tmp;
    // else if (k == 23)
    //   Op6_121212 = tmp;
  }

  par[0] = lType;
  par[1] = lL1;
  par[2] = lL2;
  par[3] = lL3;
  par[4] = lL4;
  par[5] = lL5;
  par[6] = lm12Sq;
  par[7] = lTanBeta;
  // m11Sq  = lm11Sq;
  // m22Sq  = lm22Sq;

  set_gen(par);
  return;
}

/**
 * Set Class Object with an CP-Conserving Point
 */
void Class_Potential_R2HDMEFTPHI6::set_gen(const std::vector<double> &par)
{

  // double *p = (double *)par;
  scale            = SMConstants.C_vev0;
  TanBeta          = par[7];
  beta             = std::atan(TanBeta);
  Type             = static_cast<int>(par[0]);
  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = std::sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = std::sqrt(C_SinBetaSquared);

  // corrected Lambdas to absorb CP-even mass shifts due to EFT
  double v1   = SMConstants.C_vev0 * C_CosBeta;
  double v1Sq = SMConstants.C_vev0 * SMConstants.C_vev0 * C_CosBetaSquared;
  double v2   = SMConstants.C_vev0 * C_SinBeta;
  double v2Sq = SMConstants.C_vev0 * SMConstants.C_vev0 * C_SinBetaSquared;

  L1tmp    = par[1];
  L2tmp    = par[2];
  L3       = par[3];
  L4tmp    = par[4];
  L5tmp    = par[5];
  m12Sqtmp = par[6];

  L1 = L1tmp +
       2 * std::pow(LambdaEFT, -0.2e1) *
           (6 * Op6_111111 * v1Sq * v1Sq + 2 * Op6_121211 * v1Sq * v2Sq +
            Op6_122111 * v1Sq * v2Sq -
            (2 * Op6_112222 + 2 * Op6_121222 + Op6_122122) * v2Sq * v2Sq) /
           (4 * v1Sq);
  L2 = L2tmp +
       2 * std::pow(LambdaEFT, -0.2e1) *
           (-2 * Op6_111122 * v1Sq * v1Sq - 2 * Op6_121211 * v1Sq * v1Sq -
            Op6_122111 * v1Sq * v1Sq + 2 * Op6_121222 * v1Sq * v2Sq +
            Op6_122122 * v1Sq * v2Sq + 6 * Op6_222222 * v2Sq * v2Sq) /
           (4 * v2Sq);
  L4 = L4tmp + std::pow(LambdaEFT, -0.2e1) *
                   (Op6_111122 * v1Sq + Op6_121211 * v1Sq + Op6_122111 * v1Sq +
                    Op6_112222 * v2Sq + Op6_121222 * v2Sq + Op6_122122 * v2Sq);
  L5 = L5tmp +
       std::pow(LambdaEFT, -0.2e1) * 0.5 *
           (2 * Op6_111122 * v1Sq + 4 * Op6_121211 * v1Sq + Op6_122111 * v1Sq +
            2 * Op6_112222 * v2Sq + 4 * Op6_121222 * v2Sq + Op6_122122 * v2Sq);
  m12Sq = m12Sqtmp + std::pow(LambdaEFT, -0.2e1) * 0.5 * v1 * v2 *
                         (2 * Op6_111122 * v1Sq + 2 * Op6_121211 * v1Sq +
                          Op6_122111 * v1Sq + 2 * Op6_112222 * v2Sq +
                          2 * Op6_121222 * v2Sq + Op6_122122 * v2Sq);

  m11Sq =
      m12Sq * TanBeta -
      SMConstants.C_vev0 * SMConstants.C_vev0 * C_SinBetaSquared *
          (L4 + L5 + L3) / 0.2e1 -
      SMConstants.C_vev0 * SMConstants.C_vev0 * C_CosBetaSquared * L1 / 0.2e1 +
      (3 * C_CosBetaSquared * C_CosBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_111111 * std::pow(LambdaEFT, -0.2e1)) /
          4. +
      (C_CosBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_111122 * std::pow(LambdaEFT, -0.2e1)) /
          2. +
      (C_SinBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_112222 * std::pow(LambdaEFT, -0.2e1)) /
          4. +
      C_CosBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
          SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
          Op6_121211 * std::pow(LambdaEFT, -0.2e1) +
      (C_SinBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_121222 * std::pow(LambdaEFT, -0.2e1)) /
          2. +
      (C_CosBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_122111 * std::pow(LambdaEFT, -0.2e1)) /
          2. +
      (C_SinBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_122122 * std::pow(LambdaEFT, -0.2e1)) /
          4;
  m22Sq =
      m12Sq * 1.0 / TanBeta -
      SMConstants.C_vev0 * SMConstants.C_vev0 * C_CosBetaSquared *
          (L4 + L5 + L3) / 0.2e1 -
      SMConstants.C_vev0 * SMConstants.C_vev0 * C_SinBetaSquared * L2 / 0.2e1 +
      (C_CosBetaSquared * C_CosBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_111122 * std::pow(LambdaEFT, -0.2e1)) /
          4. +
      (C_CosBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_112222 * std::pow(LambdaEFT, -0.2e1)) /
          2. +
      (C_CosBetaSquared * C_CosBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_121211 * std::pow(LambdaEFT, -0.2e1)) /
          2. +
      C_CosBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
          SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
          Op6_121222 * std::pow(LambdaEFT, -0.2e1) +
      (C_CosBetaSquared * C_CosBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_122111 * std::pow(LambdaEFT, -0.2e1)) /
          4. +
      (C_CosBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_122122 * std::pow(LambdaEFT, -0.2e1)) /
          2. +
      (3 * C_SinBetaSquared * C_SinBetaSquared * SMConstants.C_vev0 *
       SMConstants.C_vev0 * SMConstants.C_vev0 * SMConstants.C_vev0 *
       Op6_222222 * std::pow(LambdaEFT, -0.2e1)) /
          4.;

  double cb = 0;

  if (Type == 1 or Type == 3) // Type I 2HDM oder Lepton Specific
  {
    cb = std::sqrt(2) * SMConstants.C_MassBottom /
         (SMConstants.C_vev0 * C_SinBeta);
  }
  if (Type == 2 or Type == 4) // Type II 2HDM oder Flipped
  {
    cb = std::sqrt(2) * SMConstants.C_MassBottom /
         (SMConstants.C_vev0 * C_CosBeta);
  }
  CTempC1 = 1.0 / 48 *
            (12 * L1 + 8 * L3 + 4 * L4 +
             3 * (3 * SMConstants.C_g * SMConstants.C_g +
                  SMConstants.C_gs * SMConstants.C_gs));
  double ct =
      std::sqrt(2) * SMConstants.C_MassTop / (SMConstants.C_vev0 * C_SinBeta);
  CTempC2 = 1.0 / 48 *
            (12 * L2 + 8 * L3 + 4 * L4 +
             3 * (3 * SMConstants.C_g * SMConstants.C_g +
                  SMConstants.C_gs * SMConstants.C_gs) +
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
  vevTreeMin[1] = SMConstants.C_vev0 * C_CosBeta;
  vevTreeMin[2] = SMConstants.C_vev0 * C_SinBeta;
  vevTreeMin[3] = 0;
  vevTree.resize(NHiggs);
  vevTree = MinimizeOrderVEV(vevTreeMin);
}

void Class_Potential_R2HDMEFTPHI6::set_CT_Pot_Par(const std::vector<double> &p)
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
  dT1         = p[8];
  dT2         = p[9];
  dTCP        = p[10];
  dTCB        = p[11];
  dOp6_111111 = p[12];
  dOp6_111122 = p[13];
  dOp6_122111 = p[13];
  dOp6_121211 = p[14];
  dOp6_222222 = p[15];
  dOp6_112222 = p[16];
  dOp6_122122 = p[17];
  dOp6_121222 = p[18];

  Curvature_Higgs_CT_L1[2] = dTCB;
  Curvature_Higgs_CT_L1[4] = dT1;
  Curvature_Higgs_CT_L1[6] = dT2;
  Curvature_Higgs_CT_L1[7] = dTCP;

  Curvature_Higgs_CT_L2[0][0] = dm11Sq;
  Curvature_Higgs_CT_L2[0][2] = -dm12Sq;
  Curvature_Higgs_CT_L2[1][1] = dm11Sq;
  Curvature_Higgs_CT_L2[1][3] = -dm12Sq;
  Curvature_Higgs_CT_L2[2][0] = -dm12Sq;
  Curvature_Higgs_CT_L2[2][2] = dm22Sq;
  Curvature_Higgs_CT_L2[3][1] = -dm12Sq;
  Curvature_Higgs_CT_L2[3][3] = dm22Sq;
  Curvature_Higgs_CT_L2[4][4] = dm11Sq;
  Curvature_Higgs_CT_L2[4][6] = -dm12Sq;
  Curvature_Higgs_CT_L2[5][5] = dm11Sq;
  Curvature_Higgs_CT_L2[5][7] = -dm12Sq;
  Curvature_Higgs_CT_L2[6][4] = -dm12Sq;
  Curvature_Higgs_CT_L2[6][6] = dm22Sq;
  Curvature_Higgs_CT_L2[7][5] = -dm12Sq;
  Curvature_Higgs_CT_L2[7][7] = dm22Sq;

  Curvature_Higgs_CT_L4[0][0][0][0] = 3 * dL1;
  Curvature_Higgs_CT_L4[0][0][1][1] = dL1;
  Curvature_Higgs_CT_L4[0][0][2][2] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[0][0][3][3] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[0][0][4][4] = dL1;
  Curvature_Higgs_CT_L4[0][0][5][5] = dL1;
  Curvature_Higgs_CT_L4[0][0][6][6] = dL3;
  Curvature_Higgs_CT_L4[0][0][7][7] = dL3;
  Curvature_Higgs_CT_L4[0][1][2][3] = dL5;
  Curvature_Higgs_CT_L4[0][2][4][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][5][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][4][7] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][5][6] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][1][1] = 3 * dL1;
  Curvature_Higgs_CT_L4[1][1][2][2] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[1][1][3][3] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[1][1][4][4] = dL1;
  Curvature_Higgs_CT_L4[1][1][5][5] = dL1;
  Curvature_Higgs_CT_L4[1][1][6][6] = dL3;
  Curvature_Higgs_CT_L4[1][1][7][7] = dL3;
  Curvature_Higgs_CT_L4[1][2][4][7] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][5][6] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][4][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][5][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][2][2] = 3 * dL2;
  Curvature_Higgs_CT_L4[2][2][3][3] = dL2;
  Curvature_Higgs_CT_L4[2][2][4][4] = dL3;
  Curvature_Higgs_CT_L4[2][2][5][5] = dL3;
  Curvature_Higgs_CT_L4[2][2][6][6] = dL2;
  Curvature_Higgs_CT_L4[2][2][7][7] = dL2;
  Curvature_Higgs_CT_L4[3][3][3][3] = 3 * dL2;
  Curvature_Higgs_CT_L4[3][3][4][4] = dL3;
  Curvature_Higgs_CT_L4[3][3][5][5] = dL3;
  Curvature_Higgs_CT_L4[3][3][6][6] = dL2;
  Curvature_Higgs_CT_L4[3][3][7][7] = dL2;
  Curvature_Higgs_CT_L4[4][4][4][4] = 3 * dL1;
  Curvature_Higgs_CT_L4[4][4][5][5] = dL1;
  Curvature_Higgs_CT_L4[4][4][6][6] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[4][4][7][7] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[4][5][6][7] = dL5;
  Curvature_Higgs_CT_L4[5][5][5][5] = 3 * dL1;
  Curvature_Higgs_CT_L4[5][5][6][6] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[5][5][7][7] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[6][6][6][6] = 3 * dL2;
  Curvature_Higgs_CT_L4[6][6][7][7] = dL2;
  Curvature_Higgs_CT_L4[7][7][7][7] = 3 * dL2;

  sym4Dim(Curvature_Higgs_CT_L4, NHiggs, NHiggs, NHiggs, NHiggs);

  Curvature_Higgs_CT_L6[0][0][0][0][0][0] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][0][0][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][0][0][2][2] =
      (-6 * dOp6_111122 - 6 * dOp6_122111 - 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][0][0][3][3] =
      (-6 * dOp6_111122 - 6 * dOp6_122111 + 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][0][0][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][0][0][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][0][0][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][0][0][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][0][1][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][0][0][2][4][6] =
      (-3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][0][2][5][7] =
      (-3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][0][3][4][7] =
      (-3 * dOp6_122111 + 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][0][3][5][6] =
      (3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][1][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][1][1][2][2] =
      (-2 * dOp6_111122 - 2 * dOp6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][1][1][3][3] =
      (-2 * dOp6_111122 - 2 * dOp6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][1][1][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][1][1][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][1][1][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][1][1][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][1][2][4][7] =
      (dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][2][5][6] =
      (-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][3][4][6] =
      (-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][3][5][7] =
      (-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][2][2][2][2] =
      (-12 * dOp6_121222 - 6 * dOp6_122122 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][3][3] =
      (-2 * dOp6_122122 - 2 * dOp6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][4][4] =
      (-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][5][5] =
      (-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][6][6] =
      (-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][7][7] =
      (-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][3][3] =
      (-6 * dOp6_122122 - 6 * dOp6_112222 + 12 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][4][4] =
      (-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][5][5] =
      (-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][6][6] =
      (-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][7][7] =
      (-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][4][4][6][6] =
      (-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][4][4][7][7] =
      (-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][0][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][5][5][6][6] =
      (-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][5][5][7][7] =
      (-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[0][0][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[0][0][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[0][1][1][1][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][1][1][2][4][6] =
      (-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][1][1][2][5][7] =
      (-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][1][1][3][4][7] =
      (-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][1][1][3][5][6] =
      (dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][1][2][2][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[0][1][2][3][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[0][1][2][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][1][2][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][1][2][3][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[0][1][2][3][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[0][2][2][2][4][6] =
      (-3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][2][2][5][7] =
      (-3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][2][3][4][7] =
      (2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][2][3][5][6] =
      (-2 * dOp6_121222 + dOp6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][3][3][4][6] =
      (-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][3][3][5][7] =
      (-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][4][4][6] =
      (-3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][4][5][7] =
      (-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][5][5][6] =
      (-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][6][6][6] =
      (-3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][6][7][7] =
      (-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][5][5][5][7] =
      (-3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][5][6][6][7] =
      (-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][5][7][7][7] =
      (-3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][3][3][3][4][7] =
      (-3 * dOp6_122122 + 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][3][3][3][5][6] =
      (3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][4][4][7] =
      (-3 * dOp6_122111 + 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][4][5][6] =
      (dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][5][5][7] =
      (-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][6][6][7] =
      (-dOp6_122122 + 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][7][7][7] =
      (-3 * dOp6_122122 + 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][3][5][5][5][6] =
      (3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][5][6][6][6] =
      (3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][5][6][7][7] =
      (dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][1][1][1] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][1][1][2][2] =
      (-6 * dOp6_111122 - 6 * dOp6_122111 + 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][1][1][3][3] =
      (-6 * dOp6_111122 - 6 * dOp6_122111 - 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][1][1][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][1][1][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][1][1][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[1][1][1][1][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[1][1][1][2][4][7] =
      (3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][2][5][6] =
      (-3 * dOp6_122111 + 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][3][4][6] =
      (-3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][3][5][7] =
      (-3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][1][2][2][2][2] =
      (-6 * dOp6_122122 + 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][3][3] =
      (-2 * dOp6_122122 - 2 * dOp6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][4][4] =
      (-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][5][5] =
      (-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][6][6] =
      (-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][7][7] =
      (-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][3][3] =
      (-6 * dOp6_122122 - 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][4][4] =
      (-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][5][5] =
      (-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][6][6] =
      (-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][7][7] =
      (-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][4][4][6][6] =
      (-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][4][4][7][7] =
      (-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[1][1][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][5][5][6][6] =
      (-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][5][5][7][7] =
      (-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[1][1][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[1][1][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[1][2][2][2][4][7] =
      (3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][2][5][6] =
      (-3 * dOp6_122122 + 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][3][4][6] =
      (-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][3][5][7] =
      (-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][3][3][4][7] =
      (dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][3][3][5][6] =
      (-dOp6_122122 + 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][4][4][7] =
      (3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][4][5][6] =
      (-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][5][5][7] =
      (dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][6][6][7] =
      (dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][7][7][7] =
      (3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][5][5][5][6] =
      (-3 * dOp6_122111 + 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][5][6][6][6] =
      (-3 * dOp6_122122 + 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][5][6][7][7] =
      (-dOp6_122122 + 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][3][3][4][6] =
      (-3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][3][3][5][7] =
      (-3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][4][4][6] =
      (-3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][4][5][7] =
      (-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][5][5][6] =
      (-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][6][6][6] =
      (-3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][6][7][7] =
      (-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][5][5][5][7] =
      (-3 * dOp6_122111 - 6 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][5][6][6][7] =
      (-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][5][7][7][7] =
      (-3 * dOp6_122122 - 6 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[2][2][2][2][2][2] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][2][2][3][3] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][2][2][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[2][2][2][2][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[2][2][2][2][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][2][2][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][3][3][3][3] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][3][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[2][2][3][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[2][2][3][3][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][3][3][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[2][2][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[2][2][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[2][2][4][4][6][6] =
      (-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][4][4][7][7] =
      (-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[2][2][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[2][2][5][5][6][6] =
      (-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][5][5][7][7] =
      (-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
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
  Curvature_Higgs_CT_L6[3][3][3][3][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[3][3][3][3][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][3][3][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[3][3][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[3][3][4][4][6][6] =
      (-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][4][4][7][7] =
      (-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[3][3][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[3][3][5][5][6][6] =
      (-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][5][5][7][7] =
      (-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][6][6][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][6][6][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[3][3][7][7][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_222222;
  Curvature_Higgs_CT_L6[4][4][4][4][4][4] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[4][4][4][4][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[4][4][4][4][6][6] =
      (-6 * dOp6_111122 - 6 * dOp6_122111 - 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][4][4][7][7] =
      (-6 * dOp6_111122 - 6 * dOp6_122111 + 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][4][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[4][4][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[4][4][5][5][6][6] =
      (-2 * dOp6_111122 - 2 * dOp6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][5][5][7][7] =
      (-2 * dOp6_111122 - 2 * dOp6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][6][6][6][6] =
      (-6 * dOp6_122122 - 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][6][6][7][7] =
      (-2 * dOp6_122122 - 2 * dOp6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][7][7][7][7] =
      (-6 * dOp6_122122 + 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][5][5][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[4][5][6][6][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[4][5][6][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[5][5][5][5][5][5] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[5][5][5][5][6][6] =
      (-6 * dOp6_111122 - 6 * dOp6_122111 + 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][5][5][7][7] =
      (-6 * dOp6_111122 - 6 * dOp6_122111 - 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][6][6][6][6] =
      (-6 * dOp6_122122 + 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][6][6][7][7] =
      (-2 * dOp6_122122 - 2 * dOp6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][7][7][7][7] =
      (-6 * dOp6_122122 - 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
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
void Class_Potential_R2HDMEFTPHI6::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "scale = " << scale << "\n";

  ss << "The parameters are : \n";
  ss << "Model = " << Model << "\n";
  ss << "v1 = " << SMConstants.C_vev0 * C_CosBeta << "\n";
  ss << "v2 = " << SMConstants.C_vev0 * C_SinBeta << "\n";
  ss << "Type = " << Type << "\n";

  ss << "beta = " << beta << "\n";
  ss << "tan(beta) = " << TanBeta << "\n";
  ss << "L1 = " << L1tmp << " ---> " << L1 << "\n";
  ss << "L2 = " << L2tmp << " ---> " << L2 << "\n";
  ss << "L3 = " << L3 << "\n";
  ss << "L4 = " << L4tmp << " ---> " << L4 << "\n";
  ss << "Re(L5) = " << L5tmp << " ---> " << L5 << "\n";
  ss << "Re(m_12^2) = " << m12Sqtmp << " ---> " << m12Sq << "\n";
  ss << "m_{11}^2 = " << m11Sq << "\n";
  ss << "m_{22}^2 = " << m22Sq << "\n";

  ss << "The counterterms are :\n";

  ss << "dL1 := " << dL1 << ";\n";
  ss << "dL2 := " << dL2 << ";\n";
  ss << "dL3 := " << dL3 << ";\n";
  ss << "dL4 := " << dL4 << ";\n";
  ss << "dL5 := " << dL5 << ";\n";
  ss << "dm11Sq := " << dm11Sq << ";\n";
  ss << "dm22Sq := " << dm22Sq << ";\n";
  ss << "dm12Sq := " << dm12Sq << ";\n";

  ss << "dT1 := " << dT1 << ";\n";
  ss << "dT2 := " << dT2 << ";\n";
  ss << "dTCP := " << dTCP << ";\n";
  ss << "dTCB:= " << dTCB << ";\n";

  ss << "dOp6_111111 := " << dOp6_111111 << ";\n";
  ss << "dOp6_111122 := " << dOp6_111122 << ";\n";
  ss << "dOp6_122111 := " << dOp6_122111 << ";\n";
  ss << "dOp6_121211 := " << dOp6_121211 << ";\n";
  ss << "dOp6_222222 := " << dOp6_222222 << ";\n";
  ss << "dOp6_112222 := " << dOp6_112222 << ";\n";
  ss << "dOp6_122122 := " << dOp6_122122 << ";\n";
  ss << "dOp6_121222 := " << dOp6_121222 << ";\n";

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
     << "\n";
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

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms in the 2HDM
 */
std::vector<double> Class_Potential_R2HDMEFTPHI6::calc_CT() const
{
  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsDone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  double v1 = SMConstants.C_vev0 * C_CosBeta;
  double v2 = SMConstants.C_vev0 * C_SinBeta;

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

  double t = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0, s6 = 0, s7 = 0,
         s8 = 0; // Values of dL4 and dOp6_ijklmn

  // dm11Sq
  parCT.push_back(
      -((-2 * t * v1 * v2 * v2 + 5 * HesseWeinberg(0, 0) * v1 +
         HesseWeinberg(1, 3) * v2 - HesseWeinberg(4, 6) * v2 -
         HesseWeinberg(4, 4) * v1 - 2 * HesseWeinberg(5, 5) * v1) /
        v1) /
          0.2e1 -
      0.3e1 / 0.4e1 * std::pow(LambdaEFT, -2) * std::pow(v1, 4) * s1 -
      std::pow(LambdaEFT, -2) * v1 * v1 * s2 * v2 * v2 / 2 -
      std::pow(LambdaEFT, -2) * v1 * v1 * s3 * v2 * v2 -
      std::pow(LambdaEFT, -2) * v1 * v1 * s4 * v2 * v2 -
      std::pow(LambdaEFT, -2) * s6 * std::pow(v2, 4) / 4 -
      0.3e1 / 0.4e1 * std::pow(LambdaEFT, -2) * s7 * std::pow(v2, 4) -
      std::pow(LambdaEFT, -2) * s8 * std::pow(v2, 4) / 2);
  // dm22Sq
  parCT.push_back(
      ((2 * t * v1 * v1 * v2 * v2 + HesseWeinberg(6, 6) * v2 * v2 -
        2 * HesseWeinberg(0, 0) * v1 * v1 - HesseWeinberg(1, 3) * v1 * v2 -
        3 * HesseWeinberg(3, 3) * v2 * v2 + HesseWeinberg(4, 6) * v1 * v2 +
        2 * v1 * v1 * HesseWeinberg(5, 5)) *
       std::pow(v2, (-2))) /
          0.2e1 -
      std::pow(LambdaEFT, -2) * s2 * std::pow(v1, 4) / 4 -
      0.3e1 / 0.4e1 * std::pow(LambdaEFT, -2) * s3 * std::pow(v1, 4) -
      std::pow(LambdaEFT, -2) * s4 * std::pow(v1, 4) / 2 -
      0.3e1 / 0.4e1 * std::pow(LambdaEFT, -2) * std::pow(v2, 4) * s5 -
      std::pow(LambdaEFT, -2) * v2 * v2 * s6 * v1 * v1 / 2 -
      std::pow(LambdaEFT, -2) * v2 * v2 * s7 * v1 * v1 -
      std::pow(LambdaEFT, -2) * v2 * v2 * s8 * v1 * v1);
  // dm12Sq
  parCT.push_back(-(-t * v1 * v2 * v2 + HesseWeinberg(0, 0) * v1 -
                    HesseWeinberg(1, 3) * v2 - HesseWeinberg(5, 5) * v1) /
                      v2 -
                  v2 * std::pow(v1, 3) * std::pow(LambdaEFT, -2) * s3 / 2 -
                  std::pow(v2, 3) * v1 * std::pow(LambdaEFT, -2) * s7 / 2);
  // dL1
  parCT.push_back((-t * v2 * v2 + 2 * HesseWeinberg(0, 0) -
                   HesseWeinberg(4, 4) - HesseWeinberg(5, 5)) *
                      std::pow(v1, -0.2e1) +
                  3 * std::pow(LambdaEFT, -2) * v1 * v1 * s1 +
                  std::pow(LambdaEFT, -2) * v2 * v2 * s2 +
                  0.3e1 / 0.2e1 * std::pow(LambdaEFT, -2) * v2 * v2 * s3 +
                  2 * std::pow(LambdaEFT, -2) * v2 * v2 * s4 +
                  std::pow(LambdaEFT, -2) * std::pow(v1, -2) * s7 *
                      std::pow(v2, 4) / 2);
  // dL2
  parCT.push_back(
      -(t * v1 * v1 * v2 * v2 + HesseWeinberg(6, 6) * v2 * v2 -
        HesseWeinberg(0, 0) * v1 * v1 - HesseWeinberg(3, 3) * v2 * v2 +
        v1 * v1 * HesseWeinberg(5, 5)) *
          std::pow(v2, -0.4e1) +
      std::pow(v2, -2) * std::pow(v1, 4) * std::pow(LambdaEFT, -2) * s3 / 2 +
      0.3e1 / 0.2e1 * v1 * v1 * std::pow(LambdaEFT, -2) * s7 +
      v1 * v1 * std::pow(LambdaEFT, -2) * s6 +
      2 * v1 * v1 * std::pow(LambdaEFT, -2) * s8 +
      3 * v2 * v2 * std::pow(LambdaEFT, -2) * s5);
  // dL3
  parCT.push_back((-t * v1 * v2 * v2 + HesseWeinberg(0, 0) * v1 +
                   HesseWeinberg(1, 3) * v2 - HesseWeinberg(4, 6) * v2 -
                   HesseWeinberg(5, 5) * v1) /
                      v1 * std::pow(v2, -0.2e1) +
                  std::pow(LambdaEFT, -2) * v1 * v1 * s2 +
                  std::pow(LambdaEFT, -2) * s3 * v1 * v1 +
                  std::pow(LambdaEFT, -2) * s4 * v1 * v1 +
                  std::pow(LambdaEFT, -2) * v2 * v2 * s6 +
                  std::pow(LambdaEFT, -2) * v2 * v2 * s7 +
                  std::pow(LambdaEFT, -2) * v2 * v2 * s8);
  // dL4
  parCT.push_back(t);
  // dL5
  parCT.push_back(
      -(-t * v2 * v2 + 2 * HesseWeinberg(0, 0) - 2 * HesseWeinberg(5, 5)) *
          std::pow(v2, (-2)) +
      -std::pow(LambdaEFT, -2) * s3 * v1 * v1 / 2 +
      std::pow(LambdaEFT, -2) * s4 * v1 * v1 -
      std::pow(LambdaEFT, -2) * v2 * v2 * s7 / 2 +
      std::pow(LambdaEFT, -2) * v2 * v2 * s8);

  // dT1
  double tmp =
      HesseWeinberg(1, 3) * v2 + HesseWeinberg(0, 0) * v1 - NablaWeinberg(4);
  if (std::abs(tmp) < 1e-9) tmp = 0;
  parCT.push_back(tmp);
  std::cout << std::cout.precision(20)
            << "HesseWeinberg(1, 3) = " << HesseWeinberg(1, 3) << std::endl;
  std::cout << "HesseWeinberg(0, 0) = " << HesseWeinberg(0, 0) << std::endl;
  std::cout << "NablaWeinberg(4) = " << NablaWeinberg(4) << std::endl;
  std::cout << "dT1 = "
            << HesseWeinberg(1, 3) * v2 + HesseWeinberg(0, 0) * v1 -
                   NablaWeinberg(4)
            << std::endl;
  // dT2
  tmp = HesseWeinberg(1, 3) * v1 + HesseWeinberg(3, 3) * v2 - NablaWeinberg(6);
  if (std::abs(tmp) < 1e-9) tmp = 0;
  parCT.push_back(tmp);
  std::cout << std::cout.precision(20)
            << "HesseWeinberg(1, 3) = " << HesseWeinberg(1, 3) << std::endl;
  std::cout << "HesseWeinberg(3, 3) = " << HesseWeinberg(3, 3) << std::endl;
  std::cout << "NablaWeinberg(6) = " << NablaWeinberg(6) << std::endl;
  std::cout << "dT2 = "
            << HesseWeinberg(1, 3) * v1 + HesseWeinberg(3, 3) * v2 -
                   NablaWeinberg(6)
            << std::endl;
  // dTCP
  tmp = -NablaWeinberg(7);
  if (std::abs(tmp) < 1e-9) tmp = 0;
  parCT.push_back(tmp);
  // dTCB
  tmp = -NablaWeinberg(2);
  if (std::abs(tmp) < 1e-9) tmp = 0;
  parCT.push_back(tmp);
  // dOp6_111111
  parCT.push_back(s1);
  // dOp6_111122
  parCT.push_back(s2);
  // dOp6_122111
  parCT.push_back(s3);
  // dOp6_121211
  parCT.push_back(s4);
  // dOp6_222222
  parCT.push_back(s5);
  // dOp6_112222
  parCT.push_back(s6);
  // dOp6_122122
  parCT.push_back(s7);
  // dOp6_121222
  parCT.push_back(s8);

  //	double Identities[5];
  //	Identities[0] = HesseWeinberg(0, 0) - HesseWeinberg(1, 1);
  //	Identities[1] = -HesseWeinberg(3, 3) + HesseWeinberg(2, 2);
  //	Identities[2] = (HesseWeinberg(1, 1) * v1 - HesseWeinberg(5, 5) * v1 +
  // HesseWeinberg(1, 3) * v2 - HesseWeinberg(5, 7) * v2) / v2; 	Identities[3]
  // = -HesseWeinberg(0, 2) + HesseWeinberg(1, 3); 	Identities[4] = -1 / v2 *
  //(HesseWeinberg(5, 7) * v1 + HesseWeinberg(7, 7) * v2 - HesseWeinberg(1, 3) *
  // v1 - HesseWeinberg(3, 3) * v2);

  return parCT;
}

void Class_Potential_R2HDMEFTPHI6::AdjustRotationMatrix()
{
}

/**
 * Calculates the corrections to the Triple higgs couplings in the mass basis.
 *
 * Use the vector TripleHiggsCorrectionsCWPhysical to save your couplings and
 * set the nTripleCouplings to the number of couplings you want as output.
 */
void Class_Potential_R2HDMEFTPHI6::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

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

void Class_Potential_R2HDMEFTPHI6::SetCurvatureArrays()
{
  initVectors();

  HiggsVev[4] = SMConstants.C_vev0 * C_CosBeta;
  HiggsVev[6] = SMConstants.C_vev0 * C_SinBeta;

  Curvature_Higgs_L2[0][0] = m11Sq;
  Curvature_Higgs_L2[0][2] = -m12Sq;
  Curvature_Higgs_L2[1][1] = m11Sq;
  Curvature_Higgs_L2[1][3] = -m12Sq;
  Curvature_Higgs_L2[2][0] = -m12Sq;
  Curvature_Higgs_L2[2][2] = m22Sq;
  Curvature_Higgs_L2[3][1] = -m12Sq;
  Curvature_Higgs_L2[3][3] = m22Sq;
  Curvature_Higgs_L2[4][4] = m11Sq;
  Curvature_Higgs_L2[4][6] = -m12Sq;
  Curvature_Higgs_L2[5][5] = m11Sq;
  Curvature_Higgs_L2[5][7] = -m12Sq;
  Curvature_Higgs_L2[6][4] = -m12Sq;
  Curvature_Higgs_L2[6][6] = m22Sq;
  Curvature_Higgs_L2[7][5] = -m12Sq;
  Curvature_Higgs_L2[7][7] = m22Sq;

  Curvature_Higgs_L4[0][0][0][0] = 3 * L1;
  Curvature_Higgs_L4[0][0][1][1] = L1;
  Curvature_Higgs_L4[0][0][2][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][0][4][4] = L1;
  Curvature_Higgs_L4[0][0][5][5] = L1;
  Curvature_Higgs_L4[0][0][6][6] = L3;
  Curvature_Higgs_L4[0][0][7][7] = L3;
  Curvature_Higgs_L4[0][1][2][3] = L5;
  Curvature_Higgs_L4[0][2][4][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][2][5][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][3][4][7] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][6] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][1][1][1] = 3 * L1;
  Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][1][3][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][1][4][4] = L1;
  Curvature_Higgs_L4[1][1][5][5] = L1;
  Curvature_Higgs_L4[1][1][6][6] = L3;
  Curvature_Higgs_L4[1][1][7][7] = L3;
  Curvature_Higgs_L4[1][2][4][7] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][6] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[1][3][4][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][3][5][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][2][2][2] = 3 * L2;
  Curvature_Higgs_L4[2][2][3][3] = L2;
  Curvature_Higgs_L4[2][2][4][4] = L3;
  Curvature_Higgs_L4[2][2][5][5] = L3;
  Curvature_Higgs_L4[2][2][6][6] = L2;
  Curvature_Higgs_L4[2][2][7][7] = L2;
  Curvature_Higgs_L4[3][3][3][3] = 3 * L2;
  Curvature_Higgs_L4[3][3][4][4] = L3;
  Curvature_Higgs_L4[3][3][5][5] = L3;
  Curvature_Higgs_L4[3][3][6][6] = L2;
  Curvature_Higgs_L4[3][3][7][7] = L2;
  Curvature_Higgs_L4[4][4][4][4] = 3 * L1;
  Curvature_Higgs_L4[4][4][5][5] = L1;
  Curvature_Higgs_L4[4][4][6][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][5][6][7] = L5;
  Curvature_Higgs_L4[5][5][5][5] = 3 * L1;
  Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][5][7][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][6][6][6] = 3 * L2;
  Curvature_Higgs_L4[6][6][7][7] = L2;
  Curvature_Higgs_L4[7][7][7][7] = 3 * L2;

  sym4Dim(Curvature_Higgs_L4, NHiggs, NHiggs, NHiggs, NHiggs);

  Curvature_Higgs_L6[0][0][0][0][0][0] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][0][0][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][0][0][2][2] =
      (-6 * Op6_111122 - 6 * Op6_122111 - 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][0][0][3][3] =
      (-6 * Op6_111122 - 6 * Op6_122111 + 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][0][0][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][0][0][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][0][0][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][0][0][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][0][1][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][0][0][2][4][6] =
      (-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][0][2][5][7] =
      (-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][0][3][4][7] =
      (-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][0][3][5][6] =
      (3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][1][1][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][1][1][2][2] =
      (-2 * Op6_111122 - 2 * Op6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][1][1][3][3] =
      (-2 * Op6_111122 - 2 * Op6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][1][1][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][1][1][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][1][1][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][1][1][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][1][2][4][7] =
      (Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][1][2][5][6] =
      (-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][1][3][4][6] =
      (-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][1][3][5][7] =
      (-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][2][2][2][2] =
      (-6 * Op6_112222 - 6 * Op6_122122 - 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][3][3] =
      (-2 * Op6_112222 - 2 * Op6_122122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][4][4] =
      (-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][5][5] =
      (-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][6][6] =
      (-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][7][7] =
      (-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][3][3] =
      (-6 * Op6_122122 - 6 * Op6_112222 + 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][4][4] =
      (-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][5][5] =
      (-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][6][6] =
      (-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][7][7] =
      (-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][4][4][6][6] =
      (-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][4][4][7][7] =
      (-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][0][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][5][5][6][6] =
      (-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][5][5][7][7] =
      (-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[0][0][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[0][0][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[0][1][1][1][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][1][1][2][4][6] =
      (-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][1][1][2][5][7] =
      (-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][1][1][3][4][7] =
      (-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][1][1][3][5][6] =
      (Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][1][2][2][2][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[0][1][2][3][3][3] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[0][1][2][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][1][2][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][1][2][3][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[0][1][2][3][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[0][2][2][2][4][6] =
      (-6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][2][2][5][7] =
      (-6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][2][3][4][7] =
      (2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][2][3][5][6] =
      (-2 * Op6_121222 + Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][3][3][4][6] =
      (-2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][3][3][5][7] =
      (-2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][4][4][4][6] =
      (-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][4][4][5][7] =
      (-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][4][5][5][6] =
      (-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][4][6][6][6] =
      (-6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][4][6][7][7] =
      (-2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][5][5][5][7] =
      (-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][5][6][6][7] =
      (-2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][5][7][7][7] =
      (-6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][3][3][4][7] =
      (6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][3][3][5][6] =
      (-6 * Op6_121222 + 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][4][4][4][7] =
      (-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][4][4][5][6] =
      (Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][4][5][5][7] =
      (-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][4][6][6][7] =
      (2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][4][7][7][7] =
      (6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][5][5][5][6] =
      (3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][5][6][6][6] =
      (-6 * Op6_121222 + 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][5][6][7][7] =
      (-2 * Op6_121222 + Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][1][1][1][1][1] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][1][1][2][2] =
      (-6 * Op6_111122 - 6 * Op6_122111 + 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][1][1][3][3] =
      (-6 * Op6_111122 - 6 * Op6_122111 - 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][1][1][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][1][1][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][1][1][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[1][1][1][1][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[1][1][1][2][4][7] =
      (3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][1][1][2][5][6] =
      (-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][1][1][3][4][6] =
      (-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][1][1][3][5][7] =
      (-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][1][2][2][2][2] =
      (-6 * Op6_122122 - 6 * Op6_112222 + 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][3][3] =
      (-2 * Op6_122122 - 2 * Op6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][4][4] =
      (-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][5][5] =
      (-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][6][6] =
      (-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][7][7] =
      (-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][3][3] =
      (-6 * Op6_112222 - 6 * Op6_122122 - 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][4][4] =
      (-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][5][5] =
      (-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][6][6] =
      (-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][7][7] =
      (-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][4][4][6][6] =
      (-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][4][4][7][7] =
      (-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[1][1][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][5][5][6][6] =
      (-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][5][5][7][7] =
      (-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[1][1][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[1][1][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[1][2][2][2][4][7] =
      (3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][2][2][5][6] =
      (-3 * Op6_122122 + 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][2][3][4][6] =
      (-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][2][3][5][7] =
      (-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][3][3][4][7] =
      (-2 * Op6_121222 + Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][3][3][5][6] =
      (-Op6_122122 + 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][4][4][4][7] =
      (3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][4][4][5][6] =
      (-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][4][5][5][7] =
      (Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][4][6][6][7] =
      (-2 * Op6_121222 + Op6_122122) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][4][7][7][7] =
      (3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][5][5][5][6] =
      (-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][5][6][6][6] =
      (-3 * Op6_122122 + 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][5][6][7][7] =
      (-Op6_122122 + 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][3][3][4][6] =
      (-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][3][3][5][7] =
      (-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][4][4][4][6] =
      (-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][4][4][5][7] =
      (-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][4][5][5][6] =
      (-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][4][6][6][6] =
      (-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][4][6][7][7] =
      (-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][5][5][5][7] =
      (-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][5][6][6][7] =
      (-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][3][5][7][7][7] =
      (-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[2][2][2][2][2][2] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][2][2][3][3] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][2][2][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[2][2][2][2][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[2][2][2][2][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][2][2][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][3][3][3][3] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][3][3][4][4] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[2][2][3][3][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[2][2][3][3][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][3][3][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[2][2][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[2][2][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[2][2][4][4][6][6] =
      (-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][4][4][7][7] =
      (-2 * Op6_112222 - Op6_122122 + 2 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[2][2][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[2][2][5][5][6][6] =
      (-2 * Op6_112222 - Op6_122122 + 2 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][5][5][7][7] =
      (-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
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
  Curvature_Higgs_L6[3][3][3][3][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[3][3][3][3][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][3][3][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][4][4][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[3][3][4][4][5][5] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[3][3][4][4][6][6] =
      (-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][4][4][7][7] =
      (-2 * Op6_112222 - Op6_122122 + 2 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[3][3][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[3][3][5][5][6][6] =
      (-2 * Op6_112222 - Op6_122122 + 2 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][5][5][7][7] =
      (-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][6][6][6][6] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][6][6][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[3][3][7][7][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[4][4][4][4][4][4] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[4][4][4][4][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[4][4][4][4][6][6] =
      (-6 * Op6_111122 - 6 * Op6_122111 - 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][4][4][7][7] =
      (-6 * Op6_111122 - 6 * Op6_122111 + 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][4][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[4][4][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[4][4][5][5][6][6] =
      (-2 * Op6_111122 - 2 * Op6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][5][5][7][7] =
      (-2 * Op6_111122 - 2 * Op6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][6][6][6][6] =
      (-6 * Op6_112222 - 6 * Op6_122122 - 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][6][6][7][7] =
      (-2 * Op6_112222 - 2 * Op6_122122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][7][7][7][7] =
      (-6 * Op6_122122 - 6 * Op6_112222 + 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][5][5][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[4][5][6][6][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[4][5][6][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[5][5][5][5][5][5] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[5][5][5][5][6][6] =
      (-6 * Op6_111122 - 6 * Op6_122111 + 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][5][5][7][7] =
      (-6 * Op6_111122 - 6 * Op6_122111 - 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][6][6][6][6] =
      (-6 * Op6_122122 - 6 * Op6_112222 + 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][6][6][7][7] =
      (-2 * Op6_112222 - 2 * Op6_122122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][7][7][7][7] =
      (-6 * Op6_112222 - 6 * Op6_122122 - 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[6][6][6][6][6][6] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[6][6][6][6][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[6][6][7][7][7][7] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;
  Curvature_Higgs_L6[7][7][7][7][7][7] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_222222;

  sym6Dim(Curvature_Higgs_L6, NHiggs, NHiggs, NHiggs, NHiggs, NHiggs, NHiggs);

  Curvature_Gauge_G2H2[0][0][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;

  Curvature_Gauge_G2H2[0][3][0][4] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][5] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][6] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][7] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][4][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][5][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][6][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][7][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;

  Curvature_Gauge_G2H2[1][1][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;

  Curvature_Gauge_G2H2[1][3][0][5] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][7] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][4][1] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][5][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][6][3] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][7][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;

  Curvature_Gauge_G2H2[2][2][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;

  Curvature_Gauge_G2H2[2][3][0][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][4][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][5][5] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][6][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][7][7] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][4] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][5] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][6] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][7] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][4][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][5][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][6][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][7][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][5] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][7] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][4][1] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][5][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][6][3] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][7][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][4][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][5][5] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][6][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][7][7] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][0][0] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][1][1] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][2][2] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][3][3] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][4][4] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][5][5] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][6][6] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][7][7] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;

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

  double v1 = SMConstants.C_vev0 * C_CosBeta;
  double v2 = SMConstants.C_vev0 * C_SinBeta;
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

  YIJR2(0, 9)  = -std::conj(V11) * SMConstants.C_MassUp / v2;
  YIJR2(0, 10) = -std::conj(V12) * SMConstants.C_MassUp / v2;
  YIJR2(0, 11) = -std::conj(V13) * SMConstants.C_MassUp / v2;

  YIJR2(1, 9)  = -std::conj(V21) * SMConstants.C_MassCharm / v2;
  YIJR2(1, 10) = -std::conj(V22) * SMConstants.C_MassCharm / v2;
  YIJR2(1, 11) = -std::conj(V23) * SMConstants.C_MassCharm / v2;

  YIJR2(2, 9)  = -std::conj(V31) * SMConstants.C_MassTop / v2;
  YIJR2(2, 10) = -std::conj(V32) * SMConstants.C_MassTop / v2;
  YIJR2(2, 11) = -std::conj(V33) * SMConstants.C_MassTop / v2;

  YIJS2(0, 6) = SMConstants.C_MassUp / v2;
  YIJS2(1, 7) = SMConstants.C_MassCharm / v2;
  YIJS2(2, 8) = SMConstants.C_MassTop / v2;

  YIJSD(3, 9)  = SMConstants.C_MassDown / vD;
  YIJSD(4, 10) = SMConstants.C_MassStrange / vD;
  YIJSD(5, 11) = SMConstants.C_MassBottom / vD;

  YIJRD(3, 6) = V11 * SMConstants.C_MassDown / vD;
  YIJRD(3, 7) = V21 * SMConstants.C_MassDown / vD;
  YIJRD(3, 8) = V31 * SMConstants.C_MassDown / vD;
  YIJRD(4, 6) = V12 * SMConstants.C_MassStrange / vD;
  YIJRD(4, 7) = V22 * SMConstants.C_MassStrange / vD;
  YIJRD(4, 8) = V32 * SMConstants.C_MassStrange / vD;
  YIJRD(5, 6) = V13 * SMConstants.C_MassBottom / vD;
  YIJRD(5, 7) = V23 * SMConstants.C_MassBottom / vD;
  YIJRD(5, 8) = V33 * SMConstants.C_MassBottom / vD;

  YIJRL(1, 6) = SMConstants.C_MassElectron / vL;
  YIJRL(3, 7) = SMConstants.C_MassMu / vL;
  YIJRL(5, 8) = SMConstants.C_MassTau / vL;

  YIJSL(0, 1) = SMConstants.C_MassElectron / vL;
  YIJSL(2, 3) = SMConstants.C_MassMu / vL;
  YIJSL(4, 5) = SMConstants.C_MassTau / vL;

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

  SetCurvatureDone = true;
}

double Class_Potential_R2HDMEFTPHI6::SymFac_Higgs_OneLoop(
    const int &i,
    const int &j,
    const std::vector<double> &point) const
{
  (void)i;
  (void)j;
  (void)point;
  return 0;
}

double Class_Potential_R2HDMEFTPHI6::SymFac_Higgs_TwoLoop(const int &i,
                                                          const int &j) const
{
  if (i == 0 and j == 0) // rho1rho1
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-2 * (6 * Op6_111111 + 2 * Op6_111122 + Op6_112222 + Op6_122111) -
            Op6_122122) /
           48.;
  }
  else if ((i == 0 and j == 2) or (i == 2 and j == 0)) // rho1rho2
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-4 * Op6_111112 - 3 * Op6_112212 -
            2 * (Op6_121221 + 2 * Op6_122222)) /
           48.;
  }
  else if (i == 2 and j == 2) // rho2rho2
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-2 * Op6_111122 - 4 * Op6_112222 - Op6_122111 -
            2 * (Op6_122122 + 6 * Op6_222222)) /
           48.;
  }
  else if (i == 1 and j == 1) // eta1eta1
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-2 * (6 * Op6_111111 + 2 * Op6_111122 + Op6_112222 + Op6_122111) -
            Op6_122122) /
           48.;
  }
  else if ((i == 1 and j == 3) or (i == 3 and j == 1)) // eta1eta2
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-4 * Op6_111112 - 3 * Op6_112212 -
            2 * (Op6_121221 + 2 * Op6_122222)) /
           48.;
  }
  else if (i == 3 and j == 3) // eta2eta2
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-2 * Op6_111122 - 4 * Op6_112222 - Op6_122111 -
            2 * (Op6_122122 + 6 * Op6_222222)) /
           48.;
  }
  else if (i == 4 and j == 4) // zeta1zeta1
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-2 * (6 * Op6_111111 + 2 * Op6_111122 + Op6_112222 + Op6_122111) -
            Op6_122122) /
           48.;
  }
  else if ((i == 4 and j == 6) or (i == 6 and j == 4)) // zeta1zeta2
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-4 * Op6_111112 - 3 * Op6_112212 -
            2 * (Op6_121221 + 2 * Op6_122222)) /
           48.;
  }
  else if (i == 6 and j == 6) // zeta2zeta2
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-2 * Op6_111122 - 4 * Op6_112222 - Op6_122111 -
            2 * (Op6_122122 + 6 * Op6_222222)) /
           48.;
  }
  else if (i == 5 and j == 5) // psi1psi1
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-2 * (6 * Op6_111111 + 2 * Op6_111122 + Op6_112222 + Op6_122111) -
            Op6_122122) /
           48.;
  }
  else if ((i == 5 and j == 7) or (i == 7 and j == 5)) // psi1psi2
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-4 * Op6_111112 - 3 * Op6_112212 -
            2 * (Op6_121221 + 2 * Op6_122222)) /
           48.;
  }
  else if (i == 7 and j == 7)
  {
    return std::pow(LambdaEFT, -0.2e1) *
           (-2 * Op6_111122 - 4 * Op6_112222 - Op6_122111 -
            2 * (Op6_122122 + 6 * Op6_222222)) /
           48.;
  }
  else
  {
    return 0;
  }
}

bool Class_Potential_R2HDMEFTPHI6::CalculateDebyeSimplified()
{
  // not implemented
  return false;
}

bool Class_Potential_R2HDMEFTPHI6::CalculateDebyeGaugeSimplified()
{
  // not implemented
  return false;
}

double Class_Potential_R2HDMEFTPHI6::VTreeSimplified(
    const std::vector<double> &v) const
{
  (void)v;
  double res = 0;
  // not implemented
  return res;
}

double Class_Potential_R2HDMEFTPHI6::VCounterSimplified(
    const std::vector<double> &v) const
{
  (void)v;
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  // not implemented
  return res;
}

void Class_Potential_R2HDMEFTPHI6::PerformVCTShift()
{
  // not implemented
}

void Class_Potential_R2HDMEFTPHI6::Debugging(const std::vector<double> &input,
                                             std::vector<double> &output) const
{
  (void)input;
  // not implemented
  (void)output;
}

} // namespace Models
} // namespace BSMPT
