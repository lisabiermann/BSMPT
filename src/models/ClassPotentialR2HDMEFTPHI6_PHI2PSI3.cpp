// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialR2HDMEFTPHI6_PHI2PSI3.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>

namespace BSMPT
{
namespace Models
{

double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Op6_111111;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Op6_111122;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Op6_122111;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Op6_121211;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Op6_222222;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Op6_112222;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Op6_122122;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Op6_121222;

double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OL_1b12b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OL_1b21b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OL_2b11b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OL_2b22b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OQu_1b12b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OQu_1b21b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OQu_2b11b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OQu_2b22b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OQd_1b12b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OQd_1b21b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OQd_2b11b;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::OQd_2b22b;

double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::L1tmp;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::L2tmp;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::L4tmp;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::L5tmp;
double BSMPT::Models::Class_Potential_R2HDMEFTPHI6_PHI2PSI3::m12Sqtmp;

Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Class_Potential_R2HDMEFTPHI6_PHI2PSI3()
{
  // TODO Auto-generated constructor stub
  Model         = ModelID::ModelIDs::R2HDMEFTPHI6_PHI2PSI3;
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

  // Set UseTwoLoopThermalMass to include two-loop thermal mass corrections
  // propto T^4
  SetUseTwoLoopThermalMass(true);
}

Class_Potential_R2HDMEFTPHI6_PHI2PSI3::~Class_Potential_R2HDMEFTPHI6_PHI2PSI3()
{
  // TODO Auto-generated destructor stub
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given inputfile
 */
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6_PHI2PSI3::addLegendCT() const
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
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6_PHI2PSI3::addLegendTemp() const
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
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6_PHI2PSI3::addLegendVEV() const
{
  std::vector<std::string> labels;
  labels.push_back("omega_CB");
  labels.push_back("omega_1");
  labels.push_back("omega_2");
  labels.push_back("omega_CP");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * higgs couplings. Use this to complement the legend of the given inputfile
 *
 */
std::vector<std::string>
Class_Potential_R2HDMEFTPHI6_PHI2PSI3::addLegendTripleCouplings() const
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

void Class_Potential_R2HDMEFTPHI6_PHI2PSI3::ReadAndSet(
    const std::string &linestr,
    std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  // old format
  //   for (int k = 1; k <= 8; k++)
  //   {
  //     ss >> tmp;
  //     if (k == 1)
  //       Type = tmp;
  //     else if (k == 2)
  //       L1 = tmp;
  //     else if (k == 3)
  //       L2 = tmp;
  //     else if (k == 4)
  //       L3 = tmp;
  //     else if (k == 5)
  //       L4 = tmp;
  //     else if (k == 6)
  //       L5 = tmp;
  //     else if (k == 7)
  //       m12Sq = tmp;
  //     else if (k == 8)
  //       TanBeta = tmp;
  //   }

  // new format
  for (int k = 1; k <= 16; k++)
  {
    ss >> tmp;
    if (k == 16)
      Type = tmp;
    else if (k == 9)
      L1 = tmp;
    else if (k == 10)
      L2 = tmp;
    else if (k == 11)
      L3 = tmp;
    else if (k == 12)
      L4 = tmp;
    else if (k == 13)
      L5 = tmp;
    else if (k == 7)
      m12Sq = tmp;
    else if (k == 6)
      TanBeta = tmp;
  }

  //	double sa = std::sin(alpha);
  //	double ca = std::cos(alpha);
  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = std::sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = std::sqrt(C_SinBetaSquared);
  //	L1 =1.0 / (C_vev0 * C_vev0 * C_CosBeta * C_CosBeta)* (ca * ca * MH *
  // MH +
  // sa * sa * Mh * Mh- m12Sq * C_SinBeta / C_CosBeta); 	L2 =1.0 / (C_vev0
  // * C_vev0 * C_SinBeta * C_SinBeta)* (sa * sa * MH * MH + ca * ca * Mh *
  // Mh- m12Sq * C_CosBeta / C_SinBeta); 	L3 = 2 * MHP * MHP / (C_vev0 *
  // C_vev0)+ sa * ca * (MH * MH - Mh * Mh)/ (C_vev0 * C_vev0
  // *C_CosBeta*C_SinBeta )- m12Sq / (C_vev0 * C_vev0 * C_SinBeta *
  // C_CosBeta); 	L4 = (MA * MA - 2 * MHP * MHP) / (C_vev0 * C_vev0)+ m12Sq /
  // (C_vev0 * C_vev0 * C_SinBeta * C_CosBeta); 	L5 = m12Sq / (C_vev0 *
  // C_vev0 * C_SinBeta * C_CosBeta) - MA
  //* MA / (C_vev0 * C_vev0);

  par[6] = TanBeta;
  par[4] = L5;
  par[0] = L1;
  par[1] = L2;
  par[2] = L3;
  par[3] = L4;
  par[5] = m12Sq;
  par[7] = Type;

  set_gen(par);
  return;
}

/**
 * Set Class Object with an CP-Conserving Point
 */
void Class_Potential_R2HDMEFTPHI6_PHI2PSI3::set_gen(
    const std::vector<double> &par)
{

  // double *p = (double *)par;
  scale = C_vev0;
  //	scale=C_MassZ;
  double L1store    = par[0];
  double L2store    = par[1];
  L3                = par[2];
  double L4store    = par[3];
  double L5store    = par[4];
  double m12Sqstore = par[5];
  TanBeta           = par[6];
  beta              = std::atan(TanBeta);
  Type              = static_cast<int>(par[7]);
  C_CosBetaSquared  = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta         = std::sqrt(C_CosBetaSquared);
  C_SinBetaSquared  = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta         = std::sqrt(C_SinBetaSquared);

  // corrected Lambdas to absorb CP-even mass shifts due to EFT

  double v1   = C_vev0 * C_CosBeta;
  double v1Sq = C_vev0 * C_vev0 * C_CosBetaSquared;
  double v2   = C_vev0 * C_SinBeta;
  double v2Sq = C_vev0 * C_vev0 * C_SinBetaSquared;

  L1 = L1store +
       2 * (double)std::pow((double)LambdaEFT, (double)(-2)) *
           (6 * Op6_111111 * v1Sq * v1Sq + 2 * Op6_121211 * v1Sq * v2Sq +
            Op6_122111 * v1Sq * v2Sq -
            (2 * Op6_112222 + 2 * Op6_121222 + Op6_122122) * v2Sq * v2Sq) /
           (4 * v1Sq);
  L2 = L2store +
       2 * (double)std::pow((double)LambdaEFT, (double)(-2)) *
           (-2 * Op6_111122 * v1Sq * v1Sq - 2 * Op6_121211 * v1Sq * v1Sq -
            Op6_122111 * v1Sq * v1Sq + 2 * Op6_121222 * v1Sq * v2Sq +
            Op6_122122 * v1Sq * v2Sq + 6 * Op6_222222 * v2Sq * v2Sq) /
           (4 * v2Sq);
  L4 =
      L4store + (double)std::pow((double)LambdaEFT, (double)(-2)) *
                    (Op6_111122 * v1Sq + Op6_121211 * v1Sq + Op6_122111 * v1Sq +
                     Op6_112222 * v2Sq + Op6_121222 * v2Sq + Op6_122122 * v2Sq);
  L5 = L5store +
       (double)std::pow((double)LambdaEFT, (double)(-2)) * 0.5 *
           (2 * Op6_111122 * v1Sq + 4 * Op6_121211 * v1Sq + Op6_122111 * v1Sq +
            2 * Op6_112222 * v2Sq + 4 * Op6_121222 * v2Sq + Op6_122122 * v2Sq);
  m12Sq =
      m12Sqstore +
      (double)std::pow((double)LambdaEFT, (double)(-2)) * 0.5 * v1 * v2 *
          (2 * Op6_111122 * v1Sq + 2 * Op6_121211 * v1Sq + Op6_122111 * v1Sq +
           2 * Op6_112222 * v2Sq + 2 * Op6_121222 * v2Sq + Op6_122122 * v2Sq);

  L1tmp    = L1;
  L2tmp    = L2;
  L4tmp    = L4;
  L5tmp    = L5;
  m12Sqtmp = m12Sq;

  m11Sq =
      m12Sq * TanBeta -
      C_vev0 * C_vev0 * C_SinBetaSquared * (L4 + L5 + L3) / 0.2e1 -
      C_vev0 * C_vev0 * C_CosBetaSquared * L1 / 0.2e1 +
      (3 * C_CosBetaSquared * C_CosBetaSquared * C_vev0 * C_vev0 * C_vev0 *
       C_vev0 * Op6_111111 *
       (double)std::pow((double)LambdaEFT, (double)(-2))) /
          4. +
      (C_CosBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_111122 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          2. +
      (C_SinBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_112222 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          4. +
      C_CosBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
          Op6_121211 * (double)std::pow((double)LambdaEFT, (double)(-2)) +
      (C_SinBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_121222 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          2. +
      (C_CosBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_122111 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          2. +
      (C_SinBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_122122 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          4;
  m22Sq =
      m12Sq * 1.0 / TanBeta -
      C_vev0 * C_vev0 * C_CosBetaSquared * (L4 + L5 + L3) / 0.2e1 -
      C_vev0 * C_vev0 * C_SinBetaSquared * L2 / 0.2e1 +
      (C_CosBetaSquared * C_CosBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_111122 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          4. +
      (C_CosBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_112222 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          2. +
      (C_CosBetaSquared * C_CosBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_121211 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          2. +
      C_CosBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
          Op6_121222 * (double)std::pow((double)LambdaEFT, (double)(-2)) +
      (C_CosBetaSquared * C_CosBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_122111 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          4. +
      (C_CosBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 * C_vev0 *
       Op6_122122 * (double)std::pow((double)LambdaEFT, (double)(-2))) /
          2. +
      (3 * C_SinBetaSquared * C_SinBetaSquared * C_vev0 * C_vev0 * C_vev0 *
       C_vev0 * Op6_222222 *
       (double)std::pow((double)LambdaEFT, (double)(-2))) /
          4.;

  //	double ML5 = 2*m12Sq/(C_vev0*C_vev0*C_SinBeta*C_CosBeta);
  //	double TripleHiggs =
  //-3.0/(C_vev0*std::sin(2*beta))*(Mh*Mh*(2*std::cos(alpha+beta)+std::sin(2*alpha)*std::sin(beta-alpha))
  //- std::cos(alpha+beta)*
  //(double)
  // std::pow((double)std::cos(beta-alpha),(double)(2))*C_vev0*C_vev0*ML5 );

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

void Class_Potential_R2HDMEFTPHI6_PHI2PSI3::set_CT_Pot_Par(
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
  dT1         = p[8];
  dT2         = p[9];
  dTCP        = p[10];
  dTCB        = p[11];
  dOp6_111111 = p[12];

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
      (double)(-6 * dOp6_111122 - 6 * dOp6_122111 - 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][0][0][3][3] =
      (double)(-6 * dOp6_111122 - 6 * dOp6_122111 + 12 * dOp6_121211) *
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
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][0][2][5][7] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][0][3][4][7] =
      (double)(-3 * dOp6_122111 + 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][0][3][5][6] =
      (double)(3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][1][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][1][1][2][2] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][1][1][3][3] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][1][1][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][1][1][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][1][1][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][1][1][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[0][0][1][2][4][7] =
      (double)(dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][2][5][6] =
      (double)(-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][3][4][6] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][1][3][5][7] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][0][2][2][2][2] =
      (double)(-12 * dOp6_121222 - 6 * dOp6_122122 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][3][3] =
      (double)(-2 * dOp6_122122 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][4][4] =
      (double)(-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][5][5] =
      (double)(-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][6][6] =
      (double)(-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][2][2][7][7] =
      (double)(-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][3][3] =
      (double)(-6 * dOp6_122122 - 6 * dOp6_112222 + 12 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][4][4] =
      (double)(-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][5][5] =
      (double)(-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][6][6] =
      (double)(-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][3][3][7][7] =
      (double)(-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][4][4][6][6] =
      (double)(-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][4][4][7][7] =
      (double)(-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[0][0][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[0][0][5][5][6][6] =
      (double)(-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[0][0][5][5][7][7] =
      (double)(-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
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
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][1][1][2][5][7] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][1][1][3][4][7] =
      (double)(-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][1][1][3][5][6] =
      (double)(dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
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
      (double)(-3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][2][2][5][7] =
      (double)(-3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][2][3][4][7] =
      (double)(2 * dOp6_121222 - dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][2][3][5][6] =
      (double)(-2 * dOp6_121222 + dOp6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][3][3][4][6] =
      (double)(-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][3][3][5][7] =
      (double)(-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][4][4][6] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][4][5][7] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][5][5][6] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][6][6][6] =
      (double)(-3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][4][6][7][7] =
      (double)(-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][5][5][5][7] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][2][5][6][6][7] =
      (double)(-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][2][5][7][7][7] =
      (double)(-3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
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
  Curvature_Higgs_CT_L6[0][3][4][5][5][7] =
      (double)(-dOp6_122111 + 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][6][6][7] =
      (double)(-dOp6_122122 + 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[0][3][4][7][7][7] =
      (double)(-3 * dOp6_122122 + 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][5][5][5][6] =
      (double)(3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][5][6][6][6] =
      (double)(3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[0][3][5][6][7][7] =
      (double)(dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][1][1][1] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][1][1][2][2] =
      (double)(-6 * dOp6_111122 - 6 * dOp6_122111 + 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][1][1][3][3] =
      (double)(-6 * dOp6_111122 - 6 * dOp6_122111 - 12 * dOp6_121211) *
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
      (double)(3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][2][5][6] =
      (double)(-3 * dOp6_122111 + 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][3][4][6] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][1][3][5][7] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][1][2][2][2][2] =
      (double)(-6 * dOp6_122122 + 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][3][3] =
      (double)(-2 * dOp6_122122 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][4][4] =
      (double)(-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][5][5] =
      (double)(-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][6][6] =
      (double)(-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][2][2][7][7] =
      (double)(-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][3][3] =
      (double)(-6 * dOp6_122122 - 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][4][4] =
      (double)(-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][5][5] =
      (double)(-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][6][6] =
      (double)(-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][3][3][7][7] =
      (double)(-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][4][4][6][6] =
      (double)(-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][4][4][7][7] =
      (double)(-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[1][1][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[1][1][5][5][6][6] =
      (double)(-dOp6_122111 + 2 * dOp6_121211 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][5][5][7][7] =
      (double)(-2 * dOp6_121211 - dOp6_122111 - 2 * dOp6_111122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[1][1][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[1][1][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[1][1][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_112222;
  Curvature_Higgs_CT_L6[1][2][2][2][4][7] =
      (double)(3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][2][5][6] =
      (double)(-3 * dOp6_122122 + 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][3][4][6] =
      (double)(-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][2][3][5][7] =
      (double)(-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
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
  Curvature_Higgs_CT_L6[1][2][4][5][5][7] =
      (double)(dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][6][6][7] =
      (double)(dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][2][4][7][7][7] =
      (double)(3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][5][5][5][6] =
      (double)(-3 * dOp6_122111 + 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][5][6][6][6] =
      (double)(-3 * dOp6_122122 + 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][2][5][6][7][7] =
      (double)(-dOp6_122122 + 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][3][3][4][6] =
      (double)(-3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][3][3][5][7] =
      (double)(-3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][4][4][6] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][4][5][7] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][5][5][6] =
      (double)(-dOp6_122111 - 2 * dOp6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][6][6][6] =
      (double)(-3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][4][6][7][7] =
      (double)(-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][5][5][5][7] =
      (double)(-3 * dOp6_122111 - 6 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_CT_L6[1][3][5][6][6][7] =
      (double)(-dOp6_122122 - 2 * dOp6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_CT_L6[1][3][5][7][7][7] =
      (double)(-3 * dOp6_122122 - 6 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1) / 0.2e1;
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
      (double)(-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][4][4][7][7] =
      (double)(-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[2][2][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[2][2][5][5][6][6] =
      (double)(-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[2][2][5][5][7][7] =
      (double)(-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
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
      (double)(-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][4][4][7][7] =
      (double)(-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121222;
  Curvature_Higgs_CT_L6[3][3][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_111122;
  Curvature_Higgs_CT_L6[3][3][5][5][6][6] =
      (double)(-dOp6_122122 + 2 * dOp6_121222 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[3][3][5][5][7][7] =
      (double)(-dOp6_122122 - 2 * dOp6_112222 - 2 * dOp6_121222) *
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
      (double)(-6 * dOp6_111122 - 6 * dOp6_122111 - 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][4][4][7][7] =
      (double)(-6 * dOp6_111122 - 6 * dOp6_122111 + 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][4][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * dOp6_121211;
  Curvature_Higgs_CT_L6[4][4][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * dOp6_111111;
  Curvature_Higgs_CT_L6[4][4][5][5][6][6] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][5][5][7][7] =
      (double)(-2 * dOp6_111122 - 2 * dOp6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][6][6][6][6] =
      (double)(-6 * dOp6_122122 - 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][6][6][7][7] =
      (double)(-2 * dOp6_122122 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[4][4][7][7][7][7] =
      (double)(-6 * dOp6_122122 + 12 * dOp6_121222 - 6 * dOp6_112222) *
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
      (double)(-6 * dOp6_111122 - 6 * dOp6_122111 + 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][5][5][7][7] =
      (double)(-6 * dOp6_111122 - 6 * dOp6_122111 - 12 * dOp6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][6][6][6][6] =
      (double)(-6 * dOp6_122122 + 12 * dOp6_121222 - 6 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][6][6][7][7] =
      (double)(-2 * dOp6_122122 - 2 * dOp6_112222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_CT_L6[5][5][7][7][7][7] =
      (double)(-6 * dOp6_122122 - 12 * dOp6_121222 - 6 * dOp6_112222) *
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
void Class_Potential_R2HDMEFTPHI6_PHI2PSI3::write() const
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

  ss << "beta = " << beta << std::endl;
  ss << "tan(beta) = " << TanBeta << std::endl;
  ss << "L1 = " << L1 << std::endl;
  ss << "L2 = " << L2 << std::endl;
  ss << "L3 = " << L3 << std::endl;
  ss << "L4 = " << L4 << std::endl;
  ss << "Re(L5) = " << L5 << std::endl;
  ss << "Re(m_12^2) = " << m12Sq << std::endl;
  ss << "m_{11}^2 = " << m11Sq << std::endl;
  ss << "m_{22}^2 = " << m22Sq << std::endl;

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
  ss << "dTCB:= " << dTCB << ";" << std::endl;

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
  ss << std::abs(HiggsRot(posN[1], 6)) << " zeta_2" << std::endl;

  if (UseTwoLoopThermalMass)
  {
    ss << "Dim-6 two-loop corrections to thermal masses are taken into "
          "account.\n";
  }
  else
  {
    ss << "Note that NO Dim-6 two-loop corrections to thermal masses are taken "
          "into "
          "account!\n";
  }
  if (UseTensorSymFac)
  {
    ss << "Usage of combined c-factor * tensor structure!\n";
  }

  ss << "The dim-6 operator are set to \nOp6_111111 = " << Op6_111111
     << std::endl;
  ss << "Op6_111122 = " << Op6_111122 << "\n";
  ss << "Op6_122111 = " << Op6_122111 << "\n";
  ss << "Op6_121211 = " << Op6_121211 << "\n";
  ss << "Op6_222222 = " << Op6_222222 << "\n";
  ss << "Op6_112222 = " << Op6_112222 << "\n";
  ss << "Op6_122122 = " << Op6_122122 << "\n";
  ss << "Op6_121222 = " << Op6_121222 << "\n";

  ss << "OL_1b12b = " << OL_1b12b << "\n";
  ss << "OL_1b21b = " << OL_1b21b << "\n";
  ss << "OL_2b11b = " << OL_2b11b << "\n";
  ss << "OL_2b22b = " << OL_2b22b << "\n";

  ss << "OL_1b12b = " << OQu_1b12b << "\n";
  ss << "OL_1b12b = " << OQu_1b21b << "\n";
  ss << "OL_1b12b = " << OQu_2b11b << "\n";
  ss << "OL_1b12b = " << OQu_2b22b << "\n";

  ss << "OL_1b12b = " << OQd_1b12b << "\n";
  ss << "OL_1b12b = " << OQd_1b21b << "\n";
  ss << "OL_1b12b = " << OQd_2b11b << "\n";
  ss << "OL_1b12b = " << OQd_2b22b << std::endl;

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms in the 2HDM
 */
std::vector<double> Class_Potential_R2HDMEFTPHI6_PHI2PSI3::calc_CT() const
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

  double t = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0, s6 = 0, s7 = 0,
         s8 = 0; // Values of dL4 and dOp6_ijklmn

  // dm11Sq
  parCT.push_back(
      -(double)((-2 * t * v1 * v2 * v2 + 5 * HesseWeinberg(0, 0) * v1 +
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
      (double)((2 * t * v1 * v1 * v2 * v2 + HesseWeinberg(6, 6) * v2 * v2 -
                2 * HesseWeinberg(0, 0) * v1 * v1 -
                HesseWeinberg(1, 3) * v1 * v2 -
                3 * HesseWeinberg(3, 3) * v2 * v2 +
                HesseWeinberg(4, 6) * v1 * v2 +
                2 * v1 * v1 * HesseWeinberg(5, 5)) *
               (double)std::pow((double)v2, (double)(-2))) /
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
  parCT.push_back((double)(-t * v2 * v2 + 2 * HesseWeinberg(0, 0) -
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
          (double)std::pow((double)v2, (double)(-2)) +
      -std::pow(LambdaEFT, -2) * s3 * v1 * v1 / 2 +
      std::pow(LambdaEFT, -2) * s4 * v1 * v1 -
      std::pow(LambdaEFT, -2) * v2 * v2 * s7 / 2 +
      std::pow(LambdaEFT, -2) * v2 * v2 * s8);

  // dT1
  double tmp =
      HesseWeinberg(1, 3) * v2 + HesseWeinberg(0, 0) * v1 - NablaWeinberg(4);
  if (std::abs(tmp) < 1e-9) tmp = 0;
  parCT.push_back(tmp);
  // dT2
  tmp = HesseWeinberg(1, 3) * v1 + HesseWeinberg(3, 3) * v2 - NablaWeinberg(6);
  if (std::abs(tmp) < 1e-9) tmp = 0;
  parCT.push_back(tmp);
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

/**
 * Calculates the corrections to the Triple higgs couplings in the mass basis.
 *
 * Use the vector TripleHiggsCorrectionsCWPhysical to save your couplings and
 * set the nTripleCouplings to the number of couplings you want as output.
 */
void Class_Potential_R2HDMEFTPHI6_PHI2PSI3::TripleHiggsCouplings()
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

void Class_Potential_R2HDMEFTPHI6_PHI2PSI3::SetCurvatureArrays()
{
  initVectors();

  HiggsVev[4] = C_vev0 * C_CosBeta;
  HiggsVev[6] = C_vev0 * C_SinBeta;

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
      (double)(-6 * Op6_111122 - 6 * Op6_122111 - 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][0][0][3][3] =
      (double)(-6 * Op6_111122 - 6 * Op6_122111 + 12 * Op6_121211) *
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
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][0][2][5][7] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][0][3][4][7] =
      (double)(-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][0][3][5][6] =
      (double)(3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][1][1][1][1] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][1][1][2][2] =
      (double)(-2 * Op6_111122 - 2 * Op6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][1][1][3][3] =
      (double)(-2 * Op6_111122 - 2 * Op6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][1][1][4][4] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][1][1][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][1][1][6][6] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][1][1][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[0][0][1][2][4][7] = (double)(Op6_122111 - 2 * Op6_121211) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][0][1][2][5][6] =
      (double)(-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][1][3][4][6] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][1][3][5][7] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][0][2][2][2][2] =
      (double)(-6 * Op6_112222 - 6 * Op6_122122 - 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][3][3] =
      (double)(-2 * Op6_112222 - 2 * Op6_122122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][4][4] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][5][5] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][6][6] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][2][2][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][3][3] =
      (double)(-6 * Op6_122122 - 6 * Op6_112222 + 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][4][4] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][5][5] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][6][6] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][3][3][7][7] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][4][4][6][6] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][4][4][7][7] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[0][0][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[0][0][5][5][6][6] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[0][0][5][5][7][7] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
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
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][1][1][2][5][7] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][1][1][3][4][7] =
      (double)(-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][1][1][3][5][6] = (double)(Op6_122111 - 2 * Op6_121211) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
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
      (double)(-6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][2][2][5][7] =
      (double)(-6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][2][3][4][7] = (double)(2 * Op6_121222 - Op6_122122) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][2][2][3][5][6] =
      (double)(-2 * Op6_121222 + Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][3][3][4][6] =
      (double)(-2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][3][3][5][7] =
      (double)(-2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][4][4][4][6] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][4][4][5][7] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][4][5][5][6] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][4][6][6][6] =
      (double)(-6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][4][6][7][7] =
      (double)(-2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][5][5][5][7] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][5][6][6][7] =
      (double)(-2 * Op6_121222 - Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][2][5][7][7][7] =
      (double)(-6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][3][3][4][7] =
      (double)(6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][3][3][5][6] =
      (double)(-6 * Op6_121222 + 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][4][4][4][7] =
      (double)(-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][4][4][5][6] = (double)(Op6_122111 - 2 * Op6_121211) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][4][5][5][7] =
      (double)(-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][4][6][6][7] = (double)(2 * Op6_121222 - Op6_122122) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[0][3][4][7][7][7] =
      (double)(6 * Op6_121222 - 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][5][5][5][6] =
      (double)(3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][5][6][6][6] =
      (double)(-6 * Op6_121222 + 3 * Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[0][3][5][6][7][7] =
      (double)(-2 * Op6_121222 + Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][1][1][1][1][1] =
      -0.90e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][1][1][2][2] =
      (double)(-6 * Op6_111122 - 6 * Op6_122111 + 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][1][1][3][3] =
      (double)(-6 * Op6_111122 - 6 * Op6_122111 - 12 * Op6_121211) *
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
      (double)(3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][1][1][2][5][6] =
      (double)(-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][1][1][3][4][6] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][1][1][3][5][7] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][1][2][2][2][2] =
      (double)(-6 * Op6_122122 - 6 * Op6_112222 + 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][3][3] =
      (double)(-2 * Op6_122122 - 2 * Op6_112222) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][4][4] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][5][5] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][6][6] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][2][2][7][7] =
      (double)(-2 * Op6_112222 + 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][3][3] =
      (double)(-6 * Op6_112222 - 6 * Op6_122122 - 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][4][4] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][5][5] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][6][6] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][3][3][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][4][4][4][4] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][4][4][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][4][4][6][6] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][4][4][7][7] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[1][1][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[1][1][5][5][6][6] =
      (double)(-2 * Op6_111122 + 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][5][5][7][7] =
      (double)(-2 * Op6_111122 - 2 * Op6_121211 - Op6_122111) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[1][1][6][6][6][6] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[1][1][6][6][7][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[1][1][7][7][7][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_112222;
  Curvature_Higgs_L6[1][2][2][2][4][7] =
      (double)(3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][2][2][5][6] =
      (double)(-3 * Op6_122122 + 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][2][3][4][6] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][2][3][5][7] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][3][3][4][7] =
      (double)(-2 * Op6_121222 + Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][3][3][5][6] =
      (double)(-Op6_122122 + 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][4][4][4][7] =
      (double)(3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][4][4][5][6] =
      (double)(-Op6_122111 + 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][4][5][5][7] = (double)(Op6_122111 - 2 * Op6_121211) *
                                         std::pow(LambdaEFT, -0.2e1) / 0.2e1;
  Curvature_Higgs_L6[1][2][4][6][6][7] =
      (double)(-2 * Op6_121222 + Op6_122122) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][4][7][7][7] =
      (double)(3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][5][5][5][6] =
      (double)(-3 * Op6_122111 + 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][5][6][6][6] =
      (double)(-3 * Op6_122122 + 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][2][5][6][7][7] =
      (double)(-Op6_122122 + 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][3][3][4][6] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][3][3][5][7] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][4][4][4][6] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][4][4][5][7] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][4][5][5][6] =
      (double)(-Op6_122111 - 2 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][4][6][6][6] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][4][6][7][7] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][5][5][5][7] =
      (double)(-3 * Op6_122111 - 6 * Op6_121211) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][5][6][6][7] =
      (double)(-Op6_122122 - 2 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
  Curvature_Higgs_L6[1][3][5][7][7][7] =
      (double)(-3 * Op6_122122 - 6 * Op6_121222) * std::pow(LambdaEFT, -0.2e1) /
      0.2e1;
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
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][4][4][7][7] =
      (double)(-2 * Op6_112222 - Op6_122122 + 2 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[2][2][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[2][2][5][5][6][6] =
      (double)(-2 * Op6_112222 - Op6_122122 + 2 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[2][2][5][5][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
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
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][4][4][7][7] =
      (double)(-2 * Op6_112222 - Op6_122122 + 2 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][4][5][6][7] =
      -0.2e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121222;
  Curvature_Higgs_L6[3][3][5][5][5][5] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_111122;
  Curvature_Higgs_L6[3][3][5][5][6][6] =
      (double)(-2 * Op6_112222 - Op6_122122 + 2 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[3][3][5][5][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_121222 - Op6_122122) *
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
      (double)(-6 * Op6_111122 - 6 * Op6_122111 - 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][4][4][7][7] =
      (double)(-6 * Op6_111122 - 6 * Op6_122111 + 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][4][5][6][7] =
      -0.6e1 * std::pow(LambdaEFT, -0.2e1) * Op6_121211;
  Curvature_Higgs_L6[4][4][5][5][5][5] =
      -0.18e2 * std::pow(LambdaEFT, -0.2e1) * Op6_111111;
  Curvature_Higgs_L6[4][4][5][5][6][6] =
      (double)(-2 * Op6_111122 - 2 * Op6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][5][5][7][7] =
      (double)(-2 * Op6_111122 - 2 * Op6_122111) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][6][6][6][6] =
      (double)(-6 * Op6_112222 - 6 * Op6_122122 - 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][6][6][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_122122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[4][4][7][7][7][7] =
      (double)(-6 * Op6_122122 - 6 * Op6_112222 + 12 * Op6_121222) *
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
      (double)(-6 * Op6_111122 - 6 * Op6_122111 + 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][5][5][7][7] =
      (double)(-6 * Op6_111122 - 6 * Op6_122111 - 12 * Op6_121211) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][6][6][6][6] =
      (double)(-6 * Op6_122122 - 6 * Op6_112222 + 12 * Op6_121222) *
      std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][6][6][7][7] =
      (double)(-2 * Op6_112222 - 2 * Op6_122122) * std::pow(LambdaEFT, -0.2e1);
  Curvature_Higgs_L6[5][5][7][7][7][7] =
      (double)(-6 * Op6_112222 - 6 * Op6_122122 - 12 * Op6_121222) *
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

  Curvature_Gauge_G2H2[0][0][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][4][4] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][5][5] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][6][6] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][7][7] = C_g * C_g / 0.2e1;

  Curvature_Gauge_G2H2[0][3][0][4] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][5] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][6] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][7] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][4][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][5][1] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][6][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][7][3] = C_g * C_gs / 0.2e1;

  Curvature_Gauge_G2H2[1][1][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][4][4] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][5][5] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][6][6] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][7][7] = C_g * C_g / 0.2e1;

  Curvature_Gauge_G2H2[1][3][0][5] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][4] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][7] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][6] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][4][1] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][5][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][6][3] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][7][2] = C_g * C_gs / 0.2e1;

  Curvature_Gauge_G2H2[2][2][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][4][4] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][5][5] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][6][6] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][7][7] = C_g * C_g / 0.2e1;

  Curvature_Gauge_G2H2[2][3][0][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][4][4] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][5][5] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][6][6] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][7][7] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][4] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][5] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][6] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][7] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][4][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][5][1] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][6][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][7][3] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][5] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][4] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][7] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][6] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][4][1] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][5][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][6][3] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][7][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][4][4] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][5][5] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][6][6] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][7][7] = -C_g * C_gs / 0.2e1;
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

  Curvature_Lepton_F2H3[0][1][0][0][6] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][0][0][7] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][0][2][4] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[0][1][0][2][5] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[0][1][0][3][4] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[0][1][0][3][5] =
      -std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[0][1][1][1][6] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][1][1][7] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][1][2][4] =
      -std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[0][1][1][2][5] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[0][1][1][3][4] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[0][1][1][3][5] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[0][1][2][2][6] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][2][2][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][3][3][6] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][3][3][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][4][4][6] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[0][1][4][4][7] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[0][1][4][5][6] =
      II * std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][4][5][7] = std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][5][5][6] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[0][1][5][5][7] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[0][1][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[0][1][6][6][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][6][7][7] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[0][1][7][7][7] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[1][6][0][0][2] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[1][6][0][0][3] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[1][6][0][1][2] =
      II * std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][0][1][3] = std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][0][4][6] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[1][6][0][4][7] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[1][6][0][5][6] =
      -std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[1][6][0][5][7] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[1][6][1][1][2] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[1][6][1][1][3] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[1][6][1][4][6] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[1][6][1][4][7] =
      -std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[1][6][1][5][6] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[1][6][1][5][7] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[1][6][2][2][2] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[1][6][2][2][3] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][2][3][3] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][2][4][4] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][2][5][5] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][2][6][6] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][2][7][7] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][3][3][3] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[1][6][3][4][4] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][3][5][5] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][3][6][6] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[1][6][3][7][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][0][0][6] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][0][0][7] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][0][2][4] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[2][3][0][2][5] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[2][3][0][3][4] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[2][3][0][3][5] =
      -std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[2][3][1][1][6] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][1][1][7] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][1][2][4] =
      -std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[2][3][1][2][5] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[2][3][1][3][4] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[2][3][1][3][5] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[2][3][2][2][6] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][2][2][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][3][3][6] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][3][3][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][4][4][6] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[2][3][4][4][7] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[2][3][4][5][6] =
      II * std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][4][5][7] = std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][5][5][6] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[2][3][5][5][7] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[2][3][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[2][3][6][6][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][6][7][7] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[2][3][7][7][7] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[3][7][0][0][2] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[3][7][0][0][3] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[3][7][0][1][2] =
      II * std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][0][1][3] = std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][0][4][6] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[3][7][0][4][7] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[3][7][0][5][6] =
      -std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[3][7][0][5][7] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[3][7][1][1][2] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[3][7][1][1][3] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[3][7][1][4][6] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[3][7][1][4][7] =
      -std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[3][7][1][5][6] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[3][7][1][5][7] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[3][7][2][2][2] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[3][7][2][2][3] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][2][3][3] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][2][4][4] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][2][5][5] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][2][6][6] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][2][7][7] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][3][3][3] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[3][7][3][4][4] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][3][5][5] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][3][6][6] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[3][7][3][7][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][0][0][6] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][0][0][7] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][0][2][4] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][2][5] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][3][4] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][0][3][5] =
      -std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][1][6] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][1][1][7] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][1][2][4] =
      -std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][2][5] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][3][4] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][1][3][5] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[4][5][2][2][6] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][2][2][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][3][3][6] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][3][3][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][4][6] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][4][7] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][5][6] =
      II * std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][4][5][7] = std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][5][5][6] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][5][5][7] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[4][5][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[4][5][6][6][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][6][7][7] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[4][5][7][7][7] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[5][8][0][0][2] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][0][3] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][1][2] =
      II * std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][1][3] = std::sqrt(0.2e1) * OL_1b21b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][0][4][6] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][0][4][7] =
      std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][0][5][6] =
      -std::sqrt(0.2e1) * II * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][0][5][7] =
      std::sqrt(0.2e1) * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][1][1][2] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][1][3] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b + OL_2b11b) / 0.2e1;
  Curvature_Lepton_F2H3[5][8][1][4][6] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][1][4][7] =
      -std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][1][5][6] =
      std::sqrt(0.2e1) * (OL_1b12b - OL_1b21b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][1][5][7] =
      std::sqrt(0.2e1) * II * (OL_1b21b + OL_1b12b) / 0.4e1;
  Curvature_Lepton_F2H3[5][8][2][2][2] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[5][8][2][2][3] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][3][3] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][4][4] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][5][5] = std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][6][6] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][2][7][7] = std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][3][3][3] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OL_2b22b;
  Curvature_Lepton_F2H3[5][8][3][4][4] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][3][5][5] =
      II * std::sqrt(0.2e1) * OL_2b11b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][3][6][6] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;
  Curvature_Lepton_F2H3[5][8][3][7][7] =
      II * std::sqrt(0.2e1) * OL_2b22b / 0.2e1;

  sym5Dim(Curvature_Lepton_F2H3, NLepton, NLepton, NHiggs, NHiggs, NHiggs);

  Curvature_Quark_F2H3[0][6][0][0][6] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][0][0][7] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][0][2][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][0][2][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][0][3][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][0][3][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][0][4][2] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][0][4][3] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][0][5][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][0][5][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][0][6][0] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][0][7][0] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][1][1][6] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][1][1][7] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][1][2][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][1][2][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][1][3][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][1][3][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][1][4][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][1][4][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][1][5][2] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][1][5][3] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][1][6][1] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][1][7][1] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][2][0][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][2][0][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][2][1][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][2][1][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][2][2][6] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][2][2][7] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][2][4][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][2][4][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][2][5][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][2][5][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][2][6][2] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][2][7][2] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][3][0][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][3][0][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][3][1][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][3][1][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][3][3][6] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][3][3][7] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][3][4][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][3][4][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][3][5][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][3][5][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][3][6][3] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][3][7][3] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][4][0][2] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][4][0][3] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][4][1][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][4][1][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][4][2][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][4][2][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][4][3][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][4][3][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][4][4][6] =
      std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][4][4][7] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][4][5][6] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][4][5][7] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][4][6][4] =
      std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][4][6][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][4][7][4] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][4][7][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][5][0][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][5][0][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][5][1][2] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][5][1][3] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][5][2][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][5][2][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][5][3][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][5][3][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][6][5][4][6] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][5][4][7] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][5][5][6] =
      -std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][5][5][7] =
      std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][5][6][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][5][6][5] =
      -std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][5][7][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][5][7][5] =
      std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][0][0] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][1][1] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][2][2] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][3][3] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][4][4] =
      std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][4][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][5][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][5][5] =
      -std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[0][6][6][6][7] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][7][6] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][6][7][7] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][0][0] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][1][1] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][2][2] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][3][3] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][4][4] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][4][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][5][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][5][5] =
      std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][6][6] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][6][7] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][7][6] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][6][7][7][7] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[0][9][0][0][2] =
      -std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][0][0][3] =
      std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][0][1][2] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][0][1][3] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][0][2][0] =
      -std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][0][2][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][0][3][0] =
      std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][0][3][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][0][4][6] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][0][4][7] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][0][5][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][0][5][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][0][6][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][0][6][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][0][7][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][0][7][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][1][0][2] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][1][0][3] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][1][1][2] =
      std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][1][1][3] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][1][2][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][1][2][1] =
      std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][1][3][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][1][3][1] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][1][4][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][1][4][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][1][5][6] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][1][5][7] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][1][6][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][1][6][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][1][7][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][1][7][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][2][0][0] =
      -std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][0][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][1][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][1][1] =
      std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][2][2] =
      -0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[0][9][2][2][3] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][3][2] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][3][3] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][4][4] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][5][5] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][6][6] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][2][7][7] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][0][0] =
      std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][0][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][1][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][1][1] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][2][2] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][2][3] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][3][2] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][3][3] =
      -0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[0][9][3][4][4] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][5][5] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][6][6] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][3][7][7] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][4][0][6] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][4][0][7] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][4][1][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][4][1][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][4][2][4] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][4][3][4] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][4][4][2] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][4][4][3] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][4][6][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][4][6][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][4][7][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][4][7][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][5][0][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][5][0][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][5][1][6] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][5][1][7] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][5][2][5] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][5][3][5] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][5][5][2] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][5][5][3] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[0][9][5][6][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][5][6][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][5][7][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][5][7][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][6][0][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][6][0][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][6][1][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][6][1][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][6][2][6] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][6][3][6] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][6][4][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][6][4][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][6][5][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][6][5][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][6][6][2] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][6][6][3] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][7][0][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][7][0][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][7][1][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][7][1][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][7][2][7] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][7][3][7] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][7][4][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][7][4][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][7][5][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][7][5][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[0][9][7][7][2] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[0][9][7][7][3] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][0][0][6] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][0][0][7] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][0][2][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][0][2][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][0][3][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][0][3][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][0][4][2] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][0][4][3] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][0][5][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][0][5][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][0][6][0] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][0][7][0] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][1][1][6] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][1][1][7] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][1][2][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][1][2][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][1][3][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][1][3][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][1][4][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][1][4][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][1][5][2] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][1][5][3] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][1][6][1] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][1][7][1] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][2][0][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][2][0][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][2][1][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][2][1][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][2][2][6] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][2][2][7] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][2][4][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][2][4][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][2][5][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][2][5][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][2][6][2] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][2][7][2] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][3][0][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][3][0][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][3][1][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][3][1][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][3][3][6] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][3][3][7] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][3][4][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][3][4][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][3][5][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][3][5][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][3][6][3] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][3][7][3] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][4][0][2] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][4][0][3] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][4][1][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][4][1][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][4][2][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][4][2][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][4][3][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][4][3][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][4][4][6] =
      std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][4][4][7] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][4][5][6] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][4][5][7] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][4][6][4] =
      std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][4][6][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][4][7][4] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][4][7][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][5][0][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][5][0][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][5][1][2] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][5][1][3] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][5][2][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][5][2][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][5][3][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][5][3][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][7][5][4][6] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][5][4][7] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][5][5][6] =
      -std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][5][5][7] =
      std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][5][6][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][5][6][5] =
      -std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][5][7][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][5][7][5] =
      std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][0][0] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][1][1] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][2][2] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][3][3] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][4][4] =
      std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][4][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][5][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][5][5] =
      -std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[1][7][6][6][7] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][7][6] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][6][7][7] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][0][0] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][1][1] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][2][2] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][3][3] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][4][4] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][4][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][5][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][5][5] =
      std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][6][6] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][6][7] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][7][6] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][7][7][7][7] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[1][10][0][0][2] =
      -std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][0][0][3] =
      std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][0][1][2] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][0][1][3] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][0][2][0] =
      -std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][0][2][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][0][3][0] =
      std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][0][3][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][0][4][6] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][0][4][7] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][0][5][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][0][5][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][0][6][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][0][6][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][0][7][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][0][7][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][1][0][2] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][1][0][3] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][1][1][2] =
      std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][1][1][3] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][1][2][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][1][2][1] =
      std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][1][3][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][1][3][1] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][1][4][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][1][4][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][1][5][6] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][1][5][7] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][1][6][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][1][6][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][1][7][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][1][7][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][2][0][0] =
      -std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][0][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][1][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][1][1] =
      std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][2][2] =
      -0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[1][10][2][2][3] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][3][2] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][3][3] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][4][4] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][5][5] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][6][6] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][2][7][7] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][0][0] =
      std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][0][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][1][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][1][1] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][2][2] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][2][3] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][3][2] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][3][3] =
      -0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[1][10][3][4][4] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][5][5] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][6][6] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][3][7][7] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][4][0][6] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][4][0][7] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][4][1][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][4][1][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][4][2][4] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][4][3][4] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][4][4][2] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][4][4][3] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][4][6][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][4][6][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][4][7][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][4][7][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][5][0][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][5][0][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][5][1][6] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][5][1][7] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][5][2][5] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][5][3][5] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][5][5][2] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][5][5][3] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[1][10][5][6][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][5][6][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][5][7][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][5][7][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][6][0][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][6][0][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][6][1][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][6][1][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][6][2][6] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][6][3][6] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][6][4][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][6][4][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][6][5][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][6][5][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][6][6][2] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][6][6][3] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][7][0][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][7][0][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][7][1][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][7][1][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][7][2][7] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][7][3][7] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][7][4][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][7][4][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][7][5][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][7][5][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[1][10][7][7][2] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[1][10][7][7][3] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][0][0][6] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][0][0][7] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][0][2][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][2][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][3][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][3][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][4][2] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][4][3] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][5][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][5][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][0][6][0] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][0][7][0] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][1][1][6] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][1][1][7] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][1][2][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][2][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][3][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][3][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][4][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][4][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][5][2] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][5][3] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][1][6][1] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][1][7][1] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][2][0][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][2][0][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][2][1][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][2][1][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][2][2][6] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][2][2][7] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][2][4][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][2][4][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][2][5][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][2][5][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][2][6][2] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][2][7][2] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][3][0][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][3][0][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][3][1][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][3][1][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][3][3][6] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][3][3][7] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][3][4][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][3][4][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][3][5][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][3][5][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][3][6][3] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][3][7][3] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][0][2] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][4][0][3] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][4][1][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][4][1][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][4][2][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][4][2][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][4][3][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][4][3][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][4][4][6] =
      std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][4][7] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][5][6] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][5][7] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][6][4] =
      std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][6][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][7][4] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][4][7][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][0][2] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][5][0][3] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][5][1][2] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][5][1][3] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][5][2][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][5][2][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][5][3][0] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][5][3][1] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][8][5][4][6] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][4][7] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][5][6] =
      -std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][5][7] =
      std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][6][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][6][5] =
      -std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][7][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][5][7][5] =
      std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][0][0] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][1][1] = std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][2][2] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][3][3] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][4][4] =
      std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][4][5] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][5][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][5][5] =
      -std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[2][8][6][6][7] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][7][6] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][6][7][7] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][0][0] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][1][1] =
      II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][2][2] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][3][3] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][4][4] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][4][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][5][4] = std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][5][5] =
      std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][6][6] =
      II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][6][7] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][7][6] = std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][8][7][7][7] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[2][11][0][0][2] =
      -std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][0][3] =
      std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][1][2] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][1][3] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][2][0] =
      -std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][2][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][3][0] =
      std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][3][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][0][4][6] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][4][7] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][5][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][5][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][6][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][6][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][7][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][0][7][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][0][2] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][0][3] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][1][2] =
      std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][1][3] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][2][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][2][1] =
      std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][3][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][3][1] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][1][4][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][4][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][5][6] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][5][7] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][6][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][6][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][7][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][1][7][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][2][0][0] =
      -std::sqrt(0.2e1) * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][0][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][1][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][1][1] =
      std::sqrt(0.2e1) * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][2][2] =
      -0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[2][11][2][2][3] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][3][2] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][3][3] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][4][4] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][5][5] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][6][6] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][2][7][7] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][0][0] =
      std::sqrt(0.2e1) * II * (OQu_1b21b - OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][0][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][1][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][1][1] =
      -std::sqrt(0.2e1) * II * (OQu_1b21b + OQu_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][2][2] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][2][3] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][3][2] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][3][3] =
      -0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQu_2b22b;
  Curvature_Quark_F2H3[2][11][3][4][4] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][5][5] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][6][6] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][3][7][7] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][4][0][6] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][4][0][7] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][4][1][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][4][1][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][4][2][4] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][4][3][4] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][4][4][2] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][4][4][3] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][4][6][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][4][6][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][4][7][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][4][7][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][5][0][6] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][5][0][7] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][5][1][6] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][5][1][7] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][5][2][5] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][5][3][5] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][5][5][2] = -std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][5][5][3] =
      -II * std::sqrt(0.2e1) * OQu_2b11b / 0.2e1;
  Curvature_Quark_F2H3[2][11][5][6][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][5][6][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][5][7][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][5][7][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][6][0][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][6][0][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][6][1][4] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][6][1][5] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][6][2][6] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][6][3][6] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][6][4][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][6][4][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][6][5][0] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][6][5][1] = std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][6][6][2] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][6][6][3] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][7][0][4] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][7][0][5] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][7][1][4] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][7][1][5] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][7][2][7] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][7][3][7] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][7][4][0] =
      II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][7][4][1] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][7][5][0] = -std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][7][5][1] =
      -II * std::sqrt(0.2e1) * OQu_1b21b / 0.4e1;
  Curvature_Quark_F2H3[2][11][7][7][2] = -std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[2][11][7][7][3] =
      -II * std::sqrt(0.2e1) * OQu_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][0][0][2] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][0][0][3] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][0][1][2] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][0][1][3] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][0][2][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][0][2][1] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][0][3][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][0][3][1] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][0][4][6] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][0][4][7] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][0][5][6] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][0][5][7] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][0][6][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][0][6][5] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][0][7][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][0][7][5] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][1][0][2] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][1][0][3] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][1][1][2] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][1][1][3] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][1][2][0] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][1][2][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][1][3][0] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][1][3][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][1][4][6] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][1][4][7] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][1][5][6] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][1][5][7] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][1][6][4] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][1][6][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][1][7][4] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][1][7][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][2][0][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][0][1] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][1][0] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][1][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][2][2] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[3][6][2][2][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][3][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][3][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][4][4] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][5][5] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][6][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][2][7][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][0][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][0][1] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][1][0] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][1][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][2][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][2][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][3][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][3][3] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[3][6][3][4][4] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][5][5] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][6][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][3][7][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][4][0][6] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][4][0][7] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][4][1][6] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][4][1][7] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][4][2][4] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][4][3][4] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][4][4][2] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][4][4][3] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][4][6][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][4][6][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][4][7][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][4][7][1] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][5][0][6] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][5][0][7] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][5][1][6] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][5][1][7] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][5][2][5] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][5][3][5] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][5][5][2] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][5][5][3] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][6][5][6][0] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][5][6][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][5][7][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][5][7][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][6][0][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][6][0][5] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][6][1][4] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][6][1][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][6][2][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][6][3][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][6][4][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][6][4][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][6][5][0] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][6][5][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][6][6][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][6][6][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][7][0][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][7][0][5] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][7][1][4] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][7][1][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][7][2][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][7][3][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][7][4][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][7][4][1] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][7][5][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][7][5][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][6][7][7][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][6][7][7][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][0][0][6] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][0][0][7] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][0][2][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][0][2][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][0][3][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][0][3][5] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][0][4][2] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][0][4][3] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][0][5][2] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][0][5][3] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][0][6][0] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][0][7][0] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][1][1][6] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][1][1][7] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][1][2][4] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][1][2][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][1][3][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][1][3][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][1][4][2] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][1][4][3] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][1][5][2] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][1][5][3] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][1][6][1] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][1][7][1] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][2][0][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][2][0][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][2][1][4] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][2][1][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][2][2][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][2][2][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][2][4][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][2][4][1] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][2][5][0] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][2][5][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][2][6][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][2][7][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][3][0][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][3][0][5] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][3][1][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][3][1][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][3][3][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][3][3][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][3][4][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][3][4][1] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][3][5][0] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][3][5][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][3][6][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][3][7][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][4][0][2] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][4][0][3] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][4][1][2] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][4][1][3] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][4][2][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][4][2][1] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][4][3][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][4][3][1] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][4][4][6] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][4][4][7] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][4][5][6] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][4][5][7] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][4][6][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][4][6][5] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][4][7][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][4][7][5] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][5][0][2] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][5][0][3] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][5][1][2] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][5][1][3] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][5][2][0] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][5][2][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][5][3][0] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][5][3][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[3][9][5][4][6] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][5][4][7] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][5][5][6] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][5][5][7] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][5][6][4] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][5][6][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][5][7][4] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][5][7][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][0][0] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][1][1] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][2][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][3][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][4][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][4][5] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][5][4] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][5][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[3][9][6][6][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][7][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][6][7][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][0][0] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][1][1] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][2][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][3][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][4][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][4][5] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][5][4] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][5][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][6][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][6][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][7][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[3][9][7][7][7] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[4][7][0][0][2] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][0][0][3] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][0][1][2] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][0][1][3] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][0][2][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][0][2][1] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][0][3][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][0][3][1] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][0][4][6] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][0][4][7] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][0][5][6] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][0][5][7] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][0][6][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][0][6][5] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][0][7][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][0][7][5] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][1][0][2] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][1][0][3] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][1][1][2] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][1][1][3] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][1][2][0] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][1][2][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][1][3][0] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][1][3][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][1][4][6] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][1][4][7] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][1][5][6] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][1][5][7] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][1][6][4] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][1][6][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][1][7][4] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][1][7][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][2][0][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][0][1] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][1][0] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][1][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][2][2] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[4][7][2][2][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][3][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][3][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][4][4] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][5][5] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][6][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][2][7][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][0][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][0][1] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][1][0] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][1][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][2][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][2][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][3][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][3][3] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[4][7][3][4][4] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][5][5] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][6][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][3][7][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][4][0][6] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][4][0][7] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][4][1][6] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][4][1][7] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][4][2][4] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][4][3][4] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][4][4][2] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][4][4][3] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][4][6][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][4][6][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][4][7][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][4][7][1] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][5][0][6] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][5][0][7] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][5][1][6] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][5][1][7] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][5][2][5] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][5][3][5] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][5][5][2] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][5][5][3] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][7][5][6][0] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][5][6][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][5][7][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][5][7][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][6][0][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][6][0][5] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][6][1][4] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][6][1][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][6][2][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][6][3][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][6][4][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][6][4][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][6][5][0] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][6][5][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][6][6][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][6][6][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][7][0][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][7][0][5] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][7][1][4] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][7][1][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][7][2][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][7][3][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][7][4][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][7][4][1] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][7][5][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][7][5][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][7][7][7][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][7][7][7][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][0][0][6] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][0][0][7] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][0][2][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][0][2][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][0][3][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][0][3][5] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][0][4][2] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][0][4][3] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][0][5][2] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][0][5][3] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][0][6][0] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][0][7][0] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][1][1][6] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][1][1][7] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][1][2][4] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][1][2][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][1][3][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][1][3][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][1][4][2] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][1][4][3] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][1][5][2] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][1][5][3] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][1][6][1] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][1][7][1] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][2][0][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][2][0][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][2][1][4] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][2][1][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][2][2][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][2][2][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][2][4][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][2][4][1] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][2][5][0] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][2][5][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][2][6][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][2][7][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][3][0][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][3][0][5] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][3][1][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][3][1][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][3][3][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][3][3][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][3][4][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][3][4][1] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][3][5][0] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][3][5][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][3][6][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][3][7][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][4][0][2] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][4][0][3] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][4][1][2] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][4][1][3] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][4][2][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][4][2][1] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][4][3][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][4][3][1] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][4][4][6] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][4][4][7] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][4][5][6] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][4][5][7] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][4][6][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][4][6][5] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][4][7][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][4][7][5] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][5][0][2] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][5][0][3] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][5][1][2] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][5][1][3] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][5][2][0] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][5][2][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][5][3][0] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][5][3][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[4][10][5][4][6] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][5][4][7] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][5][5][6] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][5][5][7] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][5][6][4] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][5][6][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][5][7][4] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][5][7][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][0][0] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][1][1] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][2][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][3][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][4][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][4][5] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][5][4] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][5][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[4][10][6][6][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][7][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][6][7][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][0][0] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][1][1] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][2][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][3][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][4][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][4][5] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][5][4] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][5][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][6][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][6][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][7][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[4][10][7][7][7] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[5][8][0][0][2] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][0][3] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][1][2] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][1][3] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][2][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][2][1] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][3][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][3][1] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][0][4][6] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][4][7] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][5][6] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][5][7] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][6][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][6][5] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][7][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][0][7][5] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][0][2] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][0][3] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][1][2] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][1][3] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][2][0] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][2][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][3][0] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][3][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][1][4][6] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][4][7] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][5][6] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][5][7] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][6][4] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][6][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][7][4] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][1][7][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][2][0][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][0][1] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][1][0] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][1][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][2][2] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[5][8][2][2][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][3][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][3][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][4][4] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][5][5] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][6][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][2][7][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][0][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][0][1] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][1][0] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][1][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][2][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][2][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][3][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][3][3] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[5][8][3][4][4] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][5][5] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][6][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][3][7][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][4][0][6] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][4][0][7] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][4][1][6] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][4][1][7] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][4][2][4] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][4][3][4] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][4][4][2] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][4][4][3] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][4][6][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][4][6][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][4][7][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][4][7][1] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][5][0][6] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][5][0][7] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][5][1][6] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][5][1][7] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][5][2][5] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][5][3][5] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][5][5][2] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][5][5][3] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][8][5][6][0] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][5][6][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][5][7][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][5][7][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][6][0][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][6][0][5] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][6][1][4] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][6][1][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][6][2][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][6][3][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][6][4][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][6][4][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][6][5][0] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][6][5][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][6][6][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][6][6][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][7][0][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][7][0][5] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][7][1][4] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][7][1][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][7][2][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][7][3][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][7][4][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][7][4][1] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][7][5][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][7][5][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][8][7][7][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][8][7][7][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][0][0][6] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][0][0][7] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][0][2][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][2][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][3][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][3][5] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][4][2] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][4][3] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][5][2] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][5][3] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][0][6][0] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][0][7][0] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][1][1][6] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][1][1][7] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][1][2][4] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][2][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][3][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][3][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][4][2] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][4][3] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][5][2] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][5][3] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][1][6][1] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][1][7][1] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][2][0][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][2][0][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][2][1][4] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][2][1][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][2][2][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][2][2][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][2][4][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][2][4][1] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][2][5][0] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][2][5][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][2][6][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][2][7][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][3][0][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][3][0][5] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][3][1][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][3][1][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][3][3][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][3][3][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][3][4][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][3][4][1] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][3][5][0] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][3][5][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][3][6][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][3][7][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][0][2] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][4][0][3] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][4][1][2] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][4][1][3] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][4][2][0] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][4][2][1] =
      -std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][4][3][0] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][4][3][1] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][4][4][6] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][4][7] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][5][6] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][5][7] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][6][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][6][5] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][7][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][4][7][5] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][0][2] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][5][0][3] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][5][1][2] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][5][1][3] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][5][2][0] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][5][2][1] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][5][3][0] =
      -std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][5][3][1] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b) / 0.4e1;
  Curvature_Quark_F2H3[5][11][5][4][6] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][4][7] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][5][6] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][5][7] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][6][4] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][6][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][7][4] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][5][7][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][0][0] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][1][1] = std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][2][2] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][3][3] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][4][4] =
      std::sqrt(0.2e1) * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][4][5] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][5][4] =
      II * std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][5][5] =
      std::sqrt(0.2e1) * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][6][6] =
      0.3e1 / 0.2e1 * std::sqrt(0.2e1) * OQd_2b22b;
  Curvature_Quark_F2H3[5][11][6][6][7] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][7][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][6][7][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][0][0] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][1][1] =
      II * std::sqrt(0.2e1) * OQd_2b11b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][2][2] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][3][3] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][4][4] =
      std::sqrt(0.2e1) * II * (OQd_1b12b - OQd_1b21b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][4][5] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][5][4] = std::sqrt(0.2e1) * OQd_1b21b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][5][5] =
      std::sqrt(0.2e1) * II * (OQd_1b21b + OQd_1b12b + OQd_2b11b) / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][6][6] =
      II * std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][6][7] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][7][6] = std::sqrt(0.2e1) * OQd_2b22b / 0.2e1;
  Curvature_Quark_F2H3[5][11][7][7][7] =
      0.3e1 / 0.2e1 * II * std::sqrt(0.2e1) * OQd_2b22b;

  sym5Dim(Curvature_Quark_F2H3, NQuarks, NQuarks, NHiggs, NHiggs, NHiggs);

  SetUseTensorSymFac(true); // true for whole SymFac*Tensor input

  SymFac_Higgs[0][0] =
      std::pow(LambdaEFT, -0.2e1) *
      (-2 * (6 * Op6_111111 + 2 * Op6_111122 + Op6_112222 + Op6_122111) -
       Op6_122122) /
      48.; // rho1rho1
  SymFac_Higgs[0][2] =
      std::pow(LambdaEFT, -0.2e1) *
      (-4 * Op6_111112 - 3 * Op6_112212 - 2 * (Op6_121221 + 2 * Op6_122222)) /
      48.; // rho1rho2
  SymFac_Higgs[1][1] =
      std::pow(LambdaEFT, -0.2e1) *
      (-2 * (6 * Op6_111111 + 2 * Op6_111122 + Op6_112222 + Op6_122111) -
       Op6_122122) /
      48.; // eta1eta1
  SymFac_Higgs[1][3] =
      std::pow(LambdaEFT, -0.2e1) *
      (-4 * Op6_111112 - 3 * Op6_112212 - 2 * (Op6_121221 + 2 * Op6_122222)) /
      48.; // eta1eta2
  SymFac_Higgs[2][0] =
      std::pow(LambdaEFT, -0.2e1) *
      (-4 * Op6_111112 - 3 * Op6_112212 - 2 * (Op6_121221 + 2 * Op6_122222)) /
      48.; // rho2rho1
  SymFac_Higgs[2][2] = std::pow(LambdaEFT, -0.2e1) *
                       (-2 * Op6_111122 - 4 * Op6_112222 - Op6_122111 -
                        2 * (Op6_122122 + 6 * Op6_222222)) /
                       48.; // rho2rho2
  SymFac_Higgs[3][1] =
      std::pow(LambdaEFT, -0.2e1) *
      (-4 * Op6_111112 - 3 * Op6_112212 - 2 * (Op6_121221 + 2 * Op6_122222)) /
      48.; // eta2eta1
  SymFac_Higgs[3][3] = std::pow(LambdaEFT, -0.2e1) *
                       (-2 * Op6_111122 - 4 * Op6_112222 - Op6_122111 -
                        2 * (Op6_122122 + 6 * Op6_222222)) /
                       48.; // eta2eta2

  SymFac_Higgs[4][4] =
      std::pow(LambdaEFT, -0.2e1) *
      (-2 * (6 * Op6_111111 + 2 * Op6_111122 + Op6_112222 + Op6_122111) -
       Op6_122122) /
      48.; // zeta1zeta1
  SymFac_Higgs[4][6] =
      std::pow(LambdaEFT, -0.2e1) *
      (-4 * Op6_111112 - 3 * Op6_112212 - 2 * (Op6_121221 + 2 * Op6_122222)) /
      48.; // zeta1zeta2
  SymFac_Higgs[5][5] =
      std::pow(LambdaEFT, -0.2e1) *
      (-2 * (6 * Op6_111111 + 2 * Op6_111122 + Op6_112222 + Op6_122111) -
       Op6_122122) /
      48.; // psi1psi1
  SymFac_Higgs[5][7] =
      std::pow(LambdaEFT, -0.2e1) *
      (-4 * Op6_111112 - 3 * Op6_112212 - 2 * (Op6_121221 + 2 * Op6_122222)) /
      48.; // psi1psi2
  SymFac_Higgs[6][4] =
      std::pow(LambdaEFT, -0.2e1) *
      (-4 * Op6_111112 - 3 * Op6_112212 - 2 * (Op6_121221 + 2 * Op6_122222)) /
      48.; // zeta2zeta1
  SymFac_Higgs[6][6] = std::pow(LambdaEFT, -0.2e1) *
                       (-2 * Op6_111122 - 4 * Op6_112222 - Op6_122111 -
                        2 * (Op6_122122 + 6 * Op6_222222)) /
                       48.; // zeta2zeta2
  SymFac_Higgs[7][5] =
      std::pow(LambdaEFT, -0.2e1) *
      (-4 * Op6_111112 - 3 * Op6_112212 - 2 * (Op6_121221 + 2 * Op6_122222)) /
      48.; // psi2psi1
  SymFac_Higgs[7][7] = std::pow(LambdaEFT, -0.2e1) *
                       (-2 * Op6_111122 - 4 * Op6_112222 - Op6_122111 -
                        2 * (Op6_122122 + 6 * Op6_222222)) /
                       48.; // psi2psi2

  // SymFac_Gauge is independent of Op6

  SetCurvatureDone = true;
}

bool Class_Potential_R2HDMEFTPHI6_PHI2PSI3::CalculateDebyeSimplified()
{
  // not implemented
  return false;
}

bool Class_Potential_R2HDMEFTPHI6_PHI2PSI3::CalculateDebyeGaugeSimplified()
{
  // not implemented
  return false;
}

double Class_Potential_R2HDMEFTPHI6_PHI2PSI3::VTreeSimplified(
    const std::vector<double> &v) const
{
  (void)v;
  double res = 0;
  // not implemented
  return res;
}

double Class_Potential_R2HDMEFTPHI6_PHI2PSI3::VCounterSimplified(
    const std::vector<double> &v) const
{
  (void)v;
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  // not implemented
  return res;
}

void Class_Potential_R2HDMEFTPHI6_PHI2PSI3::Debugging(
    const std::vector<double> &input,
    std::vector<double> &output) const
{
  (void)input;
  // not implemented
  (void)output;
}

} // namespace Models
} // namespace BSMPT
