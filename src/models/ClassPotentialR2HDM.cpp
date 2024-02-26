// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialR2HDM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>

namespace BSMPT
{
namespace Models
{

Class_Potential_R2HDM::Class_Potential_R2HDM(const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  // TODO Auto-generated constructor stub
  Model         = ModelID::ModelIDs::R2HDM;
  NNeutralHiggs = 4;
  NChargedHiggs = 4;

  NHiggs  = NNeutralHiggs + NChargedHiggs;
  NGauge  = 4;
  NLepton = 9;
  NQuarks = 12;

  nPar   = 8;
  nParCT = 11;

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

Class_Potential_R2HDM::~Class_Potential_R2HDM()
{
  // TODO Auto-generated destructor stub
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given inputfile
 */
std::vector<std::string> Class_Potential_R2HDM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("Dm11sq");
  labels.push_back("Dm22sq");
  labels.push_back("Dm12sq");
  labels.push_back("DL1");
  labels.push_back("DL2");
  labels.push_back("DL3");
  labels.push_back("DL4");
  labels.push_back("DL5");
  labels.push_back("DT1");
  labels.push_back("DT2");
  labels.push_back("DT3");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * inputfile
 */
std::vector<std::string> Class_Potential_R2HDM::addLegendTemp() const
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
std::vector<std::string> Class_Potential_R2HDM::addLegendVEV() const
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
std::vector<std::string> Class_Potential_R2HDM::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;

  particles.push_back("G^+");
  particles.push_back("G^-");
  particles.push_back("H^+");
  particles.push_back("H^-");
  particles.push_back("G^0");
  particles.push_back("A");
  particles.push_back("H");
  particles.push_back("h");

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

void Class_Potential_R2HDM::ReadAndSet(const std::string &linestr,
                                       std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 8; k++)
  {
    ss >> tmp;
    if (k == 1)
      Type = tmp;
    else if (k == 2)
      L1 = tmp;
    else if (k == 3)
      L2 = tmp;
    else if (k == 4)
      L3 = tmp;
    else if (k == 5)
      L4 = tmp;
    else if (k == 6)
      RL5 = tmp;
    else if (k == 7)
      RealMMix = tmp;
    else if (k == 8)
      TanBeta = tmp;
  }

  //	double sa = std::sin(alpha);
  //	double ca = std::cos(alpha);
  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = sqrt(C_SinBetaSquared);
  //	L1 =1.0 / (SMConstants.C_vev0 * SMConstants.C_vev0 * C_CosBeta *
  // C_CosBeta)* (ca * ca * MH * MH +
  // sa * sa * Mh * Mh- RealMMix * C_SinBeta / C_CosBeta); 	L2 =1.0 /
  // (SMConstants.C_vev0 * SMConstants.C_vev0 * C_SinBeta * C_SinBeta)* (sa * sa
  // * MH * MH + ca * ca * Mh * Mh- RealMMix * C_CosBeta / C_SinBeta); 	L3 = 2
  // * MHP * MHP / (SMConstants.C_vev0 * SMConstants.C_vev0)+ sa * ca * (MH * MH
  // - Mh * Mh)/ (SMConstants.C_vev0 * SMConstants.C_vev0 *C_CosBeta*C_SinBeta
  // )- RealMMix / (SMConstants.C_vev0 * SMConstants.C_vev0 * C_SinBeta *
  // C_CosBeta); 	L4 = (MA * MA - 2 * MHP * MHP) / (SMConstants.C_vev0 *
  // SMConstants.C_vev0)+ RealMMix / (SMConstants.C_vev0 * SMConstants.C_vev0 *
  // C_SinBeta * C_CosBeta); 	RL5 = RealMMix / (SMConstants.C_vev0 *
  // SMConstants.C_vev0 * C_SinBeta * C_CosBeta) - MA
  //* MA / (SMConstants.C_vev0 * SMConstants.C_vev0);

  par[6] = TanBeta;
  par[4] = RL5;
  par[0] = L1;
  par[1] = L2;
  par[2] = L3;
  par[3] = L4;
  par[5] = RealMMix;
  par[7] = Type;

  set_gen(par);
  return;
}

/**
 * Set Class Object with an CP-Conserving Point
 */
void Class_Potential_R2HDM::set_gen(const std::vector<double> &par)
{

  // double *p = (double *)par;
  scale = SMConstants.C_vev0;
  //	scale=C_MassZ;
  L1               = par[0];
  L2               = par[1];
  L3               = par[2];
  L4               = par[3];
  RL5              = par[4];
  RealMMix         = par[5];
  TanBeta          = par[6];
  beta             = std::atan(TanBeta);
  Type             = static_cast<int>(par[7]);
  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = sqrt(C_SinBetaSquared);

  if (TanBeta < 0) // for beta in 4th quadrant
  {
    C_SinBeta *= -1;
  }

  u1 = RealMMix * TanBeta -
       SMConstants.C_vev0 * SMConstants.C_vev0 * C_SinBetaSquared *
           (L4 + RL5 + L3) / 0.2e1 -
       SMConstants.C_vev0 * SMConstants.C_vev0 * C_CosBetaSquared * L1 / 0.2e1;
  u2 = RealMMix * 1.0 / TanBeta -
       SMConstants.C_vev0 * SMConstants.C_vev0 * C_CosBetaSquared *
           (L4 + RL5 + L3) / 0.2e1 -
       SMConstants.C_vev0 * SMConstants.C_vev0 * C_SinBetaSquared * L2 / 0.2e1;

  //	double ML5 =
  // 2*RealMMix/(SMConstants.C_vev0*SMConstants.C_vev0*C_SinBeta*C_CosBeta);
  //	double TripleHiggs =
  //-3.0/(SMConstants.C_vev0*std::sin(2*beta))*(Mh*Mh*(2*std::cos(alpha+beta)+std::sin(2*alpha)*std::sin(beta-alpha))
  //-
  // std::cos(alpha+beta)*std::pow(std::cos(beta-alpha),2)*SMConstants.C_vev0*SMConstants.C_vev0*ML5
  //);

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

void Class_Potential_R2HDM::set_CT_Pot_Par(const std::vector<double> &p)
{
  //	double *p = (double *)par;

  Du1CT  = p[0];
  Du2CT  = p[1];
  DRu3CT = p[2];
  DL1CT  = p[3];
  DL2CT  = p[4];
  DL3CT  = p[5];
  DL4CT  = p[6];
  DRL5CT = p[7];

  DT1 = p[8];
  DT2 = p[9];
  DT3 = p[10];

  double DIL5CT = 0;
  double DIu3CT = 0;

  Curvature_Higgs_CT_L1[0] = 0;
  Curvature_Higgs_CT_L1[1] = 0;
  Curvature_Higgs_CT_L1[2] = 0;
  Curvature_Higgs_CT_L1[3] = DTCharged;
  Curvature_Higgs_CT_L1[4] = DT1;
  Curvature_Higgs_CT_L1[5] = 0;
  Curvature_Higgs_CT_L1[6] = DT2;
  Curvature_Higgs_CT_L1[7] = DT3;

  Curvature_Higgs_CT_L2[0][0] = Du1CT;
  Curvature_Higgs_CT_L2[0][1] = 0;
  Curvature_Higgs_CT_L2[0][2] = -DRu3CT;
  Curvature_Higgs_CT_L2[0][3] = DIu3CT;
  Curvature_Higgs_CT_L2[0][4] = 0;
  Curvature_Higgs_CT_L2[0][5] = 0;
  Curvature_Higgs_CT_L2[0][6] = 0;
  Curvature_Higgs_CT_L2[0][7] = 0;
  Curvature_Higgs_CT_L2[1][0] = 0;
  Curvature_Higgs_CT_L2[1][1] = Du1CT;
  Curvature_Higgs_CT_L2[1][2] = -DIu3CT;
  Curvature_Higgs_CT_L2[1][3] = -DRu3CT;
  Curvature_Higgs_CT_L2[1][4] = 0;
  Curvature_Higgs_CT_L2[1][5] = 0;
  Curvature_Higgs_CT_L2[1][6] = 0;
  Curvature_Higgs_CT_L2[1][7] = 0;
  Curvature_Higgs_CT_L2[2][0] = -DRu3CT;
  Curvature_Higgs_CT_L2[2][1] = -DIu3CT;
  Curvature_Higgs_CT_L2[2][2] = Du2CT;
  Curvature_Higgs_CT_L2[2][3] = 0;
  Curvature_Higgs_CT_L2[2][4] = 0;
  Curvature_Higgs_CT_L2[2][5] = 0;
  Curvature_Higgs_CT_L2[2][6] = 0;
  Curvature_Higgs_CT_L2[2][7] = 0;
  Curvature_Higgs_CT_L2[3][0] = DIu3CT;
  Curvature_Higgs_CT_L2[3][1] = -DRu3CT;
  Curvature_Higgs_CT_L2[3][2] = 0;
  Curvature_Higgs_CT_L2[3][3] = Du2CT;
  Curvature_Higgs_CT_L2[3][4] = 0;
  Curvature_Higgs_CT_L2[3][5] = 0;
  Curvature_Higgs_CT_L2[3][6] = 0;
  Curvature_Higgs_CT_L2[3][7] = 0;
  Curvature_Higgs_CT_L2[4][0] = 0;
  Curvature_Higgs_CT_L2[4][1] = 0;
  Curvature_Higgs_CT_L2[4][2] = 0;
  Curvature_Higgs_CT_L2[4][3] = 0;
  Curvature_Higgs_CT_L2[4][4] = Du1CT;
  Curvature_Higgs_CT_L2[4][5] = 0;
  Curvature_Higgs_CT_L2[4][6] = -DRu3CT;
  Curvature_Higgs_CT_L2[4][7] = DIu3CT;
  Curvature_Higgs_CT_L2[5][0] = 0;
  Curvature_Higgs_CT_L2[5][1] = 0;
  Curvature_Higgs_CT_L2[5][2] = 0;
  Curvature_Higgs_CT_L2[5][3] = 0;
  Curvature_Higgs_CT_L2[5][4] = 0;
  Curvature_Higgs_CT_L2[5][5] = Du1CT;
  Curvature_Higgs_CT_L2[5][6] = -DIu3CT;
  Curvature_Higgs_CT_L2[5][7] = -DRu3CT;
  Curvature_Higgs_CT_L2[6][0] = 0;
  Curvature_Higgs_CT_L2[6][1] = 0;
  Curvature_Higgs_CT_L2[6][2] = 0;
  Curvature_Higgs_CT_L2[6][3] = 0;
  Curvature_Higgs_CT_L2[6][4] = -DRu3CT;
  Curvature_Higgs_CT_L2[6][5] = -DIu3CT;
  Curvature_Higgs_CT_L2[6][6] = Du2CT;
  Curvature_Higgs_CT_L2[6][7] = 0;
  Curvature_Higgs_CT_L2[7][0] = 0;
  Curvature_Higgs_CT_L2[7][1] = 0;
  Curvature_Higgs_CT_L2[7][2] = 0;
  Curvature_Higgs_CT_L2[7][3] = 0;
  Curvature_Higgs_CT_L2[7][4] = DIu3CT;
  Curvature_Higgs_CT_L2[7][5] = -DRu3CT;
  Curvature_Higgs_CT_L2[7][6] = 0;
  Curvature_Higgs_CT_L2[7][7] = Du2CT;

  Curvature_Higgs_CT_L4[0][0][0][0] = 3 * DL1CT;
  Curvature_Higgs_CT_L4[0][0][1][1] = DL1CT;
  Curvature_Higgs_CT_L4[0][0][2][2] = DL3CT + DL4CT + DRL5CT;
  Curvature_Higgs_CT_L4[0][0][2][3] = -DIL5CT;
  Curvature_Higgs_CT_L4[0][0][3][3] = DL3CT + DL4CT - DRL5CT;
  Curvature_Higgs_CT_L4[0][0][4][4] = DL1CT;
  Curvature_Higgs_CT_L4[0][0][5][5] = DL1CT;
  Curvature_Higgs_CT_L4[0][0][6][6] = DL3CT;
  Curvature_Higgs_CT_L4[0][0][7][7] = DL3CT;
  Curvature_Higgs_CT_L4[0][1][2][2] = DIL5CT;
  Curvature_Higgs_CT_L4[0][1][2][3] = DRL5CT;
  Curvature_Higgs_CT_L4[0][1][3][3] = -DIL5CT;
  Curvature_Higgs_CT_L4[0][2][4][6] = DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][4][7] = -DIL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][5][6] = DIL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][5][7] = DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][4][6] = -DIL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][4][7] = DL4CT / 0.2e1 - DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][5][6] = -DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][5][7] = -DIL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][1][1] = 3 * DL1CT;
  Curvature_Higgs_CT_L4[1][1][2][2] = DL3CT + DL4CT - DRL5CT;
  Curvature_Higgs_CT_L4[1][1][2][3] = DIL5CT;
  Curvature_Higgs_CT_L4[1][1][3][3] = DL3CT + DL4CT + DRL5CT;
  Curvature_Higgs_CT_L4[1][1][4][4] = DL1CT;
  Curvature_Higgs_CT_L4[1][1][5][5] = DL1CT;
  Curvature_Higgs_CT_L4[1][1][6][6] = DL3CT;
  Curvature_Higgs_CT_L4[1][1][7][7] = DL3CT;
  Curvature_Higgs_CT_L4[1][2][4][6] = DIL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][4][7] = -DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][5][6] = DL4CT / 0.2e1 - DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][5][7] = DIL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][4][6] = DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][4][7] = -DIL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][5][6] = DIL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][5][7] = DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][2][2] = 3 * DL2CT;
  Curvature_Higgs_CT_L4[2][2][3][3] = DL2CT;
  Curvature_Higgs_CT_L4[2][2][4][4] = DL3CT;
  Curvature_Higgs_CT_L4[2][2][5][5] = DL3CT;
  Curvature_Higgs_CT_L4[2][2][6][6] = DL2CT;
  Curvature_Higgs_CT_L4[2][2][7][7] = DL2CT;
  Curvature_Higgs_CT_L4[3][3][3][3] = 3 * DL2CT;
  Curvature_Higgs_CT_L4[3][3][4][4] = DL3CT;
  Curvature_Higgs_CT_L4[3][3][5][5] = DL3CT;
  Curvature_Higgs_CT_L4[3][3][6][6] = DL2CT;
  Curvature_Higgs_CT_L4[3][3][7][7] = DL2CT;
  Curvature_Higgs_CT_L4[4][4][4][4] = 3 * DL1CT;
  Curvature_Higgs_CT_L4[4][4][5][5] = DL1CT;
  Curvature_Higgs_CT_L4[4][4][6][6] = DL3CT + DL4CT + DRL5CT;
  Curvature_Higgs_CT_L4[4][4][6][7] = -DIL5CT;
  Curvature_Higgs_CT_L4[4][4][7][7] = DL3CT + DL4CT - DRL5CT;
  Curvature_Higgs_CT_L4[4][5][6][6] = DIL5CT;
  Curvature_Higgs_CT_L4[4][5][6][7] = DRL5CT;
  Curvature_Higgs_CT_L4[4][5][7][7] = -DIL5CT;
  Curvature_Higgs_CT_L4[5][5][5][5] = 3 * DL1CT;
  Curvature_Higgs_CT_L4[5][5][6][6] = DL3CT + DL4CT - DRL5CT;
  Curvature_Higgs_CT_L4[5][5][6][7] = DIL5CT;
  Curvature_Higgs_CT_L4[5][5][7][7] = DL3CT + DL4CT + DRL5CT;
  Curvature_Higgs_CT_L4[6][6][6][6] = 3 * DL2CT;
  Curvature_Higgs_CT_L4[6][6][7][7] = DL2CT;
  Curvature_Higgs_CT_L4[7][7][7][7] = 3 * DL2CT;

  for (std::size_t k1 = 0; k1 < NHiggs; k1++)
  {
    for (std::size_t k2 = k1; k2 < NHiggs; k2++)
    {
      for (std::size_t k3 = k2; k3 < NHiggs; k3++)
      {
        for (std::size_t k4 = k3; k4 < NHiggs; k4++)
        {
          Curvature_Higgs_CT_L4[k2][k3][k4][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k4][k1][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k1][k2][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k1][k3][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k2][k1][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k4][k2][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k3][k4][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k2][k1][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k3][k2][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k4][k3][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k1][k4][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k2][k3][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k4][k2][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k1][k4][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k3][k1][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k3][k2][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k1][k3][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k4][k1][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k2][k4][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k2][k4][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k1][k2][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k3][k1][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k4][k3][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
        }
      }
    }
  }

  return;
}

/**
 * Console-Output of all Parameters
 */
void Class_Potential_R2HDM::write() const
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
  ss << "Lambda1 = " << L1 << "\n";
  ss << "Lambda2 = " << L2 << "\n";
  ss << "Lambda3 = " << L3 << "\n";
  ss << "Lambda4 = " << L4 << "\n";
  ss << "Re(Lambda5) = " << RL5 << "\n";
  ss << "Re(m_12^2) = " << RealMMix << "\n";
  ss << "m_{11}^2 = " << u1 << "\n";
  ss << "m_{22}^2 = " << u2 << "\n";

  ss << "The counterterms are :\n";

  ss << "DL1 := " << DL1CT << ";\n";
  ss << "DL2 := " << DL2CT << ";\n";
  ss << "DL3 := " << DL3CT << ";\n";
  ss << "DL4 := " << DL4CT << ";\n";
  ss << "DRL5 := " << DRL5CT << ";\n";
  ss << "Du1 := " << Du1CT << ";\n";
  ss << "Du2 := " << Du2CT << ";\n";
  ss << "DRu3 := " << DRu3CT << ";\n";

  ss << "DT1 := " << DT1 << ";\n";
  ss << "DT2 := " << DT2 << ";\n";
  ss << "DT3:= " << DT3 << ";\n";

  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  ss << "The mass spectrum is given by :\n"
     << "m_{G^+} = " << std::sqrt(HiggsMasses[posGp]) << " GeV \n"
     << "m_{G^-} = " << std::sqrt(HiggsMasses[posGm]) << " GeV \n"
     << "m_{H^+} = " << std::sqrt(HiggsMasses[posHp]) << " GeV \n"
     << "m_{H^-} = " << std::sqrt(HiggsMasses[posHm]) << " GeV \n"
     << "m_{G^0} = " << std::sqrt(HiggsMasses[posG0]) << " GeV \n"
     << "m_A = " << std::sqrt(HiggsMasses[posA]) << " GeV \n"
     << "m_H = " << std::sqrt(HiggsMasses[posH]) << " GeV \n"
     << "m_h = " << std::sqrt(HiggsMasses[posh]) << " GeV \n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms in the 2HDM
 */
std::vector<double> Class_Potential_R2HDM::calc_CT() const
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

  double freepar = 0; // Value of DL4CT

  // Du1CT
  parCT.push_back(
      -(double)((-2 * freepar * v1 * v2 * v2 + 5 * HesseWeinberg(0, 0) * v1 +
                 HesseWeinberg(1, 3) * v2 - HesseWeinberg(4, 6) * v2 -
                 HesseWeinberg(4, 4) * v1 - 2 * HesseWeinberg(5, 5) * v1) /
                v1) /
      0.2e1);
  // Du2CT
  parCT.push_back((double)((2 * freepar * v1 * v1 * v2 * v2 +
                            HesseWeinberg(6, 6) * v2 * v2 -
                            2 * HesseWeinberg(0, 0) * v1 * v1 -
                            HesseWeinberg(1, 3) * v1 * v2 -
                            3 * HesseWeinberg(3, 3) * v2 * v2 +
                            HesseWeinberg(4, 6) * v1 * v2 +
                            2 * v1 * v1 * HesseWeinberg(5, 5)) *
                           (double)pow((double)v2, (double)(-2))) /
                  0.2e1);
  // DRu3CT
  parCT.push_back(-(-freepar * v1 * v2 * v2 + HesseWeinberg(0, 0) * v1 -
                    HesseWeinberg(1, 3) * v2 - HesseWeinberg(5, 5) * v1) /
                  v2);
  // DL1CT
  parCT.push_back((double)(-freepar * v2 * v2 + 2 * HesseWeinberg(0, 0) -
                           HesseWeinberg(4, 4) - HesseWeinberg(5, 5)) *
                  pow(v1, -0.2e1));
  // DL2CT
  parCT.push_back(
      -(freepar * v1 * v1 * v2 * v2 + HesseWeinberg(6, 6) * v2 * v2 -
        HesseWeinberg(0, 0) * v1 * v1 - HesseWeinberg(3, 3) * v2 * v2 +
        v1 * v1 * HesseWeinberg(5, 5)) *
      pow(v2, -0.4e1));
  // DL3CT
  parCT.push_back((-freepar * v1 * v2 * v2 + HesseWeinberg(0, 0) * v1 +
                   HesseWeinberg(1, 3) * v2 - HesseWeinberg(4, 6) * v2 -
                   HesseWeinberg(5, 5) * v1) /
                  v1 * pow(v2, -0.2e1));
  // DL4CT
  parCT.push_back(freepar);
  // DRL5CT
  parCT.push_back(-(-freepar * v2 * v2 + 2 * HesseWeinberg(0, 0) -
                    2 * HesseWeinberg(5, 5)) *
                  (double)pow((double)v2, (double)(-2)));

  // DT1
  double tmp =
      HesseWeinberg(1, 3) * v2 + HesseWeinberg(0, 0) * v1 - NablaWeinberg(4);
  if (std::abs(tmp) < 1e-9) tmp = 0;
  parCT.push_back(tmp);
  // DT2
  tmp = HesseWeinberg(1, 3) * v1 + HesseWeinberg(3, 3) * v2 - NablaWeinberg(6);
  if (std::abs(tmp) < 1e-9) tmp = 0;
  parCT.push_back(tmp);
  // DT3
  tmp = -(-v1 * v1 * HesseWeinberg(4, 5) - HesseWeinberg(4, 7) * v1 * v2 +
          NablaWeinberg(7) * v2) /
        v2;
  if (std::abs(tmp) < 1e-9) tmp = 0;
  parCT.push_back(tmp);

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
 * Ensures the correct rotation matrix convention
 */
void Class_Potential_R2HDM::AdjustRotationMatrix()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

  if (!CheckRotationMatrix()) // Check whether generically generated rotation
                              // matrix is proper rotation matrix
  {
    throw std::runtime_error("Error in rotation matrix.");
  }

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  // std::cout << HiggsRot << "\n";

  // initialize position indices (new initialization for each point in multiline
  // files)
  posGp = -1, posGm = -1, posHp = -1, posHm = -1, posG0 = -1, posA = -1,
  posH = -1, posh = -1;

  // interaction basis
  // rho1, eta1, rho2, eta2, zeta1, psi1, zeta2, psi2
  int pos_rho1 = 0, pos_eta1 = 1, pos_rho2 = 2, pos_eta2 = 3, pos_zeta1 = 4,
      pos_psi1 = 5, pos_zeta2 = 6, pos_psi2 = 7;

  // higgsbasis = {rho1, eta1, rho2, eta2, zeta1, psi1, zeta2, psi2}
  for (std::size_t i = 0; i < NHiggs;
       i++) // mass base index i corresponds to mass vector sorted in ascending
            // mass
  {
    // charged submatrices
    if (std::abs(HiggsRot(i, pos_rho1)) + std::abs(HiggsRot(i, pos_rho2)) >
        ZeroThreshold)
    {
      if (posGm == -1)
      {
        posGm = i;
      }
      else
      {
        posHp = i;
      }
    }
    if (std::abs(HiggsRot(i, pos_eta1)) + std::abs(HiggsRot(i, pos_eta2)) >
        ZeroThreshold)
    {
      if (posGp == -1)
      {
        posGp = i;
      }
      else
      {
        posHm = i;
      }
    }
    if (std::abs(HiggsRot(i, pos_zeta1)) + std::abs(HiggsRot(i, pos_zeta2)) >
        ZeroThreshold) // use that mh < mH
    {
      if (posh == -1)
      {
        posh = i;
      }
      else
      {
        posH = i;
      }
    }
    if (std::abs(HiggsRot(i, pos_psi1)) + std::abs(HiggsRot(i, pos_psi2)) >
        ZeroThreshold) // use that 0 = mG0 < mA
    {
      if (posG0 == -1)
      {
        posG0 = i;
      }
      else
      {
        posA = i;
      }
    }
  }

  MatrixXd HiggsRotFixed(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotFixed.row(i) = HiggsRot.row(i);
  }

  // charged submatrix
  if (HiggsRotFixed(posGp, pos_eta1) < 0) // Gp eta1 (+ cos(beta))
  {
    HiggsRotFixed.row(posGp) *= -1;
  }
  if (HiggsRotFixed(posGm, pos_rho1) < 0) // Gm rho1 (+ cos(beta))
  {
    HiggsRotFixed.row(posGm) *= -1;
  }
  if (HiggsRotFixed(posHp, pos_rho2) < 0) // Hp rho2 (+ cos(beta))
  {
    HiggsRotFixed.row(posHp) *= -1;
  }
  if (HiggsRotFixed(posHm, pos_eta2) < 0) // Hm eta2 (+ cos(beta))
  {
    HiggsRotFixed.row(posHm) *= -1;
  }

  // check neutral, CP-odd submatrix
  if (HiggsRotFixed(posG0, pos_psi1) < 0) // G0 psi1 (+ cos(beta))
  {
    HiggsRotFixed.row(posG0) *= -1; // G0
  }
  if (HiggsRotFixed(posA, pos_psi2) < 0) // A psi2 (+ cos(beta))
  {
    HiggsRotFixed.row(posA) *= -1; // A
  }

  // // check neutral, CP-even submatrix
  if (HiggsRotFixed(posH, pos_zeta1) < 0) // H zeta1 (+ cos(alpha))
  {
    // if negative, rotate H
    HiggsRotFixed.row(posH) *= -1; // H
  }
  if (HiggsRotFixed(posh, pos_zeta2) < 0) // h zeta2 (+ cos(alpha))
  {
    // if negative, rotate h
    HiggsRotFixed.row(posh) *= -1; // h
  }

  alpha = std::asin(HiggsRotFixed(posH, pos_zeta2)); // H zeta2 (+ sin(alpha))

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRotationMatrixEnsuredConvention[i][j] = HiggsRotFixed(i, j);
    }
  }
}

/**
 * Calculates the corrections to the Triple higgs couplings in the mass basis.
 *
 * Use the vector TripleHiggsCorrectionsCWPhysical to save your couplings and
 * set the nTripleCouplings to the number of couplings you want as output.
 */
void Class_Potential_R2HDM::TripleHiggsCouplings()
{
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

  std::vector<double> HiggsOrder(NHiggs);

  HiggsOrder[0] = posGp;
  HiggsOrder[1] = posGm;
  HiggsOrder[2] = posHp;
  HiggsOrder[3] = posHm;
  HiggsOrder[4] = posG0;
  HiggsOrder[5] = posA;
  HiggsOrder[6] = posH;
  HiggsOrder[7] = posh;

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrixEnsuredConvention[i][j];
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

  CheckTrilinearTreeCouplings();
}

void Class_Potential_R2HDM::CheckTrilinearTreeCouplings()
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);
  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  // std::cout << HiggsMasses << "\n";

  // std::cout << posGp << "," << posGm << "," << posHp << "," << posHm << ","
  //           << posG0 << "," << posA << "," << posH << "," << posh <<
  //           "\n";

  double bma    = beta - alpha;
  double cbma   = std::cos(bma);
  double sbma   = std::sin(bma);
  double cot2b  = 1. / std::tan(2 * beta);
  double mhsq   = HiggsMasses[posh];
  double mHsq   = HiggsMasses[posH];
  double mAsq   = HiggsMasses[posA];
  double mHpmsq = HiggsMasses[posHp];

  double mbarsq = RealMMix / (C_CosBeta * C_SinBeta);

  double tree_hhh = -3 / scale *
                    (sbma * sbma * sbma * mhsq +
                     sbma * cbma * cbma * (3 * mhsq - 2 * mbarsq) +
                     2 * cbma * cbma * cbma * cot2b * (mhsq - mbarsq));
  double tree_hhH = cbma / scale *
                    (sbma * sbma * (2 * mhsq + mHsq - 4 * mbarsq) +
                     2 * sbma * cbma * cot2b * (2 * mhsq + mHsq - 3 * mbarsq) -
                     cbma * cbma * (2 * mhsq + mHsq - 2 * mbarsq));
  double tree_hHH = -sbma / scale *
                    (sbma * sbma * (mhsq + 2 * mHsq - 2 * mbarsq) +
                     2 * sbma * cbma * cot2b * (mhsq + 2 * mHsq - 3 * mbarsq) -
                     cbma * cbma * (mhsq + 2 * mHsq - 4 * mbarsq));
  double tree_HHH = -3 / scale *
                    (cbma * cbma * cbma * mHsq +
                     sbma * sbma * cbma * (3 * mHsq - 2 * mbarsq) -
                     2 * sbma * sbma * sbma * cot2b * (mHsq - mbarsq));
  double tree_hAA = -1. / scale *
                    (sbma * (2 * mAsq + mhsq - 2 * mbarsq) +
                     2 * cbma * cot2b * (mhsq - mbarsq));
  double tree_HAA = -1. / scale *
                    (cbma * (2 * mAsq + mHsq - 2 * mbarsq) -
                     2 * sbma * cot2b * (mHsq - mbarsq));
  double tree_hHpHm = -1. / scale *
                      (sbma * (mhsq + 2 * mHpmsq - 2 * mbarsq) +
                       2 * cbma * cot2b * (mhsq - mbarsq));
  double tree_HHpHm = -1. / scale *
                      (cbma * (mHsq + 2 * mHpmsq - 2 * mbarsq) -
                       2 * sbma * cot2b * (mHsq - mbarsq));
  double tree_hAG0    = cbma / scale * (mAsq - mhsq);
  double tree_HAG0    = -sbma / scale * (mAsq - mHsq);
  double tree_hHpmGmp = -cbma / scale * (mhsq - mHpmsq);
  double tree_HHpmGmp = sbma / scale * (mHsq - mHpmsq);
  double tree_AHpGm   = (mHpmsq - mAsq) / scale;
  double tree_AHmGp   = -tree_AHpGm;
  double tree_hG0G0   = -sbma / scale * mhsq;
  double tree_HG0G0   = -cbma / scale * mHsq;
  double tree_hGpmGmp = -sbma / scale * mhsq;
  double tree_HGpmGmp = -cbma / scale * mHsq;

  // HiggsOrder: 0:Gp, 1:Gm, 2:Hp, 3:Hm, 4:G0, 5:A, 6:H, 7:h

  ss << "------ Checking Trilinear Tree Couplings ------\n";
  ss << "The mass spectrum is given by :\n"
     << "m_{G^+} = " << std::sqrt(HiggsMasses[posGp]) << " GeV \n"
     << "m_{G^-} = " << std::sqrt(HiggsMasses[posGm]) << " GeV \n"
     << "m_{H^+} = " << std::sqrt(HiggsMasses[posHp]) << " GeV \n"
     << "m_{H^-} = " << std::sqrt(HiggsMasses[posHm]) << " GeV \n"
     << "m_{G^0} = " << std::sqrt(HiggsMasses[posG0]) << " GeV \n"
     << "m_A = " << std::sqrt(HiggsMasses[posA]) << " GeV \n"
     << "m_H = " << std::sqrt(HiggsMasses[posH]) << " GeV \n"
     << "m_h = " << std::sqrt(HiggsMasses[posh]) << " GeV \n\n";
  ss << "coupling: true_value | BSMPT_value\n";
  int count = 0;

  // charge-breaking trilinears (check if the vanish)
  double GpGpG0 = TripleHiggsCorrectionsTreePhysical[0][0][4];
  double GmGmG0 = TripleHiggsCorrectionsTreePhysical[0][0][5];
  double HpHpG0 = TripleHiggsCorrectionsTreePhysical[0][0][6];
  double HmHmG0 = TripleHiggsCorrectionsTreePhysical[0][0][7];

  if (not almost_the_same(0, GpGpG0))
  {
    count += 1;
    ss << "GpGpG0: 0 | " << GpGpG0 << "\n";
  }
  if (not almost_the_same(0, GmGmG0))
  {
    count += 1;
    ss << "GmGmG0 : 0 | " << GmGmG0 << "\n";
  }
  if (not almost_the_same(0, HpHpG0))
  {
    count += 1;
    ss << "HpHpG0 : 0 | " << HpHpG0 << "\n";
  }
  if (not almost_the_same(0, HmHmG0))
  {
    count += 1;
    ss << "HmHmG0 : 0 | " << HmHmG0 << "\n";
  }

  double GpGpA = TripleHiggsCorrectionsTreePhysical[1][1][4];
  double GmGmA = TripleHiggsCorrectionsTreePhysical[1][1][5];
  double HpHpA = TripleHiggsCorrectionsTreePhysical[1][1][6];
  double HmHmA = TripleHiggsCorrectionsTreePhysical[1][1][7];

  if (not almost_the_same(0, GpGpA))
  {
    count += 1;
    ss << "GpGpA : 0 | " << GpGpA << "\n";
  }
  if (not almost_the_same(0, GmGmA))
  {
    count += 1;
    ss << "GmGmA : 0 | " << GmGmA << "\n";
  }
  if (not almost_the_same(0, HpHpA))
  {
    count += 1;
    ss << "HpHpA : 0 | " << HpHpA << "\n";
  }
  if (not almost_the_same(0, HmHmA))
  {
    count += 1;
    ss << "HmHmA : 0 | " << HmHmA << "\n";
  }

  double GpGpH = TripleHiggsCorrectionsTreePhysical[2][2][4];
  double GmGmH = TripleHiggsCorrectionsTreePhysical[2][2][5];
  double HpHpH = TripleHiggsCorrectionsTreePhysical[2][2][6];
  double HmHmH = TripleHiggsCorrectionsTreePhysical[2][2][7];

  if (not almost_the_same(0, GpGpH))
  {
    count += 1;
    ss << "GpGpH : 0 | " << GpGpH << "\n";
  }
  if (not almost_the_same(0, GmGmH))
  {
    count += 1;
    ss << "GmGmH : 0 | " << GmGmH << "\n";
  }
  if (not almost_the_same(0, HpHpH))
  {
    count += 1;
    ss << "HpHpH : 0 | " << HpHpH << "\n";
  }
  if (not almost_the_same(0, HmHmH))
  {
    count += 1;
    ss << "HmHmH : 0 | " << HmHmH << "\n";
  }

  double GpGph = TripleHiggsCorrectionsTreePhysical[3][3][4];
  double GmGmh = TripleHiggsCorrectionsTreePhysical[3][3][5];
  double HpHph = TripleHiggsCorrectionsTreePhysical[3][3][6];
  double HmHmh = TripleHiggsCorrectionsTreePhysical[3][3][7];

  if (not almost_the_same(0, GpGph))
  {
    count += 1;
    ss << "GpGph : 0 | " << GpGph << "\n";
  }
  if (not almost_the_same(0, GmGmh))
  {
    count += 1;
    ss << "GmGmh : 0 | " << GmGmh << "\n";
  }
  if (not almost_the_same(0, HpHph))
  {
    count += 1;
    ss << "HpHph : 0 | " << HpHph << "\n";
  }
  if (not almost_the_same(0, HmHmh))
  {
    count += 1;
    ss << "HmHmh : 0 | " << HmHmh << "\n";
  }

  // hhh
  if (not almost_the_same(tree_hhh,
                          -TripleHiggsCorrectionsTreePhysical[7][7][7]))
  {
    count += 1;
    ss << "hhh: " << tree_hhh << " | "
       << -TripleHiggsCorrectionsTreePhysical[7][7][7] << "\n";
  }
  // hhH
  if (not almost_the_same(tree_hhH,
                          -TripleHiggsCorrectionsTreePhysical[7][7][6]))
  {
    count += 1;
    ss << "hhH: " << tree_hhH << " | "
       << -TripleHiggsCorrectionsTreePhysical[7][7][6] << "\n";
    ;
  }
  // hHH
  if (not almost_the_same(tree_hHH,
                          -TripleHiggsCorrectionsTreePhysical[7][6][6]))
  {
    count += 1;
    ss << "hHH: " << tree_hHH << " | "
       << -TripleHiggsCorrectionsTreePhysical[7][6][6] << "\n";
  }
  // HHH
  if (not almost_the_same(tree_HHH,
                          -TripleHiggsCorrectionsTreePhysical[6][6][6]))
  {
    count += 1;
    ss << "HHH: " << tree_HHH << " | "
       << -TripleHiggsCorrectionsTreePhysical[6][6][6] << "\n";
  }
  // hAA
  if (not almost_the_same(tree_hAA,
                          -TripleHiggsCorrectionsTreePhysical[7][5][5]))
  {
    count += 1;
    ss << "hAA: " << tree_hAA << " | "
       << -TripleHiggsCorrectionsTreePhysical[7][5][5] << "\n";
  }
  // HAA
  if (not almost_the_same(tree_HAA,
                          -TripleHiggsCorrectionsTreePhysical[6][5][5]))
  {
    count += 1;
    ss << "HAA: " << tree_HAA << " | "
       << -TripleHiggsCorrectionsTreePhysical[6][5][5] << "\n";
  }
  // hHpHm
  if (not almost_the_same(tree_hHpHm,
                          -TripleHiggsCorrectionsTreePhysical[7][2][3]))
  {
    count += 1;
    ss << "hHpHm: " << tree_hHpHm << " | "
       << -TripleHiggsCorrectionsTreePhysical[7][2][3] << "\n";
  }
  // HHpHm
  if (not almost_the_same(tree_HHpHm,
                          -TripleHiggsCorrectionsTreePhysical[6][2][3]))
  {
    count += 1;
    ss << "HHpHm: " << tree_HHpHm << " | "
       << -TripleHiggsCorrectionsTreePhysical[6][2][3] << "\n";
  }
  // hAG0
  if (not almost_the_same(tree_hAG0,
                          -TripleHiggsCorrectionsTreePhysical[7][5][4]))
  {
    count += 1;
    ss << "hAG0: " << tree_hAG0 << " | "
       << -TripleHiggsCorrectionsTreePhysical[7][5][4] << "\n";
  }
  // HAG0
  if (not almost_the_same(tree_HAG0,
                          -TripleHiggsCorrectionsTreePhysical[6][5][4]))
  {
    count += 1;
    ss << "HAG0: " << tree_HAG0 << " | "
       << -TripleHiggsCorrectionsTreePhysical[6][5][4] << "\n";
  }
  // hHpmGmp
  if (not almost_the_same(tree_hHpmGmp,
                          -TripleHiggsCorrectionsTreePhysical[7][2][1]) or
      not almost_the_same(tree_hHpmGmp,
                          -TripleHiggsCorrectionsTreePhysical[7][3][0]))
  {
    count += 1;
    ss << "hHpmGmp: " << tree_hHpmGmp << " | "
       << -TripleHiggsCorrectionsTreePhysical[7][2][1] << ","
       << -TripleHiggsCorrectionsTreePhysical[7][3][0] << "\n";
  }
  // HHpmGmp
  if (not almost_the_same(tree_HHpmGmp,
                          -TripleHiggsCorrectionsTreePhysical[6][2][1]) or
      not almost_the_same(tree_HHpmGmp,
                          -TripleHiggsCorrectionsTreePhysical[6][3][0]))
  {
    count += 1;
    ss << "HHpmGmp: " << tree_HHpmGmp << " | "
       << -TripleHiggsCorrectionsTreePhysical[6][2][1] << ","
       << -TripleHiggsCorrectionsTreePhysical[6][3][0] << "\n";
  }
  // AHpGm
  if (not almost_the_same(tree_AHpGm,
                          -TripleHiggsCorrectionsTreePhysical[5][2][1]))
  {
    count += 1;
    ss << "AHpGm: " << tree_AHpGm << " | "
       << -TripleHiggsCorrectionsTreePhysical[5][2][1] << "\n";
  }
  // AHmGp
  if (not almost_the_same(tree_AHmGp,
                          -TripleHiggsCorrectionsTreePhysical[5][3][0]))
  {
    count += 1;
    ss << "AHpGm: " << tree_AHmGp << " | "
       << -TripleHiggsCorrectionsTreePhysical[5][3][0] << "\n";
  }
  // hG0G0
  if (not almost_the_same(tree_hG0G0,
                          -TripleHiggsCorrectionsTreePhysical[7][4][4]))
  {
    count += 1;
    ss << "hG0G0: " << tree_hG0G0 << " | "
       << -TripleHiggsCorrectionsTreePhysical[7][4][4] << "\n";
  }
  // HG0G0
  if (not almost_the_same(tree_HG0G0,
                          -TripleHiggsCorrectionsTreePhysical[6][4][4]))
  {
    count += 1;
    ss << "HG0G0: " << tree_HG0G0 << " | "
       << -TripleHiggsCorrectionsTreePhysical[6][4][4] << "\n";
  }
  // hGpmGmp
  if (not almost_the_same(tree_hGpmGmp,
                          -TripleHiggsCorrectionsTreePhysical[7][0][1]) or
      not almost_the_same(tree_hGpmGmp,
                          -TripleHiggsCorrectionsTreePhysical[7][1][0]))
  {
    count += 1;
    ss << "hGpmGmp: " << tree_hGpmGmp << " | "
       << -TripleHiggsCorrectionsTreePhysical[7][0][1] << ","
       << -TripleHiggsCorrectionsTreePhysical[7][1][0] << "\n";
  }
  // HGpmGmp
  if (not almost_the_same(tree_HGpmGmp,
                          -TripleHiggsCorrectionsTreePhysical[6][0][1]) or
      not almost_the_same(tree_HGpmGmp,
                          -TripleHiggsCorrectionsTreePhysical[6][1][0]))
  {
    count += 1;
    ss << "HGpmGmp: " << tree_HGpmGmp << " | "
       << -TripleHiggsCorrectionsTreePhysical[6][0][1] << ","
       << -TripleHiggsCorrectionsTreePhysical[6][1][0] << "\n";
  }

  if (count == 0)
  {
    ss << "All trilinears agree. Success!\n";
  }
  else
  {
    ss << count << " trilinears differ. Error!\n";
  }

  Logger::Write(LoggingLevel::Default, ss.str());
}

void Class_Potential_R2HDM::SetCurvatureArrays()
{
  initVectors();

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    Curvature_Higgs_L1[i] = 0;
    HiggsVev[i]           = 0;
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
        Curvature_Higgs_L3[i][j][k] = 0;
    }
  }

  HiggsVev[4] = SMConstants.C_vev0 * C_CosBeta;
  HiggsVev[6] = SMConstants.C_vev0 * C_SinBeta;

  Curvature_Higgs_L2[0][0] = u1;
  Curvature_Higgs_L2[0][1] = 0;
  Curvature_Higgs_L2[0][2] = -RealMMix;
  Curvature_Higgs_L2[0][3] = 0;
  Curvature_Higgs_L2[0][4] = 0;
  Curvature_Higgs_L2[0][5] = 0;
  Curvature_Higgs_L2[0][6] = 0;
  Curvature_Higgs_L2[0][7] = 0;
  Curvature_Higgs_L2[1][0] = 0;
  Curvature_Higgs_L2[1][1] = u1;
  Curvature_Higgs_L2[1][2] = 0;
  Curvature_Higgs_L2[1][3] = -RealMMix;
  Curvature_Higgs_L2[1][4] = 0;
  Curvature_Higgs_L2[1][5] = 0;
  Curvature_Higgs_L2[1][6] = 0;
  Curvature_Higgs_L2[1][7] = 0;
  Curvature_Higgs_L2[2][0] = -RealMMix;
  Curvature_Higgs_L2[2][1] = 0;
  Curvature_Higgs_L2[2][2] = u2;
  Curvature_Higgs_L2[2][3] = 0;
  Curvature_Higgs_L2[2][4] = 0;
  Curvature_Higgs_L2[2][5] = 0;
  Curvature_Higgs_L2[2][6] = 0;
  Curvature_Higgs_L2[2][7] = 0;
  Curvature_Higgs_L2[3][0] = 0;
  Curvature_Higgs_L2[3][1] = -RealMMix;
  Curvature_Higgs_L2[3][2] = 0;
  Curvature_Higgs_L2[3][3] = u2;
  Curvature_Higgs_L2[3][4] = 0;
  Curvature_Higgs_L2[3][5] = 0;
  Curvature_Higgs_L2[3][6] = 0;
  Curvature_Higgs_L2[3][7] = 0;
  Curvature_Higgs_L2[4][0] = 0;
  Curvature_Higgs_L2[4][1] = 0;
  Curvature_Higgs_L2[4][2] = 0;
  Curvature_Higgs_L2[4][3] = 0;
  Curvature_Higgs_L2[4][4] = u1;
  Curvature_Higgs_L2[4][5] = 0;
  Curvature_Higgs_L2[4][6] = -RealMMix;
  Curvature_Higgs_L2[4][7] = 0;
  Curvature_Higgs_L2[5][0] = 0;
  Curvature_Higgs_L2[5][1] = 0;
  Curvature_Higgs_L2[5][2] = 0;
  Curvature_Higgs_L2[5][3] = 0;
  Curvature_Higgs_L2[5][4] = 0;
  Curvature_Higgs_L2[5][5] = u1;
  Curvature_Higgs_L2[5][6] = 0;
  Curvature_Higgs_L2[5][7] = -RealMMix;
  Curvature_Higgs_L2[6][0] = 0;
  Curvature_Higgs_L2[6][1] = 0;
  Curvature_Higgs_L2[6][2] = 0;
  Curvature_Higgs_L2[6][3] = 0;
  Curvature_Higgs_L2[6][4] = -RealMMix;
  Curvature_Higgs_L2[6][5] = 0;
  Curvature_Higgs_L2[6][6] = u2;
  Curvature_Higgs_L2[6][7] = 0;
  Curvature_Higgs_L2[7][0] = 0;
  Curvature_Higgs_L2[7][1] = 0;
  Curvature_Higgs_L2[7][2] = 0;
  Curvature_Higgs_L2[7][3] = 0;
  Curvature_Higgs_L2[7][4] = 0;
  Curvature_Higgs_L2[7][5] = -RealMMix;
  Curvature_Higgs_L2[7][6] = 0;
  Curvature_Higgs_L2[7][7] = u2;

  Curvature_Higgs_L4[0][0][0][0] = 3 * L1;
  Curvature_Higgs_L4[0][0][1][1] = L1;
  Curvature_Higgs_L4[0][0][2][2] = L3 + L4 + RL5;
  Curvature_Higgs_L4[0][0][2][3] = 0;
  Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - RL5;
  Curvature_Higgs_L4[0][0][4][4] = L1;
  Curvature_Higgs_L4[0][0][5][5] = L1;
  Curvature_Higgs_L4[0][0][6][6] = L3;
  Curvature_Higgs_L4[0][0][7][7] = L3;
  Curvature_Higgs_L4[0][1][2][2] = 0;
  Curvature_Higgs_L4[0][1][2][3] = RL5;
  Curvature_Higgs_L4[0][1][3][3] = 0;
  Curvature_Higgs_L4[0][2][4][6] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[0][2][4][7] = 0;
  Curvature_Higgs_L4[0][2][5][6] = 0;
  Curvature_Higgs_L4[0][2][5][7] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][4][6] = 0;
  Curvature_Higgs_L4[0][3][4][7] = L4 / 0.2e1 - RL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][6] = -L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][7] = 0;
  Curvature_Higgs_L4[1][1][1][1] = 3 * L1;
  Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - RL5;
  Curvature_Higgs_L4[1][1][2][3] = 0;
  Curvature_Higgs_L4[1][1][3][3] = L3 + L4 + RL5;
  Curvature_Higgs_L4[1][1][4][4] = L1;
  Curvature_Higgs_L4[1][1][5][5] = L1;
  Curvature_Higgs_L4[1][1][6][6] = L3;
  Curvature_Higgs_L4[1][1][7][7] = L3;
  Curvature_Higgs_L4[1][2][4][6] = 0;
  Curvature_Higgs_L4[1][2][4][7] = -L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][6] = L4 / 0.2e1 - RL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][7] = 0;
  Curvature_Higgs_L4[1][3][4][6] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[1][3][4][7] = 0;
  Curvature_Higgs_L4[1][3][5][6] = 0;
  Curvature_Higgs_L4[1][3][5][7] = L4 / 0.2e1 + RL5 / 0.2e1;
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
  Curvature_Higgs_L4[4][4][6][6] = L3 + L4 + RL5;
  Curvature_Higgs_L4[4][4][6][7] = 0;
  Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - RL5;
  Curvature_Higgs_L4[4][5][6][6] = 0;
  Curvature_Higgs_L4[4][5][6][7] = RL5;
  Curvature_Higgs_L4[4][5][7][7] = 0;
  Curvature_Higgs_L4[5][5][5][5] = 3 * L1;
  Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - RL5;
  Curvature_Higgs_L4[5][5][6][7] = 0;
  Curvature_Higgs_L4[5][5][7][7] = L3 + L4 + RL5;
  Curvature_Higgs_L4[6][6][6][6] = 3 * L2;
  Curvature_Higgs_L4[6][6][7][7] = L2;
  Curvature_Higgs_L4[7][7][7][7] = 3 * L2;

  for (std::size_t k1 = 0; k1 < NHiggs; k1++)
  {
    for (std::size_t k2 = k1; k2 < NHiggs; k2++)
    {
      for (std::size_t k3 = k2; k3 < NHiggs; k3++)
      {
        for (std::size_t k4 = k3; k4 < NHiggs; k4++)
        {
          Curvature_Higgs_L4[k2][k3][k4][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k4][k1][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k1][k2][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k1][k3][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k2][k1][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k4][k2][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k3][k4][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k2][k1][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k3][k2][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k4][k3][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k1][k4][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k2][k3][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k4][k2][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k1][k4][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k3][k1][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k3][k2][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k1][k3][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k4][k1][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k2][k4][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k2][k4][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k1][k2][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k3][k1][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k4][k3][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
        }
      }
    }
  }

  for (std::size_t a = 0; a < NGauge; a++)
  {
    for (std::size_t b = 0; b < NGauge; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        for (std::size_t j = 0; j < NHiggs; j++)
          Curvature_Gauge_G2H2[a][b][i][j] = 0;
      }
    }
  }

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

bool Class_Potential_R2HDM::CalculateDebyeSimplified()
{
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

  DebyeHiggs[0][0] = CTempC1;
  DebyeHiggs[1][1] = CTempC1;
  DebyeHiggs[2][2] = CTempC2;
  DebyeHiggs[3][3] = CTempC2;
  DebyeHiggs[4][4] = CTempC1;
  DebyeHiggs[5][5] = CTempC1;
  DebyeHiggs[6][6] = CTempC2;
  DebyeHiggs[7][7] = CTempC2;

  return true;
}

bool Class_Potential_R2HDM::CalculateDebyeGaugeSimplified()
{
  DebyeGauge[0][0] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[1][1] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[2][2] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[3][3] = 2 * SMConstants.C_gs * SMConstants.C_gs;

  return true;
}

double
Class_Potential_R2HDM::VTreeSimplified(const std::vector<double> &v) const
{
  (void)v;
  double res = 0;

  return res;
}

double
Class_Potential_R2HDM::VCounterSimplified(const std::vector<double> &v) const
{
  (void)v;
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  return res;
}

void Class_Potential_R2HDM::Debugging(const std::vector<double> &input,
                                      std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
