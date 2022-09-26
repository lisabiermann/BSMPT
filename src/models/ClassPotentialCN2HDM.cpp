// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialCN2HDM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
namespace Models
{
Class_Potential_CN2HDM::Class_Potential_CN2HDM()
{
  Model = ModelID::ModelIDs::CN2HDM;

  NNeutralHiggs = 6; // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 4; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 15; // number of parameters in the tree-Level Lagrangian
  nParCT = 27; // number of parameters in the counterterm potential

  nVEV = 6; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  // defines which scalar field gets which VEV
  VevOrder[0] = 2; // wCB
  VevOrder[1] = 4; // w1
  VevOrder[2] = 6; // w2
  VevOrder[3] = 7; // wCP
  VevOrder[4] = 8; // wS
  VevOrder[5] = 9; // wDM

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_CN2HDM::~Class_Potential_CN2HDM()
{
  // TODO Auto-generated destructor stub
}

/**
 * string with chronological order of the counterterms.
 */
std::vector<std::string> Class_Potential_CN2HDM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dm11Sq");
  labels.push_back("dm22Sq");
  labels.push_back("dRem12Sq");
  labels.push_back("dImm12Sq");
  labels.push_back("dmsSq");
  labels.push_back("dRebs");
  labels.push_back("dImbs");
  labels.push_back("dL1");
  labels.push_back("dL2");
  labels.push_back("dL3");
  labels.push_back("dL4");
  labels.push_back("dReL5");
  labels.push_back("dImL5");
  labels.push_back("dImL7");
  labels.push_back("dL8");
  labels.push_back("dL9");
  labels.push_back("dImL10");
  labels.push_back("dImL11");
  labels.push_back("dImL12");
  labels.push_back("dReL13");
  labels.push_back("dLS");
  labels.push_back("dTCB");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dTCP");
  labels.push_back("dTS");
  labels.push_back("dTDM");
  return labels;
}

/**
 * string with chronological order of the VEVs and critical temperature
 */
std::vector<std::string> Class_Potential_CN2HDM::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c"); // critical temperature
  labels.push_back("v_c"); // critical vev
  labels.push_back("xi_c");
  labels.push_back("wCB");
  labels.push_back("w1");
  labels.push_back("w2");
  labels.push_back("wCP");
  labels.push_back("wS");
  labels.push_back("wDM");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 *
 */
std::vector<std::string>
Class_Potential_CN2HDM::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);
  particles[0] = "G^+";
  particles[1] = "G^-";
  particles[2] = "H^+";
  particles[3] = "H^-";
  particles[4] = "G^0";
  particles[5] = "h1";
  particles[6] = "h2";
  particles[7] = "h3";
  particles[8] = "h4";
  particles[9] = "hDM";

  for (size_t i = 0; i < NHiggs; i++)
  {
    for (size_t j = i; j < NHiggs; j++)
    {
      for (size_t k = j; k < NHiggs; k++)
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
 * string with chronological order of the VEVs
 */
std::vector<std::string> Class_Potential_CN2HDM::addLegendVEV() const
{
  std::vector<std::string> labels;
  labels.push_back("wCB");
  labels.push_back("w1");
  labels.push_back("w2");
  labels.push_back("wCP");
  labels.push_back("wS");
  labels.push_back("wDM");
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_CN2HDM::ReadAndSet(const std::string &linestr,
                                        std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  double lType = 0, lL1 = 0, lL2 = 0, lL3 = 0, lL4 = 0, lReL5 = 0, lImL5 = 0,
         lLS = 0, lL8 = 0, lL9 = 0, lRem12Sq = 0, lmRebs = 0, lTanBeta = 0,
         lvs = 0;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 14; k++)
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
      lReL5 = tmp;
    else if (k == 7)
      lImL5 = tmp;
    else if (k == 8)
      lLS = tmp;
    else if (k == 9)
      lL8 = tmp;
    else if (k == 10)
      lL9 = tmp;
    else if (k == 11)
      lRem12Sq = tmp;
    else if (k == 12)
      lmRebs = tmp;
    else if (k == 13)
      lTanBeta = tmp;
    else if (k == 14)
      lvs = tmp;
  }
  par[0]  = lType;
  par[1]  = lL1;
  par[2]  = lL2;
  par[3]  = lL3;
  par[4]  = lL4;
  par[5]  = lReL5;
  par[6]  = lImL5;
  par[7]  = lLS;
  par[8]  = lL8;
  par[9]  = lL9;
  par[10] = lRem12Sq;
  par[11] = -lmRebs; // account for different prefactor convention in input
  par[12] = lTanBeta;
  par[13] = lvs;

  set_gen(par); // set model parameters
  return;
}

void Class_Potential_CN2HDM::set_gen(const std::vector<double> &par)
{

  Type    = par[0];
  vs      = par[13];
  TanBeta = par[12];
  Rem12Sq = par[10];
  Rebs    = par[11];
  L1      = par[1];
  L2      = par[2];
  L3      = par[3];
  L4      = par[4];
  ReL5    = par[5];
  ImL5    = par[6];
  L8      = par[8];
  L9      = par[9];
  LS      = par[7];

  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = std::sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = std::sqrt(C_SinBetaSquared);

  v1 = C_vev0 * C_CosBeta;
  v2 = C_vev0 * C_SinBeta;

  // define parameters fixed by tadpole conditions
  m11Sq = -v1 * v1 * L1 / 2.0 - v2 * v2 / 2.0 * (L3 + L4 + ReL5) -
          vs * vs / 4.0 * L8 + Rem12Sq * v2 / v1;
  m22Sq = -v2 * v2 * L2 / 2.0 - v1 * v1 / 2.0 * (L3 + L4 + ReL5) -
          vs * vs / 4.0 * L9 + Rem12Sq * v1 / v2;
  msSq = -v1 * v1 / 2.0 * L8 - v2 * v2 / 2.0 * L9 - vs * vs / 4.0 * LS - Rebs;
  Imm12Sq = ImL5 * v1 * v2 / 2.0;

  // tree-level minimum
  vevTreeMin.resize(nVEV);
  vevTreeMin[0] = 0;
  vevTreeMin[1] = v1;
  vevTreeMin[2] = v2;
  vevTreeMin[3] = 0;
  vevTreeMin[4] = vs;
  vevTreeMin[5] = 0;
  vevTree.resize(NHiggs);
  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * counterterm parameters and curvature tensors
 */
void Class_Potential_CN2HDM::set_CT_Pot_Par(const std::vector<double> &par)
{
  dm11Sq   = par[0];
  dm22Sq   = par[1];
  dRem12Sq = par[2];
  dImm12Sq = par[3];
  dmsSq    = par[4];
  dRebs    = par[5];
  dImbs    = par[6];
  dL1      = par[7];
  dL2      = par[8];
  dL3      = par[9];
  dL4      = par[10];
  dReL5    = par[11];
  dImL5    = par[12];
  dImL7    = par[13];
  dL8      = par[14];
  dL9      = par[15];
  dImL10   = par[16];
  dImL11   = par[17];
  dImL12   = par[18];
  dReL13   = par[19];
  dLS      = par[20];
  dTCB     = par[21];
  dT1      = par[22];
  dT2      = par[23];
  dTCP     = par[24];
  dTS      = par[25];
  dTDM     = par[26];

  Curvature_Higgs_CT_L1[2] = dTCB;
  Curvature_Higgs_CT_L1[4] = dT1;
  Curvature_Higgs_CT_L1[6] = dT2;
  Curvature_Higgs_CT_L1[7] = dTCP;
  Curvature_Higgs_CT_L1[8] = dTS;
  Curvature_Higgs_CT_L1[9] = dTDM;

  Curvature_Higgs_CT_L2[0][0] = dm11Sq;
  Curvature_Higgs_CT_L2[0][2] = -dRem12Sq;
  Curvature_Higgs_CT_L2[0][3] = dImm12Sq;
  Curvature_Higgs_CT_L2[1][1] = dm11Sq;
  Curvature_Higgs_CT_L2[1][2] = -dImm12Sq;
  Curvature_Higgs_CT_L2[1][3] = -dRem12Sq;
  Curvature_Higgs_CT_L2[2][2] = dm22Sq;
  Curvature_Higgs_CT_L2[3][3] = dm22Sq;
  Curvature_Higgs_CT_L2[4][4] = dm11Sq;
  Curvature_Higgs_CT_L2[4][6] = -dRem12Sq;
  Curvature_Higgs_CT_L2[4][7] = dImm12Sq;
  Curvature_Higgs_CT_L2[5][5] = dm11Sq;
  Curvature_Higgs_CT_L2[5][6] = -dImm12Sq;
  Curvature_Higgs_CT_L2[5][7] = -dRem12Sq;
  Curvature_Higgs_CT_L2[6][6] = dm22Sq;
  Curvature_Higgs_CT_L2[7][7] = dm22Sq;
  Curvature_Higgs_CT_L2[8][8] = dmsSq / 0.2e1 + dRebs / 0.2e1;
  Curvature_Higgs_CT_L2[8][9] = -dImbs / 0.2e1;
  Curvature_Higgs_CT_L2[9][9] = dmsSq / 0.2e1 - dRebs / 0.2e1;

  sym2Dim(Curvature_Higgs_CT_L2, NHiggs, NHiggs);

  Curvature_Higgs_CT_L4[0][0][0][0] = 3 * dL1;
  Curvature_Higgs_CT_L4[0][0][1][1] = dL1;
  Curvature_Higgs_CT_L4[0][0][2][2] = dL3 + dL4 + dReL5;
  Curvature_Higgs_CT_L4[0][0][2][3] = -dImL5;
  Curvature_Higgs_CT_L4[0][0][3][3] = dL3 + dL4 - dReL5;
  Curvature_Higgs_CT_L4[0][0][4][4] = dL1;
  Curvature_Higgs_CT_L4[0][0][5][5] = dL1;
  Curvature_Higgs_CT_L4[0][0][6][6] = dL3;
  Curvature_Higgs_CT_L4[0][0][7][7] = dL3;
  Curvature_Higgs_CT_L4[0][0][8][8] = dL8 / 0.2e1;
  Curvature_Higgs_CT_L4[0][0][8][9] = -dImL11 / 0.2e1;
  Curvature_Higgs_CT_L4[0][0][9][9] = dL8 / 0.2e1;
  Curvature_Higgs_CT_L4[0][1][2][2] = dImL5;
  Curvature_Higgs_CT_L4[0][1][2][3] = dReL5;
  Curvature_Higgs_CT_L4[0][1][3][3] = -dImL5;
  Curvature_Higgs_CT_L4[0][2][2][3] = -dImL7;
  Curvature_Higgs_CT_L4[0][2][4][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][4][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][5][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][5][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][8][8] = dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][9][9] = -dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][3][3] = -3 * dImL7;
  Curvature_Higgs_CT_L4[0][3][4][6] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][4][7] = dL4 / 0.2e1 - dReL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][5][6] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][5][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][6][6] = -dImL7;
  Curvature_Higgs_CT_L4[0][3][7][7] = -dImL7;
  Curvature_Higgs_CT_L4[0][3][8][8] = -dImL10;
  Curvature_Higgs_CT_L4[0][3][8][9] = -dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][9][9] = -dImL10;
  Curvature_Higgs_CT_L4[1][1][1][1] = 3 * dL1;
  Curvature_Higgs_CT_L4[1][1][2][2] = dL3 + dL4 - dReL5;
  Curvature_Higgs_CT_L4[1][1][2][3] = dImL5;
  Curvature_Higgs_CT_L4[1][1][3][3] = dL3 + dL4 + dReL5;
  Curvature_Higgs_CT_L4[1][1][4][4] = dL1;
  Curvature_Higgs_CT_L4[1][1][5][5] = dL1;
  Curvature_Higgs_CT_L4[1][1][6][6] = dL3;
  Curvature_Higgs_CT_L4[1][1][7][7] = dL3;
  Curvature_Higgs_CT_L4[1][1][8][8] = dL8 / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][8][9] = -dImL11 / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][9][9] = dL8 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][2][2] = 3 * dImL7;
  Curvature_Higgs_CT_L4[1][2][3][3] = dImL7;
  Curvature_Higgs_CT_L4[1][2][4][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][4][7] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][5][6] = dL4 / 0.2e1 - dReL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][5][7] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][6][6] = dImL7;
  Curvature_Higgs_CT_L4[1][2][7][7] = dImL7;
  Curvature_Higgs_CT_L4[1][2][8][8] = dImL10;
  Curvature_Higgs_CT_L4[1][2][8][9] = dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][9][9] = dImL10;
  Curvature_Higgs_CT_L4[1][3][4][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][4][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][5][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][5][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][8][8] = dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][9][9] = -dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][2][2] = 3 * dL2;
  Curvature_Higgs_CT_L4[2][2][3][3] = dL2;
  Curvature_Higgs_CT_L4[2][2][4][4] = dL3;
  Curvature_Higgs_CT_L4[2][2][4][7] = -dImL7;
  Curvature_Higgs_CT_L4[2][2][5][5] = dL3;
  Curvature_Higgs_CT_L4[2][2][5][6] = dImL7;
  Curvature_Higgs_CT_L4[2][2][6][6] = dL2;
  Curvature_Higgs_CT_L4[2][2][7][7] = dL2;
  Curvature_Higgs_CT_L4[2][2][8][8] = dL9 / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][8][9] = -dImL12 / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][9][9] = dL9 / 0.2e1;
  Curvature_Higgs_CT_L4[3][3][3][3] = 3 * dL2;
  Curvature_Higgs_CT_L4[3][3][4][4] = dL3;
  Curvature_Higgs_CT_L4[3][3][4][7] = -dImL7;
  Curvature_Higgs_CT_L4[3][3][5][5] = dL3;
  Curvature_Higgs_CT_L4[3][3][5][6] = dImL7;
  Curvature_Higgs_CT_L4[3][3][6][6] = dL2;
  Curvature_Higgs_CT_L4[3][3][7][7] = dL2;
  Curvature_Higgs_CT_L4[3][3][8][8] = dL9 / 0.2e1;
  Curvature_Higgs_CT_L4[3][3][8][9] = -dImL12 / 0.2e1;
  Curvature_Higgs_CT_L4[3][3][9][9] = dL9 / 0.2e1;
  Curvature_Higgs_CT_L4[4][4][4][4] = 3 * dL1;
  Curvature_Higgs_CT_L4[4][4][5][5] = dL1;
  Curvature_Higgs_CT_L4[4][4][6][6] = dL3 + dL4 + dReL5;
  Curvature_Higgs_CT_L4[4][4][6][7] = -dImL5;
  Curvature_Higgs_CT_L4[4][4][7][7] = dL3 + dL4 - dReL5;
  Curvature_Higgs_CT_L4[4][4][8][8] = dL8 / 0.2e1;
  Curvature_Higgs_CT_L4[4][4][8][9] = -dImL11 / 0.2e1;
  Curvature_Higgs_CT_L4[4][4][9][9] = dL8 / 0.2e1;
  Curvature_Higgs_CT_L4[4][5][6][6] = dImL5;
  Curvature_Higgs_CT_L4[4][5][6][7] = dReL5;
  Curvature_Higgs_CT_L4[4][5][7][7] = -dImL5;
  Curvature_Higgs_CT_L4[4][6][6][7] = -dImL7;
  Curvature_Higgs_CT_L4[4][6][8][8] = dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[4][6][9][9] = -dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][7][7] = -3 * dImL7;
  Curvature_Higgs_CT_L4[4][7][8][8] = -dImL10;
  Curvature_Higgs_CT_L4[4][7][8][9] = -dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][9][9] = -dImL10;
  Curvature_Higgs_CT_L4[5][5][5][5] = 3 * dL1;
  Curvature_Higgs_CT_L4[5][5][6][6] = dL3 + dL4 - dReL5;
  Curvature_Higgs_CT_L4[5][5][6][7] = dImL5;
  Curvature_Higgs_CT_L4[5][5][7][7] = dL3 + dL4 + dReL5;
  Curvature_Higgs_CT_L4[5][5][8][8] = dL8 / 0.2e1;
  Curvature_Higgs_CT_L4[5][5][8][9] = -dImL11 / 0.2e1;
  Curvature_Higgs_CT_L4[5][5][9][9] = dL8 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][6][6] = 3 * dImL7;
  Curvature_Higgs_CT_L4[5][6][7][7] = dImL7;
  Curvature_Higgs_CT_L4[5][6][8][8] = dImL10;
  Curvature_Higgs_CT_L4[5][6][8][9] = dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][9][9] = dImL10;
  Curvature_Higgs_CT_L4[5][7][8][8] = dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[5][7][9][9] = -dReL13 / 0.2e1;
  Curvature_Higgs_CT_L4[6][6][6][6] = 3 * dL2;
  Curvature_Higgs_CT_L4[6][6][7][7] = dL2;
  Curvature_Higgs_CT_L4[6][6][8][8] = dL9 / 0.2e1;
  Curvature_Higgs_CT_L4[6][6][8][9] = -dImL12 / 0.2e1;
  Curvature_Higgs_CT_L4[6][6][9][9] = dL9 / 0.2e1;
  Curvature_Higgs_CT_L4[7][7][7][7] = 3 * dL2;
  Curvature_Higgs_CT_L4[7][7][8][8] = dL9 / 0.2e1;
  Curvature_Higgs_CT_L4[7][7][8][9] = -dImL12 / 0.2e1;
  Curvature_Higgs_CT_L4[7][7][9][9] = dL9 / 0.2e1;
  Curvature_Higgs_CT_L4[8][8][8][8] = 0.3e1 / 0.4e1 * dLS;
  Curvature_Higgs_CT_L4[8][8][9][9] = dLS / 0.4e1;
  Curvature_Higgs_CT_L4[9][9][9][9] = 0.3e1 / 0.4e1 * dLS;

  sym4Dim(Curvature_Higgs_CT_L4, NHiggs, NHiggs, NHiggs, NHiggs);
}

/**
 * console output of all parameters
 */
void Class_Potential_CN2HDM::write() const
{
  std::stringstream ss;
  ss.precision(std::numeric_limits<double>::max_digits10);

  ss << "scale = " << scale << "\n";
  ss << "The parameters are :  \n";
  ss << "Model = " << Model << "\n";
  ss << "Type: " << Type << "\n"
     << "m11Sq = " << m11Sq << "\n"
     << "m22Sq = " << m22Sq << "\n"
     << "Rem12Sq = " << Rem12Sq << "\n"
     << "Imm12Sq = " << Imm12Sq << "\n"
     << "msSq = " << msSq << "\n"
     << "Rebs = " << Rebs << "\n"
     << "L1 = " << L1 << "\n"
     << "L2 = " << L2 << "\n"
     << "L3 = " << L3 << "\n"
     << "L4 = " << L4 << "\n"
     << "ReL5 = " << ReL5 << "\n"
     << "ImL5 = " << ImL5 << "\n"
     << "L8 = " << L8 << "\n"
     << "L9 = " << L9 << "\n"
     << "LS = " << LS << "\n";

  ss << "The counterterms are :\n";
  ss << "dm11Sq = " << dm11Sq << "\n"
     << "dm22Sq = " << dm22Sq << "\n"
     << "dRem12Sq = " << dRem12Sq << "\n"
     << "dImm12Sq = " << dImm12Sq << "\n"
     << "dmsSq = " << dmsSq << "\n"
     << "dRebs = " << dRebs << "\n"
     << "dImbs = " << dImbs << "\n"
     << "dL1 = " << dL1 << "\n"
     << "dL2 = " << dL2 << "\n"
     << "dL3 = " << dL3 << "\n"
     << "dL4 = " << dL4 << "\n"
     << "dReL5 = " << dReL5 << "\n"
     << "dImL5 = " << dImL5 << "\n"
     << "dImL7 = " << dImL7 << "\n"
     << "dL8 = " << dL8 << "\n"
     << "dL9 = " << dL9 << "\n"
     << "dImL10 = " << dImL10 << "\n"
     << "dImL11 = " << dImL11 << "\n"
     << "dImL12 = " << dImL12 << "\n"
     << "dReL13 = " << dReL13 << "\n"
     << "dLS = " << dLS << "\n"
     << "dTCB = " << dTCB << "\n"
     << "dT1 = " << dT1 << "\n"
     << "dT2 = " << dT2 << "\n"
     << "dTCP = " << dTCP << "\n"
     << "dTS = " << dTDM << "\n";

  ss << "The scale is given by mu = " << C_vev0 << " GeV\n";
  ss << "The tree-level VEVs are given by: \n";
  ss << "TanBeta = " << TanBeta << "\n";
  ss << "v1 		= " << v1 << " GeV\n";
  ss << "v2 		= " << v2 << " GeV\n";
  ss << "vs 		= " << vs << " GeV\n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * calculation of finite counterterms
 */
std::vector<double> Class_Potential_CN2HDM::calc_CT() const
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
  for (size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  // chosen counterterm scheme sets all t-parameters to zero
  double t1 = 0, t2 = 0;
  parCT.push_back(
      (2 * t1 * v1 * v1 * v2 * v2 - 4 * HesseWeinberg(3, 3) * v2 * v2 +
       HesseWeinberg(4, 4) * v1 * v1 + HesseWeinberg(4, 6) * v1 * v2 +
       HesseWeinberg(4, 8) * v1 * vs - 3 * v1 * v1 * HesseWeinberg(5, 5) -
       HesseWeinberg(5, 7) * v1 * v2 + vs * v2 * HesseWeinberg(7, 9) +
       4 * HesseWeinberg(7, 7) * v2 * v2) *
      std::pow(v1, -2) / 2); // dm11Sq
  parCT.push_back((2 * t1 * v1 * v1 * v2 - 4 * HesseWeinberg(3, 3) * v2 +
                   HesseWeinberg(4, 6) * v1 - HesseWeinberg(5, 7) * v1 +
                   vs * HesseWeinberg(7, 9) + HesseWeinberg(7, 7) * v2 +
                   HesseWeinberg(6, 6) * v2 + HesseWeinberg(6, 8) * vs) /
                  v2 / 2); // dm22Sq
  parCT.push_back((2 * t1 * v1 * v1 * v2 - 4 * HesseWeinberg(3, 3) * v2 +
                   2 * HesseWeinberg(5, 7) * v1 + vs * HesseWeinberg(7, 9) +
                   4 * HesseWeinberg(7, 7) * v2) /
                  v1 / 2); // dRem12Sq
  parCT.push_back((HesseWeinberg(4, 7) * v1 + 3 * HesseWeinberg(5, 6) * v1 +
                   4 * HesseWeinberg(6, 7) * v2 + vs * HesseWeinberg(7, 8)) /
                  v1 / 2); // dImm12Sq
  parCT.push_back((2 * HesseWeinberg(4, 8) * v1 + 2 * HesseWeinberg(6, 8) * v2 +
                   4 * HesseWeinberg(7, 9) * v2 -
                   vs * (t2 * vs * vs + 2 * HesseWeinberg(8, 8) +
                         2 * HesseWeinberg(9, 9))) /
                  vs / 2); // dmsSq
  parCT.push_back((-8 * HesseWeinberg(7, 9) * v2 -
                   vs * (t2 * vs * vs + 4 * HesseWeinberg(8, 8) -
                         4 * HesseWeinberg(9, 9))) /
                  vs / 4); // dRebs
  parCT.push_back((-HesseWeinberg(4, 9) * v1 - HesseWeinberg(6, 9) * v2 +
                   2 * HesseWeinberg(8, 9) * vs) /
                  vs); // dImbs
  parCT.push_back((2 * HesseWeinberg(3, 3) * v2 * v2 -
                   HesseWeinberg(4, 4) * v1 * v1 +
                   v1 * v1 * HesseWeinberg(5, 5) -
                   v2 * v2 * (t1 * v1 * v1 + 2 * HesseWeinberg(7, 7))) *
                  std::pow(v1, -4)); // dL1
  parCT.push_back((-t1 * v1 * v1 + 2 * HesseWeinberg(3, 3) -
                   HesseWeinberg(6, 6) - HesseWeinberg(7, 7)) *
                  std::pow(v2, -2)); // dL2
  parCT.push_back((-t1 * v1 * v2 - HesseWeinberg(4, 6) + HesseWeinberg(5, 7)) /
                  v2 / v1); // dL3
  parCT.push_back(t1);      // dL4
  parCT.push_back(
      (t1 * v1 * v1 - 2 * HesseWeinberg(3, 3) + 2 * HesseWeinberg(7, 7)) *
      std::pow(v1, -2)); // dReL5
  parCT.push_back((2 * HesseWeinberg(4, 7) * v1 + 2 * HesseWeinberg(5, 6) * v1 +
                   2 * HesseWeinberg(6, 7) * v2) *
                  std::pow(v1, -2) / v2); // dImL5
  parCT.push_back((-HesseWeinberg(4, 7) - HesseWeinberg(5, 6)) *
                  std::pow(v2, -2)); // dImL7
  parCT.push_back(
      (-2 * HesseWeinberg(4, 8) * v1 - 2 * HesseWeinberg(7, 9) * v2) *
      std::pow(v1, -2) / vs); // dL8
  parCT.push_back((-2 * HesseWeinberg(6, 8) - 2 * HesseWeinberg(7, 9)) / vs /
                  v2);                                // dL9
  parCT.push_back(HesseWeinberg(7, 8) / vs / v1);     // dImL10
  parCT.push_back(2 * HesseWeinberg(4, 9) / vs / v1); // dImL11
  parCT.push_back(2 * HesseWeinberg(6, 9) / vs / v2); // dImL12
  parCT.push_back(2 * HesseWeinberg(7, 9) / vs / v1); // dReL13
  parCT.push_back(t2);                                // dLS
  parCT.push_back(-NablaWeinberg(2));                 // dTCB
  parCT.push_back(HesseWeinberg(5, 5) * v1 + HesseWeinberg(5, 7) * v2 -
                  NablaWeinberg(4)); // dT1
  parCT.push_back(HesseWeinberg(5, 7) * v1 + HesseWeinberg(7, 7) * v2 -
                  NablaWeinberg(6)); // dT2
  parCT.push_back(-HesseWeinberg(5, 6) * v1 - HesseWeinberg(6, 7) * v2 -
                  NablaWeinberg(7)); // dTCP
  parCT.push_back(t2 * std::pow(vs, 3) / 4 + HesseWeinberg(8, 8) * vs -
                  NablaWeinberg(8));                            // dTS
  parCT.push_back(HesseWeinberg(8, 9) * vs - NablaWeinberg(9)); // dTDM

  return parCT;
}

void Class_Potential_CN2HDM::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  std::vector<double> HiggsOrder(NHiggs);

  // todo: set mass order correctly!
  for (size_t i = 0; i < NHiggs; i++)
  {
    HiggsOrder[i] = i;
  }

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  for (size_t i = 0; i < NHiggs; i++)
  {
    for (size_t j = 0; j < NHiggs; j++)
    {
      for (size_t k = 0; k < NHiggs; k++)
      {
        GaugeBasis[i][j][k] =
            TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
      }
    }
  }

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (size_t i = 0; i < NHiggs; i++)
  {
    for (size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  MatrixXd HiggsRotSort(NHiggs, NHiggs);

  for (size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
  }

  TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
  TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
  TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
  for (size_t i = 0; i < NHiggs; i++)
  {
    TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
    for (size_t j = 0; j < NHiggs; j++)
    {
      TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
    }
  }

  for (size_t i = 0; i < NHiggs; i++)
  {
    for (size_t j = 0; j < NHiggs; j++)
    {
      for (size_t k = 0; k < NHiggs; k++)
      {
        TripleHiggsCorrectionsCWPhysical[i][j][k]   = 0;
        TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
        TripleHiggsCorrectionsCTPhysical[i][j][k]   = 0;
        for (size_t l = 0; l < NHiggs; l++)
        {
          for (size_t m = 0; m < NHiggs; m++)
          {
            for (size_t n = 0; n < NHiggs; n++)
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

void Class_Potential_CN2HDM::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  Curvature_Higgs_L2[0][0] = m11Sq;
  Curvature_Higgs_L2[0][2] = -Rem12Sq;
  Curvature_Higgs_L2[0][3] = Imm12Sq;
  Curvature_Higgs_L2[1][1] = m11Sq;
  Curvature_Higgs_L2[1][2] = -Imm12Sq;
  Curvature_Higgs_L2[1][3] = -Rem12Sq;
  Curvature_Higgs_L2[2][2] = m22Sq;
  Curvature_Higgs_L2[3][3] = m22Sq;
  Curvature_Higgs_L2[4][4] = m11Sq;
  Curvature_Higgs_L2[4][6] = -Rem12Sq;
  Curvature_Higgs_L2[4][7] = Imm12Sq;
  Curvature_Higgs_L2[5][5] = m11Sq;
  Curvature_Higgs_L2[5][6] = -Imm12Sq;
  Curvature_Higgs_L2[5][7] = -Rem12Sq;
  Curvature_Higgs_L2[6][6] = m22Sq;
  Curvature_Higgs_L2[7][7] = m22Sq;
  Curvature_Higgs_L2[8][8] = msSq / 0.2e1 + Rebs / 0.2e1;
  Curvature_Higgs_L2[9][9] = msSq / 0.2e1 - Rebs / 0.2e1;

  sym2Dim(Curvature_Higgs_L2, NHiggs, NHiggs);

  Curvature_Higgs_L4[0][0][0][0] = 3 * L1;
  Curvature_Higgs_L4[0][0][1][1] = L1;
  Curvature_Higgs_L4[0][0][2][2] = L3 + L4 + ReL5;
  Curvature_Higgs_L4[0][0][2][3] = -ImL5;
  Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - ReL5;
  Curvature_Higgs_L4[0][0][4][4] = L1;
  Curvature_Higgs_L4[0][0][5][5] = L1;
  Curvature_Higgs_L4[0][0][6][6] = L3;
  Curvature_Higgs_L4[0][0][7][7] = L3;
  Curvature_Higgs_L4[0][0][8][8] = L8 / 0.2e1;
  Curvature_Higgs_L4[0][0][9][9] = L8 / 0.2e1;
  Curvature_Higgs_L4[0][1][2][2] = ImL5;
  Curvature_Higgs_L4[0][1][2][3] = ReL5;
  Curvature_Higgs_L4[0][1][3][3] = -ImL5;
  Curvature_Higgs_L4[0][2][4][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
  Curvature_Higgs_L4[0][2][4][7] = -ImL5 / 0.2e1;
  Curvature_Higgs_L4[0][2][5][6] = ImL5 / 0.2e1;
  Curvature_Higgs_L4[0][2][5][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][4][6] = -ImL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][4][7] = L4 / 0.2e1 - ReL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][6] = -L4 / 0.2e1 + ReL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][7] = -ImL5 / 0.2e1;
  Curvature_Higgs_L4[1][1][1][1] = 3 * L1;
  Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - ReL5;
  Curvature_Higgs_L4[1][1][2][3] = ImL5;
  Curvature_Higgs_L4[1][1][3][3] = L3 + L4 + ReL5;
  Curvature_Higgs_L4[1][1][4][4] = L1;
  Curvature_Higgs_L4[1][1][5][5] = L1;
  Curvature_Higgs_L4[1][1][6][6] = L3;
  Curvature_Higgs_L4[1][1][7][7] = L3;
  Curvature_Higgs_L4[1][1][8][8] = L8 / 0.2e1;
  Curvature_Higgs_L4[1][1][9][9] = L8 / 0.2e1;
  Curvature_Higgs_L4[1][2][4][6] = ImL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][4][7] = -L4 / 0.2e1 + ReL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][6] = L4 / 0.2e1 - ReL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][7] = ImL5 / 0.2e1;
  Curvature_Higgs_L4[1][3][4][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
  Curvature_Higgs_L4[1][3][4][7] = -ImL5 / 0.2e1;
  Curvature_Higgs_L4[1][3][5][6] = ImL5 / 0.2e1;
  Curvature_Higgs_L4[1][3][5][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
  Curvature_Higgs_L4[2][2][2][2] = 3 * L2;
  Curvature_Higgs_L4[2][2][3][3] = L2;
  Curvature_Higgs_L4[2][2][4][4] = L3;
  Curvature_Higgs_L4[2][2][5][5] = L3;
  Curvature_Higgs_L4[2][2][6][6] = L2;
  Curvature_Higgs_L4[2][2][7][7] = L2;
  Curvature_Higgs_L4[2][2][8][8] = L9 / 0.2e1;
  Curvature_Higgs_L4[2][2][9][9] = L9 / 0.2e1;
  Curvature_Higgs_L4[3][3][3][3] = 3 * L2;
  Curvature_Higgs_L4[3][3][4][4] = L3;
  Curvature_Higgs_L4[3][3][5][5] = L3;
  Curvature_Higgs_L4[3][3][6][6] = L2;
  Curvature_Higgs_L4[3][3][7][7] = L2;
  Curvature_Higgs_L4[3][3][8][8] = L9 / 0.2e1;
  Curvature_Higgs_L4[3][3][9][9] = L9 / 0.2e1;
  Curvature_Higgs_L4[4][4][4][4] = 3 * L1;
  Curvature_Higgs_L4[4][4][5][5] = L1;
  Curvature_Higgs_L4[4][4][6][6] = L3 + L4 + ReL5;
  Curvature_Higgs_L4[4][4][6][7] = -ImL5;
  Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - ReL5;
  Curvature_Higgs_L4[4][4][8][8] = L8 / 0.2e1;
  Curvature_Higgs_L4[4][4][9][9] = L8 / 0.2e1;
  Curvature_Higgs_L4[4][5][6][6] = ImL5;
  Curvature_Higgs_L4[4][5][6][7] = ReL5;
  Curvature_Higgs_L4[4][5][7][7] = -ImL5;
  Curvature_Higgs_L4[5][5][5][5] = 3 * L1;
  Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - ReL5;
  Curvature_Higgs_L4[5][5][6][7] = ImL5;
  Curvature_Higgs_L4[5][5][7][7] = L3 + L4 + ReL5;
  Curvature_Higgs_L4[5][5][8][8] = L8 / 0.2e1;
  Curvature_Higgs_L4[5][5][9][9] = L8 / 0.2e1;
  Curvature_Higgs_L4[6][6][6][6] = 3 * L2;
  Curvature_Higgs_L4[6][6][7][7] = L2;
  Curvature_Higgs_L4[6][6][8][8] = L9 / 0.2e1;
  Curvature_Higgs_L4[6][6][9][9] = L9 / 0.2e1;
  Curvature_Higgs_L4[7][7][7][7] = 3 * L2;
  Curvature_Higgs_L4[7][7][8][8] = L9 / 0.2e1;
  Curvature_Higgs_L4[7][7][9][9] = L9 / 0.2e1;
  Curvature_Higgs_L4[8][8][8][8] = 0.3e1 / 0.4e1 * LS;
  Curvature_Higgs_L4[8][8][9][9] = LS / 0.4e1;
  Curvature_Higgs_L4[9][9][9][9] = 0.3e1 / 0.4e1 * LS;

  sym4Dim(Curvature_Higgs_L4, NHiggs, NHiggs, NHiggs, NHiggs);

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

  // different types of CN2HDM
  // vL: leptons, vu: u-type quarks, vd: d-type quarks
  // default: Type 1
  double vL = v2;
  double vd = v2;
  double vu = v2;

  if (Type == 2)
  {
    vL = v1;
    vd = v1;
  }
  else if (Type == 3)
    vd = v1;
  else if (Type == 4)
    vL = v1;

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

  Curvature_Lepton_F2H1[0][1][4] = 0.1e1 / vL * C_MassElectron;
  Curvature_Lepton_F2H1[0][1][5] = II / vL * C_MassElectron;
  Curvature_Lepton_F2H1[1][0][4] = 0.1e1 / vL * C_MassElectron;
  Curvature_Lepton_F2H1[1][0][5] = II / vL * C_MassElectron;
  Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / vL * C_MassElectron;
  Curvature_Lepton_F2H1[1][6][1] = II / vL * C_MassElectron;
  Curvature_Lepton_F2H1[2][3][4] = 0.1e1 / vL * C_MassMu;
  Curvature_Lepton_F2H1[2][3][5] = II / vL * C_MassMu;
  Curvature_Lepton_F2H1[3][2][4] = 0.1e1 / vL * C_MassMu;
  Curvature_Lepton_F2H1[3][2][5] = II / vL * C_MassMu;
  Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / vL * C_MassMu;
  Curvature_Lepton_F2H1[3][7][1] = II / vL * C_MassMu;
  Curvature_Lepton_F2H1[4][5][4] = 0.1e1 / vL * C_MassTau;
  Curvature_Lepton_F2H1[4][5][5] = II / vL * C_MassTau;
  Curvature_Lepton_F2H1[5][4][4] = 0.1e1 / vL * C_MassTau;
  Curvature_Lepton_F2H1[5][4][5] = II / vL * C_MassTau;
  Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / vL * C_MassTau;
  Curvature_Lepton_F2H1[5][8][1] = II / vL * C_MassTau;
  Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / vL * C_MassElectron;
  Curvature_Lepton_F2H1[6][1][1] = II / vL * C_MassElectron;
  Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / vL * C_MassMu;
  Curvature_Lepton_F2H1[7][3][1] = II / vL * C_MassMu;
  Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / vL * C_MassTau;
  Curvature_Lepton_F2H1[8][5][1] = II / vL * C_MassTau;

  Curvature_Quark_F2H1[0][1][6]   = 0.1e1 / vu * C_MassUp;
  Curvature_Quark_F2H1[0][1][7]   = -II / vu * C_MassUp;
  Curvature_Quark_F2H1[0][3][0]   = 0.1e1 / vd * C_MassDown * V11;
  Curvature_Quark_F2H1[0][3][1]   = II / vd * C_MassDown * V11;
  Curvature_Quark_F2H1[0][7][0]   = 0.1e1 / vd * C_MassStrange * V12;
  Curvature_Quark_F2H1[0][7][1]   = II / vd * C_MassStrange * V12;
  Curvature_Quark_F2H1[0][11][0]  = 0.1e1 / vd * C_MassBottom * V13;
  Curvature_Quark_F2H1[0][11][1]  = II / vd * C_MassBottom * V13;
  Curvature_Quark_F2H1[1][0][6]   = 0.1e1 / vu * C_MassUp;
  Curvature_Quark_F2H1[1][0][7]   = -II / vu * C_MassUp;
  Curvature_Quark_F2H1[1][2][2]   = -0.1e1 / vu * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[1][2][3]   = II / vu * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[1][6][2]   = -0.1e1 / vu * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[1][6][3]   = II / vu * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[1][10][2]  = -0.1e1 / vu * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[1][10][3]  = II / vu * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[2][1][2]   = -0.1e1 / vu * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[2][1][3]   = II / vu * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[2][3][4]   = 0.1e1 / vd * C_MassDown;
  Curvature_Quark_F2H1[2][3][5]   = II / vd * C_MassDown;
  Curvature_Quark_F2H1[2][5][2]   = -0.1e1 / vu * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[2][5][3]   = II / vu * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[2][9][2]   = -0.1e1 / vu * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][9][3]   = II / vu * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[3][0][0]   = 0.1e1 / vd * C_MassDown * V11;
  Curvature_Quark_F2H1[3][0][1]   = II / vd * C_MassDown * V11;
  Curvature_Quark_F2H1[3][2][4]   = 0.1e1 / vd * C_MassDown;
  Curvature_Quark_F2H1[3][2][5]   = II / vd * C_MassDown;
  Curvature_Quark_F2H1[3][4][0]   = V21 / vd * C_MassDown;
  Curvature_Quark_F2H1[3][4][1]   = II * V21 / vd * C_MassDown;
  Curvature_Quark_F2H1[3][8][0]   = 0.1e1 / vd * C_MassDown * V31;
  Curvature_Quark_F2H1[3][8][1]   = II / vd * C_MassDown * V31;
  Curvature_Quark_F2H1[4][3][0]   = V21 / vd * C_MassDown;
  Curvature_Quark_F2H1[4][3][1]   = II * V21 / vd * C_MassDown;
  Curvature_Quark_F2H1[4][5][6]   = 0.1e1 / vu * C_MassCharm;
  Curvature_Quark_F2H1[4][5][7]   = -II / vu * C_MassCharm;
  Curvature_Quark_F2H1[4][7][0]   = V22 / vd * C_MassStrange;
  Curvature_Quark_F2H1[4][7][1]   = II * V22 / vd * C_MassStrange;
  Curvature_Quark_F2H1[4][11][0]  = 0.1e1 / vd * C_MassBottom * V23;
  Curvature_Quark_F2H1[4][11][1]  = II / vd * C_MassBottom * V23;
  Curvature_Quark_F2H1[5][2][2]   = -0.1e1 / vu * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[5][2][3]   = II / vu * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[5][4][6]   = 0.1e1 / vu * C_MassCharm;
  Curvature_Quark_F2H1[5][4][7]   = -II / vu * C_MassCharm;
  Curvature_Quark_F2H1[5][6][2]   = -0.1e1 / vu * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[5][6][3]   = II / vu * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[5][10][2]  = -0.1e1 / vu * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[5][10][3]  = II / vu * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[6][1][2]   = -0.1e1 / vu * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[6][1][3]   = II / vu * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[6][5][2]   = -0.1e1 / vu * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[6][5][3]   = II / vu * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[6][7][4]   = 0.1e1 / vd * C_MassStrange;
  Curvature_Quark_F2H1[6][7][5]   = II / vd * C_MassStrange;
  Curvature_Quark_F2H1[6][9][2]   = -0.1e1 / vu * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[6][9][3]   = II / vu * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[7][0][0]   = 0.1e1 / vd * C_MassStrange * V12;
  Curvature_Quark_F2H1[7][0][1]   = II / vd * C_MassStrange * V12;
  Curvature_Quark_F2H1[7][4][0]   = V22 / vd * C_MassStrange;
  Curvature_Quark_F2H1[7][4][1]   = II * V22 / vd * C_MassStrange;
  Curvature_Quark_F2H1[7][6][4]   = 0.1e1 / vd * C_MassStrange;
  Curvature_Quark_F2H1[7][6][5]   = II / vd * C_MassStrange;
  Curvature_Quark_F2H1[7][8][0]   = 0.1e1 / vd * C_MassStrange * V32;
  Curvature_Quark_F2H1[7][8][1]   = II / vd * C_MassStrange * V32;
  Curvature_Quark_F2H1[8][3][0]   = 0.1e1 / vd * C_MassDown * V31;
  Curvature_Quark_F2H1[8][3][1]   = II / vd * C_MassDown * V31;
  Curvature_Quark_F2H1[8][7][0]   = 0.1e1 / vd * C_MassStrange * V32;
  Curvature_Quark_F2H1[8][7][1]   = II / vd * C_MassStrange * V32;
  Curvature_Quark_F2H1[8][9][6]   = 0.1e1 / vu * C_MassTop;
  Curvature_Quark_F2H1[8][9][7]   = -II / vu * C_MassTop;
  Curvature_Quark_F2H1[8][11][0]  = 0.1e1 / vd * C_MassBottom * V33;
  Curvature_Quark_F2H1[8][11][1]  = II / vd * C_MassBottom * V33;
  Curvature_Quark_F2H1[9][2][2]   = -0.1e1 / vu * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][2][3]   = II / vu * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][6][2]   = -0.1e1 / vu * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[9][6][3]   = II / vu * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[9][8][6]   = 0.1e1 / vu * C_MassTop;
  Curvature_Quark_F2H1[9][8][7]   = -II / vu * C_MassTop;
  Curvature_Quark_F2H1[9][10][2]  = -0.1e1 / vu * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[9][10][3]  = II / vu * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[10][1][2]  = -0.1e1 / vu * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[10][1][3]  = II / vu * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[10][5][2]  = -0.1e1 / vu * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[10][5][3]  = II / vu * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[10][9][2]  = -0.1e1 / vu * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[10][9][3]  = II / vu * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[10][11][4] = 0.1e1 / vd * C_MassBottom;
  Curvature_Quark_F2H1[10][11][5] = II / vd * C_MassBottom;
  Curvature_Quark_F2H1[11][0][0]  = V13 / vd * C_MassBottom;
  Curvature_Quark_F2H1[11][0][1]  = II / vd * C_MassBottom * V13;
  Curvature_Quark_F2H1[11][4][0]  = V23 / vd * C_MassBottom;
  Curvature_Quark_F2H1[11][4][1]  = II / vd * C_MassBottom * V23;
  Curvature_Quark_F2H1[11][8][0]  = V33 / vd * C_MassBottom;
  Curvature_Quark_F2H1[11][8][1]  = II / vd * C_MassBottom * V33;
  Curvature_Quark_F2H1[11][10][4] = 0.1e1 / vd * C_MassBottom;
  Curvature_Quark_F2H1[11][10][5] = II / vd * C_MassBottom;
}

bool Class_Potential_CN2HDM::CalculateDebyeSimplified()
{
  return false;
}

bool Class_Potential_CN2HDM::CalculateDebyeGaugeSimplified()
{
  return false;
}

double
Class_Potential_CN2HDM::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;

  double res = 0;

  double wcb, w1, w2, wcp, ws, wdm;
  wcb = v[2];
  w1  = v[4];
  w2  = v[6];
  wcp = v[7];
  ws  = v[8];
  wdm = v[9];

  double C11 = 0.5 * std::pow(w1, 2);
  double C22 = 0.5 * (std::pow(wcb, 2) + std::pow(w2, 2) + std::pow(wcp, 2));
  double CSS = 0.5 * (std::pow(wdm, 2) + std::pow(ws, 2));

  res += m11Sq * C11 + m22Sq * C22 + w2;
  res += w1 * (Imm12Sq * wcp - Rem12Sq * w2);
  res += 0.5 * msSq * CSS;
  res += 0.25 * Rebs * std::pow(ws, 2) - 0.25 * Rebs * std::pow(wdm, 2);
  res += 0.5 * L1 * std::pow(C11, 2);
  res += 0.5 * L2 * std::pow(C22, 2);
  res += L3 * C11 * C22;
  res += 0.5 * L4 * C11 * (std::pow(w2, 2) + std::pow(wcp, 2));
  res += 0.5 * C11 *
         ((std::pow(w2, 2) - std::pow(wcp, 2)) * ReL5 - 2 * wcp * w2 * ImL5);
  res += 0.5 * L8 * C11 * CSS;
  res += 0.5 * L9 * C22 * CSS;
  res += 0.125 * LS * std::pow(CSS, 2);

  return res;
}

double
Class_Potential_CN2HDM::VCounterSimplified(const std::vector<double> &v) const
{
  if (not UseVCounterSimplified) return 0;
  double res = 0;

  double wcb, w1, w2, wcp, ws, wdm;
  wcb = v[2];
  w1  = v[4];
  w2  = v[6];
  wcp = v[7];
  ws  = v[8];
  wdm = v[9];

  double C11 = std::pow(w1, 2);
  double C22 = std::pow(wcb, 2) + std::pow(w2, 2) + std::pow(wcp, 2);
  double CSS = std::pow(ws, 2) + std::pow(wdm, 2);

  res += C11 * dm11Sq;
  res += C22 * dm22Sq;
  res += -w1 * w2 * dRem12Sq;
  res += w1 * wcp * dImm12Sq;
  res += 0.5 * CSS * dmsSq;
  res += 0.25 * (std::pow(ws, 2) - std::pow(wdm, 2)) * dRebs;
  res += -0.5 * wdm * ws * dImbs;
  res += 0.5 * std::pow(C11, 2) * dL1;
  res += 0.5 * std::pow(C22, 2) * dL2;
  res += C11 * C22 * dL3;
  res += 0.5 * C11 * (std::pow(w2, 2) + std::pow(wcp, 2)) * dL4;
  res += 0.5 * C11 * (std::pow(w2, 2) - std::pow(wcp, 2)) * dReL5;
  res += -0.5 * w2 * wcp * C11 * dImL5;
  res += -w1 * wcp * C22 * dImL7;
  res += 0.5 * C11 * CSS * dL8;
  res += 0.5 * C22 * CSS * dL9;
  res += -w1 * wcp * CSS * dImL10;
  res += -0.5 * C11 * ws * wdm * dImL11;
  res += -0.5 * C22 * ws * wdm * dImL12;
  res += (0.25 * w1 * w2 * (std::pow(ws, 2) - std::pow(wdm, 2)) -
          2 * w1 * wcp * ws * wdm) *
         dReL13;
  res += 0.125 * std::pow(CSS, 2) * dLS;
  res += wcb * dTCB;
  res += w1 * dT1;
  res += w2 * dT2;
  res += wcp * dTCP;
  res += ws * dTS;
  res += wdm * dTDM;

  return res;
}

void Class_Potential_CN2HDM::Debugging(const std::vector<double> &input,
                                       std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
