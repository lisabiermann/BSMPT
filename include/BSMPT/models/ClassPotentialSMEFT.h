// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Model file for the SMEFT
 */

#pragma once

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT
{
namespace Models
{

/**
 * @brief The Class_Potential_SMEFT class
 * Implementation of the SMEFT
 */

class Class_Potential_SMEFT : public Class_Potential_Origin
{
public:
  Class_Potential_SMEFT(const ISMConstants &smConstants);
  virtual ~Class_Potential_SMEFT();

  double muSq, lambda;
  double Ouphi;

  double dmuSq, dlambda, dT1, dT2, dT3, dT4;
  double dZtL, dZtR, dZbL;
  double dZtLv2, dZtRv2, dZbLv2;
  double dZtLyv2, dZtRyv2, dZbLyv2;
  double dZtLv3, dZtRv3, dZbLv3;
  double dZtLy1v3, dZtRy1v3, dZbLy1v3;
  double dZtLy2v3, dZtRy2v3, dZbLy2v3;

  double LambdaEFT = 1000; // EFT scale

  double v0;

  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;
  std::vector<std::string> addLegendEFT() const override;
  std::vector<double> getParamsEFT() const override;

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void PerformVCTShift() override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;

  double SymFac_Higgs_OneLoop(const int &i,
                              const int &j,
                              const std::vector<double> &point) const override;
  double SymFac_Higgs_TwoLoop(const int &i, const int &j) const override;
};

} // namespace Models
} // namespace BSMPT
