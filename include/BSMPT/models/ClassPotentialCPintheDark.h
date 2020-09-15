/*
 * ClassPotentialCPintheDark.h
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

		This program is free software: you can redistribute it and/or modify
		it under the terms of the GNU General Public License as published by
		the Free Software Foundation, either version 3 of the License, or
		(at your option) any later version.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.

		You should have received a copy of the GNU General Public License
		along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
  * @file
  */

#pragma once

#include <string>                               // for string
#include <vector>                               // for vector

#include <BSMPT/models/ClassPotentialOrigin.h>
namespace BSMPT{
namespace Models{

/**
 * @brief The Class_Potential_CPintheDark class
 * Template for implementing a new model
 */
class Class_Potential_CPintheDark : public Class_Potential_Origin
{
public:
  Class_Potential_CPintheDark ();
  virtual
  ~Class_Potential_CPintheDark ();


  // Add here your parameters for the Lagrangian as well as for the counterterm potential
  // Add here your variables in which you will save the Debye correction factors

  // parameters of scalar potential
  double m11s, m22s, mSs, ReA, ImA, L1, L2, L3, L4, L5, L6, L7, L8;
  // counterterms
  double dm11s, dm22s, dmSs, dReA, dImA, dL1, dL2, dL3, dL4, dL5, dL6, dL7, dL8, dTCB, dT1, dT2, dTCP, dTS;
  // vev
  double v1;

  void ReadAndSet(const std::string& linestr, std::vector<double>& par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double>& par) override;
  void set_CT_Pot_Par(const std::vector<double>& par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double>& v) const override;
  double VCounterSimplified(const std::vector<double>& v) const override;
  void Debugging(const std::vector<double>& input, std::vector<double>& output) const override;
};

}
}

