// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

using Approx = Catch::Approx;

#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>

namespace
{

const std::vector<double> example_point_4D{/* lambdaS = */ 3,
                                           /* lambdaHS = */ 0.55,
                                           /* vs = */ 0,
                                           /* mS = */ 74,
                                           /* etakS = */ 0};

const std::vector<double> example_point_6D{/* lambdaS = */ 3,
                                           /* lambdaHS = */ 0.55,
                                           /* vs = */ 0,
                                           /* mS = */ 74,
                                           /* etakS = */ 30};

const std::vector<double> point = {1, 2, 3, 4, 5};

} // namespace

/**
 * @test Check if internal momentum dependent dim6 alpha factor modifies
 * one-loop potential like expected for 4D case
 */
TEST_CASE("Test AlphaVec 4D", "[eft]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer_4D =
      ModelID::FChoose(ModelID::ModelIDs::RXSMEFT, SMConstants);
  modelPointer_4D->initModel(example_point_4D);

  auto alpha_vec = modelPointer_4D->GetAlphaVec(point);

  REQUIRE(alpha_vec.at(0) == 1);
  REQUIRE(alpha_vec.at(1) == 1);
  REQUIRE(alpha_vec.at(2) == 1);
  REQUIRE(alpha_vec.at(3) == 1);
  REQUIRE(alpha_vec.at(4) == 1);
}

/**
 * @test Check if internal momentum dependent dim6 alpha factor modifies
 * one-loop potential like expected for 6D case
 */
TEST_CASE("Test AlphaVec 6D", "[eft]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer_6D =
      ModelID::FChoose(ModelID::ModelIDs::RXSMEFT, SMConstants);
  modelPointer_6D->initModel(example_point_6D);

  auto alpha_vec = modelPointer_6D->GetAlphaVec(point);

  REQUIRE(alpha_vec.at(0) == 1);
  REQUIRE(alpha_vec.at(1) == 1);
  REQUIRE(alpha_vec.at(2) == 1);
  REQUIRE(alpha_vec.at(3) == 1);
  REQUIRE(alpha_vec.at(4) == 1 - 30 * point.at(4) * point.at(4) / 1000 / 1000);
}

/**
 * @brief Check internal momentum dependent dim6 modification of
 * zero-temperature part of boson function
 */
TEST_CASE("Test VCW 6D", "[eft]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer_6D =
      ModelID::FChoose(ModelID::ModelIDs::RXSMEFT, SMConstants);
  modelPointer_6D->initModel(example_point_6D);

  REQUIRE(std::pow(2, -4) * modelPointer_6D->boson(10, 0, 2, 0, 1) ==
          modelPointer_6D->boson(10, 0, 2, 0, 2));
}

/**
 * @brief Check internal momentum dependent dim6 modification of
 * finite-temperature part of boson function
 */
TEST_CASE("Test VT 6D", "[eft]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer_6D =
      ModelID::FChoose(ModelID::ModelIDs::RXSMEFT, SMConstants);
  modelPointer_6D->initModel(example_point_6D);

  double MassSquared = 10;
  double Temp        = 1;
  double alpha       = 2;
  double cb          = 2;
  double diff        = 0;

  double Ratio   = MassSquared / std::pow(alpha * Temp, 2);
  double VCWterm = modelPointer_6D->CWTerm(std::abs(MassSquared), cb, diff);
  double VTterm  = std::pow(alpha * Temp, 4) / (2 * std::pow(M_PI, 2)) *
                  ThermalFunctions::JbosonInterpolated(Ratio);

  REQUIRE(std::pow(alpha, -4) * (VCWterm + VTterm) ==
          modelPointer_6D->boson(MassSquared, Temp, cb, diff, alpha));
}