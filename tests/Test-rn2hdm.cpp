#include "catch.hpp"

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/ModelTestfunctions.h>
#include <BSMPT/utility.h>

const std::vector<double> example_point_RN2HDM{/* lambda_1 = */ 0.300812,
                                               /* lambda_2 = */ 0.321809,
                                               /* lambda_3 = */ -0.133425,
                                               /* lambda_4 = */ 4.11105,
                                               /* lambda_5 = */ -3.84178,
                                               /* lambda_6 = */ 9.46329,
                                               /* lambda_7 = */ -0.750455,
                                               /* lambda_8 = */ 0.743982,
                                               /* tan(beta) = */ 5.91129,
                                               /* v_s = */ 293.035,
                                               /* m_{12}^2 = */ 4842.28,
                                               /* Yukawa Type = */ 1};

constexpr auto Model = BSMPT::ModelID::ModelIDs::RN2HDM;

TEST_CASE("Checking NLOVEV for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(Model);
  modelPointer->initModel(example_point_RN2HDM);
  std::vector<double> Check;
  auto sol = Minimizer::Minimize_gen_all(modelPointer,
                                         0,
                                         Check,
                                         modelPointer->get_vevTreeMin(),
                                         Minimizer::WhichMinimizerDefault);
  for (std::size_t i{0}; i < sol.size(); ++i)
  {
    auto expected = std::abs(modelPointer->get_vevTreeMin(i));
    auto res      = std::abs(sol.at(i));
    REQUIRE(std::abs(res - expected) <= 1e-4);
  }
}

TEST_CASE("Checking EWPT for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(Model);
  modelPointer->initModel(example_point_RN2HDM);
  std::vector<double> Check;
  auto EWPT = Minimizer::PTFinder_gen_all(
      modelPointer, 0, 300, Minimizer::WhichMinimizerDefault);
  const double omega_c_expected = 180.5917335676395;
  const double Tc_expected      = 120.7305908203125;
  const std::vector<double> min_expected{
      0, 0, -32.70827526931041, -177.6050195289305, -297.0418903961274};

  REQUIRE(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS);
  std::cout << "Passed" << std::endl;
  REQUIRE(std::abs(omega_c_expected - EWPT.vc) / omega_c_expected <= 1e-2);
  std::cout << "Passed" << std::endl;
  REQUIRE(std::abs(Tc_expected - EWPT.Tc) / Tc_expected <= 1e-2);
  std::cout << "Passed" << std::endl;
  for (std::size_t i{0}; i < EWPT.EWMinimum.size(); ++i)
  {
    auto res      = std::abs(EWPT.EWMinimum.at(i));
    auto expected = std::abs(min_expected.at(i));
    if (expected != 0)
    {
      std::cout << std::abs(res - expected) / expected << std::endl;
      REQUIRE(std::abs(res - expected) / expected <= 1e-2);
    }
    else
    {
      REQUIRE(res <= 1e-2);
    }
  }
}

TEST_CASE("Checking number of CT parameters for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  modelPointer->initModel(example_point_RN2HDM);
  auto result = ModelTests::CheckNumberOfCTParameters(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of VEV labels for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  modelPointer->initModel(example_point_RN2HDM);
  auto result = ModelTests::CheckNumberOfVEVLabels(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE(
    "Checking number of labels for temperature dependend results for N2HDM",
    "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  modelPointer->initModel(example_point_RN2HDM);
  auto result = ModelTests::CheckLegendTemp(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of triple Higgs couplings for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  modelPointer->initModel(example_point_RN2HDM);
  auto result = ModelTests::CheckNumberOfTripleCouplings(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking Gauge Boson masses for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  modelPointer->initModel(example_point_RN2HDM);
  auto result = ModelTests::CheckGaugeBosonMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking fermion and quark masses masses for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  modelPointer->initModel(example_point_RN2HDM);
  auto result = ModelTests::CheckFermionicMasses(*modelPointer);
  REQUIRE(result.first == ModelTests::TestResults::Pass);
  REQUIRE(result.second == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level minimum for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  modelPointer->initModel(example_point_RN2HDM);
  auto result = ModelTests::CheckTreeLevelMin(*modelPointer,
                                              Minimizer::WhichMinimizerDefault);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level tadpoles for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  modelPointer->initModel(example_point_RN2HDM);
  auto result = ModelTests::CheckTadpoleRelations(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking NLO masses matching tree level masses for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  modelPointer->initModel(example_point_RN2HDM);
  auto result = ModelTests::CheckNLOMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking VTreeSimplified for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  if (modelPointer->UseVTreeSimplified)
  {
    modelPointer->initModel(example_point_RN2HDM);
    auto result = ModelTests::CheckVTreeSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking VCounterSimplified for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::RN2HDM);
  if (modelPointer->UseVCounterSimplified)
  {
    modelPointer->initModel(example_point_RN2HDM);
    auto result = ModelTests::CheckVCounterSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking CKM Unitarity for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  REQUIRE(ModelTests::CheckCKMUnitarity() == ModelTests::TestResults::Pass);
}