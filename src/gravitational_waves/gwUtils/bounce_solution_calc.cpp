// Copyright (C) 2023 Lisa Biermann, João Viana
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file bounce solution calculator class
 */

#include <BSMPT/gravitational_waves/gwUtils/bounce_solution_calc.h>

namespace BSMPT
{
namespace GravitationalWaves
{

/**
 * @brief Default constructor
 *
 */
BounceSolCalc::BounceSolCalc()
{
}

/**
 * @brief Constructor used for testing
 */
BounceSolCalc::BounceSolCalc(
    std::shared_ptr<Class_Potential_Origin> &pointer_in)
{
  modelPointer = pointer_in;
}

/**
 * @brief Constructor
 */
BounceSolCalc::BounceSolCalc(
    std::shared_ptr<Class_Potential_Origin> &pointer_in,
    CoexPhases &coex_in,
    int indexTrueCandidatePhase_in,
    double UserDefined_vwall_in,
    double UserDefined_epstol_in,
    int MaxPathIntegrations_in)
{
  coex                    = coex_in;
  indexTrueCandidatePhase = indexTrueCandidatePhase_in;
  Tc                      = coex.CalculateTc(indexTrueCandidatePhase);
  Tm                      = coex.T_low;
  modelPointer            = pointer_in;
  UserDefined_vwall       = UserDefined_vwall_in;
  epstol                  = UserDefined_epstol_in;
  MaxPathIntegrations     = MaxPathIntegrations_in;
  this->CalcGstarPureRad(); // initialize degrees of freedom for purely
                            // radiative universe

  if (Tc > 0)
    BounceSolCalc::GWInitialScan();
  else
    failed = true;
}

void BounceSolCalc::GWInitialScan()
{
  if (Tc < 0)
  {
    // Transition is never viable
    return;
  }

  double last_action = -1;
  double dT          = (Tc - coex.T_low) / 7.;
  std::vector<double> TrueVacuum, FalseVacuum, last_TrueVacuum,
      last_FalseVacuum;
  std::vector<std::vector<double>> last_path, path;

  for (double T = Tc - dT; T >= coex.T_low + dT; T -= dT)
  {
    BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed,
                         "T = " + std::to_string(T));

    TrueVacuum  = coex.vec_phases[indexTrueCandidatePhase].Get(T).point;
    FalseVacuum = coex.vec_phases[0].Get(T).point;
    std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
    {
      // Potential wrapper
      return coex.modelPointer->VEff(coex.modelPointer->MinimizeOrderVEV(vev),
                                     T);
    };
    if (last_action < 0)
    {
      path = {TrueVacuum, FalseVacuum};
    }
    else
    {
      path = coex.MinTracer.WarpPath(last_path,
                                     last_TrueVacuum,
                                     last_FalseVacuum,
                                     TrueVacuum,
                                     FalseVacuum);
    }
    BSMPT::GravitationalWaves::BounceActionInt bc(
        path, TrueVacuum, FalseVacuum, V, T, MaxPathIntegrations);
    bc.CalculateAction();

    last_path        = bc.path;
    last_TrueVacuum  = bc.TrueVacuum;
    last_FalseVacuum = bc.FalseVacuum;

    // Comment this is you want dumb paths!!
    last_action = bc.action;
    if (bc.action / T > 0)
    {
      SolutionList.insert(
          std::upper_bound(SolutionList.begin(),
                           SolutionList.end(),
                           bc,
                           [](BSMPT::GravitationalWaves::BounceActionInt a,
                              BSMPT::GravitationalWaves::BounceActionInt b)
                           { return a.T < b.T; }),
          bc);
      // SolutionList.push_back(bc);
    }

    if (bc.action / T < 40 and bc.action > 0) break;
  }
  GWSecondaryScan();
}

void BounceSolCalc::CalculateActionAt(double T, bool smart)
{
  // Action outside allowed range
  if (T < Tm or T > Tc) return;
  BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed,
                       " T = " + std::to_string(T));
  // Find the closest solution to our goal temperature
  if (SolutionList.size() > 0)
  {
    auto it =
        std::min_element(SolutionList.begin(),
                         SolutionList.end(),
                         [T](BSMPT::GravitationalWaves::BounceActionInt a,
                             BSMPT::GravitationalWaves::BounceActionInt b)
                         { return std::abs(T - a.T) < std::abs(T - b.T); });
    BSMPT::GravitationalWaves::BounceActionInt Nearest_bc = *it;

    if (abs(Nearest_bc.T - T) < 0.001) return;

    std::vector<double> TrueVacuum =
        coex.vec_phases[indexTrueCandidatePhase].Get(T).point;
    std::vector<double> FalseVacuum = coex.vec_phases[0].Get(T).point;
    std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
    {
      // Potential wrapper
      return coex.modelPointer->VEff(coex.modelPointer->MinimizeOrderVEV(vev),
                                     T);
    };
    std::vector<std::vector<double>> path;

    if (smart)
      path = coex.MinTracer.WarpPath(Nearest_bc.path,
                                     Nearest_bc.TrueVacuum,
                                     Nearest_bc.FalseVacuum,
                                     TrueVacuum,
                                     FalseVacuum);
    else
      path = {TrueVacuum, FalseVacuum};

    BSMPT::GravitationalWaves::BounceActionInt bc(
        path, TrueVacuum, FalseVacuum, V, T, MaxPathIntegrations);
    bc.CalculateAction();
    if (bc.action / T > 0)
    {
      SolutionList.insert(
          std::upper_bound(SolutionList.begin(),
                           SolutionList.end(),
                           bc,
                           [](BSMPT::GravitationalWaves::BounceActionInt a,
                              BSMPT::GravitationalWaves::BounceActionInt b)
                           { return a.T < b.T; }),
          bc);
    }
  }
  else
  {
    std::vector<double> TrueVacuum =
        coex.vec_phases[indexTrueCandidatePhase].Get(T).point;
    std::vector<double> FalseVacuum = coex.vec_phases[0].Get(T).point;
    std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
    {
      // Potential wrapper
      return coex.modelPointer->VEff(coex.modelPointer->MinimizeOrderVEV(vev),
                                     T);
    };
    std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

    BSMPT::GravitationalWaves::BounceActionInt bc(
        path, TrueVacuum, FalseVacuum, V, T, MaxPathIntegrations);
    bc.CalculateAction();
    if (bc.action / T > 0)
    {
      SolutionList.push_back(bc);
    }
  }
}

void BounceSolCalc::GWSecondaryScan()
{
  if (SolutionList.size() == 0)
  {
    BSMPT::Logger::Write(
        BSMPT::LoggingLevel::GWDetailed,
        "No solution was found during the initial scan.\n Abort!\n");
    // Failed to find the solution
    failed = true;
    return;
  }
  else if (SolutionList.size() == 1)
  {
    BSMPT::Logger::Write(
        BSMPT::LoggingLevel::GWDetailed,
        "Only one solution was found. Searching near what was found\n");
    CalculateActionAt(SolutionList[0].T - 0.01);
    if (SolutionList.size() == 1)
    {
      // Failed to find the solution
      SetBounceSol();
      failed = true;
      return;
    };
  }

  GWScanTowardsLowAction();
  GWScanTowardsHighAction();

  std::size_t NumOfSol = SolutionList.size();

  for (int i = 0; i < 2; i++)
  {
    std::vector<double> nextTList;

    for (std::size_t sol = 0; sol < SolutionList.size() - 1; sol++)
    {
      if (SolutionList[sol].action / SolutionList[sol].T > 200) continue;
      if (SolutionList[sol + 1].action / SolutionList[sol + 1].T < 50) continue;

      double NumOfT =
          ceil((SolutionList[sol + 1].action / SolutionList[sol + 1].T -
                SolutionList[sol].action / SolutionList[sol].T) /
               20);

      double dT = (SolutionList[sol + 1].T - SolutionList[sol].T) / NumOfT;
      for (int it_T = 1; it_T < NumOfT; it_T++)
      {
        nextTList.push_back(SolutionList[sol].T + it_T * dT);
      }
    }

    for (auto j : nextTList)
    {
      CalculateActionAt(j);
    }

    if (NumOfSol == SolutionList.size()) break;
    NumOfSol = SolutionList.size();
  }
  GWScanTowardsLowAction();
  GWScanTowardsHighAction();
  SetBounceSol();
}

void BounceSolCalc::GWScanTowardsHighAction()
{
  for (int i = 0;
       i <=
       1. + (200 - SolutionList.back().action / SolutionList.back().T) / 5.;
       i++)
  {
    if (SolutionList.back().action / SolutionList.back().T > 200) break;

    // Try to calculate until s3/T = 200
    double t1            = SolutionList[SolutionList.size() - 2].T;
    double t2            = SolutionList[SolutionList.size() - 1].T;
    double s1            = SolutionList[SolutionList.size() - 2].action / t1;
    double s2            = SolutionList[SolutionList.size() - 1].action / t2;
    double goal          = s2 + 5;
    std::size_t NumOfSol = SolutionList.size();

    double T = ((s1 - goal) * t2 - (s2 - goal) * t1) / (s1 - s2);

    // Action is not monotonic
    if (T < t2) return;

    if (T > this->Tc) // Action is flat
    {
      if (t2 + (t2 - t1) > this->Tc) return; // Close already
      CalculateActionAt(t2 + (t2 - t1));
      i--;
    }
    else
    {
      CalculateActionAt(T);
      if (NumOfSol == SolutionList.size())
        CalculateActionAt(
            T, false); // No solution was found. Calculate using dumb method
    }
    if (NumOfSol == SolutionList.size()) return; // No solution was found. Abort
  }
}

void BounceSolCalc::GWScanTowardsLowAction()
{
  for (int i = 0;
       i <=
       1. + (SolutionList.front().action / SolutionList.front().T - 50) / 10.;
       i++)
  {
    if (SolutionList.front().action / SolutionList.front().T < 50) break;
    // Try to calculate it at s3/T = 300
    double t1            = SolutionList[0].T;
    double t2            = SolutionList[1].T;
    double s1            = SolutionList[0].action / t1;
    double s2            = SolutionList[1].action / t2;
    double goal          = s1 - 10;
    std::size_t NumOfSol = SolutionList.size();

    double T = ((s1 - goal) * t2 - (s2 - goal) * t1) / (s1 - s2);

    // Action is not monotonic
    if (T > t1) return;

    if (T < this->Tm) // Action is flat
    {
      if (t1 - (t2 - t1) < this->Tm) return; // Close already
      CalculateActionAt(t1 - (t2 - t1));
      i--;
    }
    else
    {
      CalculateActionAt(T);
      if (NumOfSol == SolutionList.size())
        CalculateActionAt(
            T, false); // No solution was found. Calculate using dumb method
    }
    if (NumOfSol == SolutionList.size()) return; // No solution was found. Abort
  }
}

void BounceSolCalc::SetBounceSol()
{

  std::vector<double> list_T, list_S3;
  std::stringstream ss;
  ss << "------------ Solution list ------------\n";
  for (auto sol : SolutionList)
  {
    ss << std::setprecision(10) << "{" << sol.T << ",\t" << sol.action << ",\t"
       << sol.action / sol.T << "},\n";
    list_T.push_back(sol.T);
    list_S3.push_back(sol.action);
  }
  ss << "---------------------------------------\n";
  BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed, ss.str());
  if (SolutionList.size() < 4)
  {
    failed = true;
    BSMPT::Logger::Write(
        BSMPT::LoggingLevel::GWDetailed,
        "There are not enough points to calculate the path. Abort.\n");
    return; // Not enough points to calculate path
  }

  S3ofT_spline.set_boundary(
      tk::spline::not_a_knot, 0.0, tk::spline::not_a_knot, 0.0);
  S3ofT_spline.set_points(list_T, list_S3);
}

/**
 * @brief Get the bubble wall velocity
 * @return vb
 */
double BounceSolCalc::GetWallVelocity()
{
  return vwall;
}

/**
 * @brief Get epsturb
 * @return epsturb
 */
double BounceSolCalc::GetEpsTurb()
{
  return epstol;
}

/**
 * @brief SetGstar Set gstar
 */
void BounceSolCalc::SetGstar(const double &gstar_in)
{
  gstar = gstar_in;
}
/**
 * @brief GetGstar Get gstar
 */
double BounceSolCalc::GetGstar()
{
  return gstar;
}
/**
 * @brief SetCriticalTemp Set critical temperature
 */
void BounceSolCalc::SetCriticalTemp(const double &T_in)
{
  Tc = T_in;
}
/**
 * @brief GetCriticalTemp Get critical temperature
 */
double BounceSolCalc::GetCriticalTemp()
{
  return Tc;
}

/**
 * @brief SetStoredTemp Set stored temperature
 */
void BounceSolCalc::SetStoredTemp(const double &T_in)
{
  store_Temp = T_in;
}
/**
 * @brief GetStoredTemp Get stored temperature
 */
double BounceSolCalc::GetStoredTemp()
{
  return store_Temp;
}

/**
 * @brief GetNucleationTemp Get nucleation temperature via exact method
 */
double BounceSolCalc::GetNucleationTemp()
{
  return Tnucl;
}

/**
 * @brief GetNucleationTempApprox Get nucleation temperature via approximate
 * method
 */
double BounceSolCalc::GetNucleationTempApprox()
{
  return Tnucl_approx;
}

/**
 * @brief GetPercolationTemp Get percolation temperature
 */
double BounceSolCalc::GetPercolationTemp()
{
  return Tperc;
}

/**
 * @brief GetPTStrength Get PT strength alpha
 */
double BounceSolCalc::GetPTStrength()
{
  return alpha;
}

/**
 * @brief Storage of the tunneling rate per volume of the transition from
 * false to true vacuum
 * @param Temp temperature
 */
double BounceSolCalc::Tunneling_Rate(const double &Temp)
{
  if (Temp < SolutionList.front().T or Temp > SolutionList.back().T)
    return 0; // Never extrapolate
  double Shat3 = GetBounceSol(Temp);
  if (Shat3 < 0)
  {
    BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed,
                         "Action spline became negative somewhere. T = " +
                             std::to_string(Temp));
    return 1e100; // Spline is unstable or ridiculous action. Return values that
                  // kills the tunnekinug rate
  }

  double amp =
      std::pow(Temp, 4) * std::pow((Shat3 / (2 * M_PI * Temp)), 3. / 2);
  return amp * std::exp(-Shat3 / Temp);
}
/**
 * @brief Storage of the temperature-dependent Hubble rate
 * @param Temp temperature
 */
double BounceSolCalc::Hubble_Rate(const double &Temp)
{
  return M_PI * std::sqrt(this->GetGstar()) / std::sqrt(90.) *
         std::pow(Temp, 2) / MPl; // radiation dominated universe
}

/**
 * @brief CalcGstar Calculate the number of effective degrees of freedom as a
 * function of temperature, higgs masses are evaluated at global minimum
 * @param Temp temperature
 */
void BounceSolCalc::CalcGstar(const double &Temp)
{
  std::size_t NHiggs = this->modelPointer->get_NHiggs();

  const std::vector<double> HiggsMassesSquared =
      this->modelPointer->HiggsMassesSquared(
          this->modelPointer->MinimizeOrderVEV(
              this->modelPointer->get_vevTreeMin()),
          Temp); // todo: calculate HiggsMassesSquared not in global, but in
                 // minimum in which we are still in (e.g. local one until
                 // T < Tperc)

  std::size_t HiggsAboveThreshold = 0;

  for (double mass : HiggsMassesSquared)
  {
    if (mass > Temp)
    {
      HiggsAboveThreshold++;
    }
  }

  double gb   = 8 * 2 + 3 * 3 + 2 + NHiggs - HiggsAboveThreshold;
  double gf   = 6 * 3 * 2 * 2 + 3 * 2 * 2 + 3 * 2;
  double geff = gb + 7. / 8 * gf;

  this->SetGstar(geff);
}

/**
 * @brief CalcGstarPureRad Calculate the number of effective degrees of
 * freedom assuming a purely radiative universe
 */
void BounceSolCalc::CalcGstarPureRad()
{
  std::size_t NHiggs = this->modelPointer->get_NHiggs();

  double gb   = 8 * 2 + 3 * 3 + 2 + NHiggs;
  double gf   = 6 * 3 * 2 * 2 + 3 * 2 * 2 + 3 * 2;
  double geff = gb + 7. / 8 * gf;

  this->SetGstar(geff);
}

/**
 * @brief Calculation of nucleation temperature
 */
void BounceSolCalc::CalculateNucleationTemp()
{
  if (failed) // Not enough points to calculate path
  {
    BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed,
                         "There are not enough points to calculate the path, "
                         "nucleation temperature calculation failed.\n");
    Tnucl = -1;
    return;
  }
  double T_up   = -1;
  double T_down = -1;
  double T_middle;

  for (auto sol = SolutionList.rbegin(); sol != SolutionList.rend(); sol++)
  {
    // Catches the first interval with the nucleation temperature
    if (T_up == -1 and
        Tunneling_Rate(sol->T) / std::pow(Hubble_Rate(sol->T), 4) < 1)
      T_up = sol->T;

    if (T_down == -1 and
        Tunneling_Rate(sol->T) / std::pow(Hubble_Rate(sol->T), 4) > 1)
      T_down = sol->T;

    if (T_up > 0 and T_down > 0) break;
  }

  if (T_up > 0 and T_down > 0)
  {
    // There is a Tn to be calculated! Use bisection method
    for (int i = 0; i < 50; i++)
    {
      T_middle = (T_up + T_down) / 2;

      if (Tunneling_Rate(T_middle) / std::pow(Hubble_Rate(T_middle), 4) < 1)
      {
        T_up = T_middle;
      }
      else
      {
        T_down = T_middle;
      }
      if (std::abs(T_up / T_down - 1) < 1e-10)
      {
        Tnucl               = T_middle;
        nucleation_temp_set = true;
        return;
      }
    }
  }
  std::stringstream ss;
  ss << "Nucleation temperature calculation failed\n";

  if (T_up < 0) ss << "Tunneling rate/H never is less than 1.\n";
  if (T_down < 0) ss << "Tunneling rate/H never is more than 1.\n";

  BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed, ss.str());

  // Unsuccessfull
  Tnucl = -1;
  return;
}

/**
 * @brief Approximate calculation of nucleation temperature
 */
void BounceSolCalc::CalculateNucleationTempApprox()
{
  if (failed) // Not enough points to calculate path
  {
    BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed,
                         "There are not enough points to calculate the path, "
                         "nucleation temperature calculation failed.\n");
    Tnucl_approx = -1;
    return;
  }

  double T_up   = -1;
  double T_down = -1;
  double T_middle;

  for (auto sol = SolutionList.rbegin(); sol != SolutionList.rend(); sol++)
  {
    // Catches the first interval with the nucleation temperature
    if (T_up == -1 and GetBounceSol(sol->T) / sol->T < 140) T_up = sol->T;

    if (T_down == -1 and GetBounceSol(sol->T) / sol->T > 140) T_down = sol->T;

    if (T_up > 0 and T_down > 0) break;
  }

  if (T_up > 0 and T_down > 0)
  {
    // There is a Tn to be calculated! Use bisection method
    for (int i = 0; i < 50; i++)
    {
      T_middle = (T_up + T_down) / 2;

      if (GetBounceSol(T_middle) / T_middle < 140)
      {
        T_up = T_middle;
      }
      else
      {
        T_down = T_middle;
      }
      if (std::abs(T_up / T_down - 1) < 1e-10)
      {
        Tnucl_approx = T_middle;
        return;
      }
    }
  }
  std::stringstream ss;
  ss << "Approximate nucleation temperature calculation failed\n";

  if (T_up < 0) ss << "S3(T)/T never is less than 140.\n";
  if (T_down < 0) ss << "S3(T)/T never is more than 140.\n";

  BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed, ss.str());

  // Unsuccessfull
  Tnucl_approx = -1;
  return;
}

/**
 * @brief Calculation of percolation temperature
 */
void BounceSolCalc::CalculatePercolationTemp()
{
  if (failed or Tnucl == -1)
  {
    Tperc = -1;
    return; // Not enough points to calculate path or Tn failed
  }

  double T_up             = -1;
  double T_down           = -1;
  double prefac           = 4. * M_PI / 3. * std::pow(vwall, 3);
  double Int_at_perc_temp = 0.34;
  double T_middle;
  for (auto sol = SolutionList.rbegin(); sol != SolutionList.rend(); sol++)
  {
    // Catches the first interval with the percolation temperature
    this->SetStoredTemp(sol->T);
    double IatT_solT = prefac * Nintegrate_Outer(*this).result;

    if (T_up == -1 and IatT_solT < Int_at_perc_temp) T_up = sol->T;
    if (T_down == -1 and IatT_solT > Int_at_perc_temp) T_down = sol->T;
    if (T_up > 0 and T_down > 0) break;
  }

  std::stringstream ss;

  ss << "T_up = " << T_up << "\n";
  ss << "T_down = " << T_down << "\n";

  if (T_up > 0 and T_down > 0)
  {
    T_middle = (T_up + T_down) / 2.;
    this->SetStoredTemp(T_middle); // update temp storage for inner integral
    double IatT = prefac * Nintegrate_Outer(*this).result;

    while (std::abs(T_up / T_down - 1) > 1e-10)
    {
      T_middle = (T_up + T_down) / 2.;
      this->SetStoredTemp(T_middle); // update temp storage for inner integral
      IatT = prefac * Nintegrate_Outer(*this).result;

      ss << "I ( T = " << T_middle << " ) = " << IatT << "\n";

      if (IatT < Int_at_perc_temp)
      {
        T_up = T_middle;
      }
      else
      {
        T_down = T_middle;
      }
    }

    Tperc = T_middle;

    if (Tperc > 0 and percolation_temp_set == false)
    {
      // Try to calculate action at Tp
      // CalculateActionAt(Tperc);
      for (std::size_t i = 0; i < SolutionList.size() - 1; i++)
      {
        if (Tperc > SolutionList[i].T and Tperc < SolutionList[i + 1].T)
          CalculateActionAt((SolutionList[i].T + SolutionList[i + 1].T) / 2);
      }
      percolation_temp_set = true;
      SetBounceSol();
      CalculatePercolationTemp();
    }

    BSMPT::Logger::Write(BSMPT::LoggingLevel::Debug, ss.str());
  }
  else // unsuccessful
  {
    std::stringstream err;
    err << "Percolation temperature calculation failed\n";

    if (T_up == -1) err << "I(T) is never less than 0.34.\n";
    if (T_down == -1) err << "I(T) is never more than 0.34.\n";

    BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed, err.str());

    Tperc = -1;
  }

  return;
}

double inner_integrand(double Temp, void *params)
{
  struct BounceSolCalc &obj = *static_cast<BounceSolCalc *>(params);
  double func               = 1. / obj.Hubble_Rate(Temp);
  return func;
}

double outer_integrand(double Temp, void *params)
{
  struct BounceSolCalc &obj = *static_cast<BounceSolCalc *>(params);

  double func = obj.Tunneling_Rate(Temp) /
                (std::pow(Temp, 4) * obj.Hubble_Rate(Temp)) *
                std::pow(Nintegrate_Inner(obj, Temp).result, 3);
  return func;
}

/**
 * @brief Nintegrate_Inner Numerical integration of inner integral over
 * inverse Hubble rate for the percolation temperature calculation
 * @param obj Class reference to pass all needed parameters
 * @param T upper integration boundary
 * @return Numerical value of integral and absolute error
 */
struct resultErrorPair Nintegrate_Inner(BounceSolCalc &obj,
                                        const double &Tprime)
{
  double abs_err = obj.AbsErr;
  double rel_err = obj.RelErr;

  std::size_t workspace_size = 1000;
  gsl_integration_workspace *w =
      gsl_integration_workspace_alloc(workspace_size);
  gsl_function F;
  F.function = &inner_integrand;
  F.params   = static_cast<void *>(&obj);

  struct resultErrorPair res;

  gsl_integration_qags(&F,
                       obj.GetStoredTemp(),
                       Tprime,
                       abs_err,
                       rel_err,
                       workspace_size,
                       w,
                       &res.result,
                       &res.error);

  gsl_integration_workspace_free(w);

  return res;
}

/**
 * @brief Nintegrate_Outer Numerical integration of outer integral for the
 * percolation temperature calculation
 * @param obj Class reference to pass all needed parameters
 * @param T lower integration boundary
 * @return Numerical value of integral and absolute error
 */
struct resultErrorPair Nintegrate_Outer(BounceSolCalc &obj)
{
  double abs_err = obj.AbsErr;
  double rel_err = obj.RelErr;

  std::size_t workspace_size = 1000;
  gsl_integration_workspace *w =
      gsl_integration_workspace_alloc(workspace_size);
  gsl_function F;
  F.function = &outer_integrand;
  F.params   = static_cast<void *>(&obj);

  struct resultErrorPair res;

  gsl_integration_qags(&F,
                       obj.GetStoredTemp(),
                       obj.GetCriticalTemp(),
                       abs_err,
                       rel_err,
                       workspace_size,
                       w,
                       &res.result,
                       &res.error);

  gsl_integration_workspace_free(w);

  return res;
}

/**
 * @brief Calculate phase transition strength alpha at percolation temperature
 * and vwall
 * @param false_min initial, false minimum
 * @param true_min final, true minimum
 */
void BounceSolCalc::CalculatePTStrength()
{
  if (not percolation_temp_set)
  {
    throw std::runtime_error("Phase transition strength cannot be calculated "
                             "because percolation temperature is not yet set.");
  }

  if (UserDefined_vwall >= 0)
    vwall = UserDefined_vwall;
  else
    vwall = 0.95; // Initial guess

  double old_alpha; // To keep track of the convergence

  for (int c = 0; c < 20; c++)
  {
    // Use recurvie method to find solution of
    // alpha = alpha(T_*)
    // T_* =  T_*(alpha, v_wall)
    // v_wall = v_wall(alpha, T_*)
    // Should converge quickly, if fails use default value of v_wall = 0.95
    old_alpha = alpha;
    CalculatePercolationTemp();
    Minimum true_min =
        coex.vec_phases[indexTrueCandidatePhase].Get(GetPercolationTemp());
    Minimum false_min = coex.vec_phases[0].Get(GetPercolationTemp());

    double Vi =
        false_min
            .potential; // potential at false vacuum and percolation temperature
    double Vf =
        true_min
            .potential; // potential at true vacuum and percolation temperature
    double dTVi = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(false_min.point),
        Tperc,
        -1); // temperature-derivative at false vacuum
    double dTVf = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(true_min.point),
        Tperc,
        -1); // temperature-derivative at true vacuum const

    double rho_gam =
        this->GetGstar() * std::pow(M_PI, 2) / 30 * std::pow(Tperc, 4);
    alpha = 1 / rho_gam * (Vi - Vf - Tperc / 4. * (dTVi - dTVf));
    CalculateWallVelocity(false_min, true_min);
    if (abs(alpha / old_alpha - 1) < 1e-7) return; // Found a solution
  }
  if (UserDefined_vwall > 0)
  {
    // Calculation failed, abort.
    BSMPT::Logger::Write(BSMPT::LoggingLevel::GWDetailed,
                         "Strenght of the transitions, alpha, could not be "
                         "calculated. vwall = " +
                             std::to_string(UserDefined_vwall));
    return;
  }
  // We could not find the solution for the system. use default value of .95
  // instead
  BSMPT::Logger::Write(
      BSMPT::LoggingLevel::GWDetailed,
      "v_wall could not be calculated. Using 0.95 for v_wall.");
  UserDefined_vwall = 0.95;
  // Call this function recursively
  CalculatePTStrength();
  return;
}

/**
 * @brief Calculate wall velocity
 * @param false_min initial, false minimum
 * @param true_min final, true minimum
 */
void BounceSolCalc::CalculateWallVelocity(const Minimum &false_min,
                                          const Minimum &true_min)
{
  if (not percolation_temp_set)
  {
    throw std::runtime_error("Wall velocity cannot be calculated "
                             "because percolation temperature is not yet set.");
  }

  // User defined
  if (UserDefined_vwall >= 0) vwall = UserDefined_vwall;
  if (UserDefined_vwall == -1)
  {
    // vwall https://arxiv.org/abs/2210.16305
    double Vi =
        false_min
            .potential; // potential at false vacuum and percolation temperature
    double Vf =
        true_min
            .potential; // potential at true vacuum and percolation temperature
    double rho_gam =
        this->GetGstar() * std::pow(M_PI, 2) / 30 * std::pow(Tperc, 4);

    double v_ChapmanJouget =
        1. / (1 + alpha) *
        (Csound + std::sqrt(std::pow(alpha, 2) + 2. / 3 * alpha));

    // Candidate wall velocity
    vwall = std::sqrt((Vi - Vf) / (alpha * rho_gam));
    // If cancidate is bigger than chapman jouget velocity, v = 1
    if (vwall > v_ChapmanJouget) vwall = 1;
  }
  if (UserDefined_vwall == -2)
  {
    // Upper bound implementation https://arxiv.org/abs/2305.02357
    double v_ChapmanJouget =
        1. / (1 + alpha) *
        (Csound + std::sqrt(std::pow(alpha, 2) + 2. / 3 * alpha));

    double dTVi = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(false_min.point),
        Tperc,
        -1); // temperature-derivative at false vacuum
    double dTVf = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(true_min.point),
        Tperc,
        -1); // temperature-derivative at true vacuum const

    double psi = dTVf / dTVi;
    double a   = 0.2233;
    double b   = 1.704;
    double p   = -3.433;

    vwall = std::pow(
        pow(abs((3 * alpha + psi - 1) / (2 * (2 - 3 * psi + std::pow(psi, 3)))),
            p / 2) +
            std::pow(
                abs(v_ChapmanJouget * (1 - a * std::pow(1 - psi, b) / alpha)),
                p / 2),
        1 / p);
  }
}

double BounceSolCalc::GetBounceSol(const double &Temp)
{
  return S3ofT_spline(Temp);
}

/**
 * @brief action_ratio friend to define input of numerical derivative in
 * calculation of inverse time scale
 */
double action_ratio(double Temp, void *params)
{
  struct BounceSolCalc &obj = *static_cast<BounceSolCalc *>(params);
  double func               = obj.GetBounceSol(Temp) / Temp;
  return func;
}

/**
 * @brief Nderive_BounceRatio Numerical derivative for the inverse time scale
 * calculation
 * @param obj Class reference to pass all needed parameters
 * @return Numerical value of derivative and absolute error
 */
struct resultErrorPair Nderive_BounceRatio(BounceSolCalc &obj)
{
  double step_size = 1e-8;

  gsl_function F;
  F.function = &action_ratio;
  F.params   = static_cast<void *>(&obj);

  struct resultErrorPair res;

  gsl_deriv_central(
      &F, obj.GetPercolationTemp(), step_size, &res.result, &res.error);

  return res;
}

/**
 * @brief Get inverse time scale of phase transition
 */
void BounceSolCalc::CalculateInvTimeScale()
{
  if (not percolation_temp_set)
  {
    throw std::runtime_error("Phase transition strength cannot be calculated "
                             "because percolation temperature is not yet set.");
  }

  struct resultErrorPair res = Nderive_BounceRatio(*this);
  this->betaH                = this->GetPercolationTemp() * res.result;
}

/**
 * @brief Get inverse time scale of phase transition
 */
double BounceSolCalc::GetInvTimeScale()
{
  if (this->betaH == -1) CalculateInvTimeScale();
  return this->betaH;
}

} // namespace GravitationalWaves
} // namespace BSMPT
