// Copyright (C) 2023 Lisa Biermann, João Viana
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#ifndef INCLUDE_BSMPT_GravitationalWaves_GWUtils_BOUNCESOL_CALC_H_
#define INCLUDE_BSMPT_GravitationalWaves_GWUtils_BOUNCESOL_CALC_H_

/**
 * @file bounce solution calculator class
 */

#include <BSMPT/gravitational_waves/gwUtils/bounce_solution/ActionCalculation.h>
#include <BSMPT/gravitational_waves/gwUtils/bounce_solution/spline.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <algorithm>             // std::max
#include <gsl/gsl_deriv.h>       // numerical derivative
#include <gsl/gsl_integration.h> // numerical integration

namespace BSMPT
{
namespace GravitationalWaves
{

/**
 * struct to store result and error
 */
struct resultErrorPair
{
  double result;
  double error;
};

class BounceSolCalc
{
protected:
  /**
   * @brief modelPointer for the used parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

public:
  /**
   * @brief epsilon of turbulence efficiency factor
   */
  double epstol = 0.1;
  /**
   * @brief wall velocity
   */
  double vwall = 0.95;

  /**
   * @brief set to true if nucleation temperature is set
   */
  bool nucleation_temp_set = false;

  /**
   * @brief set to true if percolation temperature is set
   */
  bool percolation_temp_set = false;

  /**
   * @brief critical temperature/highest temperature when transition can occur
   */
  double Tc;

  /**
   * @brief lowest temperature when a transition can occur
   */
  double Tm;

  /**
   * @brief nucleation temperature
   */
  double Tnucl;

  /**
   * @brief approximate nucleation temperature
   */
  double Tnucl_approx;

  /**
   * @brief percolation temperature
   */
  double Tperc;

  /**
   * @brief stored temperature
   */
  double store_Temp;

  /**
   * @brief PT strength
   */
  double alpha = -1;

  /**
   * @brief Inverse time scale \f$ \frac{\beta}{H} \f$
   *
   */
  double betaH = -1;

  /**
   * @brief number of effective degrees of freedom
   */
  double gstar;

  /**
   * @brief index of the true vacuum phase candidate in the coex list
   *
   */
  int indexTrueCandidatePhase;

  /**
   * @brief spline used to interpolate the action as a function of the
   * temperature
   *
   */
  tk::spline S3ofT_spline;

  /**
   * @brief Set of BounceActionInt objects with valid solutions.
   *
   */
  std::vector<BSMPT::GravitationalWaves::BounceActionInt> SolutionList;

  /**
   * @brief Storage of the tunneling rate per volume of the transition from
   * false to true vacuum
   * @param Temp temperature
   */
  double Tunneling_Rate(const double &Temp);
  /**
   * @brief Storage of the temperature-dependent Hubble rate
   * @param Temp temperature
   */
  double Hubble_Rate(const double &Temp);

  /**
   * @brief inner_integrand friend to define inner integrand of percolation
   * temperature integral
   */
  friend double inner_integrand(double var, void *params);

  /**
   * @brief outer_integrand friend to define outer integrand of percolation
   * temperature integral
   */
  friend double outer_integrand(double var, void *params);

  /** Calculate euclidian action at temperature T*/
  double GetBounceSol(const double &Temp);

  /**
   * @brief action_ratio friend to define input of numerical derivative in
   * calculation of inverse time scale
   */
  friend double action_ratio(double var, void *params);

public:
  /**
   * @brief AbsErr absolute error for numerical integration
   */
  const double AbsErr = 0;

  /**
   * @brief RelErr relative error for numerical integration
   */
  const double RelErr = 1e-6;

  /**
   * @brief Set of coexisting phases for the corresponding gravitational waves.
   * We assume only 1 gravitational wave per coex phase.
   *
   */
  CoexPhases coex;

  /**
   * @brief Keeps track if bounce solver worked or failed
   *
   */
  bool failed = false;

  /**
   * @brief \f$ v_{\text{wall}}\f$ defined by the user as an input parameter.
   *  If \f$ v_{\text{wall}}\f = -1$ then we use the approximation coming from
   * https://arxiv.org/abs/2210.16305
   * If \f$ v_{\text{wall}}\f = -2$ then we use the upper bound from
   * https://arxiv.org/abs/2305.02357
   *
   */
  double UserDefined_vwall;

  /**
   * @brief Number of integration of the bounce
   *
   */
  int MaxPathIntegrations;

  /**
   * @brief Default constructor
   *
   */
  BounceSolCalc();

  /**
   * @brief Construct a new Bounce Sol Calc object. Used for testing
   * @param pointer_in model pointer
   */
  BounceSolCalc(std::shared_ptr<Class_Potential_Origin> &pointer_in);

  /**
   * @brief Construct a new Bounce Sol Calc object. This class takes as input a
   * set of coexisting phases and the corresponding true vacuum candidate
   *
   * @param pointer_in model pointer
   * @param coex_in set of coexisting phases
   * @param indexTrueCandidatePhase_in index of the true vacuum candidadte in
   * the coex phases list
   * @param UserDefined_vwall_in is the input value for v_wall. If = -1$ then we
   * use the approximation coming from https://arxiv.org/abs/2210.16305. If =
   * -2$ then we use the upper bound from https://arxiv.org/abs/2305.02357
   * @param UserDefined_epstol_in is the input value for epstol. If [0..1] set
   * to value, for -1 we use the upper bound from
   * https://arxiv.org/abs/1704.05871
   * @param MaxPathIntegrations_in max number of path integrations
   */
  BounceSolCalc(std::shared_ptr<Class_Potential_Origin> &pointer_in,
                CoexPhases &coex_in,
                int indexTrueCandidatePhase_in,
                double UserDefined_vwall_in,
                double UserDefined_epstol_in,
                int MaxPathIntegrations_in);

  /**
   * @brief Initially we have no idea where the transition can occur, therefore
   * we scan the complete temperature range
   *
   */
  void GWInitialScan();

  /**
   * @brief Calculate the euclidian action of the transition from index 0 of
   * CoexPhases coex to indexTrueCandidatePhase.
   *
   * @param T temperature
   */

  void CalculateActionAt(double T, bool smart = true);

  /**
   * @brief If solution were found by the GWInitialScan() then we scan
   * temperature range in the vicinity such that we are get a enough sample to
   * then do the extrapolation.
   *
   */
  void GWSecondaryScan();

  /**
   * @brief Do linear extrapolations to calculate action at higher temperatures
   *
   */
  void GWScanTowardsHighAction();

  /**
   * @brief Do linear extrapolations to calculate action at lower temperatures
   *
   */
  void GWScanTowardsLowAction();

  /**
   * @brief Set the Bounce Sol object
   *
   */
  void SetBounceSol();

  /**
   * @brief Get the bubble wall velocity
   * @return vb
   */
  double GetWallVelocity();

  /**
   * @brief Get epsturb
   * @return epsturb
   */
  double GetEpsTurb();

  /**
   * @brief SetGstar Set gstar
   */
  void SetGstar(const double &gstar_in);

  /**
   * @brief GetGstar Get gstar
   */
  double GetGstar();

  /**
   * @brief SetCriticalTemp Set critical temperature
   */
  void SetCriticalTemp(const double &T_in);

  /**
   * @brief GetCriticalTemp Get critical temperature
   */
  double GetCriticalTemp();

  /**
   * @brief SetStoredTemp Set stored temperature
   */
  void SetStoredTemp(const double &T_in);
  /**
   * @brief GetStoredTemp Get stored temperature
   */
  double GetStoredTemp();

  /**
   * @brief GetNucleationTemp Get nucleation temperature via exact method
   */
  double GetNucleationTemp();

  /**
   * @brief GetNucleationTempApprox Get nucleation temperature via approximate
   * method
   */
  double GetNucleationTempApprox();

  /**
   * @brief GetPercolationTemp Get percolation temperature
   */
  double GetPercolationTemp();

  /**
   * @brief GetPTStrength Get PT strength alpha
   */
  double GetPTStrength();

  /**
   * @brief CalcGstar Calculate the number of effective degrees of freedom as a
   * function of temperature, higgs masses are evaluated at global minimum
   * @param Temp temperature
   */
  void CalcGstar(const double &Temp);

  /**
   * @brief CalcGstarPureRad Calculate the number of effective degrees of
   * freedom assuming a purely radiative universe
   */
  void CalcGstarPureRad();

  /**
   * @brief Calculation of nucleation temperature
   */
  void CalculateNucleationTemp();

  /**
   * @brief Approximate calculation of nucleation temperature
   */
  void CalculateNucleationTempApprox();

  /**
   * @brief Calculation of percolation temperature
   */
  void CalculatePercolationTemp();

  /**
   * @brief Calculate phase transition strength alpha at percolation temperature
   */
  void CalculatePTStrength();

  /**
   * @brief Calculate wall velocity
   * @param false_min initial, false minimum
   * @param true_min final, true minimum
   */
  void CalculateWallVelocity(const Minimum &false_min, const Minimum &true_min);

  /**
   * @brief Calculate inverse time scale of phase transition
   */
  void CalculateInvTimeScale();

  /**
   * @brief Get inverse time scale of phase transition
   */
  double GetInvTimeScale();
};

/**
 * @brief Nintegrate_Inner Numerical integration of inner integral over inverse
 * Hubble rate for the percolation temperature calculation
 * @param obj Class reference to pass all needed parameters
 * @param T upper integration boundary
 * @return Numerical value of integral and absolute error
 */
struct resultErrorPair Nintegrate_Inner(BounceSolCalc &obj,
                                        const double &Tprime);

/**
 * @brief Nintegrate_Outer Numerical integration of outer integral for the
 * percolation temperature calculation
 * @param obj Class reference to pass all needed parameters
 * @return Numerical value of integral and absolute error
 */
struct resultErrorPair Nintegrate_Outer(BounceSolCalc &obj);

/**
 * @brief Nderive_BounceRatio Numerical derivative for the inverse time scale
 * calculation
 * @param obj Class reference to pass all needed parameters
 * @return Numerical value of derivative and absolute error
 */
struct resultErrorPair Nderive_BounceRatio(BounceSolCalc &obj);

/**
 * @brief speed of sound
 */
const double Csound = 0.5773502691896258; // 1/sqrt(3)

/**
 * @brief reduced Planck mass = MPl / (8 Pi)
 */
const double MPl = 2.4e18;

} // namespace GravitationalWaves
} // namespace BSMPT

#endif /* INCLUDE_BSMPT_GravitationalWaves_GWUtils_BOUNCESOL_CALC_H_ */