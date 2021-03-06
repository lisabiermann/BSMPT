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
 * Implementation of the model CP in the Dark as given in the manual.
 *
 * 
 * 
 * * \f$-\mathcal{L}_S = m_{11}^2|\Phi_1|^2+m_{22}^2|\Phi_2|^2+\frac{m_S^2}{2}\Phi_S^2+(\Re(A)+i\Im(A))\Phi_1^\dagger\Phi_2\Phi_S+(\Re(A)-i \Im(A))\Phi_S \Phi_2^\dagger \Phi_1 +\frac{\lambda_1}{2}|\Phi_1|^4+\frac{\lambda_2}{2}|\Phi_2|^4+\lambda_3|\Phi_1|^2|\Phi_2|^2+\lambda_4|\Phi_1^\dagger\Phi_2|^2+\frac{\lambda_5}{2}[(\Phi_1^\dagger\Phi_2)^2+(\Phi_2^\dagger\Phi_1)^2]+\frac{\lambda_6}{8}\Phi_S^4+\frac{\lambda_7}{2}|\Phi_1|^2\Phi_S^2+\frac{\lambda_8}{2}|\Phi_2|^2\Phi_S^2\f$
 * * \f$-\mathcal{L}_{Y,\text{lepton}} = y_\tau\bar{L}_{\tau}\Phi_1\tau_R+y_\mu\bar{L}_{\mu}\Phi_1\mu_R+y_e\bar{L}_{e}\Phi_1 e_R\\ \f$
 * with \f$\bar{L}_\tau = \overline{\begin{pmatrix}\tau_L\\\nu_{\tau_L}\end{pmatrix}} =\begin{pmatrix}\nu_{\tau_L}&\tau_L\end{pmatrix}\f$ for the third lepton generation (first and second work the same way)
 * * \f$ \mathcal{L}_{Y,\text{quark}} = \bar{D}_R\text{diag}(y_d,y_s,y_b)V^\dagger \Phi_1^{+*} U_L+\bar{D}_R\text{diag}(y_d,y_s,y_b)D_L\Phi_1^{0*}-\bar{U}_R\text{diag}(y_u,y_c,y_t)\Phi_1^{+}V D_L+\bar{U}_R\text{diag}(y_u,y_c,y_t)\Phi_1^0 U_L\\ \f$
 * with \f$ U = \begin{pmatrix}u\\c\\t\end{pmatrix} \f$ and \f$ D = \begin{pmatrix}d\\s\\b\end{pmatrix} \f$. \f$ V\f$ denotes the Cabibbo-Kobayashi-Maskawa matrix.
 * * \f$\mathcal{L}_G = (D^\mu\Phi_1)^\dagger(D_\mu\Phi_1)+(D^\mu\Phi_2)^\dagger(D_\mu\Phi_2)\\ \f$
 * using \f$ D_\mu=-i\frac{C_g}{2}\vec{\sigma}\vec{W}_\mu-i\frac{C_{gs}}{2}\sigma_0 B_\mu \f$
 *
 * with
 *
 * * \f$ \Phi_1 = \begin{pmatrix} \Phi_1^+ \\  \Phi_1^0 \end{pmatrix} =  \frac{1}{\sqrt{2}} \begin{pmatrix} \rho_1+i\eta_1\\  \zeta_1+\omega_1+i\Psi_1 \end{pmatrix}\,, \f$
 * * \f$ \Phi_2 = \begin{pmatrix} \Phi_2^+ \\  \Phi_2^0 \end{pmatrix} =  \frac{1}{\sqrt{2}} \begin{pmatrix} \rho_2 + \omega_\text{CB} + i \eta_2 \\ \zeta_2 + \omega_2 + i (\Psi_2 + \omega_\text{CP}) \end{pmatrix}\,, \f$
 * * \f$ \Phi_S = \zeta_S+\omega_S \f$
 *
 * using the configurations
 *
 * * \f$\text{higgsbase} = \{\rho_1,\eta_1,\rho_2,\eta_2,\zeta_1,\Psi_1,\zeta_2,\Psi_2,\zeta_S\} \,, \f$
 * * \f$\text{higgsvevFiniteTemp} = \{0,0,\omega_\text{CB},0,\omega_1,0,\omega_2,\omega_\text{CP},\omega_S\} \,, \f$
 * * \f$\text{higgsvevZeroTemp} = \{0,0,0,0,v_1,0,0,0,0\} \f$
 * 
 * and 
 *
 * * \f$ \langle \Phi_1 \rangle\big|_{T=0} = \frac{1}{\sqrt{2}} \begin{pmatrix}0\\v_1\end{pmatrix}\,, \f$
 * * \f$ \langle \Phi_2 \rangle\big|_{T=0} = \begin{pmatrix}0\\0\end{pmatrix}\,, \f$
 * * \f$ \langle \Phi_S \rangle\big|_{T=0} = 0\,. \f$
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

