// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "R2HDM.h" 
Compare_R2HDM::Compare_R2HDM()
{
  std::size_t NHiggs = 8;
  CheckTripleTree = Matrix3D{NHiggs, Matrix2D{NHiggs,   std::vector<double>(NHiggs, 0)}};
  CheckTripleCW =   Matrix3D{NHiggs, Matrix2D{NHiggs,   std::vector<double>(NHiggs, 0)}};
  CheckTripleCT =   Matrix3D{NHiggs, Matrix2D{NHiggs,   std::vector<double>(NHiggs, 0)}};
  EWPTPerSetting[4].Tc = 108.838;
  EWPTPerSetting[4].vc = 238.589;
  EWPTPerSetting[4].EWMinimum.push_back(0);
  EWPTPerSetting[4].EWMinimum.push_back(-54.1076);
  EWPTPerSetting[4].EWMinimum.push_back(-232.373);
  EWPTPerSetting[4].EWMinimum.push_back(0);
  EWPTPerSetting[1].Tc = 155.283;
  EWPTPerSetting[1].vc = 190.738;
  EWPTPerSetting[1].EWMinimum.push_back(0);
  EWPTPerSetting[1].EWMinimum.push_back(54.2467);
  EWPTPerSetting[1].EWMinimum.push_back(182.861);
  EWPTPerSetting[1].EWMinimum.push_back(0);
  EWPTPerSetting[5].Tc = 155.283;
  EWPTPerSetting[5].vc = 190.738;
  EWPTPerSetting[5].EWMinimum.push_back(0);
  EWPTPerSetting[5].EWMinimum.push_back(54.2467);
  EWPTPerSetting[5].EWMinimum.push_back(182.861);
  EWPTPerSetting[5].EWMinimum.push_back(1.81177e-05);
  EWPTPerSetting[2].Tc = 155.283;
  EWPTPerSetting[2].vc = 190.738;
  EWPTPerSetting[2].EWMinimum.push_back(1.66586e-05);
  EWPTPerSetting[2].EWMinimum.push_back(-54.2468);
  EWPTPerSetting[2].EWMinimum.push_back(-182.861);
  EWPTPerSetting[2].EWMinimum.push_back(0);
  EWPTPerSetting[6].Tc = 155.283;
  EWPTPerSetting[6].vc = 190.738;
  EWPTPerSetting[6].EWMinimum.push_back(2.05901e-05);
  EWPTPerSetting[6].EWMinimum.push_back(54.2467);
  EWPTPerSetting[6].EWMinimum.push_back(182.861);
  EWPTPerSetting[6].EWMinimum.push_back(-1.29244e-05);
  EWPTPerSetting[3].Tc = 155.283;
  EWPTPerSetting[3].vc = 190.738;
  EWPTPerSetting[3].EWMinimum.push_back(0);
  EWPTPerSetting[3].EWMinimum.push_back(54.2467);
  EWPTPerSetting[3].EWMinimum.push_back(182.861);
  EWPTPerSetting[3].EWMinimum.push_back(0);
  EWPTPerSetting[7].Tc = 155.283;
  EWPTPerSetting[7].vc = 190.738;
  EWPTPerSetting[7].EWMinimum.push_back(0);
  EWPTPerSetting[7].EWMinimum.push_back(54.2467);
  EWPTPerSetting[7].EWMinimum.push_back(182.861);
  EWPTPerSetting[7].EWMinimum.push_back(0);
  CheckTripleTree.at(0).at(0).at(6) = -61.4198;
  CheckTripleCT.at(0).at(0).at(6) = 77.8048;
  CheckTripleCW.at(0).at(0).at(6) = -77.8048;
  CheckTripleTree.at(0).at(0).at(7) = 45.1709;
  CheckTripleCT.at(0).at(0).at(7) = -24.2999;
  CheckTripleCW.at(0).at(0).at(7) = 24.2999;
  CheckTripleTree.at(0).at(1).at(5) = -2.22045e-16;
  CheckTripleCT.at(0).at(1).at(5) = 1.38778e-17;
  CheckTripleCW.at(0).at(1).at(5) = 3.1225e-17;
  CheckTripleTree.at(0).at(2).at(6) = -174.929;
  CheckTripleCT.at(0).at(2).at(6) = 5.08295;
  CheckTripleCW.at(0).at(2).at(6) = -5.08295;
  CheckTripleTree.at(0).at(2).at(7) = -549.757;
  CheckTripleCT.at(0).at(2).at(7) = 5.48284;
  CheckTripleCW.at(0).at(2).at(7) = -5.48284;
  CheckTripleCW.at(0).at(3).at(4) = -1.11022e-16;
  CheckTripleTree.at(0).at(3).at(5) = 44.255;
  CheckTripleCT.at(0).at(3).at(5) = -3.27917;
  CheckTripleCW.at(0).at(3).at(5) = 3.27917;
  CheckTripleTree.at(0).at(4).at(1) = -4.44089e-16;
  CheckTripleCT.at(0).at(4).at(1) = 2.77556e-17;
  CheckTripleCW.at(0).at(4).at(1) = -2.77556e-17;
  CheckTripleCW.at(0).at(4).at(3) = -1.11022e-16;
  CheckTripleTree.at(0).at(5).at(1) = -2.22045e-16;
  CheckTripleCT.at(0).at(5).at(1) = 1.38778e-17;
  CheckTripleCW.at(0).at(5).at(1) = 3.1225e-17;
  CheckTripleTree.at(0).at(5).at(3) = 44.255;
  CheckTripleCT.at(0).at(5).at(3) = -3.27917;
  CheckTripleCW.at(0).at(5).at(3) = 3.27917;
  CheckTripleTree.at(0).at(6).at(0) = -61.4198;
  CheckTripleCT.at(0).at(6).at(0) = 77.8048;
  CheckTripleCW.at(0).at(6).at(0) = -77.8048;
  CheckTripleTree.at(0).at(6).at(2) = -174.929;
  CheckTripleCT.at(0).at(6).at(2) = 5.08295;
  CheckTripleCW.at(0).at(6).at(2) = -5.08295;
  CheckTripleTree.at(0).at(7).at(0) = 45.1709;
  CheckTripleCT.at(0).at(7).at(0) = -24.2999;
  CheckTripleCW.at(0).at(7).at(0) = 24.2999;
  CheckTripleTree.at(0).at(7).at(2) = -549.757;
  CheckTripleCT.at(0).at(7).at(2) = 5.48284;
  CheckTripleCW.at(0).at(7).at(2) = -5.48284;
  CheckTripleTree.at(1).at(0).at(5) = 2.22045e-16;
  CheckTripleCT.at(1).at(0).at(5) = -1.38778e-17;
  CheckTripleCW.at(1).at(0).at(5) = -3.1225e-17;
  CheckTripleTree.at(1).at(1).at(6) = -61.4198;
  CheckTripleCT.at(1).at(1).at(6) = 77.8048;
  CheckTripleCW.at(1).at(1).at(6) = -77.8048;
  CheckTripleTree.at(1).at(1).at(7) = 45.1709;
  CheckTripleCT.at(1).at(1).at(7) = -24.2999;
  CheckTripleCW.at(1).at(1).at(7) = 24.2999;
  CheckTripleCW.at(1).at(2).at(4) = 1.11022e-16;
  CheckTripleTree.at(1).at(2).at(5) = -44.255;
  CheckTripleCT.at(1).at(2).at(5) = 3.27917;
  CheckTripleCW.at(1).at(2).at(5) = -3.27917;
  CheckTripleTree.at(1).at(3).at(6) = -174.929;
  CheckTripleCT.at(1).at(3).at(6) = 5.08295;
  CheckTripleCW.at(1).at(3).at(6) = -5.08295;
  CheckTripleTree.at(1).at(3).at(7) = -549.757;
  CheckTripleCT.at(1).at(3).at(7) = 5.48284;
  CheckTripleCW.at(1).at(3).at(7) = -5.48284;
  CheckTripleTree.at(1).at(4).at(0) = 4.44089e-16;
  CheckTripleCT.at(1).at(4).at(0) = -2.77556e-17;
  CheckTripleCW.at(1).at(4).at(0) = 2.77556e-17;
  CheckTripleCW.at(1).at(4).at(2) = 1.11022e-16;
  CheckTripleTree.at(1).at(5).at(0) = 2.22045e-16;
  CheckTripleCT.at(1).at(5).at(0) = -1.38778e-17;
  CheckTripleCW.at(1).at(5).at(0) = -3.1225e-17;
  CheckTripleTree.at(1).at(5).at(2) = -44.255;
  CheckTripleCT.at(1).at(5).at(2) = 3.27917;
  CheckTripleCW.at(1).at(5).at(2) = -3.27917;
  CheckTripleTree.at(1).at(6).at(1) = -61.4198;
  CheckTripleCT.at(1).at(6).at(1) = 77.8048;
  CheckTripleCW.at(1).at(6).at(1) = -77.8048;
  CheckTripleTree.at(1).at(6).at(3) = -174.929;
  CheckTripleCT.at(1).at(6).at(3) = 5.08295;
  CheckTripleCW.at(1).at(6).at(3) = -5.08295;
  CheckTripleTree.at(1).at(7).at(1) = 45.1709;
  CheckTripleCT.at(1).at(7).at(1) = -24.2999;
  CheckTripleCW.at(1).at(7).at(1) = 24.2999;
  CheckTripleTree.at(1).at(7).at(3) = -549.757;
  CheckTripleCT.at(1).at(7).at(3) = 5.48284;
  CheckTripleCW.at(1).at(7).at(3) = -5.48284;
  CheckTripleTree.at(2).at(0).at(6) = -174.929;
  CheckTripleCT.at(2).at(0).at(6) = 5.08295;
  CheckTripleCW.at(2).at(0).at(6) = -5.08295;
  CheckTripleTree.at(2).at(0).at(7) = -549.757;
  CheckTripleCT.at(2).at(0).at(7) = 5.48284;
  CheckTripleCW.at(2).at(0).at(7) = -5.48284;
  CheckTripleTree.at(2).at(1).at(4) = 5.55112e-17;
  CheckTripleCT.at(2).at(1).at(4) = 3.46945e-18;
  CheckTripleCW.at(2).at(1).at(4) = 1.14492e-16;
  CheckTripleTree.at(2).at(1).at(5) = -44.255;
  CheckTripleCT.at(2).at(1).at(5) = 3.27917;
  CheckTripleCW.at(2).at(1).at(5) = -3.27917;
  CheckTripleTree.at(2).at(2).at(6) = -1307.14;
  CheckTripleCT.at(2).at(2).at(6) = 97.6374;
  CheckTripleCW.at(2).at(2).at(6) = -156.454;
  CheckTripleTree.at(2).at(2).at(7) = 449.738;
  CheckTripleCT.at(2).at(2).at(7) = 15.6005;
  CheckTripleCW.at(2).at(2).at(7) = 79.7108;
  CheckTripleTree.at(2).at(3).at(5) = 2.22045e-16;
  CheckTripleCT.at(2).at(3).at(5) = -1.38778e-17;
  CheckTripleCW.at(2).at(3).at(5) = -3.1225e-17;
  CheckTripleCT.at(2).at(4).at(1) = 6.93889e-18;
  CheckTripleCW.at(2).at(4).at(1) = 1.11022e-16;
  CheckTripleTree.at(2).at(4).at(3) = -4.44089e-16;
  CheckTripleCT.at(2).at(4).at(3) = 2.77556e-17;
  CheckTripleCW.at(2).at(4).at(3) = -2.77556e-17;
  CheckTripleTree.at(2).at(5).at(1) = -44.255;
  CheckTripleCT.at(2).at(5).at(1) = 3.27917;
  CheckTripleCW.at(2).at(5).at(1) = -3.27917;
  CheckTripleTree.at(2).at(5).at(3) = 2.77556e-16;
  CheckTripleCT.at(2).at(5).at(3) = -1.73472e-17;
  CheckTripleCW.at(2).at(5).at(3) = -2.77556e-17;
  CheckTripleTree.at(2).at(6).at(0) = -174.929;
  CheckTripleCT.at(2).at(6).at(0) = 5.08295;
  CheckTripleCW.at(2).at(6).at(0) = -5.08295;
  CheckTripleTree.at(2).at(6).at(2) = -1307.14;
  CheckTripleCT.at(2).at(6).at(2) = 97.6374;
  CheckTripleCW.at(2).at(6).at(2) = -156.454;
  CheckTripleTree.at(2).at(7).at(0) = -549.757;
  CheckTripleCT.at(2).at(7).at(0) = 5.48284;
  CheckTripleCW.at(2).at(7).at(0) = -5.48284;
  CheckTripleTree.at(2).at(7).at(2) = 449.738;
  CheckTripleCT.at(2).at(7).at(2) = 15.6005;
  CheckTripleCW.at(2).at(7).at(2) = 79.7108;
  CheckTripleTree.at(3).at(0).at(4) = -5.55112e-17;
  CheckTripleCT.at(3).at(0).at(4) = -3.46945e-18;
  CheckTripleCW.at(3).at(0).at(4) = -1.14492e-16;
  CheckTripleTree.at(3).at(0).at(5) = 44.255;
  CheckTripleCT.at(3).at(0).at(5) = -3.27917;
  CheckTripleCW.at(3).at(0).at(5) = 3.27917;
  CheckTripleTree.at(3).at(1).at(6) = -174.929;
  CheckTripleCT.at(3).at(1).at(6) = 5.08295;
  CheckTripleCW.at(3).at(1).at(6) = -5.08295;
  CheckTripleTree.at(3).at(1).at(7) = -549.757;
  CheckTripleCT.at(3).at(1).at(7) = 5.48284;
  CheckTripleCW.at(3).at(1).at(7) = -5.48284;
  CheckTripleTree.at(3).at(2).at(5) = -2.22045e-16;
  CheckTripleCT.at(3).at(2).at(5) = 1.38778e-17;
  CheckTripleCW.at(3).at(2).at(5) = 3.1225e-17;
  CheckTripleTree.at(3).at(3).at(6) = -1307.14;
  CheckTripleCT.at(3).at(3).at(6) = 97.6374;
  CheckTripleCW.at(3).at(3).at(6) = -156.454;
  CheckTripleTree.at(3).at(3).at(7) = 449.738;
  CheckTripleCT.at(3).at(3).at(7) = 15.6005;
  CheckTripleCW.at(3).at(3).at(7) = 79.7108;
  CheckTripleCT.at(3).at(4).at(0) = -6.93889e-18;
  CheckTripleCW.at(3).at(4).at(0) = -1.11022e-16;
  CheckTripleTree.at(3).at(4).at(2) = 4.44089e-16;
  CheckTripleCT.at(3).at(4).at(2) = -2.77556e-17;
  CheckTripleCW.at(3).at(4).at(2) = 2.77556e-17;
  CheckTripleTree.at(3).at(5).at(0) = 44.255;
  CheckTripleCT.at(3).at(5).at(0) = -3.27917;
  CheckTripleCW.at(3).at(5).at(0) = 3.27917;
  CheckTripleTree.at(3).at(5).at(2) = -2.77556e-16;
  CheckTripleCT.at(3).at(5).at(2) = 1.73472e-17;
  CheckTripleCW.at(3).at(5).at(2) = 2.77556e-17;
  CheckTripleTree.at(3).at(6).at(1) = -174.929;
  CheckTripleCT.at(3).at(6).at(1) = 5.08295;
  CheckTripleCW.at(3).at(6).at(1) = -5.08295;
  CheckTripleTree.at(3).at(6).at(3) = -1307.14;
  CheckTripleCT.at(3).at(6).at(3) = 97.6374;
  CheckTripleCW.at(3).at(6).at(3) = -156.454;
  CheckTripleTree.at(3).at(7).at(1) = -549.757;
  CheckTripleCT.at(3).at(7).at(1) = 5.48284;
  CheckTripleCW.at(3).at(7).at(1) = -5.48284;
  CheckTripleTree.at(3).at(7).at(3) = 449.738;
  CheckTripleCT.at(3).at(7).at(3) = 15.6005;
  CheckTripleCW.at(3).at(7).at(3) = 79.7108;
  CheckTripleTree.at(4).at(0).at(1) = -4.44089e-16;
  CheckTripleCT.at(4).at(0).at(1) = 2.77556e-17;
  CheckTripleCW.at(4).at(0).at(1) = -2.77556e-17;
  CheckTripleCW.at(4).at(0).at(3) = -1.11022e-16;
  CheckTripleTree.at(4).at(1).at(0) = 4.44089e-16;
  CheckTripleCT.at(4).at(1).at(0) = -2.77556e-17;
  CheckTripleCW.at(4).at(1).at(0) = 2.77556e-17;
  CheckTripleCW.at(4).at(1).at(2) = 1.11022e-16;
  CheckTripleTree.at(4).at(2).at(1) = 2.22045e-16;
  CheckTripleCT.at(4).at(2).at(1) = -1.38778e-17;
  CheckTripleCW.at(4).at(2).at(1) = 7.97973e-17;
  CheckTripleTree.at(4).at(2).at(3) = -4.44089e-16;
  CheckTripleCT.at(4).at(2).at(3) = 2.77556e-17;
  CheckTripleCW.at(4).at(2).at(3) = -2.77556e-17;
  CheckTripleTree.at(4).at(3).at(0) = -2.22045e-16;
  CheckTripleCT.at(4).at(3).at(0) = 1.38778e-17;
  CheckTripleCW.at(4).at(3).at(0) = -7.97973e-17;
  CheckTripleTree.at(4).at(3).at(2) = 4.44089e-16;
  CheckTripleCT.at(4).at(3).at(2) = -2.77556e-17;
  CheckTripleCW.at(4).at(3).at(2) = 2.77556e-17;
  CheckTripleTree.at(4).at(4).at(6) = -61.4198;
  CheckTripleCT.at(4).at(4).at(6) = 77.8048;
  CheckTripleCW.at(4).at(4).at(6) = -77.8048;
  CheckTripleTree.at(4).at(4).at(7) = 45.1709;
  CheckTripleCT.at(4).at(4).at(7) = -24.2999;
  CheckTripleCW.at(4).at(4).at(7) = 24.2999;
  CheckTripleTree.at(4).at(5).at(6) = 163.564;
  CheckTripleCT.at(4).at(5).at(6) = -4.24085;
  CheckTripleCW.at(4).at(5).at(6) = 4.24085;
  CheckTripleTree.at(4).at(5).at(7) = 506.986;
  CheckTripleCT.at(4).at(5).at(7) = -2.31364;
  CheckTripleCW.at(4).at(5).at(7) = 2.31364;
  CheckTripleTree.at(4).at(6).at(4) = -61.4198;
  CheckTripleCT.at(4).at(6).at(4) = 77.8048;
  CheckTripleCW.at(4).at(6).at(4) = -77.8048;
  CheckTripleTree.at(4).at(6).at(5) = 163.564;
  CheckTripleCT.at(4).at(6).at(5) = -4.24085;
  CheckTripleCW.at(4).at(6).at(5) = 4.24085;
  CheckTripleTree.at(4).at(7).at(4) = 45.1709;
  CheckTripleCT.at(4).at(7).at(4) = -24.2999;
  CheckTripleCW.at(4).at(7).at(4) = 24.2999;
  CheckTripleTree.at(4).at(7).at(5) = 506.986;
  CheckTripleCT.at(4).at(7).at(5) = -2.31364;
  CheckTripleCW.at(4).at(7).at(5) = 2.31364;
  CheckTripleTree.at(5).at(0).at(1) = 5.55112e-17;
  CheckTripleCT.at(5).at(0).at(1) = -3.46945e-18;
  CheckTripleCW.at(5).at(0).at(1) = 3.46945e-18;
  CheckTripleTree.at(5).at(0).at(3) = 44.255;
  CheckTripleCT.at(5).at(0).at(3) = -3.27917;
  CheckTripleCW.at(5).at(0).at(3) = 3.27917;
  CheckTripleTree.at(5).at(1).at(0) = -5.55112e-17;
  CheckTripleCT.at(5).at(1).at(0) = 3.46945e-18;
  CheckTripleCW.at(5).at(1).at(0) = -3.46945e-18;
  CheckTripleTree.at(5).at(1).at(2) = -44.255;
  CheckTripleCT.at(5).at(1).at(2) = 3.27917;
  CheckTripleCW.at(5).at(1).at(2) = -3.27917;
  CheckTripleTree.at(5).at(2).at(1) = -44.255;
  CheckTripleCT.at(5).at(2).at(1) = 3.27917;
  CheckTripleCW.at(5).at(2).at(1) = -3.27917;
  CheckTripleTree.at(5).at(2).at(3) = 5.55112e-17;
  CheckTripleCT.at(5).at(2).at(3) = -3.46945e-18;
  CheckTripleCW.at(5).at(2).at(3) = 3.46945e-18;
  CheckTripleTree.at(5).at(3).at(0) = 44.255;
  CheckTripleCT.at(5).at(3).at(0) = -3.27917;
  CheckTripleCW.at(5).at(3).at(0) = 3.27917;
  CheckTripleTree.at(5).at(3).at(2) = -5.55112e-17;
  CheckTripleCT.at(5).at(3).at(2) = 3.46945e-18;
  CheckTripleCW.at(5).at(3).at(2) = -3.46945e-18;
  CheckTripleTree.at(5).at(4).at(6) = 163.564;
  CheckTripleCT.at(5).at(4).at(6) = -4.24085;
  CheckTripleCW.at(5).at(4).at(6) = 4.24085;
  CheckTripleTree.at(5).at(4).at(7) = 506.986;
  CheckTripleCT.at(5).at(4).at(7) = -2.31364;
  CheckTripleCW.at(5).at(4).at(7) = 2.31364;
  CheckTripleTree.at(5).at(5).at(6) = -1221.6;
  CheckTripleCT.at(5).at(5).at(6) = 91.299;
  CheckTripleCW.at(5).at(5).at(6) = -137.072;
  CheckTripleTree.at(5).at(5).at(7) = 427.008;
  CheckTripleCT.at(5).at(5).at(7) = 17.2847;
  CheckTripleCW.at(5).at(5).at(7) = 71.1447;
  CheckTripleTree.at(5).at(6).at(4) = 163.564;
  CheckTripleCT.at(5).at(6).at(4) = -4.24085;
  CheckTripleCW.at(5).at(6).at(4) = 4.24085;
  CheckTripleTree.at(5).at(6).at(5) = -1221.6;
  CheckTripleCT.at(5).at(6).at(5) = 91.299;
  CheckTripleCW.at(5).at(6).at(5) = -137.072;
  CheckTripleTree.at(5).at(7).at(4) = 506.986;
  CheckTripleCT.at(5).at(7).at(4) = -2.31364;
  CheckTripleCW.at(5).at(7).at(4) = 2.31364;
  CheckTripleTree.at(5).at(7).at(5) = 427.008;
  CheckTripleCT.at(5).at(7).at(5) = 17.2847;
  CheckTripleCW.at(5).at(7).at(5) = 71.1447;
  CheckTripleTree.at(6).at(0).at(0) = -61.4198;
  CheckTripleCT.at(6).at(0).at(0) = 77.8048;
  CheckTripleCW.at(6).at(0).at(0) = -77.8048;
  CheckTripleTree.at(6).at(0).at(2) = -174.929;
  CheckTripleCT.at(6).at(0).at(2) = 5.08295;
  CheckTripleCW.at(6).at(0).at(2) = -5.08295;
  CheckTripleTree.at(6).at(1).at(1) = -61.4198;
  CheckTripleCT.at(6).at(1).at(1) = 77.8048;
  CheckTripleCW.at(6).at(1).at(1) = -77.8048;
  CheckTripleTree.at(6).at(1).at(3) = -174.929;
  CheckTripleCT.at(6).at(1).at(3) = 5.08295;
  CheckTripleCW.at(6).at(1).at(3) = -5.08295;
  CheckTripleTree.at(6).at(2).at(0) = -174.929;
  CheckTripleCT.at(6).at(2).at(0) = 5.08295;
  CheckTripleCW.at(6).at(2).at(0) = -5.08295;
  CheckTripleTree.at(6).at(2).at(2) = -1307.14;
  CheckTripleCT.at(6).at(2).at(2) = 97.6374;
  CheckTripleCW.at(6).at(2).at(2) = -156.454;
  CheckTripleTree.at(6).at(3).at(1) = -174.929;
  CheckTripleCT.at(6).at(3).at(1) = 5.08295;
  CheckTripleCW.at(6).at(3).at(1) = -5.08295;
  CheckTripleTree.at(6).at(3).at(3) = -1307.14;
  CheckTripleCT.at(6).at(3).at(3) = 97.6374;
  CheckTripleCW.at(6).at(3).at(3) = -156.454;
  CheckTripleTree.at(6).at(4).at(4) = -61.4198;
  CheckTripleCT.at(6).at(4).at(4) = 77.8048;
  CheckTripleCW.at(6).at(4).at(4) = -77.8048;
  CheckTripleTree.at(6).at(4).at(5) = 163.564;
  CheckTripleCT.at(6).at(4).at(5) = -4.24085;
  CheckTripleCW.at(6).at(4).at(5) = 4.24085;
  CheckTripleTree.at(6).at(5).at(4) = 163.564;
  CheckTripleCT.at(6).at(5).at(4) = -4.24085;
  CheckTripleCW.at(6).at(5).at(4) = 4.24085;
  CheckTripleTree.at(6).at(5).at(5) = -1221.6;
  CheckTripleCT.at(6).at(5).at(5) = 91.299;
  CheckTripleCW.at(6).at(5).at(5) = -137.072;
  CheckTripleTree.at(6).at(6).at(6) = -170.222;
  CheckTripleCT.at(6).at(6).at(6) = 229.769;
  CheckTripleCW.at(6).at(6).at(6) = -325.36;
  CheckTripleTree.at(6).at(6).at(7) = 30.1014;
  CheckTripleCT.at(6).at(6).at(7) = -23.3707;
  CheckTripleCW.at(6).at(6).at(7) = 65.133;
  CheckTripleTree.at(6).at(7).at(6) = 30.1014;
  CheckTripleCT.at(6).at(7).at(6) = -23.3707;
  CheckTripleCW.at(6).at(7).at(6) = 65.133;
  CheckTripleTree.at(6).at(7).at(7) = -156.506;
  CheckTripleCT.at(6).at(7).at(7) = 109.139;
  CheckTripleCW.at(6).at(7).at(7) = -96.4585;
  CheckTripleTree.at(7).at(0).at(0) = 45.1709;
  CheckTripleCT.at(7).at(0).at(0) = -24.2999;
  CheckTripleCW.at(7).at(0).at(0) = 24.2999;
  CheckTripleTree.at(7).at(0).at(2) = -549.757;
  CheckTripleCT.at(7).at(0).at(2) = 5.48284;
  CheckTripleCW.at(7).at(0).at(2) = -5.48284;
  CheckTripleTree.at(7).at(1).at(1) = 45.1709;
  CheckTripleCT.at(7).at(1).at(1) = -24.2999;
  CheckTripleCW.at(7).at(1).at(1) = 24.2999;
  CheckTripleTree.at(7).at(1).at(3) = -549.757;
  CheckTripleCT.at(7).at(1).at(3) = 5.48284;
  CheckTripleCW.at(7).at(1).at(3) = -5.48284;
  CheckTripleTree.at(7).at(2).at(0) = -549.757;
  CheckTripleCT.at(7).at(2).at(0) = 5.48284;
  CheckTripleCW.at(7).at(2).at(0) = -5.48284;
  CheckTripleTree.at(7).at(2).at(2) = 449.738;
  CheckTripleCT.at(7).at(2).at(2) = 15.6005;
  CheckTripleCW.at(7).at(2).at(2) = 79.7108;
  CheckTripleTree.at(7).at(3).at(1) = -549.757;
  CheckTripleCT.at(7).at(3).at(1) = 5.48284;
  CheckTripleCW.at(7).at(3).at(1) = -5.48284;
  CheckTripleTree.at(7).at(3).at(3) = 449.738;
  CheckTripleCT.at(7).at(3).at(3) = 15.6005;
  CheckTripleCW.at(7).at(3).at(3) = 79.7108;
  CheckTripleTree.at(7).at(4).at(4) = 45.1709;
  CheckTripleCT.at(7).at(4).at(4) = -24.2999;
  CheckTripleCW.at(7).at(4).at(4) = 24.2999;
  CheckTripleTree.at(7).at(4).at(5) = 506.986;
  CheckTripleCT.at(7).at(4).at(5) = -2.31364;
  CheckTripleCW.at(7).at(4).at(5) = 2.31364;
  CheckTripleTree.at(7).at(5).at(4) = 506.986;
  CheckTripleCT.at(7).at(5).at(4) = -2.31364;
  CheckTripleCW.at(7).at(5).at(4) = 2.31364;
  CheckTripleTree.at(7).at(5).at(5) = 427.008;
  CheckTripleCT.at(7).at(5).at(5) = 17.2847;
  CheckTripleCW.at(7).at(5).at(5) = 71.1447;
  CheckTripleTree.at(7).at(6).at(6) = 30.1014;
  CheckTripleCT.at(7).at(6).at(6) = -23.3707;
  CheckTripleCW.at(7).at(6).at(6) = 65.133;
  CheckTripleTree.at(7).at(6).at(7) = -156.506;
  CheckTripleCT.at(7).at(6).at(7) = 109.139;
  CheckTripleCW.at(7).at(6).at(7) = -96.4585;
  CheckTripleTree.at(7).at(7).at(6) = -156.506;
  CheckTripleCT.at(7).at(7).at(6) = 109.139;
  CheckTripleCW.at(7).at(7).at(6) = -96.4585;
  CheckTripleTree.at(7).at(7).at(7) = 450.506;
  CheckTripleCT.at(7).at(7).at(7) = 47.0721;
  CheckTripleCW.at(7).at(7).at(7) = 86.0667;
}