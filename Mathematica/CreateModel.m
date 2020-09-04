(* 
    CreateModel.m

    author: Lisa Biermann
    date of last modification: 04th September 2020

    A Mathematica package to calculate the tensors for the model CP in the Dark needed as input for BSMPT.

    CP in the Dark is a N2HDM (Next-to-2HDM) model with an extended scalar sector consisting of two Higgs doublets and a real singlet.
    The fermion and gauge sector are identical to the SM.

    An additional Z_2-symmetry is imposed on the model. As a result, only the first doublet couples to fermions. 
    And only this doublet acquieres a neutral non-zero vev which is assumed to be identical to the SM-Higgs-vev.
    Both Higgs doublets couple to gauge bosons.
*)

(* 
-------    Higgs potential  --------
definition of the scalar potential and calculation of scalar tensors
------------------------------------
*)

(* initialize file to save output *)
CurvatureFile = FileNameJoin[{NotebookDirectory[],"CurvatureFile.txt"}];
strm = OpenWrite[CurvatureFile, FormatType -> OutputForm];
AppendTo[$Output,strm];

higgsbase = {Gplus,Hplus,h,Gzero,rho,eta,s};
higgsvev = {0,0,v,0,0,0,0};
higgsvevFiniteTemp = {0,0,w,0,0,0,0};

higgsEVOzero = higgsbase + higgsvev;
higgsEVOfinite = higgsbase + higgsvevFiniteTemp;

nHiggs = Length[higgsbase];

par = {m11s,m22s,mSs,ReA,ImA,L1,L2,L3,L4,L5,L6,L7,L8};

(* define Higgs doublets and real singlet *)
Phi1zero = {{higgsEVOzero[[1]]} , {1/Sqrt[2]*(higgsEVOzero[[3]] + I * higgsEVOzero[[4]])}};
Phi2zero = {{higgsEVOzero[[2]]}, {1/Sqrt[2]*(higgsEVOzero[[5]] + I * higgsEVOzero[[6]])}};
PhiSzero = {{higgsEVOzero[[7]]}};

(* define field combinations FC for potential *)
$Assumptions = {Element[higgsEVOzero, Reals]};

Phi1zerodagger = Transpose[Refine[Conjugate[Phi1zero]]];
Phi2zerodagger = Transpose[Refine[Conjugate[Phi2zero]]];

FCm11s = Flatten[Phi1zerodagger.Phi1zero//Simplify]; 
FCm22s = Flatten[Phi2zerodagger.Phi2zero//Simplify];
FCmSs = Flatten[PhiSzero.PhiSzero];

FCARe = Flatten[Phi1zerodagger.Phi2zero.PhiSzero//Simplify];
FCAIm = Flatten[PhiSzero.Phi2zerodagger.Phi1zero//Simplify];

FCL1 = Flatten[FCm11s^2//Expand];
FCL2 = Flatten[FCm22s^2//Expand];
FCL3 = Flatten[FCm11s.FCm22s//Expand];
FCL4 = Flatten[Phi2zerodagger.Phi1zero.Phi1zerodagger.Phi2zero//Expand];
FCL5 = Flatten[((Phi1zerodagger.Phi2zero)^2+(Phi2zerodagger.Phi1zero)^2)//Expand];
FCL6 = Flatten[FCmSs^2];
FCL7 = Flatten[FCm11s.FCmSs//Expand];
FCL8 = Flatten[FCm22s.FCmSs//Expand];

(* define Higgs potential *)
VHiggs = par[[1]]*FCm11s+par[[2]]*FCm22s+par[[3]]/2*FCmSs+(par[[4]] + I*par[[5]])*FCARe+(par[[4]] - I*par[[5]])*FCAIm+par[[6]]/2*FCL1+par[[7]]/2*FCL2+par[[8]]*FCL3+par[[9]]*FCL4+par[[10]]/2*FCL5+par[[11]]/4*FCL6+par[[12]]/2*FCL7+par[[13]]/2*FCL8//Simplify;

(* calculate Higgs Curvature Tensors L_i,L_ij,L_ijk,L_ijkl *)
RepAllZero = {v->0,Gplus->0,Hplus->0,h->0,Gzero->0,rho->0,eta->0,s->0};

Print["----------------------------------------------------------------------"];
Print["                    Higgs potential"]
Print["----------------------------------------------------------------------"];

(* calculate L_i *)
nCount = 0;
For[i=1,i<=nHiggs,i++,
    tmp = D[VHiggs,higgsbase[[i]]]/.RepAllZero;
    If[tmp =!= {0}, 
        Print["Curvature_Higgs_L1[",i-1,"] = ", tmp[[1]]//InputForm,";"];
        nCount += 1;
    ];                                  (*print only if not zero*)
    If[i == nHiggs && nCount == 0,
        Print["Curvature_Higgs_L1[i] is ZERO for all entries."]
    ];                                  (*if whole tensor is zero*)
];

Print["----------------------------------------------------------------------"];

(* calculate L_ij *)
nCount = 0;
For[i=1,i<=nHiggs,i++,
    For[j=1,j<=nHiggs,j++,
        tmp = D[VHiggs,higgsbase[[i]],higgsbase[[j]]]/.RepAllZero;
        If[tmp =!= {0},
            Print["Curvature_Higgs_L2[",i-1,"][",j-1,"] = ",tmp[[1]]//InputForm,";"];
            nCount += 1;
        ];
    ];
    If[i == nHiggs && nCount == 0,
        Print["Curvature_Higgs_L2[i,j] is ZERO for all entries."];
    ];
];

Print["----------------------------------------------------------------------"];

(* calculate L_ijk *)
nCount = 0;
For[i=1,i<=nHiggs,i++,
    For[j=1,j<=nHiggs,j++,
        For[k=1,k<=nHiggs,k++,
            tmp = D[VHiggs,higgsbase[[i]],higgsbase[[j]],higgsbase[[k]]]/.RepAllZero;
            If[tmp =!= {0},
                Print["Curvature_Higgs_L3[",i-1,"][",j-1,"][",k-1,"] = ", tmp[[1]]//InputForm,";"];
                nCount += 1;
            ];
        ];
    ];
    If[i == nHiggs && nCount == 0,
        Print["Curvature_Higgs_L3[i,j,k] is ZERO for all entries."];
    ];
];

Print["----------------------------------------------------------------------"];

(* calculate L_ijkl *)
nCount = 0;
For[i=1,i<=nHiggs,i++,
    For[j=1,j<=nHiggs,j++,
        For[k=1,k<=nHiggs,k++,
            For[l=1,l<=nHiggs,l++,
                tmp = D[VHiggs,higgsbase[[i]],higgsbase[[j]],higgsbase[[k]],higgsbase[[l]]]/.RepAllZero//Simplify;
                If[tmp =!= {0},
                    Print["Curvature_Higgs_L4[",i-1,"][",j-1,"][",k-1,"][",l-1,"] = ", tmp[[1]]//InputForm,";"];
                    nCount += 1;
                ];
            ];
        ];
    ];
    If[i == nHiggs && nCount == 0,
        Print["Curvature_Higgs_L4[i,j,k,l] is ZERO for all entries."];
    ];
];

(* 
-------    Counterterm potential  --------
definition of the CT potential and derivation of explicit forms for counterterms from renormalization conditions
------------------------------------------
*)

Print["----------------------------------------------------------------------"];
Print["                    Counterterm potential"]
Print["----------------------------------------------------------------------"];

RepAllVEV = {Gplus->0,Hplus->0,h->0,Gzero->0,rho->0,eta->0,s->0};

CT = {dm11s,dm22s,dmSs,dReA,dImA,dL1,dL2,dL3,dL4,dL5,dL6,dL7,dL8,dT};

(* counterterm potential constructed from counterterms CT and field combinations FC*)
VCT = CT[[1]]*FCm11s+CT[[2]]*FCm22s+CT[[3]]/2*FCmSs+(CT[[4]] + I*CT[[5]])*FCARe+(CT[[4]] - I*CT[[5]])*FCAIm+CT[[6]]/2*FCL1+CT[[7]]/2*FCL2+CT[[8]]*FCL3+CT[[9]]*FCL4+CT[[10]]/2*FCL5+CT[[11]]/4*FCL6+CT[[12]]/2*FCL7+CT[[13]]/2*FCL8+CT[[14]]*(h+v)//Simplify;

(* store set of equations with one derivative *)
NCW = Table[0, {i, 1, nHiggs}];
Table[NCW[[i]]= -D[VCT,higgsbase[[i]]]/.RepAllVEV//Simplify,{i,1,nHiggs}]//MatrixForm;

(* store set of equations with two derivatives *)
HCW = Table[0,{i,1,nHiggs},{j,1,nHiggs}];
Table[HCW[[i]][[j]]= -D[VCT,higgsbase[[i]],higgsbase[[j]]]/.RepAllVEV//Simplify,{i,1,nHiggs},{j,1,nHiggs}]//MatrixForm;

(* get 10 equations for 14 counterterms *)
eq1 := NCW3 == NCW[[3]][[1]];
eq2 := HCW11 == HCW[[1]][[1]][[1]];
eq3 := HCW22 == HCW[[2]][[2]][[1]];
eq4 := HCW33 == HCW[[3]][[3]][[1]];
eq5 := HCW44 == HCW[[4]][[4]][[1]]; 
eq6 := HCW55 == HCW[[5]][[5]][[1]];
eq7 := HCW66 == HCW[[6]][[6]][[1]];
eq8 := HCW77 == HCW[[7]][[7]][[1]];
eq9 := HCW57 == HCW[[5]][[7]][[1]];
eq10 := HCW67 == HCW[[6]][[7]][[1]];

(* ReA *)
Solve[eq9,CT[[4]]]
(* ImA *)
Solve[eq10,CT[[5]]]

(* dT, dL1, dm11s *)
r1 = Solve[eq1, CT[[14]]];
r2 = Solve[eq2, CT[[1]]];
r3 = Solve[eq4 /. r2, CT[[6]]]
r4 = Solve[eq2 /. r3, CT[[1]]]
r5 = eq5 /. r3 /. r4 // Simplify  (* has to be fulfilled ??*)
r6 = r1 /. r3 /. r4 // Simplify

(* dL4, dL5, dm22s assuming dL3 = t *)
r7 = Solve[eq3 /. dL3 -> t, CT[[2]]]
r8 = Solve[eq6 /. r7 /. dL3 -> t, CT[[10]]];
r9 = Solve[eq7 /. r8 /. dL3 -> t, CT[[9]]]
r10 = r8 /. r9 // Simplify

(* dmSs assuming dL7 = u *)
r11 = Solve[eq8 /. dL7 -> u, CT[[3]]]

(* dL2, dL6, dL8 are not constrained!!! *)

(* calculate curvature Higgs CT *)
(* calculate L_i *)
nCount = 0;
For[i=1,i<=nHiggs,i++,
    tmp = D[VCT,higgsbase[[i]]]/.RepAllZero;
    If[tmp =!= {0}, 
        Print["Curvature_Higgs_CT_L1[",i-1,"] = ", tmp[[1]]//InputForm,";"];
        nCount += 1;
    ];                                  (*print only if not zero*)
    If[i == nHiggs && nCount == 0,
        Print["Curvature_Higgs_CT_L1[i] is ZERO for all entries."]
    ];                                  (*if whole tensor is zero*)
];

Print["----------------------------------------------------------------------"];

(* calculate L_ij *)
nCount = 0;
For[i=1,i<=nHiggs,i++,
    For[j=1,j<=nHiggs,j++,
        tmp = D[VCT,higgsbase[[i]],higgsbase[[j]]]/.RepAllZero;
        If[tmp =!= {0},
            Print["Curvature_Higgs_CT_L2[",i-1,"][",j-1,"] = ",tmp[[1]]//InputForm,";"];
            nCount += 1;
        ];
    ];
    If[i == nHiggs && nCount == 0,
        Print["Curvature_Higgs_CT_L2[i,j] is ZERO for all entries."];
    ];
];

Print["----------------------------------------------------------------------"];

(* calculate L_ijk *)
nCount = 0;
For[i=1,i<=nHiggs,i++,
    For[j=1,j<=nHiggs,j++,
        For[k=1,k<=nHiggs,k++,
            tmp = D[VCT,higgsbase[[i]],higgsbase[[j]],higgsbase[[k]]]/.RepAllZero;
            If[tmp =!= {0},
                Print["Curvature_Higgs_CT_L3[",i-1,"][",j-1,"][",k-1,"] = ", tmp[[1]]//InputForm,";"];
                nCount += 1;
            ];
        ];
    ];
    If[i == nHiggs && nCount == 0,
        Print["Curvature_Higgs_CT_L3[i,j,k] is ZERO for all entries."];
    ];
];

Print["----------------------------------------------------------------------"];

(* calculate L_ijkl *)
nCount = 0;
For[i=1,i<=nHiggs,i++,
    For[j=1,j<=nHiggs,j++,
        For[k=1,k<=nHiggs,k++,
            For[l=1,l<=nHiggs,l++,
                tmp = D[VCT,higgsbase[[i]],higgsbase[[j]],higgsbase[[k]],higgsbase[[l]]]/.RepAllZero//Simplify;
                If[tmp =!= {0},
                    Print["Curvature_Higgs_CT_L4[",i-1,"][",j-1,"][",k-1,"][",l-1,"] = ", tmp[[1]]//InputForm,";"];
                    nCount += 1;
                ];
            ];
        ];
    ];
    If[i == nHiggs && nCount == 0,
        Print["Curvature_Higgs_CT_L4[i,j,k,l] is ZERO for all entries."];
    ];
];

(* 
-------    Gauge interaction  --------
definition of the gauge potential and derivation of the corresponding tensor
--------------------------------------
*)

Print["----------------------------------------------------------------------"];
Print["                    Gauge interaction"]
Print["----------------------------------------------------------------------"];

gaugebase = {W1,W2,W3,B0};

nGauge = Length[gaugebase];

sigma0 = {{1,0},{0,1}};
sigma1 = {{0,1},{1,0}};
sigma2 = {{0,-I},{I,0}};
sigma3 = {{1,0},{0,-1}};

Dmu = I/2*(C_g*(sigma1*gaugebase[[1]]+sigma2*gaugebase[[2]]+sigma3*gaugebase[[3]])+C_gp*sigma0*gaugebase[[4]])//Simplify;

Dmudagger = -I/2*(C_g*(ConjugateTranspose[sigma1]*gaugebase[[1]] + ConjugateTranspose[sigma2]*gaugebase[[2]] + ConjugateTranspose[sigma3]*gaugebase[[3]]) + C_gp*sigma0*gaugebase[[4]])//Simplify;

VGauge = Flatten[Phi1zerodagger.Dmudagger.Dmu.Phi1zero + Phi2zerodagger.Dmudagger.Dmu.Phi2zero]//Simplify;

nCount = 0;
For[a=1,a<=nGauge,a++,
    For[b=1,b<=nGauge,b++,
        For[i=1,i<=nHiggs,i++,
            For[j=1,j<=nHiggs,j++,
                tmp = D[VGauge,gaugebase[[a]],gaugebase[[b]],higgsbase[[i]],higgsbase[[j]]]/.RepAllZero//Simplify;
                If[tmp =!= {0},
                    Print["Curvature_Gauge[",a-1,"][",b-1,"][",i-1,"][",j-1,"] = ",tmp[[1]]//InputForm,";"];
                    nCount += 1;
                ];
            ];
        ];
    ];
    If[a == nGauge && nCount == 0,
        Print["Curvature_Gauge[a,b,i,j] is ZERO for all entries."];
    ];
];

(* 
-------    Fermion interaction  ----------
definition of the Yukawa potential and derivation of the corresponding tensor
------------------------------------------
*)

Print["----------------------------------------------------------------------"];
Print["                    Fermion interaction"]
Print["----------------------------------------------------------------------"];

(* lepton sector *)
Ll = {{nueL, eL}, {numuL, muL}, {nutauL, tauL}};
Er = {eR,muR,tauR};

lepbase = {Ll[[1]][[2]],Er[[1]],Ll[[2]][[2]],Er[[2]],Ll[[3]][[2]],Er[[3]],Ll[[1]][[1]],Ll[[2]][[1]],Ll[[3]][[1]]};

nLep = Length[lepbase];

(* lepton yukawa couplings *)
yukLep = {yel,ymu,ytau};

(* lepton potential ( = - L_lep ) *)
VLepton = Flatten[Sum[yukLep[[i]]*Ll[[i]].Conjugate[Phi1zero]*Er[[i]],{i,1,3}]][[1]]//Simplify;

massLep = Table[D[VLepton,lepbase[[i]],lepbase[[j]]]/.RepAllVEV,{i,1,nLep},{j,1,nLep}];
MLep = Eigenvalues[massLep];
RepLepMass = Flatten[Solve[{MLep[[5]] == CMassElectron,MLep[[7]] == CMassMu,MLep[[9]] == CMassTau}, {yel, ymu, ytau}]];

nCount = 0;
For[i=1,i<=nLep,i++,
    For[j=1,j<=nLep,j++,
        For[k=1,k<=nHiggs,k++,
            tmp = D[VLepton,lepbase[[i]],lepbase[[j]],higgsbase[[k]]]/.RepLepMass//Simplify;
            If[tmp =!= 0,
                    Print["Curvature_Lepton[",i-1,"][",j-1,"][",k-1,"] = ",tmp//InputForm,";"];
                    nCount += 1;
                ];
        ];
    ];
    If[i == nLep && nCount == 0,
        Print["Curvature_Lep[I,J,k] is ZERO for all entries."];
    ];
];

Print["----------------------------------------------------------------------"];

(* quark sector *)
UL = {uL,cL,tL};
UR = {uR,cR,tR};
DL = {dL,sL,bL};
DR = {dR,sR,bR};

quarkbase = {uL,uR,dL,dR,cL,cR,sL,sR,tL,tR,bL,bR};

nQuark = Length[quarkbase];

VCKM = {{Vud,Vus,Vub},{Vcd,Vcs,Vcb},{Vtd,Vts,Vtb}};

yukd = {{yd,0,0},{0,ys,0},{0,0,yb}};
yuku = {{yu,0,0},{0,yc,0},{0,0,yt}};

VQuark = DR.yukd.ConjugateTranspose[VCKM].UL*Conjugate[Phi1zero[[1]]]+DR.yukd.DL*Conjugate[Phi1zero[[2]]]-UR.yuku.VCKM.DL*Phi1zero[[1]]+UR.yuku.UL*Phi1zero[[2]]//Simplify;

massQuark = Table[D[VQuark[[1]],quarkbase[[i]],quarkbase[[j]]]/.RepAllVEV,{i,1,nQuark},{j,1,nQuark}];

MQuark = Eigenvalues[massQuark];

RepQuarkMass = Flatten[Solve[{MQuark[[2]] == CMassBottom,MQuark[[4]] == CMassCharm,MQuark[[6]] == CMassDown,MQuark[[8]] == CMassStrange,MQuark[[10]] == CMassTop,MQuark[[12]] == CMassUp},{yb,yc,yd,ys,yt,yu}]];

nCount = 0;
For[i=1,i<=nQuark,i++,
    For[j=1,j<=nQuark,j++,
        For[k=1,k<=nHiggs,k++,
            tmp = D[VQuark,quarkbase[[i]],quarkbase[[j]],higgsbase[[k]]]/.RepQuarkMass//Simplify;
            If[tmp =!= {0},
                Print["Curvature_Quark[",i-1,"][",j-1,"][",k-1,"] = ",tmp[[1]]//InputForm,";"];
                nCount += 1;
            ];
        ];
    ];
    If[i == nQuark && nCount == 0,
        Print["Curvature_Quark[I,J,k] is ZERO for all entries."];
    ];
];

Close[strm];