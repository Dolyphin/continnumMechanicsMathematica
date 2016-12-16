(* ::Package:: *)

(*
 ---------------------------------------------------------------------------------------------------
 -
 -         CONTINUUM MECHANICS USING MATHEMATICA: FUNDAMENTALS, APPLICATIONS AND SCIENTIFIC COMPUTING
 -
 -                     A. Romano & A. Marasco
 -
 -
 -         THE PACKAGE MECHANICS.m
 -
 -             Title: MECHANICS.m
 -
 -             Authors: A. Romano, R. Lancellotta & A. Marasco
 -
 -             Copyright: Copyright 2014 by Birkauser-Springer. All rights reserved.
 -
 -             Mathematica Versions: 9.0
 -
 ----------------------------------------------------------------------------------------------------


15 giugno 2014*)



BeginPackage["MECHANICS`"]



(* = = = = = = = = = = = = = = = = = =  U S A G E  = = = = = = = = = = = = = = = = = = = *)

HelpVectorSys::usage:="Help for VectorSys"
UsageVectorSys::usage:="Use VectorSys"
VectorSys::usage:="The program VectorSys"


HelpEigenSystemAG::usage:="Help for EigenSystemAG"
UsageEigenSystemAG::usage:="Use EigenSystemAG"
EigenSystemAG::usage:="The program EigenSystemAG"


HelpOperator::usage:="Help for Operator"
UsageOperator::usage:="Use Operator"
Operator::usage:="The program Operator"


HelpDeformation::usage:="Help for Deformation"
UsageDeformation::usage:="Use Deformation"
Deformation::usage:="The program Deformation"

HelpVelocity::usage:="Help for Velocity"
UsageVelocity::usage:="Use Velocity"
Velocity::usage:="The program Velocity"

HelpLinElasticityTensor::usage:="Help for LinElasticityTensor"
UsageLinElasticityTensor::usage:="Use LinElasticityTensor"
LinElasticityTensor::usage:="The program LinElasticityTensor"

HelpPdeSysClass::usage:="Help for PdeSysClass"
UsagePdeSysClass::usage:="Use PdeSysClass"
PdeSysClass::usage:="The program PdeSysClass"

HelpPdeEqClass::usage:="Help for PdeEqClass"
UsagePdeEqClass::usage:="Use PdeEqClass"
PdeEqClass::usage:="The program PdeEqClass"

HelpWavesI::usage:="Help for WavesI"
UsageWavesI::usage:="Use WavesI"
WavesI::usage:="The program WavesI"

HelpWavesII::usage:="Help for WavesII"
UsageWavesII::usage:="Use WavesII"
WavesII::usage:="The program WavesII"


HelpPotential::usage:="Help for Potential"
UsagePotential::usage:="Use Potential"
Potential::usage:="The program Potential"

HelpWing::usage:="Help for Wing"
UsageWing::usage:="Use Wing"
Wing::usage:="The program Wing"

HelpJoukowsky::usage:="Help for Joukowsky"
UsageJoukowsky::usage:="Use Joukowsky"
Joukowsky::usage:="The program Joukowsky"

HelpJoukowskyMap::usage:="Help for JoukowskyMap"
UsageJoukowskyMap::usage:="Use JoukowskyMap"
JoukowskyMap::usage:="The program JoukowskyMap"



(* = = = = = = = = = = = = = =  B E G I N  P A C K A G E  = = = = = = = = = = = = = = = = *)


Begin["`Private`"]



(* = = = = = = = = = =  C O N T E N T S  O F  P A C K A G E  MECHANICS.m  = = = = = = = = = = *)


StylePrint["This package has been written by using the version 9.0 of Mathematica. The more recent versions could not work properly.", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
StylePrint["The user can understand the aims and the use of any program contained in this package by respectively typing HelpProgram[] o UsageProgram[], where Program is the name of each of them.", "Output", FontFamily -> "Times-Bold", FontSize -> 12, CellDingbat -> "\[FilledSquare]"];
StylePrint["List of programs contained in the package MECHANICS.m", "Output", FontFamily -> "Times-Bold", FontSize -> 12, CellDingbat -> "\[EmptyDiamond]"];
StylePrint["VectorSys:  Systems of equivalent vectors.","Output",FontFamily->"Times-Bold",CellDingbat->"\[FilledDiamond]"];
StylePrint["EigenSystemAG: Spectrum of a 2-tensor.","Output",FontFamily->"Times-Bold",CellDingbat->"\[FilledDiamond]"];
StylePrint["Operator:  Differential operators in curvilinear orthogonal coordinates.","Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["Deformation:  Analysis of finite deformations in curvilinear orthogonal coordinates.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["Velocity:  Dynamical characteristics of a velocity field in curvilinear orthogonal coordinates.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["LinElasticityTensor: Elastic potential and stress tensor in linear elasticity.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["PdeEqClass:  Second-order quasilinear PDE classification.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["PdeSysClass:  Classification of first-order quasilinear systems of PDEs.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["WavesI:  Characteristic equations and advancing velocity of waves in first-order quasilinear systems of PDEs.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["WavesII: Characteristic equations and advancing velocity of waves in quasilinear second-order PDE.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["Potential:  Complex potential and velocity. Level line of Stokes and kinetics potential.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["Wing:  Wing profile by Joukowsky's mapping.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["Joukowsky:  Streamlines around the cylinder and the wing profile.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];
StylePrint["JoukowskyMap: Geometrical analysis Joukowsky's map.", "Output", FontFamily -> "Times-Bold", CellDingbat -> "\[FilledDiamond]"];




(* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = *)


Off[General::spell];
Off[General::spell1];
Off[Solve::svars];
Off[Plot::plnr];
Off[Part::partd];
Off[Part::partw];
Off[Part::pspec];
Off[FindRoot::regex];
Off[FindRoot::cvnwt];
Off[InterpolatingFunction::dmval];
Off[General::wrsym];


(* = = = = = = = = = = = = = = =  T H E   P R O G R A M S  = = = = = = = = = = = = = = = *)


(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  VectorSys:
 -
 -                 UsageVectorSys[]
 -                 HelpVectorSys[]
 -                 VectorSys[A_, V_, P_]
 -
 -----------------------------------------------------------------------------------------
*)


UsageVectorSys[] :=
Module[{},
StylePrint["Aims of the program VectorSys", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["For a given applied vector system  \[CapitalSigma]={(\!\(A\_i\),\!\(v\_i\))\!\(\(}\_\(i = 1, \[CenterEllipsis], \
n\)\)\)  and a point P belonging to \!\(\[DoubleStruckCapitalR]\^3\), the program determines the momentum of \[CapitalSigma] with respect to P, an equivalent vector system \[CapitalSigma]' as well as the central axis, when the scalar invariant vanishes and the resultant does not vanish. Moreover, both the systems \[CapitalSigma] and \[CapitalSigma] ' are represented in the space, together with the central axis of \[CapitalSigma].", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
StylePrint["The command raw is", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["VectorSys[A, V, P],", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["where A is the list of the application points of vectors belonging to \[CapitalSigma], V is the list of the corresponding components, and P is the pole with respect to which the momentum of \[CapitalSigma] has to be evaluated.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]


HelpVectorSys[]:=
Module[{},
StylePrint["How to use the program VectorSys", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
StylePrint["For any point P, the program evaluates a vector system equivalent to an assigned one also determining the central axis.","Output", FontFamily->"Times-Plain",FontSize->12];
StylePrint["To run the program the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["A = application point list of the vectors of \[CapitalSigma] (f.e. A={{0, 1, 0}, {1, 0, 0}, {0, 1, 2}};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["V = component list of the vectors of \[CapitalSigma] (f.e. V={{1, 0, 1}, {2, 1, 0}, {3, 0, 0}};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["P = pole with respect to evaluate the moment of \[CapitalSigma] (f.e. P = {0, 0, 0};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["VectorSys[A, V, P]", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]


VectorSys[A_,V_,P_]:=

Module[{R,m,a,comp,MP,inv,e,d,arrow,arrowR,coppia,prodMR,asse,puntoAC, maxcomp,ranget,plotasse},

(*Dimension of applied vector system*)
n=Length[A];

(*Test on input data*)

Which[Dimensions[A]=!=Dimensions[V]||Length[P]=!=Dimensions[A][[2]],
      StylePrint["ERROR: The length of the application point list does not coincide with the length of the vector list, or the three lists defining elements in different spaces.","Output",FontFamily->"Times-Plain",FontSize->12];
      Goto[errorBIS],

     Dimensions[V][[2]]===Length[P]===2,
     StylePrint["ERROR: The input vectors belong to \!\(\[DoubleStruckCapitalR]\^2\). Apply again the program to a suitable system of \!\(\[DoubleStruckCapitalR]\^3\):", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
     StylePrint["Application point list","Output",FontFamily->"Times-Plain",FontSize->12];
     Print["A = ",Table[Append[A[[i]],0],{i,1,n}]];
     StylePrint["Vector list","Output",FontFamily->"Times-Plain",FontSize->12];
     Print["V = ",Table[Append[V[[i]],0],{i,1,n}]];
     StylePrint["Pole of \[CapitalSigma]","Output",FontFamily->"Times-Plain",FontSize->12];
     Print["P = ",Append[P,0],"."];
     Goto[errorBIS],

     Dimensions[A][[2]]=!=2&&Dimensions[A][[2]]=!=3,
     StylePrint["ERROR: The vectors do not belong to \!\(\[DoubleStruckCapitalR]\^3\).", "Output", FontFamily -> "Times-
     Plain", FontSize -> 12];
     Goto[errorBIS],

     Dimensions[A][[2]]===3,Goto[begin]];


Label[begin];

(*-------------EXECUTIVE SECTION--------------*)


(*Resultant*)
R=Sum[V[[i]],{i,1,n}];

(*Applied vectors at A-P*)
m[i_]:={(A[[i]]-P),V[[i]]};

(*Momentum with respect to P*)
MP=Sum[Cross[m[i][[1]],m[i][[2]]],{i,1,n}];

(*Scalar invariant*)
inv=MP.R;

(*System of applied vectors*)

arrow[i_]:={Graphics3D[{PointSize[0.03],Point[A[[i]]]}],Graphics3D[Line[{A[[i]],A[[i]]+V[[i]]}]],Graphics3D[{PointSize[0.02],Point[A[[i]]+V[[i]]]}]};
arrow=Table[arrow[i],{i,1,n}];

(*System equivalent to zero*)
If[N[R]===N[{0,0,0}]&&N[MP]===N[{0,0,0}],
   StylePrint["The applied vector system \[CapitalSigma] is equivalent to zero vector.","Output",FontFamily->
   "Times-Plain",FontSize->12];
   Print[StyleForm["System of applied vectors \[CapitalSigma].","Output",FontFamily->"Times-Plain",FontSize->12]];
   Print[Show[arrow,Axes->True,AspectRatio->Automatic,PlotRange->All]];
   Goto[end],Goto[1]];

Label[1];

(*System equivalent to the resultant applied at P, plus a couple having a momentum MP*)
  If[N[inv]=!=N[0],
     Print[StyleForm["The vector system \[CapitalSigma] is equivalent to the resultant ","Output",FontFamily->
     "Times-Plain",FontSize->12],"R = ",R,StyleForm[" applied at P, plus a couple having a momentum ",
     "Output",FontFamily->"Times-Plain",FontSize->12],"\!\(M\_P\) = ",MP,"."];
     Print[StyleForm["System of applied vectors \[CapitalSigma].","Output",FontFamily->"Times-Plain",FontSize->12]];
     Print[Show[arrow,Axes->True,AspectRatio->Automatic,PlotRange->All]];
     Goto[notvanInv],Goto[vaninv]];

Label[notvanInv];

(*Couple*)
coppia={Graphics3D[{Hue[0.01],PointSize[0.03],Point[P+{0.5,0,0}]}],
      Graphics3D[{Hue[0.01],Thickness[0.0008],
          Line[{P+{0.5,0,0},P+{0.5,0,0}+MP}]}],
      Graphics3D[{Hue[0.01],PointSize[0.02],Point[P+{0.5,0,0}+MP]}],
      Graphics3D[{Hue[0.01],PointSize[0.03],Point[P+{-0.5,0,0}]}],
      Graphics3D[{Hue[0.01],Thickness[0.0008],
          Line[{P+{-0.5,0,0},P+{-0.5,0,0}-MP}]}],
      Graphics3D[{Hue[0.01],PointSize[0.02],Point[P+{-0.5,0,0}-MP]}]};

(*Resultant*)
arrowR={Graphics3D[{Hue[0.01],PointSize[0.03],Point[P]}],
      Graphics3D[{Hue[0.01],Thickness[0.008],Line[{P,P+R}]}],
      Graphics3D[{Hue[0.01],PointSize[0.02],Point[P+R]}]};

Print[StyleForm["Equivalent vector system \[CapitalSigma]'.","Output",FontFamily->"Times-Plain",FontSize->12]];
Print[Show[{coppia,arrowR},Axes->True,AspectRatio->Automatic,PlotRange-> All]];
StylePrint["System of applied vectors \[CapitalSigma] and equivalent system \[CapitalSigma]'.","Output",FontFamily->"Times-Plain",FontSize->12];
Print[Show[{Flatten[arrow],coppia,arrowR},Axes->True,AspectRatio->Automatic,PlotRange->All]];

Goto[end];

Label[vaninv];

(*System equivalent to a resultant applied at a point of the central axis*)
If[N[R]=!=N[{0,0,0}],
   Print[StyleForm["The vector system \[CapitalSigma] is equivalent to the resultant ","Output",FontFamily->"Times-
   Plain",FontSize->12],"R = ",R,StyleForm[" applied at a point of the central axis.","Output",
   FontFamily->"Times-Plain",FontSize->12]];
   Goto[notvanR],Goto[vanR]];

  Label[notvanR];

  (*Central axis*)
  prodMR=Cross[MP,R];
  asse=t*R+prodMR/(R[[1]]*R[[1]]+R[[2]]*R[[2]]+R[[3]]*R[[3]]);

  (*Point of central axis*)
  puntoAC=asse/.t->0;

  Print[StyleForm["The central axis is the straight line of equations","Output",FontFamily->"Times-
  Plain",FontSize->12]];
  Print["x(t) = ",asse[[1]]/.t->Global`t,","];
  Print["y(t) = ",asse[[2]]/.t->Global`t,",    \[ForAll]t\[Element]\[DoubleStruckCapitalR]"];
  Print["z(t) = ",asse[[3]]/.t->Global`t,","];

  Print[StyleForm["and the system \[CapitalSigma] is equivalent to the resultant ","Output",FontFamily->"Times-
  Plain",FontSize->12],"R = ",R,StyleForm[" applied at any point of the central axis.","Output",
  FontFamily->"Times-Plain",FontSize->12]];

  (*Graphics*)
  Print[""];
  StylePrint["Central axis.","Output",FontFamily->"Times-Plain",FontSize->12];

  maxcomp=Table[A[[i]]+V[[i]],{i,1,n}];
  ranget=Max[maxcomp];

  plotasse=ParametricPlot3D[{asse[[1]],asse[[2]],asse[[3]]},{t,0,2*ranget},PlotRange->All];
  Print[plotasse];

  Print[StyleForm["System of applied vectors \[CapitalSigma].","Output",FontFamily->"Times-Plain",FontSize->12]];

  Print[Show[arrow,Axes->True,AspectRatio->Automatic,PlotRange->All]];

  (*Resultant applied at a point of central axis*)

  arrowRAC={Graphics3D[{Hue[0.01],PointSize[0.03],Point[puntoAC]}],
        Graphics3D[{Hue[0.01],Thickness[0.008],Line[{puntoAC,puntoAC+R}]}],
        Graphics3D[{Hue[0.01],PointSize[0.02],Point[puntoAC+R]}]};
  Print[""];
  Print[StyleForm["Equivalent vector system \[CapitalSigma]'.","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[Show[arrowRAC,plotasse,Axes->True,AspectRatio->Automatic,PlotRange->All]];

  StylePrint["System of applied vectors \[CapitalSigma] and equivalent system \[CapitalSigma]'.","Output",FontFamily->"Times-
  Plain",FontSize->12];
  Print[Show[{Flatten[arrow],Flatten[arrowRAC]},Axes->True,AspectRatio->Automatic,PlotRange->All]];
  Goto[end];


  Label[vanR];
  Print[StyleForm["The vector system \[CapitalSigma] is equivalent to a couple having a momentum ","Output",FontFamily->"Times-Plain",FontSize->12],"\!\(M\_P\) = ",MP,"."];

  Print[StyleForm["System of applied vectors \[CapitalSigma].","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[Show[arrow,Axes->True,AspectRatio->Automatic,PlotRange->All]];

  (*Couple*)
  coppia={Graphics3D[{Hue[0.01],PointSize[0.03],Point[P+{0.5,0,0}]}],
        Graphics3D[{Hue[0.01],Thickness[0.0008],
            Line[{P+{0.5,0,0},P+{0.5,0,0}+MP}]}],
        Graphics3D[{Hue[0.01],PointSize[0.02],Point[P+{0.5,0,0}+MP]}],
        Graphics3D[{Hue[0.01],PointSize[0.03],Point[P+{-0.5,0,0}]}],
        Graphics3D[{Hue[0.01],Thickness[0.0008],
            Line[{P+{-0.5,0,0},P+{-0.5,0,0}-MP}]}],
        Graphics3D[{Hue[0.01],PointSize[0.02],Point[P+{-0.5,0,0}-MP]}]};

  Print[StyleForm["Equivalent vector system \[CapitalSigma]'.","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[Show[{coppia},Axes->True,AspectRatio->Automatic,PlotRange->All]];
  StylePrint["System of applied vectors \[CapitalSigma] and equivalent system \[CapitalSigma]'.","Output",FontFamily->"Times-
  Plain",FontSize->12];
  Print[Show[{Flatten[arrow],coppia},Axes->True,AspectRatio->Automatic,PlotRange->All]];

  Goto[end];

  Label[end];

  StylePrint["Legend of Graphics","Output",FontFamily->"Times-Plain",FontSize->12];
  StylePrint[ "\[FilledCircle]\[LongDash]\[LongDash]\[LongDash]\[FilledSmallCircle]      Vectors of system \[CapitalSigma].","Output",FontFamily->"Times-Plain",FontSize->12];
  StylePrint[ "\*StyleBox[\"\[FilledCircle]\[LongDash]\[LongDash]\[LongDash]\\[FilledSmallCircle]\",FontColor->RGBColor[1, 0, 0]]      Vectors of system \[CapitalSigma]'.","Output",FontFamily->"Times-Plain",FontSize->12];
  StylePrint["\*StyleBox[\"\[FilledCircle]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\[KeyBar]\\[FilledSmallCircle]\",FontColor->RGBColor[1, 0, 0]]   Resultant of system \[CapitalSigma] applied on P.","Output",FontFamily->"Times-Plain",FontSize->12];



  Label[errorBIS];

  ]




(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  EigenSystemAG:
 -
 -                 UsageEigenSystemAG[]
 -                 HelpEigenSystemAG[]
 -                 EigenSystemAG[matrix_]
 -
 -----------------------------------------------------------------------------------------
*)

UsageEigenSystemAG[] :=
Module[{},
   StylePrint["Aims of the program EigenSystemAG", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
   Print[StyleForm["A 2-tensor by the matrix of its mixed components T=(\!\(T\_i\%j\)) being given, the program evaluates: ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
   Print[StyleForm["\[FilledSmallCircle] the distinct eigenvalues \!\(TraditionalForm\`\[Lambda]\_i\) with the corresponding algebraic and geometric multiplicities;", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
   Print[StyleForm["\[FilledSmallCircle] the eigenspace of each eigenvalue; ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
   Print[StyleForm["\[FilledSmallCircle] an eigenvector basis when the tensor admits a diagonal form.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
   Print[StyleForm["The command raw is", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
   Print[StyleForm["EigenSystemAG[matrix],", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
   Print[StyleForm["where matrix represents the tensor T.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
]



HelpEigenSystemAG[] :=
Module[{},
  StylePrint["How to use the program EigenSystemAG", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  StylePrint["The program evaluates the eigenvalues of a 2-tensor with their algebraic and geometric multiplicities as well as the corresponding eigenspaces.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["To run the program the following input datum has to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print[StyleForm["matrix = 2-tensor expressed by its mixed components (i.e. matrix = {{1,1,0,1},{1,0,1,0},{0,0,1,0},{0,1,0,1}};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["EigenSystemAG[matrix]", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];

 ]


EigenSystemAG[matrix_] :=
Module[{eig, eigve, eigM, molt, fSys, aut, eigind, matrixDiag},

(*Dimension of the input matrix*)
order = Dimensions[matrix];

(*Test on dimension of the input matrix*)
If[Length[order] === 2 && order[[1]] === order[[2]], Goto[begin], StylePrint["Input error: The matrix is not squared!","Output",FontSize->12,FontFamily->"Times-Plain"]; Goto[end]];

Label[begin];

(*Eigenvalues and eigenvectors: Algebraic and Geometric multiplicity*)

(*Eigenvalues*)
eig = Eigenvalues[matrix];

(*Eigenvectors*)
eigve = DeleteCases[Eigenvectors[matrix], Table[0, {i, 1, Length[matrix]}]];

eigM = Split[eig];

(*Algebraic Multiplicity*)
molt = Table[Length[eigM[[h]]], {h, 1, Length[eigM]}];

(*Eigenvectors belonging to the same eigenvalue*)
fSys[h_] := Table[If[Simplify[matrix . eigve[[i]]] === Simplify[eigM[[h,1]]*eigve[[i]]], eigve[[i]], no], {i, 1, Length[eigve]}];

(*Eigenspace*)
aut[h_] := DeleteCases[fSys[h], no];

(*Geometric Multiplicity*)
eigind = Table[Length[aut[h]], {h, 1, Length[eigM]}];
StylePrint["Algebraic and geometric multiplicity of distinct eigenvalues:","Output",FontSize->12,FontFamily->"Times-Plain"];
Print[""];
Do[Print[Subscript["\[Lambda]", i], " = ", eigM[[i,1]], ":   AlgMult = ", molt[[i]], "   GeomMult = ", eigind[[i]]];
Print[StyleForm["Eigenspace relative to ","Output",FontSize->12,FontFamily->"Times-Plain"], Subscript["\[Lambda]", i], ":  ", Subscript["V", i], " = ", aut[i], "."]; Print[],
{i, 1, Length[eigM]}];

(*Test on diagonal properties of input matrix*)
matrixDiag = DiagonalMatrix[Table[matrix[[i,i]], {i, 1, order[[1]]}]];
Which[matrix === matrixDiag,
                 StylePrint["The matrix is diagonal!","Output",FontSize->12,FontFamily->"Times-Plain"],
      matrix =!= matrixDiag && order === Dimensions[eigve],
                 StylePrint["The matrix is diagonal in the basis of its eigenvectors:","Output",FontSize->12,FontFamily->"Times-Plain"];
                 Print["\[ScriptCapitalB] = ", eigve, "."],
     matrix =!= matrixDiag && order =!= Dimensions[eigve],
                 StylePrint["The matrix cannot be reduced to a diagonal form!","Output",FontSize->12,FontFamily->"Times-Plain"];
                 Goto[end]];

Label[end];

]




(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  Operator:
 -
 -                 UsageOperator[]
 -                 HelpOperator[]
 -                 Operator[tensor_, var_, transform_, characteristic_, option_]
 -
 -----------------------------------------------------------------------------------------
*)


UsageOperator[] :=
  Module[{},
  StylePrint["Aims of the program Operator", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  Print[StyleForm["If a coordinate transformation from \!\(\*StyleBox[\"curvilinear\",\nFontWeight->\"Bold\"]\)\!\(\*StyleBox[\" \",\nFontWeight->\"Bold\"]\)coordinates to Cartesian ones is given, the program evaluates: ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["\[FilledSmallCircle] the Jacobian determinant of this transformation;", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["\[FilledSmallCircle] the holonomic bases associated to the new coordinates; ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["\[FilledSmallCircle] the metric matrix and its inverse;", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["\[FilledSmallCircle] Christoffel's symbols.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["Moreover, the program evaluates the unit basis and all the differential operators relative to scalar, vector, and tensor fields in the particular case of \!\(\*StyleBox[\"orthogonal curvilinear\", FontWeight->\"Bold\"]\) coordinates, provided that the aforesaid fields are given in the \!\(\*StyleBox[\"unit holonomic base\", FontWeight->\"Bold\"]\) or in the Cartesian frame. ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["The command raw is", "Output", FontFamily -> "Times-Plain",FontSize -> 12]];
  Print[StyleForm["Operator[tensor, var, transform, characteristic, option, simplifyoption],", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["where", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["tensor = {}, {scalar field} or components of a vector or second order tensor field in the unit base of curvilinear coordinates or in the Cartesian frame; ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["var = curvilinear or Cartesian coordinates; ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["transform =  right-hand sides of the transformation from curvilinear coordinates to Cartesian ones or of identity transformation; ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["characteristic = symbolic or numeric (see Chapter 2);", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["option = metric, operator or all (see Chapter 2);", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["simplifyoption = true or false (see Chapter 2).", "Output", FontFamily -> "Times-Plain", FontSize -> 12]]

]


HelpOperator[] :=
  Module[{},
  StylePrint["How to use the program Operator", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  StylePrint["The program determines the metric characteristics and the differential operators of scalar, vector, or second order tensor fields in a given system of curvilinear coordinates.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["To run the program, the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print[StyleForm["tensor = {}, {scalar field} or components of a vector or second order tensor field in the unit base of curvilinear coordinates or in the Cartesian frame (f.e. tensor = {F};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["var = curvilinear or Cartesian coordinates (f.e. var = {r, \[CurlyPhi]};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["transform = right-hand sides of the transformation from curvilinear coordinates to Cartesian ones or of identity transformation (f.e. transform = {r Cos[\[CurlyPhi]], r Sin[\[CurlyPhi]]};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["characteristic = symbolic or numeric (f.e. characteristic = symbolic;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["option = metric, operator or all (f.e. option = all;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["simplifyoption = true or false (f.e. option = true;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["Operator[tensor, var, transform, characteristic, option, simplifyoption]", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];

]


Operator[tensor_,var_,transform_,characteristic_,option_,simplifyoption_]:=

 Module[{n,e,u,m,PaDer,Jac,bvec,mco,met,invmet,cr,ort,cart,tensor1,f1,sostgradf,sostlap,gradf,hd,lap,tensorFin,sostgradV,
         sostZero,gradV,divVector,eLC,curlV,sostT,sostZeroT,covDer,divTensor},

     (*Test simplifyoption datum input*)
    If[simplifyoption === Global`true||simplifyoption === Global`false,
      Goto[go],
      StylePrint["ERROR: The user must give one of the choices between true or false to the input datum simplifyoption.",
      "Output",FontSize->12,FontFamily->"Times-Plain"];Goto[end]];

    Label[go];

    If[simplifyoption === Global`true,
       Unprotect[Sqrt];
       Unprotect[Times];
       Sqrt[a_/b_]*Sqrt[c_/d_] := Sqrt[(a*c)/(b*d)] // Simplify;
       Sqrt[a_/b_]*Sqrt[c_] := Sqrt[(a*c)/b] // Simplify;
       a_/(Sqrt[a_*b_]) := Sqrt[a/b] // Simplify;
       Sqrt[(a_)^(num_ /; EvenQ[num])] := a^(num/2);
       Sqrt[1/((a_)^(num_ /; EvenQ[num]))] := 1/(a^(num/2));
       Sqrt[(a_)^(num_ /; EvenQ[num])*(b_)^(numb_ /; EvenQ[numb])] := (a^(num/2))*(b^(numb/2));
       Sqrt[(a_)^(num_ /; EvenQ[num])*(b_)^(numb_ /; OddQ[numb])] := (a^(num/2))*Sqrt[b^numb];
       Sqrt[(a_)^(num_ /; EvenQ[num])*b_] := (a^(num/2))*Sqrt[b];
       Sqrt[(a_)^(num_ /; EvenQ[num])/(b_)^(numb_ /; OddQ[numb])] := (a^(num/2))/Sqrt[b^numb];
       Sqrt[(a_)^(num_ /; OddQ[num])/(b_)^(numb_ /; EvenQ[numb])] := Sqrt[a^num]/(b^(numb/2)), Goto[lecture]];

    Label[lecture];

    (*Dimension of the space*)
    n=Length[var];

    (*Order of tensor field*)
    m=Dimensions[tensor];

    (*Test on space dimension*)
    If[m=!={0}&&m=!={1}&&m[[1]]=!=n,
      StylePrint["ERROR: Component number of vector or tensor are not appropriate to the space dimension.","Output",
      FontSize->12,FontFamily->"Times-Plain"];Goto[end],Goto[test]];
    Label[test];

    If[n=!=Length[transform],
      StylePrint["ERROR: The list of the variables  and  the list of the right-hand side of the transformation functions have different lengths!","Output",
      FontSize->12,FontFamily->"Times-Plain"];Goto[end],Goto[test1]];
    Label[test1];

    (*Test of tensor input datum*)
    If[m==={0}||m==={1}||m==={2}||m==={3}||m==={2,2}||m==={3,3},Goto[OK],
      StylePrint["ERROR: The program evaluates the differential operators of scalar, vector or 2-tensor fields.",
      "Output",FontSize->12,FontFamily->"Times-Plain"];Goto[end]];

    Label[OK];

    (*Test characteristic datum input*)
    If[characteristic===Global`symbolic||characteristic===Global`numeric,
      Goto[continue],
      StylePrint["ERROR: The user must give one of the choices between symbolic or numeric to the input datum characteristic.",
      "Output",FontSize->12,FontFamily->"Times-Plain"];Goto[end]];

    Label[continue];

    (*Jacobian Matrix*)
    PaDer[i_,j_]:=D[transform[[i]],var[[j]]]//Simplify;
    Jac=Table[PaDer[i,j],{i,1,n},{j,1,n}]//Simplify;

    (*Test on coordinate transformations*)
    If[Det[Jac]=!=0,Goto[begin],
      Print[StyleForm["ERROR: The functions ","Output",
      FontSize->12,FontFamily->"Times-Plain"],transform,StyleForm[" do not define a coordinate change.","Output",
      FontSize->12,FontFamily->"Times-Plain"]];Goto[end]];


    (*---------------------EXECUTIVE SECTION------------------------*)

    Label[begin];
    (*Holonomic basis*)
    bvec[i_]:=Table[Jac[[j,i]],{j,1,n}]//Simplify;

    (*Metric coefficients and inverse ones*)
    mco[i_,j_]:=Simplify[bvec[i].bvec[j]];
    met=Table[mco[i,j],{i,1,n},{j,1,n}];
    invmet=Inverse[met]//Simplify;

    (*Christoffel's symbols*)
    cr[l_,j_,h_]:=(1/2)Sum[invmet[[l,i]](D[mco[i,j],var[[h]]]+D[mco[h,i],var[[j]]]-D[mco[j,h],var[[i]]]),{i,1,n}]//Simplify;

    (*Test of option input datum*)
    Which[option===Global`metric,Goto[metric],
          option===Global`all,Goto[metric],
          option===Global`operator,Goto[executive],
          option=!=Global`metric||option=!=Global`all||option=!=Global`operator,
                 StylePrint["ERROR: The user must give a value to input datum option.","Output",FontSize->12,
                 FontFamily->"Times-Plain"];
                 Print[StyleForm["The possible values are: metric, operator and all, to which correspond outputs relative to the geometric characteristics of curvilinear coordinates, to the differential operators, or to them both.","Output",
                 FontSize->12,FontFamily->"Times-Plain"]];
      Goto[end]];

    Label[metric];
    StylePrint["Jacobian matrix","Output",FontSize->12,FontFamily->"Times-Plain"];
    Print["J = ",MatrixForm[Jac]];
    StylePrint["Holonomic basis","Output",FontSize->12,FontFamily->"Times-Plain"];
    Do[Print[Subscript["e",i]," = ",Sum[bvec[i][[j]]*Subscript["u",j],{j,1,n}]],{i,1,n}];
    StylePrint["Metric matrix","Output",FontSize->12,FontFamily->"Times-Plain"];
    Print["G = ",MatrixForm[met]];
    StylePrint["Inverse metric matrix","Output",FontSize->12,FontFamily->"Times-Plain"];
    Print["\!\(G\^\(-1\)\) = ",MatrixForm[invmet]];
    taChr=Union[Table[cr[l,j,h],{l,1,n},{j,1,n},{h,1,n}]//Flatten];

    If[taChr=!={0},StylePrint["Christoffel's symbols","Output",FontSize->12,FontFamily->"Times-Plain"];

    ta=Do[
          If[cr[l,j,h]=!=0&&j<=h,
            Print[DisplayForm[
                SubscriptBox[
                  SubscriptBox[SuperscriptBox["\[CapitalGamma]",l],j],h]],
              " = ",cr[l,j,h]]],{l,1,n},{j,1,n},{h,1,n}],
      StylePrint[
        "Christoffel's symbols vanish since the coordinates are linear.",
        "Output",FontSize->12,FontFamily->"Times-Plain"]];

    (*Test on the orthogonality of coordinates var*)
    taUnit=Table[If[i=!=j&&mco[i,j]=!=0,mco,0],{i,1,n},{j,1,n}]//Flatten;
    If[taUnit===Table[0,{i,1,Length[taUnit]}],Goto[unit],Print[StyleForm["The coordinate system ","Output",
       FontSize->12,FontFamily->"Times-Plain"],var,StyleForm[" is not orthogonal.","Output",
       FontSize->12,FontFamily->"Times-Plain"]];Goto[end]];

    Label[unit];

    (*Unit Basis*)
    ort[i_]:=Simplify[bvec[i]/Sqrt[mco[i,i]]];
    StylePrint["Unit basis","Output",FontSize->12,FontFamily->"Times-Plain"];

    Do[Print[Subscript["a",i]," = ",Sum[ort[i][[j]]*Subscript["u",j],{j,1,n}]],{i,1,n}];
    Print[""];
    Which[option===Global`metric,Goto[end],
          option===Global`all,Goto[executive],
          option===Global`operator,Goto[executive]];
    Label[executive];

    (*Test on the orthogonality of coordinates var*)
    taUnit=Table[If[i=!=j&&mco[i,j]=!=0,mco,0],{i,1,n},{j,1,n}]//Flatten;
    If[taUnit===Table[0,{i,1,Length[taUnit]}],Goto[executiveB],
      Print[StyleForm["The coordinate system ","Output",FontSize->12,FontFamily->"Times-Plain"],var,
        StyleForm[" is not orthogonal.","Output",FontSize->12,FontFamily->"Times-Plain"]];Goto[end]];

    Label[executiveB];

    cart=Which[n===2&&var==={Global`x,Global`y},True,
               n===3&&var==={Global`x,Global`y,Global`z},True,
               n===2&&var=!={Global`x,Global`y},False,
               n===3&&var=!={Global`x,Global`y,Global`z},False];

    tensor1=
      Which[m=!={0}&&m=!={1}&&m=!={2,2}&&m=!={3,3}&&characteristic===Global`symbolic,Table[Apply[tensor[[i]],var],{i,1,n}],
            m=!={0}&&m=!={1}&&m=!={2}&&m=!={3}&&characteristic===Global`symbolic,Table[Apply[tensor[[i,j]],var],{i,1,n},{j,1,n}],
            m=!={0}&&m=!={1}&&m=!={2,2}&&m=!={3,3}&&characteristic===Global`numeric,Table[tensor[[i]],{i,1,n}],
            m=!={0}&&m=!={1}&&m=!={2}&&m=!={3}&&characteristic===Global`numeric,Table[tensor[[i,j]],{i,1,n},{j,1,n}],
            m==={0},tensor,
            m==={1}&&characteristic===Global`symbolic,{Apply[tensor[[1]],var]},
            m==={1}&&characteristic===Global`numeric,tensor];

    (*Substitution*)
    sostCart=Which[n===2,{Global`x->transform[[1]],Global`y->transform[[2]]},
                   n===3,{Global`x->transform[[1]],Global`y->transform[[2]],Global`z->transform[[3]]}];

    If[m==={1},Goto[1],Goto[2]];

    (*----------------------SCALAR FIELDS---------------------*)

    Label[1];

    (*Tensor in the new coordinates*)
    f1=Simplify[Which[n===2,tensor1/.sostCart,n===3,tensor1/.sostCart]][[1]];

    Which[Simplify[f1]=!=Simplify[tensor1[[1]]],Print[StyleForm["The scalar field ","Output",FontSize->12,FontFamily->"Times-Plain"],
                                                "F = ",tensor1[[1]],StyleForm[" in the coordinates ","Output",FontSize->12,FontFamily->"Times-Plain"],
                                                var,StyleForm["is written as ","Output",FontSize->12,FontFamily->"Times-Plain"]];
                                                Print["F = ",f1];Print[""];
                                                StylePrint["Differential operators in the unit basis.","Output",FontSize->12,FontFamily->"Times-Plain"],
          Simplify[f1]===Simplify[tensor1[[1]]]&&cart=!=True,Print[StyleForm["Differential operators in the unit basis, relative to the scalar field ","Output",FontSize->12,
                                                             FontFamily->"Times-Plain"]," F = ",tensor1[[1]]],
          Simplify[f1]===Simplify[tensor1[[1]]]&&cart===True,StylePrint["Differential operators in Cartesian coordinates.","Output",FontSize->12,FontFamily->"Times-Plain"]];

    (*Substitution for the gradient*)
    sostgradf=
      Which[n===2,{tensor[[1]]@@var->tensor[[1]],Derivative[1,0][tensor[[1]]]@@var->Subscript["\[PartialD]",
                  var[[1]]]@@tensor,Derivative[0,1][tensor[[1]]]@@var->Subscript["\[PartialD]",var[[2]]]@@tensor},
            n===3,{tensor[[1]]@@var->tensor[[1]],Derivative[1,0,0][tensor[[1]]]@@var->Subscript["\[PartialD]",
                  var[[1]]]@@tensor,Derivative[0,1,0][tensor[[1]]]@@var->Subscript["\[PartialD]",var[[2]]]@@tensor,
                  Derivative[0,0,1][tensor[[1]]]@@var->Subscript["\[PartialD]",var[[3]]]@@tensor}];

    (*Substitution for the Laplacian operator*)
    sostlap=Which[n===2,Union[sostgradf,{Derivative[2,0][tensor[[1]]]@@var->Subscript["\[PartialD]",
                        var[[1]]^2]@@tensor,Derivative[0,2][tensor[[1]]]@@var->Subscript["\[PartialD]",
                        var[[2]]^2]@@tensor}],
                  n===3,Union[sostgradf,{Derivative[2,0,0][tensor[[1]]]@@var->Subscript["\[PartialD]",
                        var[[1]]^2]@@tensor,Derivative[0,2,0][tensor[[1]]]@@var->Subscript["\[PartialD]",
                        var[[2]]^2]@@tensor,Derivative[0,0,2][tensor[[1]]]@@var->Subscript["\[PartialD]",
                        var[[3]]^2]@@tensor}]];

    (*Function gradient*)
    Print[""];
    StylePrint["Gradient","Output",FontSize->12,FontFamily->"Times-Plain"];
    gradf=Table[Sqrt[invmet[[i,i]]] D[f1,var[[i]]],{i,1,n}]//Simplify;

    Which[characteristic===Global`symbolic,Print["\[Del]F = ",MatrixForm[gradf/.sostgradf]],
          characteristic===Global`numeric,Print["\[Del]F = ",MatrixForm[gradf]]];

    (*Laplacian operator*)
    hd=Sqrt[Product[met[[i,i]],{i,1,n}]];
    lap=Simplify[Sum[(1/hd)*D[hd*(1/met[[i,i]])*D[f1,var[[i]]],var[[i]]],{i,1,n}]];
    Print[""];
    StylePrint["Laplacian operator","Output",FontSize->12,FontFamily->"Times-Plain"];

    Which[characteristic===Global`symbolic,Print["\[CapitalDelta]F = ",lap/.sostlap],
          characteristic===Global`numeric,Print["\[CapitalDelta]F = ",lap]];

    Goto[end];
    Label[2];
    If[m==={2}||m==={3},Goto[3],Goto[4]];
    Label[3];

    (*----------------------VECTOR FIELDS-------------------*)

    (*Tensor components in the unit basis*)
    tensorFin=
      Which[n===2&&Simplify[tensor1/.sostCart]===Simplify[tensor1],tensor1,
            n===3&&Simplify[tensor1/.sostCart]===Simplify[tensor1],tensor1,
            Simplify[tensor1/.sostCart]=!=Simplify[tensor1]&&n===2,Table[Simplify[Sum[Inverse[Jac][[i,j]]*Sqrt[met[[i,i]]]
                                                                   *tensor1[[j]],{j,1,n}]/.sostCart],{i,1,n}],
            Simplify[(tensor1/.sostCart)]=!=Simplify[tensor1]&&n===3,Table[Simplify[Sum[Inverse[Jac][[i,j]]*Sqrt[met[[i,i]]]*
                                                                     tensor1[[j]],{j,1,n}]/.sostCart],{i,1,n}]];

    Which[Simplify[tensorFin]=!=Simplify[tensor1],StylePrint["The vector field","Output",FontSize->12,FontFamily->"Times-Plain"];
                                                  Print["v = ",Sum[tensor1[[i]]*Subscript["u",i],{i,1,n}]];
                                                  Print[StyleForm[" in the unit basis of the coordinates ","Output",FontSize->12,FontFamily->"Times-Plain"],var,
                                                  StyleForm[" is written as","Output",FontSize->12,FontFamily->"Times-Plain"]];
                                                  Print["v = ",Sum[tensorFin[[i]]*Subscript["a",i],{i,n}]];
                                                  Print[""];
                                                  StylePrint["Differential Operators in the unit basis.","Output",FontSize->12,FontFamily->"Times-Plain"],
          Simplify[tensorFin]===Simplify[tensor1]&&cart=!=True,StylePrint["Differential operators in the unit basis.","Output",FontSize->12,FontFamily->"Times-Plain"],
          Simplify[tensorFin]===Simplify[tensor1]&&cart===True,StylePrint["Differential operators in the Cartesian coordinates.","Output",FontSize->12,FontFamily->"Times-Plain"]];


    sostgradV2[i_]:={tensor[[i]]@@var->tensor[[i]],Derivative[1,0][tensor[[i]]]@@var->Subscript["\[PartialD]",var[[1]]]@@{tensor[[i]]},
                    Derivative[0,1][tensor[[i]]]@@var->Subscript["\[PartialD]",var[[2]]]@@{tensor[[i]]}};

    sostgradV3[i_]:={tensor[[i]]@@var->tensor[[i]],Derivative[1,0,0][tensor[[i]]]@@var->Subscript["\[PartialD]",var[[1]]]@@{tensor[[i]]},
                     Derivative[0,1,0][tensor[[i]]]@@var->Subscript["\[PartialD]",var[[2]]]@@{tensor[[i]]},
                     Derivative[0,0,1][tensor[[i]]]@@var->Subscript["\[PartialD]",var[[3]]]@@{tensor[[i]]}};

    sostZero=Table[Subscript["\[PartialD]",var[[i]]][0]->0,{i,1,n}];

    sostgradV=Which[n===2,Union[sostgradV2[1],sostgradV2[2]],
                    n===3,Union[sostgradV3[1],sostgradV3[2],sostgradV3[3]]]/.sostZero;

    (*Gradient*)
    Print[""];
    StylePrint["Gradient","Output",FontSize->12,FontFamily->"Times-Plain"];

    gradV=
      Table[Simplify[(D[tensorFin[[k]]/Sqrt[met[[k,k]]],var[[h]]]+
                Sum[cr[k,l,h]*tensorFin[[l]]/Sqrt[met[[l,l]]],{l,1,n}])(Sqrt[
                  met[[k,k]]]*(1/Sqrt[met[[h,h]]]))],{h,1,n},{k,1,n}];

    Which[characteristic===Global`symbolic,Print["\[Del]v = ",MatrixForm[gradV/.sostgradV]],
          characteristic===Global`numeric,Print["\[Del]v = ",MatrixForm[gradV]]];

    (*Divergence*)
    hd=Sqrt[Product[met[[i,i]],{i,1,n}]];
    divVector=Sum[(1/hd)*D[hd*tensorFin[[i]]/Sqrt[met[[i,i]]],var[[i]]],{i,1,n}]//Simplify;
    Print[""];
    StylePrint["Divergence","Output",FontSize->12,FontFamily->"Times-Plain"];

    Which[n===2&&characteristic===Global`symbolic,Print["\[Del]\[CenterDot]v = ",divVector/.sostgradV],
          n===3&&characteristic===Global`symbolic,Print["\[Del]\[CenterDot]v = ",divVector/.sostgradV],
          characteristic===Global`numeric,Print["\[Del]\[CenterDot]v = ",divVector]];

    (*Curly*)
    (*Levi-Civita's symbol*)

    eLC[i_,h_,k_]:=Signature[{i,h,k}];
    curlV=Table[Simplify[(1/(2*hd))*(Sum[eLC[i,h,k]*(D[tensor1[[k]]*Sqrt[met[[k,k]]],var[[h]]]-
                         D[tensor1[[h]]*Sqrt[met[[h,h]]],var[[k]]]),{h,1,n},{k,1,n}])*Sqrt[met[[i,i]]]],{i,1,n}];

    Which[n===3&&characteristic===Global`symbolic,Print[""];StylePrint["Curly","Output",FontSize->12,FontFamily->"Times-Plain"];
                                                  Print["\[Del]\[Times]v = ",MatrixForm[curlV/.sostgradV]],
          n===3&&characteristic===Global`numeric,Print[""];StylePrint["Curly","Output",FontSize->12,FontFamily->"Times-Plain"];
                                                 Print["\[Del]\[Times]v = ",MatrixForm[curlV]]];
    Goto[end];

    (*----------------------2-TENSOR FIELDS--------------------*)

    Label[4];
    If[m==={2,2}||m==={3,3},Goto[tensor2],Goto[end]];
    Label[tensor2];

    (*Tensor in the new coordinates*)
    tensorFin=
      Which[n===2&&Simplify[tensor1/.sostCart]===Simplify[tensor1],tensor1,
            n===3&&Simplify[tensor1/.sostCart]===Simplify[tensor1],tensor1,
            Simplify[tensor1/.sostCart]=!=Simplify[tensor1]&&n===2,Table[Simplify[Sum[Inverse[Jac][[i,h]]*Sqrt[met[[i,i]]]*Inverse[Jac][[j,k]]*
                                                                   Sqrt[met[[j,j]]]*tensor1[[h,k]],{h,1,n},{k,1,n}]/.sostCart],{i,1,n},{j,1,n}],
            Simplify[tensor1/.sostCart]=!=Simplify[tensor1]&&n===3,Table[Simplify[Sum[Inverse[Jac][[i,h]]*Sqrt[met[[i,i]]]*Inverse[Jac][[j,k]]*
                                                                   Sqrt[met[[j,j]]]*tensor1[[h,k]],{h,1,n},{k,1,n}]/.sostCart],{i,1,n},{j,1,n}]];

    Which[Simplify[tensorFin]=!=Simplify[tensor1],StylePrint["The tensor field in Cartesian coordinates","Output",FontSize->12,FontFamily->"Times-Plain"];
                                                  Print["T = ",Sum[tensor1[[i,j]]*Subscript["u",i]\[CircleTimes]Subscript["u",j],{i,1,n},{j,1,n}]];
                                                  Print[StyleForm["in the unit basis relative to the coordinates ","Output",FontSize->12,FontFamily->"Times-Plain"],var,
                                                  StyleForm[" is written as","Output",FontSize->12,FontFamily->"Times-Plain"]];
                                                  Print["T = ",Sum[tensorFin[[i,j]]*Subscript["a",i]\[CircleTimes]Subscript["a",j],{i,1,n},{j,1,n}]];
                                                  Print[""];
                                                  StylePrint["Differential operators.","Output",FontSize->12,FontFamily->"Times-Plain"],
          Simplify[tensorFin]===Simplify[tensor1]&&cart=!=True,StylePrint["Differential operators in the unit basis, for the tensor field","Output",FontSize->12,FontFamily->"Times-Plain"];
                                                               Print["T = ",Sum[tensorFin[[i,j]]*Subscript["a",i]\[CircleTimes]Subscript["a",j],{i,n},{j,1,n}]],
          Simplify[tensorFin]===Simplify[tensor1]&&cart===True,StylePrint["Differential operators in Cartesian coordinates","Output",FontSize->12,FontFamily->"Times-Plain"]];


    sostT2[i_,j_]:={tensor[[i,j]]@@var->tensor[[i,j]],Derivative[1,0][tensor[[i,j]]]@@var->Subscript["\[PartialD]",
                    var[[1]]]@@{tensor[[i,j]]},Derivative[0,1][tensor[[i,j]]]@@var->Subscript["\[PartialD]",var[[2]]]@@{tensor[[i,j]]}};

    sostT3[i_,j_]:={tensor[[i,j]]@@var->tensor[[i,j]],Derivative[1,0,0][tensor[[i,j]]]@@var->Subscript["\[PartialD]",
                   var[[1]]]@@{tensor[[i,j]]},Derivative[0,1,0][tensor[[i,j]]]@@var->Subscript["\[PartialD]",var[[2]]]@@{tensor[[i,j]]},
                   Derivative[0,0,1][tensor[[i,j]]]@@var->Subscript["\[PartialD]",var[[3]]]@@{tensor[[i,j]]}};

    sostZeroT=Table[Subscript["\[PartialD]",var[[i]]][0]->0,{i,1,n}];

    sostT=Which[n===2,Union[sostT2[1,1],sostT2[1,2],sostT2[2,1],sostT2[2,2]],
                n===3,Union[sostT3[1,1],sostT3[1,2],sostT3[1,3],sostT3[2,1],sostT3[2,2],sostT3[2,3],sostT3[3,1],sostT3[3,2],
                      sostT3[3,3]]]/.sostZeroT;

    (*Covariant derivative*)
    covDer=Sum[Simplify[(D[tensorFin[[k,l]]/Sqrt[met[[k,k]]*met[[l,l]]],var[[h]]]+
                  Sum[cr[k,p,h]*(tensorFin[[p,l]]/Sqrt[met[[p,p]]*met[[l,l]]])+
                      cr[l,p,h]*(tensorFin[[k,p]]/Sqrt[met[[k,k]]*met[[p,p]]]),{p,1,n}])*(1/Sqrt[met[[h,h]]])*Sqrt[met[[k,k]]]*
                      Sqrt[met[[l,l]]]]*Subscript["a",h]\[CircleTimes]Subscript["a",k]\[CircleTimes]Subscript["a",l],{h,1,n},{k,1,n},{l,1,n}];
    Print[""];
    StylePrint["Covariant derivative","Output",FontSize->12,FontFamily->"Times-Plain"];

    Which[n===2&&characteristic===Global`symbolic,Print["\[Del]T = ",covDer/.sostT],
          n===3&&characteristic===Global`symbolic,Print["\[Del]T = ",covDer/.sostT],
          characteristic===Global`numeric,Print["\[Del]T = ",covDer]];

    (*Divergence*)
    divTensor=Sum[Simplify[(Sum[D[tensorFin[[h,l]]/Sqrt[met[[h,h]]*met[[l,l]]],
                      var[[h]]],{h,1,n}]+Sum[cr[h,k,h]*(tensorFin[[k,l]]/Sqrt[met[[k,k]]*met[[l,l]]])+
                      cr[l,k,h]*(tensorFin[[h,k]]/Sqrt[met[[k,k]]*met[[h,h]]]),{h,1,n},{k,1,n}])*
                      Sqrt[met[[l,l]]]]*Subscript["a",l],{l,1,n}];
    Print[""];
    StylePrint["Divergence","Output",FontSize->12,FontFamily->"Times-Plain"];

    Which[n===2&&characteristic===Global`symbolic,Print["\[Del]\[CenterDot]T = ",divTensor/.sostT],
          n===3&&characteristic===Global`symbolic,Print["\[Del]\[CenterDot]T = ",divTensor/.sostT],
          characteristic===Global`numeric,Print["\[Del]\[CenterDot]T = ",divTensor]];
    Label[end];

    If[simplifyoption===Global`true,Clear[Sqrt];Protect[Sqrt];Protect[Times],Goto[end2]];

    Label[end2];

]





(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  Deformation:
 -
 -                 UsageDeformation[]
 -                 HelpDeformation[]
 -                 Deformation[func_, var_, point_, vers1_, vers2_, transform_, option_]
 -
 -----------------------------------------------------------------------------------------
*)


UsageDeformation[]:=
  Module[{},
  StylePrint["Aims of the program Deformation","Output", FontFamily->"Times-Bold",FontSize->12];
  Print[StyleForm["A transformation from \!\(\* StyleBox[\"orthogonal curvilinear coordinates\",\nFontWeight->\"Bold\"]\) var to Cartesian ones and two or three deformations functions being given, the program evaluates: ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["\[FilledSmallCircle] the deformation gradient F;  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["\[FilledSmallCircle] the right Cauchy-Green tensor C; ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["\[FilledSmallCircle] the eigenvalues, an orthonormal eigenvector basis, and the principal invariants of C;  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["\[FilledSmallCircle] the left Cauchy-Green tensor B;","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["\[FilledSmallCircle] an orthonormal eigenvector basis of B, and the inverse tensor of B;  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["\[FilledSmallCircle] the right stretching tensor U;  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["\[FilledSmallCircle] the rotation tensor R;  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["\[FilledSmallCircle] the stretching ratios in the directions vers1 and vers2, evaluated at point;  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["\[FilledSmallCircle] the shearing angle between the directions vers1 and vers2 evaluated at point.  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["The command line is","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["Deformation[func, var, transform, point, vers1, vers2, option, simplifyoption]","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["where func denotes the right-hand side of finite deformation vector field, var the curvilinear or Cartesian coordinates, transform the right-hand side of coordinate transformation, point the coordinates of the point at which the stretching and the shearing have to be evaluated, option denotes an option to calculate the stretching and rotation tensor. Finally, simplifyoption allows the simplification of the irrational expressions. The components of vers1 and vers2 are given in the unit basis associated to the orthogonal curvilinear coordinates var.","Output",FontFamily->"Times-Plain",FontSize->12]];
]


HelpDeformation[]:=
  Module[{},
  StylePrint["How to use the program Deformation","Output",FontFamily->"Times-Bold",FontSize->12];
  StylePrint["The program analyzes the finite deformations in orthogonal curvilinear coordinates.","Output",FontFamily->"Times-Plain",FontSize->12];
  StylePrint["To run the program the following input data have to be given:", "Output",FontFamily->"Times-Plain",FontSize->12];
  Print[StyleForm["func = right-hand side of the finite deformation field (f.e. func = {\!\(X\_1\) + \!\(kX\_2\), \!\(X\_2\)};);  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["var = curvilinear or Cartesian coordinates in the reference configuration (f.e. var = {\!\(X\_1\), \!\(X\_2\)};)","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["transform =  right-hand side of the transformation from curvilinear coordinates var to Cartesian ones or the identity transformation (f.e. transform = {\!\(X\_1\), \!\(X\_2\)};)  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["point = coordinates of the point at which the stretching and the shearing between vers1 and vers2 have to be evaluated (f.e. point = {0, 0};)  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["vers1 = components of the first vector (f.e. vers1 = {1, 0};) ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["vers2 = components of the second vector (f.e. point = {0, 1};)  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["option = symbolic, numeric or null (f.e option = null;)  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["simplifyoption = true or false (f.e simplifyoption = true;)  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  Print[StyleForm["Deformation[func, var, transform, point, vers1, vers2, option, simplifyoption]","Output",FontFamily->"Times-Plain",FontSize->12]];
]




Deformation[func_, var_, transform_, point_, vers1_, vers2_, option_, simplifyoption_] :=
 Module[{dim, g, jacma, Jac, PaDer, jacobian, sost, rg, defgr, DCG, autc, eigc, mult, mod, inv, SCG,
 matdil, matrot, stretch, scorr,\[CapitalTheta]},

 (*Finite Deformations*)

 (*A new definition of Sqrt in the real field*)

 If[simplifyoption === Global`true, Unprotect[Sqrt]; Unprotect[Times]; Unprotect[Plus];
      Sqrt[a_/b_]*Sqrt[c_/d_] := Simplify[Sqrt[(a*c)/(b*d)]];
      Sqrt[a_/b_]*Sqrt[c_] := Simplify[Sqrt[(a*c)/b]];
      a_/Sqrt[a_*b_] := Simplify[Sqrt[a/b]];
      Sqrt[a_^(num_ /; EvenQ[num])] := a^(num/2);
      -Sqrt[a_^(num_ /; EvenQ[num])] := -a^(num/2);
      a_ + Sqrt[b_^(numb_ /; EvenQ[numb])] := a + b^(numb/2);
      Sqrt[1/a_^(num_ /; EvenQ[num])] := 1/a^(num/2);
      Sqrt[a_^(num_ /; EvenQ[num])*b_^(numb_ /; EvenQ[numb])] := a^(num/2)*b^(numb/2);
      Sqrt[a_^(num_ /; EvenQ[num])*b_^(numb_ /; OddQ[numb])] := a^(num/2)*Sqrt[b^numb];
      Sqrt[a_^(num_ /; EvenQ[num])*b_] := a^(num/2)*Sqrt[b];
      Sqrt[a_^(num_ /; EvenQ[num])/b_^(numb_ /; OddQ[numb])] := a^(num/2)/Sqrt[b^numb];
      Sqrt[a_^(num_ /; OddQ[num])/b_^(numb_ /; EvenQ[numb])] := Sqrt[a^num]/b^(numb/2), Goto[lecture]];

 Label[lecture];

 (*Test on input data*)
 If[Length[func] === Length[var] === Length[transform] === Length[vers1] === Length[vers2] ===
      Length[point], Goto[executive],
      StylePrint["ERROR: The input functions and/or the vectors do not have the same dimension.",
      "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Goto[end]];
 Label[executive];
 If[option === Global`symbolic || option === Global`numeric || option === Global`null, Goto[executive1],
    StylePrint["ERROR: The user must give one of the choices among symbolic, numeric, and null to the input datum option.",
    "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Goto[end]];

 Label[executive1];

 (*Dimension of the space*)
 dim = Length[func];

 (*Metric Coefficients*)
 g[i_, j_] := Simplify[Sum[D[transform[[h]], var[[i]]]*D[transform[[h]], var[[j]]], {h, 1, dim}]];

 (*Test on the orthogonality of coordinates var*)
 Do[Do[If[i=!= j && g[i, j] =!= 0, StylePrint["The curvilinear coordinate system ", var, " is not orthogonal.", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Goto[end], Goto[begin1]],
      {j, 1, dim}], {i, 1, dim}];

 Label[begin1];
 jacma[i_, j_] := D[func[[i]], var[[j]]];
 Jac = Table[Simplify[jacma[i, j]], {i, 1, dim}, {j, 1, dim}];

 (*Test on the coordinate transformation*)
 PaDer[i_, j_] := Simplify[D[transform[[i]], var[[j]]]];
 jacobian = Simplify[Table[PaDer[i, j], {i, 1, dim}, {j, 1, dim}]];

 If[Det[jacobian] =!= 0, Goto[begin], StylePrint["ERROR: The functions", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Do[Print[Subscript["x", i], "=", transform[[i]]], {i, 1, dim}];
 StylePrint["do not define a coordinate change.", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Goto[end]];

 Label[begin];

 (*Deformation gradient*)
 sost = Flatten[Table[{var[[i]] -> func[[i]]}, {i, 1, dim}]];
 rg[i_, j_] := Simplify[(Sqrt[g[i, i]] /. sost)/Sqrt[g[j, j]]];
 jacma1[i_, j_] := rg[i, j]*jacma[i, j];

 defgr = Table[Simplify[jacma1[i, j]], {i, 1, dim}, {j, 1, dim}];

 If[var =!= transform, Print[StyleForm["All tensor and vector components are relative to the basis of unit vectors
 associated with the natural basis of the curvilinear coordinates ", "Output", FontSize -> 12, FontFamily -> "Times-Plain"], var, "."]];

 Print[" "];
 StylePrint["Deformation gradient", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 Print["F = ", MatrixForm[defgr]];

 StylePrint["Right Cauchy-Green tensor", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 DCG = Transpose[defgr] . defgr;
 Print["C =", MatrixForm[DCG]];

 StylePrint["Eigenvalues of C", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 autc = Eigenvalues[DCG];
 eigc = Split[autc];
 mult = Table[Length[eigc[[h]]], {h, 1, Length[eigc]}];
 Do[Print[Subscript["\[Lambda]", i], " = ", Factor[eigc[[i,1]]], ": AlgMult = GeoMult = ", mult[[i]]],
     {i, 1, Length[eigc]}];

 StylePrint["An orthonormal basis of eigenvectors of C", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 eigcnun = Eigenvectors[DCG];
 modc[i_] := Sqrt[Sum[eigcnun[[i,j]]^2, {j, 1, dim}]];
 eigc = Table[Simplify[eigcnun[[i,j]]/modc[i]], {i, 1, dim}, {j, 1, dim}];
 Table[Print[Subscript["u", i], " = ", eigc[[i]]], {i, 1, dim}];

 (*Principal Invariants of C*)
 inv1 = Simplify[Sum[DCG[[i,i]], {i, 1, dim}]];
 If[dim === 2, Goto[2], Goto[1]];
 Label[1];
 inv2 = Simplify[Det[{{DCG[[1,1]], DCG[[1,2]]}, {DCG[[2,1]], DCG[[2,2]]}}] + Det[{{DCG[[1,1]], DCG[[1,3]]}, {DCG[[3,1]], DCG[[3,3]]}}] + Det[{{DCG[[2,2]], DCG[[2,3]]}, {DCG[[3,2]], DCG[[3,3]]}}]];

 inv3 = Simplify[Det[DCG]];
 StylePrint["Principal invariants of C", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 Print["\!\(I\_C\) = ", inv1];
 Print["\!\(II\_C\) = ", inv2];
 Print["\!\(III\_C\) = ", inv3]; Goto[3];

 Label[2];
 inv3 = Simplify[Det[DCG]];
 StylePrint["Principal invariants of C", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 Print["\!\(I\_C\) = ", inv1];
 Print["\!\(II\_C\) = ", inv3];

 Label[3];
 SCG = Simplify[Factor[defgr . Transpose[defgr]]];

 StylePrint["Left Cauchy-Green tensor", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 Print["B = ", MatrixForm[SCG]];

 StylePrint["An orthonormal basis of eigenvectors of B", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 eigbnun = Eigenvectors[SCG];
 modb[i_] := Sqrt[Sum[eigbnun[[i,j]]^2, {j, 1, dim}]];
 eigb = Table[eigbnun[[i,j]]/modb[i], {i, 1, dim}, {j, 1, dim}];
 Table[Print[Subscript["v", i], " = ", Simplify[eigb[[i]]]], {i, 1, dim}];

 StylePrint["Inverse of the left Cauchy-Green tensor", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 invb = Simplify[Factor[Inverse[SCG]]];
 Print["\!\(B\^\(-1\)\) = ", MatrixForm[invb]];

 If[option === Global`symbolic || option === Global`numeric, Goto[matRU], Goto[delta]];

 Label[matRU];

 (*Right Stretching Tensor*)
 matdil = Table[Simplify[Sum[Sqrt[autc[[i]]]*eigc[[i,h]]*eigc[[i,j]], {i, 1, dim}]], {h, 1, dim}, {j, 1, dim}];

 (*Rotation Tensor*)
 matrot = Simplify[Factor[defgr . Inverse[matdil]]];

 Which[option === Global`symbolic, StylePrint["Right stretching tensor", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Print["U = ", MatrixForm[matdil]];
      StylePrint["Rotation tensor", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Print["R = ", MatrixForm[matrot]], option === Global`numeric,
     StylePrint["Right stretching tensor", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Print["U = ", MatrixForm[Chop[N[matdil], 10^(-3)]]]; StylePrint["Rotation tensor", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Print["R = ", MatrixForm[Chop[N[matrot], 10^(-3)]]], option === Global`null, Goto[delta]];

 Label[delta];

 (*Stretching and Shearing*)
 sost1 = Flatten[Table[{var[[i]] -> point[[i]]}, {i, 1, dim}]];
 DCG0 = DCG //. sost1;

 (*Modules of vers1 and vers2 in the unit basis*)
 mod1 = Sqrt[Sum[vers1[[i]]^2, {i, 1, dim}]];
 mod2 = Sqrt[Sum[vers2[[i]]^2, {i, 1, dim}]];

 vers = {vers1, vers2};
 mod = {mod1, mod2};

 (*Stretching*)
 stretch[i_] := Sqrt[Simplify[(vers[[i]]/mod[[i]]) . DCG0 . (vers[[i]]/mod[[i]])]];

 (*Shearing between the directions vers1 and vers2*)
 scorr = (vers1/mod1) . DCG0 . (vers2/mod2)/(Sqrt[(vers1/mod1) . DCG0 . (vers1/mod1)]*Sqrt[(vers2/mod2) . DCG0 . (vers2/mod2)]);

 (*Initial and final angles between vers1 and vers2*)
 \[CapitalTheta]in = ArcCos[vers1 . vers2/(mod1*mod2)]*(180/Pi);
 \[CapitalTheta]fin = ArcCos[scorr]*(180/Pi);

 Which[vers[[1]] === vers[[2]], Print[StyleForm["Stretch ratio in the direction ", "Output", FontSize -> 12,
                                FontFamily -> "Times-Plain"], vers1/mod1];
                                Print["\[Delta]", " = ", stretch[1]],
       vers[[1]] =!= vers[[2]], Do[Print[StyleForm["Stretch ratio in the direction ", "Output", FontSize -> 12,
                                FontFamily -> "Times-Plain"], vers[[i]]/mod[[i]]];
                                Print[Subscript["\[Delta]", i], " = ", stretch[i]], {i, 1, 2}];
                                Print[StyleForm["Shear angle  between the directions ", "Output", FontSize -> 12,
                                FontFamily -> "Times-Plain"], vers1/mod1, StyleForm[" and ", "Output", FontSize -> 12,
                                FontFamily -> "Times-Plain"], vers2/mod2];
                                Print["cos \!\(\[CapitalTheta]\_12\) = ", scorr];
                                StylePrint["Initial value of the angle (in degree)", "Output", FontSize -> 12,
                                FontFamily -> "Times-Plain"];
                                Print["\!\(\[CapitalTheta]\_12\^i\) = ", \[CapitalTheta]in^Degree/.{(180/Pi)^Degree->1}];
                                StylePrint["Final value of the angle (in degree)", "Output", FontSize -> 12,
                                FontFamily -> "Times-Plain"];
                                Print["\!\(\[CapitalTheta]\_12\^f\) = ", \[CapitalTheta]fin^Degree/.{(180/Pi)^Degree->1}]];

 Label[end];

 If[simplifyoption === Global`true, Clear[Sqrt]; Protect[Sqrt], Goto[end2]];

 Label[end2];

]




(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  Velocity:
 -
 -                 UsageVelocity[]
 -                 HelpVelocity[]
 -                 Velocity[vel_, var_, transform_, characteristic_, simplifyoption_]
 -
 -----------------------------------------------------------------------------------------
*)


UsageVelocity[] :=
  Module[{},
  StylePrint["Aims of the program Velocity", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  Print[StyleForm["A coordinate transformation from orthogonal curvilinear coordinates to Cartesian ones being given together with the components of a velocity field, in the unit basis of the curvilinear coordinates, the program evaluates the acceleration, the gradient, the divergence, and the angular velocity.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  StylePrint["The command line is", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["Velocity[vel, var, transform, characteristic, simplifyoption],", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["where vel are the components of the velocity field, var the curvilinear orthogonal coordinates, transform the right-hand side of the transformation from curvilinear to Cartesian coordinates, characteristic is an option assuming the values symbolic or numeric, and finally simplifyoption allows the simplification of the irrational expressions.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]



HelpVelocity[] :=
  Module[{},
  StylePrint["How to use the program Velocity", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  StylePrint["Velocity field given, the program evaluates the acceleration, the gradient, the divergence, and the angular velocity.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["To run the program the following input data have to be given:", "Output",FontFamily->"Times-Plain",FontSize->12];
  StylePrint["vel = velocity components (f.e. vel = {\!\(v\_r\), \!\(v\_\[CurlyPhi]\), \!\(v\_z\)};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["var = orthogonal curvilinear or Cartesian coordinates (f.e. var = {r, \[CurlyPhi], z};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["transform = right-hand side terms of the transformation from curvilinear to Cartesian coordinates (f.e. transform = {r Cos[\[CurlyPhi]], r Sin[\[CurlyPhi]], z};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["characteristic = symbolic or numeric (f.e. characteristic = symbolic;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print[StyleForm["simplifyoption = true or false (f.e simplifyoption = true;)  ","Output",FontFamily->"Times-Plain",FontSize->12]];
  StylePrint["Velocity[vel, var, transform, characteristic, simplifyoption]", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]



Velocity[vel_, var_, transform_, characteristic_, simplifyoption_] :=
  Module[{dim, dimVel, matrix, Chr, speed, sost, GradV, Accel, DivV, curlV},

(*A new definition of Sqrt in the real field*)
If[simplifyoption === Global`true, Unprotect[Sqrt];
      Sqrt[alpha_*alpha_] := alpha;
      Sqrt[(alpha_/beta_)^2] := alpha/beta;
      Sqrt[(alpha_*beta_)^2] := alpha*beta;
      Sqrt[((alpha_)^2)*beta_] := alpha*Sqrt[beta];
      Sqrt[alpha_*((beta_)^2)] := Sqrt[alpha]*beta;
      Sqrt[alpha_/((beta_)^2)] := Sqrt[alpha]/beta;
      Sqrt[((alpha_)^2)/beta_] := alpha/Sqrt[beta];
      Sqrt[1/((beta_)^2)] := 1/beta;
      Unprotect[Times];
      Sqrt[alpha_/beta_]*Sqrt[gamma_/delta_] := Sqrt[(alpha*gamma)/(beta*delta)] // Simplify;
      Sqrt[alpha_/beta_]*Sqrt[gamma_] := Sqrt[(alpha*gamma)/beta] // Simplify;
      alpha_/(Sqrt[alpha_*beta_]) := Sqrt[alpha/beta] // Simplify, Goto[lecture]];

Label[lecture];

(*Test on input data*)
If[Length[vel] === Length[var] === Length[transform], Goto[executive],
      StylePrint["ERROR: The input functions and/or the vectors do not have the same dimension.",
       "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Goto[end]];

Label[executive];

If[characteristic === Global`symbolic || characteristic === Global`numeric, Goto[executive1],
      StylePrint["ERROR: The user must give one of the choices between symbolic and numeric to the input datum characteristic.",
      "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Goto[end]];

Label[executive1];

(*Dimension of the space*)
dim = Length[var];
dimVel = Dimensions[vel];

(*Jacobian of coordinate transformation*)
matrixJac = Table[D[transform[[i]], var[[j]]], {i, 1, dim}, {j, 1, dim}] // Simplify;

(*Metric Matrix*)
 matrixG = Table[Sum[matrixJac[[h, i]]*matrixJac[[h, j]], {h, 1, dim}], {i, 1, dim}, {j, 1, dim}] // Simplify;

 (*Inverse Metric Matrix*)
 matrixG1 = Inverse[matrixG] // Simplify;

 (*Christoffel's symbols*)
 Chr[i_, j_, h_] := (1/2)*Sum[matrixG1[[i, k]]*(D[matrixG[[k, j]], var[[h]]] + D[matrixG[[h, k]], var[[j]]] -
                  D[matrixG[[j, h]], var[[k]]]), {k, 1, dim}] // Simplify;

(*Velocity*)
speed = Which[(dimVel === {2} || dimVel === {3}) && characteristic === Global`symbolic,
        Table[vel[[i]] @@ Join[var, {Global`t}], {i, 1, dim}],
        (dimVel === {2} || dimVel === {3}) && characteristic === Global`numeric, Table[vel[[i]], {i, 1, dim}]];

(*Substitutions*)
sostV2[i_] := {vel[[i]] @@ Join[var, {Global`t}] -> vel[[i]], Derivative[1, 0, 0][vel[[i]]] @@ Join[var, {Global`t}] ->
          Subscript["\[PartialD]", var[[1]]] @@ {vel[[i]]}, Derivative[0, 1, 0][vel[[i]]] @@ Join[var, {Global`t}] ->
          Subscript["\[PartialD]", var[[2]]] @@ {vel[[i]]}, Derivative[0, 0, 1][vel[[i]]] @@ Join[var, {Global`t}] ->
          Subscript["\[PartialD]", Global`t] @@ {vel[[i]]}};

sostV3[i_] := {vel[[i]] @@ Join[var, {Global`t}] -> vel[[i]], Derivative[1, 0, 0, 0][vel[[i]]] @@ Join[var, {Global`t}] ->
          Subscript["\[PartialD]", var[[1]]] @@ {vel[[i]]}, Derivative[0, 1, 0, 0][vel[[i]]] @@ Join[var, {Global`t}] ->
          Subscript["\[PartialD]", var[[2]]] @@ {vel[[i]]}, Derivative[0, 0, 1, 0][vel[[i]]] @@ Join[var, {Global`t}] ->
          Subscript["\[PartialD]", var[[3]]] @@ {vel[[i]]}, Derivative[0, 0, 0, 1][vel[[i]]] @@ Join[var, {Global`t}] ->
          Subscript["\[PartialD]", Global`t] @@ {vel[[i]]}};

(*Nuovo comando introdotto per la versione 9*)
sostZero =Union[Union[
  Table[Subscript["\[PartialD]", var[[i]]][0] -> 0, {i, 1, 
    dim}], {Subscript["\[PartialD]", Global`t][0] -> 0}], {
\!\(\*SuperscriptBox[\(Subscript\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[_, _] @@ Join[var, {Global`t}] -> 0}];

sost = Which[dim === 2, Union[sostV2[1], sostV2[2]], dim === 3, Union[sostV3[1], sostV3[2], sostV3[3]]] //. sostZero;

If[var =!= transform, Print[StyleForm["The components of any vector or tensor quantity are relative to the unit basis associated with the holonomic basis
of the curvilinear coordinates ", "Output", FontSize -> 12, FontFamily -> "Times-Plain"], var, "."]];

(* Gradient of velocity *)
GradV1 = Table[Simplify[(Sqrt[matrixG[[i, i]]]/Sqrt[matrixG[[j, j]]])*(D[speed[[i]]/Sqrt[matrixG[[i, i]]],
                    var[[j]]] + Sum[Chr[i, j, h]*(speed[[h]]/Sqrt[matrixG[[h, h]]]), {h, 1, dim}])], {j, 1, dim}, {i, 1, dim}] // Simplify;

StylePrint["Spatial gradient of velocity", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];

Which[dim === 2 && characteristic === Global`symbolic, Print["\[Del]v = ", MatrixForm[(GradV1 //. sostZero)//.sost]],
      dim === 3 && characteristic === Global`symbolic, Print["\[Del]v = ", MatrixForm[(GradV1 //. sostZero)//.sost]],
      characteristic === Global`numeric, Print["\[Del]v = ", MatrixForm[GradV1]]];

(*Acceleration*)
Accel1[i_] := D[speed[[i]], Global`t] + Sum[speed[[j]]*GradV1[[i, j]], {j, 1, dim}];
Accel = Table[Simplify[(Accel1[i] //. sostZero)//.sost], {i, 1, dim}];
StylePrint["Acceleration", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
Print["a = ", MatrixForm[Accel]];

(*Divergence*)
DivV = Simplify[Sum[(GradV1[[i, i]] //. sostZero)//.sost, {i, 1, dim}]];
StylePrint["Divergence of velocity", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
Print["\[Del]\[CenterDot]v = ", DivV];

(*Levi - Civita's symbol*)
eLC[i_, h_, k_] := Signature[{i, h, k}];

curlV = Table[Simplify[(1/(2*Sqrt[Product[matrixG[[i, i]], {i, 1, dim}]]))*Sum[eLC[i, h, k]*(D[speed[[k]]*
             Sqrt[matrixG[[k, k]]], var[[h]]] - D[speed[[h]]*Sqrt[matrixG[[h, h]]], var[[k]]]), {h, 1, dim},
              {k, 1, dim}]*Sqrt[matrixG[[i, i]]]], {i, 1, dim}];

Which[dim === 3 && characteristic === Global`symbolic, StylePrint["Angular velocity", "Output", FontSize -> 12,
        FontFamily -> "Times-Plain"];Print["\[Omega] = ", MatrixForm[((1/2)*curlV //. sostZero)//.sost]],
      dim === 3 && characteristic === Global`numeric, StylePrint["Angular velocity", "Output", FontSize -> 12,
        FontFamily -> "Times-Plain"];Print["\[Omega] = ", MatrixForm[(1/2)*curlV]]];

Label[end];

If[simplifyoption === Global`true, Clear[Sqrt]; Protect[Sqrt];
      Unprotect[Times]; Clear[Times]; Protect[Times], Goto[end1]];

Label[end1];

]



(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  LinElasticityTensor:
 -
 -                 UsageLinElasticityTensor[]
 -                 HelpLinElasticityTensor[]
 -                 UsageLinElasticityTensor[class]
 -
 -----------------------------------------------------------------------------------------
*)


UsageLinElasticityTensor[] :=
Module[{},
StylePrint["Aims of the program LinElasticityTensor", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["The program applies to \*StyleBox[\"isotropic linear\", FontWeight->\"Bold\"] material or to \*StyleBox[\"anisotropic \", FontWeight->\"Bold\"] ones belonging to the \*StyleBox[\"monoclinic\", FontWeight->\"Bold\"],  \*StyleBox[\"triclinic\", FontWeight->\"Bold\"], and \*StyleBox[\"rombic\", FontWeight->\"Bold\"] crystal classes.","Output",FontFamily->"Times-Plain",FontSize->12]];
StylePrint["It determines:","Output",FontFamily->"Times-Plain",FontSize->12];
Print[StyleForm["\[FilledSmallCircle] the elastic potential \[Psi] as a function of the principal invariants of E or of the invariants of the crystal class to which the solid belongs;", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["\[FilledSmallCircle] the components of \[DoubleStruckCapitalC] in Voigt's notation; ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["\[FilledSmallCircle] the expression of \[Psi] in Voigt's notation; ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["\[FilledSmallCircle] the component of stress tensor T in Voigt's notation.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["The command raw is", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["LinElasticityTensor[class],", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["where class characterizes the solid and can assume the values isotropic, triclinic,
monoclinic, or rombic. ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
]


HelpLinElasticityTensor[] :=
Module[{},
StylePrint["How to use the program LinElasticityTensor", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
StylePrint["The program evaluates the elastic potential and the stress tensor for an isotropic or anisotropic material belonging to the monoclinic, triclinic, or rombic classes.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["To run the program, the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["class = isotropic, monoclinic, triclinic, or rombic (f.e. option = isotropic;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["LinElasticityTensor[class]", "Output", FontFamily -> "Times-Plain", FontSize -> 12]
]



LinElasticityTensor[class_] :=
 Module[{matrixC, sostsymC, Inv, n, ind, b, bfin, gb, hb, PolyFin, tensor, A, t2bis, t2bisFIN, polyquadr1, polyquadr2, polyquadrtot},

 (*Test on input datum class*)
 If[class === Global`isotropic || class === Global`triclinic || class === Global`monoclinic|| class === Global`rombic, Goto[executive],
    StylePrint["ERROR: The user must give one of the choices among isotropic, monoclinic, triclinic, and rombic to the input datum class.",
    "Output", FontSize -> 12, FontFamily -> "Times-Plain"]; Goto[end]];
 Label[executive];

 (*Principal Invariants*)
 matrixC = DeleteCases[Flatten[Table[Which[i <= j, \[Epsilon]*Subscript[C, i, j], i > j, Null], {i, 1, 3}, {j, 1, 3}]], Null];

 Inv1 = Sum[\[Epsilon]*Subscript[C, i, i], {i, 1, 3}];

 matrixC1 = Table[\[Epsilon]*Subscript[C, i, j], {i, 1, 2}, {j, 1, 2}];
 matrixC2 = Table[\[Epsilon]*Subscript[C, i, j], {i, 2, 3}, {j, 2, 3}];
 matrixC3 = Table[\[Epsilon]*Subscript[C, 2*i - 1, 2*j - 1], {i, 1, 2}, {j, 1, 2}];
 sostsymC = {Subscript[C, 2, 1] -> Subscript[C, 1, 2], Subscript[C, 3, 1] -> Subscript[C, 1, 3], Subscript[C, 3, 2] -> Subscript[C, 2, 3]};

 Inv2 = Det[matrixC1] + Det[matrixC2] + Det[matrixC3] /. sostsymC;

 StylePrint["\[FilledDiamond] \ Principal invariants in linear elasticity", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 Which[class === Global`triclinic,
       Inv = matrixC; Print["\[ScriptCapitalB] = ", {Inv /. \[Epsilon] -> 1} //. {C -> "E"}],

       class === Global`monoclinic,
       Inv = Table[Which[h <= 3, \[Epsilon]*Subscript[C, h, h], h == 4, (\[Epsilon]*Subscript[C, 1, 2])^2,
         h == 5, (\[Epsilon]*Subscript[C, 1, 3])^2, h == 6, \[Epsilon]*Subscript[C, 2, h/2], h == 7,
         \[Epsilon]*Subscript[C, 1, 3]*\[Epsilon]*Subscript[C, 1, 2]], {h, 1, 7}];
       Print["\[ScriptCapitalB] = ", Inv /. \[Epsilon] -> 1 //. {C -> "E"}],

       class === Global`rombic,
       Inv = Table[Which[h <= 3, \[Epsilon]*Subscript[C, h, h], h == 4, (\[Epsilon]*Subscript[C, 2, 3])^2,
         h == 5, (\[Epsilon]*Subscript[C, 1, 3])^2, h == 6, (\[Epsilon]*Subscript[C, 1, 2])^2], {h, 1, 6}];
       Print["\[ScriptCapitalB] = ", Inv /. \[Epsilon] -> 1 //. {C -> "E"}],

       class === Global`isotropic,
       Inv = {Inv1, Inv2}; Print["\[ScriptCapitalB] = {\!\(I\_E\)=(", Inv[[1]] /. \[Epsilon] -> 1 //. {C -> "E"},
       "), \!\(II\_E\)=(", Inv[[2]] /. \[Epsilon] -> 1 //. {C -> "E"}, ")}"]];

 (*Homogeneous Polynomial*)
 n = Length[Inv];
 ind = Join[Table[j[i], {i, 1, n}], {0}];
 funzioneIND[{b__}] := Subscript["a", b];
 coef = funzioneIND[Table[l[i], {i, 1, Length[Inv]}]];
 b[0] = coef;
 b[i_] := b[i - 1] //. {l[i] -> ind[[i]] - ind[[i + 1]]};
 bfin = Table[b[i], {i, 1, n}] //. Table[l[s] -> 0, {s, 2, Length[Inv]}];
 gb[{w__}] := Sum[bfin[[n]]*Product[Inv[[i]]^(ind[[i]] - ind[[i + 1]]), {i, 1, n}], w];
 hb = Join[{{ind[[1]], 0, 2}}, Table[{ind[[i]], 0, ind[[i - 1]]}, {i, 2, n}]];
 g[hb];

 PolyFin = gb[hb] - coef //. Table[l[s] -> 0, {s, 1, Length[Inv]}];
 PolyFin2 = Normal[Series[PolyFin, {\[Epsilon], 0, 2}]] /. {\[Epsilon] -> 1};

 StylePrint["\[FilledDiamond] Elastic potential as a function of the principal invariants", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 Print["\[Psi] = ", PolyFin2 //. {C -> "E"}];

 (*Tensor of linear elasticity*)
 tensor = Table[4*"\!\(\(\[Rho]\_*\)\)"*D[PolyFin2, Subscript[C, i, j], Subscript[C, h, k]], {i, 1, 3}, {j, 1, 3},
      {h, 1, 3}, {k, 1, 3}];

(*Voigt' s notation for C*)
 Subscript[A, i_, j_] := If[i <= 3 && j <= 3, tensor[[i,i,j,j]], Null];
 Subscript[A, i_, 4] := If[i <= 3, tensor[[i,i,2,3]], Null];
 Subscript[A, i_, 5] := If[i <= 3, tensor[[i,i,1,3]], Null];
 Subscript[A, i_, 6] := If[i <= 3, tensor[[i,i,1,2]], Null];
 Subscript[A, 4, 4] = tensor[[2,3,2,3]];
 Subscript[A, 5, 5] = tensor[[1,3,1,3]];
 Subscript[A, 6, 6] = tensor[[1,2,1,2]];
 Subscript[A, 4, 5] = tensor[[2,3,1,3]];
 Subscript[A, 4, 6] = tensor[[2,3,1,2]];
 Subscript[A, 5, 6] = tensor[[1,3,1,2]];

 t2bis = DeleteCases[Flatten[Table[Which[i <= j, Subscript["\[DoubleStruckCapitalC]", i, j] == Subscript[A, i, j],
         i > j, Null], {i, 1, 6}, {j, 1, 6}]], Null];
 t2bisFIN = DeleteCases[t2bis, _ == 0];
 StylePrint["\[FilledDiamond] Independent components of linear elasticity tensor \[DoubleStruckCapitalC]=(\!\(\[DoubleStruckCapitalC]\_ijhk\)) in Voigt's Notation", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 Do[Print[t2bisFIN[[i,1]], " = ",
      t2bisFIN[[i,2]]], {i, 1, Length[t2bisFIN]}];

 (*Tensor of Linear Elasticity*)
 tableCcomplete = Table[Subscript[C, i, j] == If[i <= j, Subscript[A, i, j],
        Subscript[A, j, i]], {i, 1, 6}, {j, 1, 6}];

 (*Free Energy on E*)
 polyCE = Sum[(1/(2*"\!\(\(\[Rho]\_*\)\)"))*tableCcomplete[[i,j,2]]*Subscript["E", i]*
       Subscript["E", j], {i, 1, 6}, {j, 1, 6}];

 StylePrint["\[FilledDiamond] Elastic potential in Voigt's Notation", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 Print["\[Psi] = ", polyCE];

 (*Stress Tensor*)
 Subscript[stressT, i_] := Sum[tableCcomplete[[i,j,2]]*Subscript["E", j], {j, 1, 6}];
 StylePrint["\[FilledDiamond] Stress tensor component", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
 Do[Print[Subscript["T", i], " = ", Subscript[stressT, i]], {i, 1, 6}];

 Label[end];

]



(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  PdeEqClass:
 -
 -                 UsagePdeEqClass[]
 -                 HelpPdeEqClass[]
 -                 PdeEqClass[eq_, var_, unk_, point_, unk0_, option_]
 -
 -----------------------------------------------------------------------------------------
*)


UsagePdeEqClass[] :=
Module[{},
StylePrint["Aims of the program PdeEqClass", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["A quasilinear second-order PDE being given, which depends on n rectilinear variables var, a point of \!\(TraditionalForm\`R\^n\), and a known solution unk0, the program determines: ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["\[FilledSmallCircle] the coefficient matrix \!\(TraditionalForm\`A\);", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["\[FilledSmallCircle] the matrix \!\(TraditionalForm\`A\_0\) obtained by evaluating \!\(TraditionalForm\`A\) at \!\(TraditionalForm\`point\) and with respect to the solution \!\(TraditionalForm\`unk0\); ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["\[FilledSmallCircle] the eigenvalues of \!\(TraditionalForm\`A\) if \!\(TraditionalForm\`option = symbolic\);", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["\[FilledSmallCircle] the eigenvalues of \!\(TraditionalForm\`A\_0\) if \!\(TraditionalForm\`option = numeric\);", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["\[FilledSmallCircle] the elliptic, hyperbolic, or parabolic character of the PDE if \!\(TraditionalForm\`option = numeric\).", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["The command line is", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["PdeEqClass[eq, var, unk, point, unk0, option, optionEig],", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["where eq is a second-order quasilinear PDE, where all the terms containing the second derivative appear on the left-hand side, while all the other terms appear on the right-hand side; var is the list of independent variables, unk is the unknown function, point are the coordinates of the point at which the PDE has to be classified, unk0 is the known solution of PDE when it is quasilinear, option is the option for the symbolic or numeric calculation of the eigenvalues of \!\(TraditionalForm\`A\) and \!\(TraditionalForm\`A\_0\), respectively; and finally optionEig controls the round-off errors of the eigenvalues of \!\(TraditionalForm\`A\_0\).", "Output", FontFamily -> "Times-Plain", FontSize -> 12]]
]


HelpPdeEqClass[] :=
Module[{},
StylePrint["How to use the program PdeEqClass", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
StylePrint["The program classifies the quasilinear second-order PDEs.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["To run the program, the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
Print[StyleForm["eq = second-order quasilinear PDE, where all the terms containing the second derivative appear on the left-hand side, while all the other terms appear on the right-hand side (f.e. eq = \!\(v\_\(x, x\) + v\_\(y, y\)\) == 0;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["var = list of independent variables (f.e. var = {x, y};) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["unk = unknown function (f.e. unk = v;) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["point = coordinates of the point at which the PDE has to be classified (f.e. point = {x, y};) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["unk0 = known solution of PDE (f.e. unk0 = \!\(v\_0\);)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["option = symbolic or numeric (f.e. option = numeric;).", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["optionEig = symbolic or numeric (f.e. optionEig = symbolic;).", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
Print[StyleForm["PdeEqClass[eq, var, unk, point, unk0, option, optionEig]", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
]

PdeEqClass[eq_,var_,unk_,point_,unk0_,option_,optionEig_]:=
  Module[{nvar,Aeq,sost,aeq,Matoeff,eig,eignew},

  (*Independent variables*)
  nvar=Length[var];

  (*Test on input data*)
  Which[Length[point]=!=nvar,StylePrint["ERROR: The input data var and point do not have the same dimension!","Output",FontSize->12,FontFamily->"Times-Plain"];
                             Goto[end],
        option=!=Global`symbolic&&option=!=Global`numeric,StylePrint["ERROR: The user must give one of the choices between numeric and symbolic to the input datum option!","Output",FontSize->12,FontFamily->"Times-Plain"];
                                                          Goto[end],
        optionEig=!=Global`symbolic&&optionEig=!=Global`numeric,StylePrint["ERROR: The user must give one of the choices between numeric and symbolic to the input datum optionEig!","Output",FontSize->12,FontFamily->"Times-Plain"];
                                                                Goto[end]];

  (*Coefficient Matrix A*)
  aeq[i_,j_]:=Which[i===j,Coefficient[eq[[1]],Subscript[unk,var[[i]],var[[j]]]],
                    i<j,(1/2)*Coefficient[eq[[1]],Subscript[unk,var[[i]],var[[j]]]],
                    i>j,(1/2)*Coefficient[eq[[1]],Subscript[unk,var[[j]],var[[i]]]]];
  Aeq=Table[aeq[i,j],{i,1,nvar},{j,1,nvar}];
  StylePrint["Coefficient matrix","Output",FontSize->12,FontFamily->"Times-Plain"];
  Print["A = ",MatrixForm[Aeq]];

  (*Coefficient matrix evaluated at point respect to unk0*)
  If[option===Global`numeric,sost=Union[Flatten[Table[{Subscript[unk,var[[i]]]->D[unk0,var[[i]]],
                                        var[[i]]->point[[i]]},{i,1,nvar}]],{unk->unk0}];
                             MatCoeff=Evaluate[Aeq/.sost],
     MatCoeff=Aeq];
  If[option===Global`numeric&&MatCoeff=!=Aeq,Print[StyleForm["Matrix \!\(A\_0\) obtained by evaluating A at point ",
          "Output",FontSize->12,FontFamily->"Times-Plain"],point,StyleForm[" and with respect to the solution ","Output",
          FontSize->12,FontFamily->"Times-Plain"],unk0];
          Print["\!\(A\_0\) = ",MatrixForm[MatCoeff]];
          StylePrint["Eigenvalues of \!\(A\_0\)","Output",FontSize->12,FontFamily->"Times-Plain"],
     StylePrint["Eigenvalues of A","Output",FontSize->12,FontFamily->"Times-Plain"]];

  (*Eigenvalues*)
  eig=Which[optionEig===Global`numeric,Chop[N[Eigenvalues[MatCoeff]]],
            optionEig===Global`symbolic,Eigenvalues[MatCoeff]];

  (*Distinct Eigenvalues*)
  eigD=Intersection[eig];
  Do[Print[Subscript["\[Lambda]",i]," = ",eigD[[i]]],{i,1,Length[eigD]}];

  If[option===Global`numeric,Goto[class],Goto[end]];
  Label[class];

  If[Cases[eig//N,_Complex]=!={},StylePrint["Round-off error: launch again the program choosing optionEig=numeric.","Output",FontSize->12,FontFamily->"Times-Plain"];
                                 Goto[end],Goto[go]];

  Label[go];

  (*Eigenvalues with their signs*)
  eignew=Table[Which[eig[[i]]>0,plus,eig[[i]]===0,0,eig[[i]]<0,minus],{i,1,nvar}];

  (*No symbolic eigenvalues*)
  Which[Complement[eignew,{plus,minus,0}]=!={}&&Intersection[eignew,{0}]=!={0},StylePrint["The equation cannot be classified since not all the eigenvalues have a definite sign!","Output",FontSize->12,FontFamily->"Times-Plain"];
                                                                               Goto[end],
        Complement[eignew,{plus,minus,0}]=!={}&&Intersection[eignew,{0}]==={0},StylePrint["The Pde is parabolic","Output",FontSize->12,FontFamily->"Times-Plain"];
                                                                               Goto[end],
        Complement[eignew,{plus,minus,0}]==={},Goto[class2]];

  Label[class2];
  TmultP=Sort[Join[Table[plus,{nvar-1}],{minus}]];
  TmultM=Sort[Join[Table[minus,{nvar-1}],{plus}]];

  (*Classification*)
  Which[Length[Intersection[eignew,{0}]]===0&&Length[Intersection[eignew,{minus}]]===0,StylePrint["The Pde is elliptic","Output",FontSize->12,FontFamily->"Times-Plain"],
        Length[Intersection[eignew,{0}]]=!=0,StylePrint["The Pde is parabolic","Output",FontSize->12,FontFamily->"Times-Plain"],
        Sort[eignew]===TmultP||Sort[eignew]===TmultM,StylePrint["The Pde is hyperbolic","Output",FontSize->12,FontFamily->"Times-Plain"]];

  Label[end];
]



(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  PdeSysClass:
 -
 -                 UsagePdeSysClass[]
 -                 HelpPdeSysClass[]
 -                 PdeSysClass[sys_, var_, unk_, point_, unk0_, option_]
 -
 -----------------------------------------------------------------------------------------
*)

UsagePdeSysClass[] :=
 Module[{},
 StylePrint["Aims of the program PdeSysClass", "Output", FontFamily -> "Times-Bold",FontSize -> 12];
 Print[StyleForm["Given a quasilinear first-order systems of PDEs in the m unknown functions depending on n (rectilinear) variables var, a point of \!\(R\^n\) and a known solution unk0, the program determines: ", "Output", FontFamily -> "Times-Plain",FontSize -> 12]];
 Print[StyleForm["\[FilledSmallCircle] the matrices \!\(A\_i\) of the coefficients of the system of PDEs;", "Output",FontFamily -> "Times-Plain",FontSize -> 12]];
 Print[StyleForm["\[FilledSmallCircle] the  matrices \!\(A\_i\^0\) obtained by evaluating \!\(A\_i\) at point and with respect to the solution unk0;", "Output", FontFamily -> "Times-Plain", FontSize->12]];
 Print[StyleForm["\[FilledSmallCircle] the characteristic equation det (\!\(IN\_1\)+(\!\(A\^1\)\!\(\()\^\(-1\)\)\)\!\(A\^2\)) = 0;","Output", FontFamily ->"Times-Plain", FontSize->12]];
 Print[StyleForm["\[FilledSmallCircle] the solutions of the characteristic equation.", "Output",FontFamily -> "Times-Plain",FontSize -> 12]];
 Print[StyleForm["Moreover, if option is equal to numeric","Output", FontFamily -> "Times-Plain",FontSize -> 12]];
 Print[StyleForm["\[FilledSmallCircle] the eigenvalues of B = (\!\(A\^1\)\!\(\()\^\(-1\)\)\)\!\(A\^2\);","Output", FontFamily -> "Times-Plain", FontSize -> 12]];
 Print[StyleForm["\[FilledSmallCircle] the algebraic and geometric multiplicities of the distinct eigenvalues of B;", "Output",FontFamily -> "Times-Plain", FontSize -> 12]];
 Print[StyleForm["\[FilledSmallCircle] the elliptic, hyperbolic, totally hyperbolic, or parabolic character of the system of PDEs.","Output", FontFamily -> "Times-Plain",FontSize -> 12]];
 Print[StyleForm["The command raw is","Output", FontFamily -> "Times-Plain",FontSize -> 12]];
 Print[StyleForm["PdeSysClass[sys, var, unk, point, unk0, option],","Output", FontFamily -> "Times-Plain",FontSize -> 12]];
 Print[StyleForm["where sys is the quasilinear first-order system of PDE, where the terms containing the first derivatives appear on the left-hand side, while the right-hand side contains the functions and the known terms; var is the list of the independent variables, unk is the list of the unknown functions, point are the coordinates of the point at which the system of PDEs has to be classified, unk0 is the known solution of the system of PDEs when it is quasilinear, and finally option is an option for the symbolic or numeric calculations of the solutions of the associated characteristic equation.", "Output",FontFamily -> "Times-Plain",FontSize -> 12]]
 ]




HelpPdeSysClass[] :=
  Module[{},
  StylePrint["How to use the program PdeSysClass", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  StylePrint["The program classifies the quasilinear first-order system of PDEs.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["To run the program, the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print[StyleForm["sys = quasilinear first-order PDE system, where the terms containing the first derivatives appear on the left-hand side, while the right-hand side contains the functions and the known terms (f.e. sys = {\!\(v\_x\) + \!\(w\_y\) == 0, \!\(v\_y\) - \!\(w\_x\) == 0};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["var = list of the independent variables (f.e. var = {x, y};) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["unk = list of the unknown functions (f.e. unk = {v, w};) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["point = coordinates of the point at which the system of PDEs has to be classified (f.e. point = {\!\(x\_0\), \!\(y\_0\)};) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["unk0 = known solution of the system of PDEs, when it is quasilinear (f.e. unk0 = {\!\(v\_0\),  \!\(w\_0\)};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["option = symbolic or numeric (f.e. option = numeric;).", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["PdeSysClass[sys, var, unk, point, unk0, option]", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
]



PdeSysClass[sys_, var_, unk_, point_, unk0_, option_] :=
  Module[{nv, mk, as, Asys, sostsys, TaTest, BB, Bfin, Mat, eqnormal,solnormal, BNew, eig, eigve, eigM, fSys, aut, molt, eigind, TaHE, TaHReal, TaHComplex, multAG},
  (*Independent variables*)
  nv = Length[var];
  (*Number of equations and scalar unknown functions*)
  mk = Length[unk];
  (*Test on input data*)
  Which[mk =!= Length[sys], StylePrint["ERROR: The number of the unknowns is not equal to the equation ones!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
                            Goto[end],
        Length[var] =!= Length[point], StylePrint["ERROR: The input data var and point do not have the same dimension!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
                                       Goto[end],
        Length[unk] =!= Length[unk0], StylePrint["ERROR: The input data unk and unk0 do not have the same dimension!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
                                      Goto[end],
        option =!= Global`symbolic && option =!= Global`numeric, StylePrint["ERROR: The user must give one of the choices between symbolic or numeric to the input datum option!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
        Goto[end]];
  (*Evaluation of matrices A^i*)
  as[i_, j_, k_] := Coefficient[sys[[j,1]], Subscript[unk[[k]], var[[i]]]];
  Asys[i_] := Table[as[i, j, k], {j, 1, mk}, {k, 1, mk}];
  StylePrint["Coefficient matrices", "Output", FontSize -> 12,FontFamily -> "Times-Plain"];

  Do[Print[Subscript["A", i], " = ",MatrixForm[Asys[i]]], {i, 1, nv}];
  (*Substitution for the classification*)
  sostsys = Union[Flatten[Table[var[[i]] -> point[[i]], {i, 1, nv}]], Flatten[Table[unk[[j]] -> unk0[[j]], {j, 1, mk}]]];
  (*Coefficient matrices evaluated at point with respect to the solution unk0*)
  TaTest = Table[If[Evaluate[Asys[i] /. sostsys] =!= Asys[i], contat = i,contat = 0], {i, 1, nv}];
  If[TaTest =!= Table[0, {nv}],Print[StyleForm["Coefficient matrices evaluated at point ", "Output", FontSize -> 12, FontFamily -> "Times-Plain"],
     point, StyleForm[" respect to the solution ", "Output", FontSize -> 12, FontFamily -> "Times-Plain"], unk0];
            Do[Print[Subscript["A"^"0", i], " = ", MatrixForm[Evaluate[Asys[i] /. sostsys]]], {i, 1, nv}], Goto[nosost]];
  Label[nosost];
  (*Test on A1*)
  If[Det[Evaluate[Asys[1] /. sostsys]] === 0,Print["ERROR: The matrix \!\(A\_1\) is singular. Try again modifying the order of the variables."];
                                             Goto[end], Goto[char]];
  Label[char];
  B2 = Inverse[Evaluate[Asys[1] /. sostsys]] . Evaluate[Asys[2] /. sostsys];
  Mat = IdentityMatrix[mk]*V[1] + B2;
  (*Characteristic Polynomial*)
  eqnormal = Det[Mat];
  StylePrint["Characteristic equation", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
  Print[Simplify[Chop[eqnormal]] //. {V[1] -> Subscript["N", 1]}, " = 0"];
  solnormal = Intersection[Simplify[Flatten[Solve[eqnormal == 0, V[1]]]]];
  StylePrint["Solutions of the characteristic equation", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
  solnormalnew = Chop[solnormal //. {V[1] -> Subscript["N", 1]}];
  Do[Print[solnormalnew[[i,1]], " = ", solnormalnew[[i,2]]],{i, 1, Length[solnormalnew]}];
  (*Classification of PDE system*)
  If[option === Global`numeric, Goto[class], Goto[end]];
  Label[class];
  BNew = Inverse[Evaluate[Asys[1] /. sostsys]] . Evaluate[Asys[2] /. sostsys];
  (*Eigenvalues and eigenvectors : Algebraic and geometric multiplicity*)
  (*Eigenvalues of BNew, i.e.the solutions N1 with the opposite sign*)
  eig = Chop[Eigenvalues[BNew]];
  eigD = Intersection[eig];
  Print[StyleForm["Eigenvalues of B = (\!\(A\_1\)\!\(\()\^\(-1\)\)\)\!\(A\_2\) = ", "Output", FontSize -> 12, FontFamily -> "Times-Plain"], MatrixForm[BNew]];
  Do[Print[Subscript["\[Lambda]", i], " = ", eigD[[i]]],{i, 1, Length[eigD]}];
  (*Not vanishing eigenvectors of BNew*)
  eigve = DeleteCases[Chop[Eigenvectors[BNew]], Table[0, {i, 1, Length[BNew]}]];
  eigM = Split[Sort[Chop[Eigenvalues[BNew]]], Chop[#1 - #2] === 0 & ];
  (*Grouping of eigenvectors corresponding to the same eigenvalue*)
  fSys[h_] := Table[If[Chop[Simplify[BNew . eigve[[i]] - eigM[[h,1]]*eigve[[i]]]] === Table[0, {i, 1, Length[BNew]}],
                       eigve[[i]], no], {i, 1, Length[eigve]}];
  (*Algebraic Multiplicity*)
  aut[h_] := DeleteCases[fSys[h], no];
  molt = Table[Length[eigM[[h]]], {h, 1, Length[eigM]}];
  (*Geometric Multiplicity*)
  eigind = Table[Length[aut[h]], {h, 1, Length[eigM]}];
  StylePrint["Algebraic and geometric multiplicity of distinct eigenvalues of B", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
  Do[Print[Subscript["\[Lambda]", i], " = ", eigM[[i,1]], ":   AlgMult = ", molt[[i]], "   GeomMult = ", eigind[[i]]],
     {i, 1, Length[eigM]}];
  (*Test on eigenvalues*)
  (*Heads of eigenvalues*)
  TaHE = Table[Head[N[Simplify[eig[[i]]]]], {i, 1, Length[eig]}];
  (*Heads Real*)
  TaHReal = Table[Real, {i, 1, Length[eig]}];
  (*Heads Complex*)
  TaHComplex = Table[Complex, {i, 1, Length[eig]}];
  If[option === Global`numeric && Complement[TaHE, TaHReal, TaHComplex] === {}, Goto[FinClass], Goto[end]];
  Label[FinClass];
  (*Comparison between the algebraic(molt) and geometric (eigind) multiplicities of eigenvalues*)
  multAG = Table[molt[[i]] > eigind[[i]], {i, 1, Length[eigM]}];
  (*Classification*)
  Which[TaHE === TaHComplex, StylePrint["The system of PDEs is elliptic","Output", FontSize -> 12, FontFamily -> "Times-Plain"],
        Length[Intersection[Flatten[eig], {0}]] === 0 && TaHE === TaHReal && Length[Intersection[eig]] === mk && molt === eigind,
        StylePrint["The system of PDEs is totally hyperbolic", "Output", FontSize -> 12, FontFamily -> "Times-Plain"],
        TaHE === TaHReal && molt === eigind, StylePrint["The system of PDEs is hyperbolic", "Output", FontSize -> 12, FontFamily -> "Times-Plain"],
        Intersection[multAG, {True}] =!= {}, StylePrint["The system of PDEs is parabolic", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]];
  Label[end];
]

(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  WavesI:
 -
 -                 UsageWavesI[]
 -                 HelpWavesI[]
 -                 WavesI[sys_, unk_, var_]
 -
 -----------------------------------------------------------------------------------------
*)

UsageWavesI[] :=
  Module[{},
  StylePrint["Aims of the program WavesI", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  Print[StyleForm["The program determines the characteristic equation of a quasilinear first-order system of PDEs as well as the advancing speeds of the characteristic surfaces.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["The command raw is", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["WavesI[sys, unk, var],", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["where sys is the quasilinear first-order system of PDEs, where the terms containing the first derivatives appear on the left-hand side, and the right-hand side contains the functions and the known terms; unk is the list of the unknown functions; and finally var is the list of the independent variables, the first one is the time.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]]
]


HelpWavesI[] :=
  Module[{},
  StylePrint["How to use the program WavesI", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  StylePrint["The program determines the characteristic equation of a quasilinear first-order system of PDEs as well as the advancing speeds of the characteristic surfaces.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["To run the program, the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print[StyleForm["sys = quasilinear first-order system of PDEs, where the terms containing the first derivatives appear on the left-hand side, and the right-hand side contains the functions and the known terms (f.e. sys = {\!\(v\_x - w\_\(\(t\)\(\ \)\) == 0\), \!\(v\_t - t\ w\_x == 0\)};", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["unk = list of the unknown functions (f.e. unk = {v, w};) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["var = list of the independent variables; the first one is the time (f.e. var = {t, x};) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["WavesI[sys, unk, var]", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
]



WavesI[sys_, unk_, var_] :=
  Module[{m, nv, a, A, d, d1, d2, sol, vel, velN},

  m = Length[unk];
  nv = Length[var];

  (*Test on input data*)
  Which[Length[unk] =!= Length[sys], StylePrint["ERROR: The number of the unknowns is not equal to the equation ones!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
                           Goto[end]];

  (*The index i  refers to the differentiation variable, j to the equation and k to the unknown*)
  a[i_, j_, k_] := Coefficient[sys[[j,1]], Subscript[unk[[k]], var[[i]]]];
  A[i_] := Table[a[i, j, k], {j, 1, m}, {k, 1, m}];
  sob2 = Join[{Subscript["f", var[[1]]] -> (-Subscript[c, n])*gradf}, Table[Subscript["f", var[[i + 1]]] -> Subscript["n", i]*gradf, {i, 1, nv - 1}]];
  syscan = Sum[MatrixForm[A[i]]*MatrixForm[Table[Subscript[unk[[j]], var[[i]]], {j, 1, m}]], {i, 1, nv}] == MatrixForm[Table[sys[[k,2]], {k, 1, m}]];
  d = Simplify[Det[Sum[A[i]*Subscript["f", var[[i]]], {i, 1, nv}]]];
  d1 = Simplify[Det[Sum[MatrixForm[A[i]]*MatrixForm[Subscript["f", var[[i]]]], {i, 1, nv}]]];
  d2 = d //. sob2;
  StylePrint["Canonical form of the system", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
  Print[syscan[[1]], " = ", syscan[[2]]];
  Print[""];
  StylePrint["Characteristic equation in a matrix form", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
  Print[d1, " = 0"];
  Print[""];
  If[d === 0, StylePrint["The characteristic equation vanishes identically!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
              Goto[end],
     StylePrint["The explicit characteristic equation", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
     Print[d, " = 0"]];
   sol = Flatten[Solve[d2 == 0, Subscript[c, n]]];
   vel = Flatten[Table[{Simplify[sol[[r,2]]] //. {Sum[Subscript["n", i]^2, {i, 1, nv - 1}] -> 1,
          Sum[-Subscript["n", i]^2, {i, 1, nv - 1}] -> -1}}, {r, 1, Length[sol]}]];
   If[vel === {}, StylePrint["shows that ordinary discontinuity waves do not exist!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
                  Goto[end],
      StylePrint["Normal speed of \[CapitalSigma](t)", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
      velN = Union[vel];
      Do[Print[Subscript["c", "n", i], " = ", velN[[i]]], {i, 1, Length[velN]}]];

  Label[end];
]



(*
 -----------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  WavesII:
 -
 -                 UsageWavesII[]
 -                 HelpWavesII[]
 -                 WavesII[sys_, unk_, var_]
 -
 -----------------------------------------------------------------------------------------
*)

UsageWavesII[] :=
  Module[{},
  StylePrint["Aims of the program WavesII", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  Print[StyleForm["The program determines the characteristic equation of a system of second-order quasilinear PDEs as well as the advancing speeds of the characteristic surfaces and the second derivative jumps. Moreover, in the 2D-case, it distinguishes between transverse and longitudinal waves.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["The command raw is", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["WavesII[sys, unk, var],", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["where sys is the second-order quasilinear system of PDEs, where all the terms containing the second derivative appear on the left-hand side, and all the other terms appear on the right-hand side; unk is the list of the unknown functions; and finally var is the list of the independent variables where the first one is the time.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]]
]


HelpWavesII[] :=
  Module[{},
  StylePrint["How to use the program WavesII", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
  StylePrint["The program determines the characteristic equation of a system of second-order quasilinear PDEs as well as the advancing speeds of the characteristic surfaces and the second derivative jumps. Moreover, in the 2D-case, it distinguishes between transverse and longitudinal waves.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["To run the program, the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print[StyleForm["sys = second-order quasilinear system of PDEs, where all the terms containing the second derivative appear on the left-hand side, and all the other terms appear on the right-hand side (f.e. sys = {\!\(u\_\(t, t\) + u\_\(y, y\) == 0\)};", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["unk = list of the unknown functions (f.e. unk = {u};) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["var = list of the independent variables where the first one is the time (f.e. var = {t, y};) ", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print[StyleForm["WavesII[sys, unk, var]", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
]


WavesII[sys_, unk_, var_] :=
  Module[{m, nv, sob2, sost, sost1, syscan, co, ma, vel, disc, eig, teig, teigeff, teigeffN, lteigeff, teigeff1, ps, pv},

  m = Length[unk];
  nv = Length[var];

  (*Test on input data*)
  Which[Length[unk] =!= Length[sys], StylePrint["ERROR: The number of the unknowns is not equal to the equation ones!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
                           Goto[end]];

  (*The index i refers to the differentiation variable, j to the equation and k to the unknown*)
  sob2 = Flatten[Join[Table[Subscript[unk[[j]], var[[1]], var[[1]]] -> Subscript["c", "n"]^2*Subscript["a", j], {j, 1, m}],
         Table[Subscript[unk[[j]], var[[1]], var[[h]]] -> (-Subscript["c", "n"])*Subscript["n", h - 1]*Subscript["a", j], {j, 1, m}, {h, 1, nv}],
         Table[Subscript[unk[[j]], var[[h]], var[[k]]] -> Subscript["n", h - 1]*Subscript["n", k - 1]*Subscript["a", j], {j, 1, m}, {h, 2, nv}, {k, 2, nv}]]];
  sost = Subscript["n", nv - 1]^2 -> 1 - Sum[Subscript["n", r]^2, {r, 1, nv - 2}];
  sost1 = {Sum[Subscript["n", r]^2, {r, 1, nv - 1}] -> 1, -Sum[Subscript["n", r]^2, {r, 1, nv - 1}] -> -1};
  syscan = sys //. sob2;
  co[i_, j_] := Coefficient[syscan[[i,1]], Subscript["a", j]];
  ma = Table[co[i, j], {i, 1, m}, {j, 1, m}];
  vel = Solve[Det[ma] == 0, Subscript["c", "n"]];
  If[(Simplify[Det[ma]] //. sost1) === 0, StylePrint["The characteristic equation vanishes identically!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
                                         Goto[fine],
     StylePrint["Characteristic equation", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
     Print[Simplify[Det[ma]] //. sost1, " = 0"]];
  If[Length[DeleteCases[vel, {}]] == 0, StylePrint["It is not possible to determine the normal speed of \[CapitalSigma](t)!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
     Goto[fine], Goto[speed]];
  Label[speed];
  Print[""];
  StylePrint["Normal speed of \[CapitalSigma](t)", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
  velN = Union[Simplify[Flatten[vel]] //. sost];
  Do[Print[Subscript["c", "n", i], " = ", velN[[i,2]]], {i, 1, Length[velN]}];
  disc[k_] := Table[(Sum[co[i, j]*Subscript["a", j], {j, 1, m}] /. vel[[k]]) == 0, {i, 1, m}];
  eig[k_] := Flatten[Solve[disc[k], Table[Subscript["a", j], {j, 1, m}]]];
  teig = Table[eig[k], {k, 1, Length[vel]}];
  teigeff = Flatten[DeleteCases[teig, {}]];
  teigeffN = Union[teigeff];
  If[Length[teigeff] == 0, Print[""]; StylePrint["The second derivative jumps have no physical  meaning!", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
     Goto[fine],
     Print[""]; StylePrint["Jump vectors", "Output", FontSize -> 12, FontFamily -> "Times-Plain"];
     Do[Print[teigeffN[[i,1]], " = ", teigeffN[[i,2]]], {i, 1, Length[teigeffN]}]; Goto[onde]];
  Label[onde];
  (*2D-Longitudinal and transverse waves*)
  lteigeff = Length[teigeff];
  teigeff1 = DeleteCases[teig, {}];
  If[Length[teigeff1[[1]]] == 1, Goto[1], Goto[fine]];
  Label[1];
  Do[If[teigeff[[1,1]] === Subscript["a", 1], vet1[k] = teigeff[[k,2]]; vet2[k] = Subscript["a", 2],
       vet1[k] = Subscript["a", 1]; vet2[k] = teigeff[[k,2]]];
       ps[k] = {Subscript["n", 1], Subscript["n", 2]} . {vet1[k], vet2[k]};
       pv[k] = Simplify[Cross[{Subscript["n", 1], Subscript["n", 2], 0}, {vet1[k], vet2[k], 0}]] //. sost1;
       Which[ps[k] =!= 0 && pv[k] === {0, 0, 0}, Print[StyleForm["The velocity ", "Output", FontSize -> 12, FontFamily -> "Times-Plain"], Subscript["c", "n", k], " = ", velN[[k,2]], StyleForm[" refers to a longitudinal wave.", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]],
       ps[k] === 0 && pv[k] =!= {0, 0, 0}, Print[StyleForm["The velocity ", "Output", FontSize -> 12, FontFamily -> "Times-Plain"], Subscript["c", "n", k], " = ", velN[[k,2]], StyleForm[" refers to a transverse wave.", "Output", FontSize -> 12, FontFamily -> "Times-Plain"]]],
       {k, 1, lteigeff}];
  Label[fine];
  Label[end];

]

(*
 --------------------------------------------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  Potential:
 -
 -                 UsagePotential[]
 -                 HelpPotential[]
 -                 Potential[F_, {a_, b_}, {c_, d_}, {val0K_, val1K_, stepK_}, {val0S_, val1S_, stepS_}, points_, option_]
 -
 --------------------------------------------------------------------------------------------------------------------------
*)


UsagePotential[]:=
Module[{},
StylePrint["Aims of the program Potential", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["The program evaluates the complex velocity, the kinetic potential, and Stokes potential, when the complex potential F = F [ z ] is given. Moreover, it plots the level curves of the above potentials.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
StylePrint["The command raw is", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["Potential[F, {a, b}, {c, d}, {val0K, val1K, stepK}, {val0S, val1S, stepS}, points, option],", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["where F is the complex potential, {a, b} and {c, d} define the window in which to plot the level curves \[Phi] = cost, between the values val0k and val1K with the step equal to stepK, and the level ones \[Psi] = cost, between the values val0S and val1S with the step equal to stepS. If option = Kinetic, the program plots only the kinetic level curves; if option = Stokes, it plots only the Stokes level curves. When option = All, both the curves are represented.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]



HelpPotential[]:=
Module[{},
StylePrint["How to use the program Potential", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["The program plots the level curves of the kinetic and Stokes potentials, when the complex potential F = F [ z ] is given.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
StylePrint["To run the program, the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["F = complex potential (f.e. F=2 z;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["{a, b} = extrema on the abscissa axis (f.e. {a, b} = {-1, 1};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["{c, d} = extrema on the ordinate axis  (f.e. {a, b} = {-1, 1};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["{val0K, val1K, stepK} = lowest and highest values of the kinetic level curves; stepK is the step to draw them (f.e. {val0K, val1K, stepK} = {-10,10,0.5};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["{val0S, val1S, stepS} = lowest and highest values of the Stokes level curves; stepS is the step to draw them; (f.e. {val0S, val1S, stepS} = {-10,10,0.5};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["points = number of plot points (f.e. points=Automatic;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["option = plot options (f.e. option=All;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["Potential[F, {a, b}, {c, d}, {val0K, val1K, stepK}, {val0S, val1S, stepS}, points, option]", "Output", FontFamily -> "Times-Plain", FontSize -> 12];

]

     
 Potential[F_, {a_, b_}, {c_, d_}, {val0K_, val1K_, stepK_}, {val0S_, val1S_, stepS_}, points_, option_] :=
 Module[{POT\[Psi], POT\[CurlyPhi], G1, G2, derG1, derG2, plpsiD, plpsi, plphi, tapsiD, tapsi, taphi},
  $Assumptions = {Global`z = UNKx + I*UNKy};
  $Assumptions = {UNKx \[Element] Reals, UNKy \[Element] Reals};
  sost = {ArcTan[r_, s_] -> ArcTan[s/r]};
  
  POT\[Psi] = ComplexExpand[Im[F], TargetFunctions -> {Im, Re}];
  POT\[CurlyPhi] = ComplexExpand[Re[F], TargetFunctions -> {Im, Re}];
 
  G1 = (FullSimplify[POT\[CurlyPhi]] //. sost) //. {UNKx -> Global`x, 
     UNKy -> Global`y};
  G2 = (FullSimplify[POT\[Psi]] //. sost) //. {UNKx -> Global`x, 
     UNKy -> Global`y};
  G = G1 + I*G2;
  
  Print[StyleForm["Complex Potential", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print["F = ", G];
  Print[""];
  derG1 = (FullSimplify[ComplexExpand[D[POT\[CurlyPhi], UNKx], TargetFunctions -> {Im, Re}]] //. sost) //. {UNKx -> Global`x, UNKy -> Global`y};
  derG2 = (FullSimplify[ComplexExpand[D[POT\[Psi], UNKx], TargetFunctions -> {Im, Re}]] //. sost) //. {UNKx -> Global`x, UNKy -> Global`y};
  derG = derG1 + I*derG2;
  Print[StyleForm["Complex velocity", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
  Print["V = ", derG //. sost];
  
  plpsiD[i_] := ContourPlot[(Evaluate[POT\[Psi] // N]) == i, {UNKx, a, b}, {UNKy, c, d}, PlotPoints -> points, ContourStyle -> {GrayLevel[0], Dashing[{0.01}]}, Axes -> True, Frame -> False];
  plpsi[i_] := ContourPlot[(Evaluate[POT\[Psi] // N]) == i, {UNKx, a, b}, {UNKy, c, d}, PlotPoints -> points, ContourStyle -> {GrayLevel[0]}, Axes -> True, Frame -> False];
  plphi[i_] := ContourPlot[(Evaluate[POT\[CurlyPhi] // N]) == i, {UNKx, a, b}, {UNKy, c, d}, PlotPoints -> points, ContourStyle -> {GrayLevel[0]}, Axes -> True, Frame -> False];
  
  tapsiD = N[Table[plpsiD[value], {value, val0S, val1S, stepS}]];
  tapsi = N[Table[plpsi[value], {value, val0S, val1S, stepS}]];
  taphi = N[Table[plphi[value], {value, val0K, val1K, stepK}]];
  Which[option === Global`Stokes, Print[Show[tapsi, AxesLabel -> {"x", "y"}, AxesOrigin -> {0, 0}]]; Print[""]; Print[StyleForm["Level line of the Stokes potential \[Psi]", 
     "Output", FontFamily -> "Times-Plain", FontSize -> 12]],
   option === Global`Kinetic, Print[Show[taphi, AxesLabel -> {"x", "y"}, AxesOrigin -> {0, 0}]]; Print[""]; Print[StyleForm["Level line of the kinetic potential \[CurlyPhi]", 
     "Output", FontFamily -> "Times-Plain", FontSize -> 12]],
   option === All, Print[Show[{tapsiD, taphi}, AxesLabel -> {"x", "y"}, AxesOrigin -> {0, 0}, AxesStyle -> Thick]]; 
   
   Print[""];
   Print[StyleForm["----- Level line of the Stokes potential \[Psi]", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
   Print[StyleForm["——— Level line of the kinetic potential \[CurlyPhi]", "Output", FontFamily -> "Times-Plain", FontSize -> 12]]];
   
   ]


(*
 --------------------------------------------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  Wing:
 -
 -                 UsageWing[]
 -                 HelpWing[]
 -                 Wing[xc_, yc_]
 -
 --------------------------------------------------------------------------------------------------------------------------
*)

UsageWing[]:=
Module[{},
StylePrint["Aims of the program Wing", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["The program draws the curve \[CapitalGamma] which represents the image under the Joukowsky map of the unit circle with a given center C \[Congruent] (xc, yc) as well as the straight line r containing the point T' \[Congruent] (2l, 0) and forming an angle \[Beta] with respect to Ox. Moreover, it supplies the parametric equations of both \[CapitalGamma] and r and the values of l and \[Beta].","Output", FontFamily->"Times-Plain",FontSize->12]];
StylePrint["The command raw is", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["Wing[xc, yc],", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["where xc and yc are the abscissa and ordinate of the center C of the unit circle, respectively.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]


HelpWing[]:=
Module[{},
StylePrint["How to use the program Wing", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["The program allows one to define the curve \[CapitalGamma], which is the image in the plane z of the unit circle, by Joukowsky's map.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
StylePrint["To run the program the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["xc = abscissa of the center of the unit circle (f.e. xc = 0.1;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["yc = ordinate of the center of the unit circle (f.e. yc = 0.2;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["Wing[xc, yc]", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]


Wing[xc_, yc_] := 
 Module[{l, beta, \[Zeta], abs\[Zeta], co\[Zeta], w, w1, w2, w1n, w2n, xT, plr, pl2, plAngleT, plC}, 
 
 l = xc + Sqrt[1 - yc^2];
 beta = ArcSin[yc/1];
 (*Joukowsky's Transformation*)
  $Assumptions = {\[Xi] \[Element] Reals, \[Eta] \[Element] Reals};
  \[Zeta] = \[Xi] + I*\[Eta];
  abs\[Zeta] = Sqrt[\[Xi]^2 + \[Eta]^2];
  co\[Zeta] = \[Xi] - I*\[Eta];
  w = \[Zeta] + l^2/abs\[Zeta]^2 co\[Zeta];
  w1 = ComplexExpand[Re[w]];
  w2 = ComplexExpand[Im[w]];
 
 (*Parametric Equations of the Circumference Image*)
  w1n = Chop[w1 /. {\[Xi] -> xc + Cos[t], \[Eta] -> yc + Sin[t]}];
  w2n = Chop[w2 /. {\[Xi] -> xc + Cos[t], \[Eta] -> yc + Sin[t]}];
  StylePrint["Parametric equations of the wing", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print["x(t) = ", w1n //. t -> Global`t];
  Print["y(t) = ", w2n //. t -> Global`t];
  (*Parametric Equations of the Straight Line including T' with a slope beta*)
  xT = Chop[w1 /. {\[Xi] -> l, \[Eta] -> 0}];
  Print[];
  StylePrint[ "Parametric equations of the straight line r containing the point T' \[Congruent] (2l, 0) and forming an angle \[Beta] with Ox-axis", 
   "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print["x(t) = ", Chop[(xT - Cos[beta]*t)] //. t -> Global`t];
  Print["y(t) = ", Chop[(Sin[beta]*t)] //. t -> Global`t];
  plr = ParametricPlot[Evaluate[{xT - Cos[beta]*t, Sin[beta]*t}], {t, 0, 10}, PlotStyle -> Dashing[{0.01, 0.01}], AspectRatio -> Automatic, 
    PlotRange -> {{-2.5*l, 2.2*l}, {-1.5, 1.5}}, PlotPoints -> 40];
  
  (*Graphics*)
  (*Wing Profile*)
  pl2 = ParametricPlot[Evaluate[{w1n, w2n}], {t, 0, 2*Pi}, AspectRatio -> Automatic, PlotRange -> All, PlotPoints -> 40];
  plAngleT = Graphics[{{Thickness[.008], RGBColor[0, 0, 1], Circle[{l/2, 0}, 2 l/3, {0, beta}]}, {Thickness[0.2], RGBColor[0, 0, 1], 
      Text["\[Beta]", {0.8, 0.05}]}, {Thickness[0.9], RGBColor[0, 0, 1], Point[{2 l, 0}]}, {Thickness[0.2], RGBColor[0, 0, 1], Text["T'", {2 l + 0.1, -0.15}]}}, 
    AspectRatio -> Automatic, Axes -> Automatic];
  plC = Graphics[{{RGBColor[1, 0, 0], Circle[{xc, yc}, 1]}, {Thickness[0.2], RGBColor[1, 0, 0], Point[{xc, yc}]}, {RGBColor[1, 0, 0], 
      Text["C", {xc + 0.1, yc + 0.1}]}, {Thickness[0.2], RGBColor[1, 0, 0], Point[{l, 0}]}, {RGBColor[1, 0, 0], Text["T", {l + 0.1, -0.15}]}}, AspectRatio -> Automatic, 
    Axes -> Automatic];
  Print[Show[plr, pl2, plAngleT, plC, AxesLabel -> {"x", "y"}, AxesStyle -> Thickness[0.01]]];
  StylePrint["——— Wing profile.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  StylePrint["----- Straight line r.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print["\[Beta] = ", beta*180/Pi, "°"];
  Print["l = ", l];
  
  ]


(*
 --------------------------------------------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  Joukowsky:
 -
 -                 UsageJoukowsky[]
 -                 HelpJoukowsky[]
 -                 Joukowsky[phi_,xc_,yc_,{a_,b_},{c_,d_},indata_,steps_,T1_,T2_]
 -
 --------------------------------------------------------------------------------------------------------------------------
*)


UsageJoukowsky[]:=
Module[{},
StylePrint["Aims of the program Joukowsky", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["The program plots the streamlines around the cylinder and the corresponding wing profile. Moreover, it allows one to change the attack angle as well as the coordinates of the center of the cylinder section.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
StylePrint["The command raw is", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["Joukowsky[\[CurlyPhi], xc, yc, {a, b}, {c, d}, indata, steps, T1, T2],", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["where \[CurlyPhi] is the attack angle, xc and yc are the coordinates of the center of the unit circle, {a, b} and {c, d} define the graphic window, indata is the set of the initial data with respect to which we wish to plot the streamlines in the interval [T1, T2], and steps denotes the number of steps of the numerical integration.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]


HelpJoukowsky[]:=
Module[{},
StylePrint["How to use the program Joukowsky", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["The program plots the streamlines around the cylinder and the corresponding wing profile. Moreover, it allows one to change the attack angle as well as the coordinates of the center of the cylinder section.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
StylePrint["To run the program the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["\[CurlyPhi] = angle of attack (f.e. phi = 10 Degree;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["xc = abscissa of the center C of the unit circle (f.e. xc = -0.2;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["yc = ordinate of the center C of the unit circle (f.e. yc = 0.1;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["{a, b}, {c, d} = definition of the graphic window in which to represent the streamlines (f.e. {a, b} = {-3, 3}; {c, d} = {-1.5, 1.6};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["indata = set of initial data for numerically integrating the differential system (f.e. indata=Join[Table[{-3,-1.4+1.4i/15},{i,0,15}],Table[{-3+2*i/5,-1.5},{i,0,5}]];)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["steps = steps of the numerical integration (f.e. steps = 10000;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["T1 = left-bound of time interval (T1, T2) (f.e. T1 = 0;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["T2 = right-bound of time interval (T1, T2) (f.e. T2 = 10;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["Joukowsky[\[CurlyPhi], xc, yc, {a, b}, {c, d}, indata, steps, T1, T2]", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]



Joukowsky[phi_, xc_, yc_, {a_, b_}, {c_, d_}, indata_, steps_, T1_, T2_] :=
  Module[{(*l, theta, beta,b1 s1, r, gamma,z,\[Zeta],re\[Zeta],im\[Zeta],c\[Zeta],rec\[Zeta],imc\[Zeta],arg\[Zeta],abs\[Zeta],sterm,resterm,imsterm,tt, JP,\[CurlyPhi],\[Psi], vex,vey, eq1,eq2,,sol1, 
pl,plo,u1,u2,absz,coz,wJ,w1,w2*)},
  l = xc + Sqrt[1 - yc^2];
  theta = Pi - ArcTan[Abs[yc/xc]];
  beta = ArcSin[yc/l]; 
b1 = Sqrt[xc^2 + yc^2];
  s1 = Cos[phi] - I*Sin[phi];
  r = Cos[theta] + I*Sin[theta];
  gamma = N[4*Pi*Sin[beta + phi]];
  z = x + I*y;
Im[x]^=0;
Re[x]^=x;
Im[y]^=0;
Re[y]^=y;
\[Zeta]=ComplexExpand[(z - b1*r)*s1];
re\[Zeta]=Re[\[Zeta]];
im\[Zeta]=Im[\[Zeta]];
c\[Zeta]=ComplexExpand[Conjugate[\[Zeta]]];
rec\[Zeta]=Re[c\[Zeta]];
imc\[Zeta]=Im[c\[Zeta]];
arg\[Zeta]=ArcTan[imc\[Zeta]/rec\[Zeta]];
abs\[Zeta]=Sqrt[(Re[\[Zeta]])^2+(Im[\[Zeta]])^2];
sterm=c\[Zeta]/abs\[Zeta]^2;
resterm=Re[sterm];
imsterm=Im[sterm];
tt=gamma/(2\[Pi]) (I Log[abs\[Zeta]]-arg\[Zeta]);


(*Joukovski Complex Potential JP=(z - b1*r)*s1 + 1/((z - b1*r)s1) + (I*gamma*Log[(z - b1*r)*s1])/(2*Pi) = \[Zeta]+Conjugate[\[Zeta]]/\[Zeta]^2+gamma/(2\[Pi])(I Log|\[Zeta]|-Arg[\[Zeta]])*)
JP=ComplexExpand[\[Zeta] + c\[Zeta]/abs\[Zeta]^2+ gamma/(2 \[Pi])*(I Log[abs\[Zeta]]-arg\[Zeta])];
\[CurlyPhi] = re\[Zeta]+rec\[Zeta]/abs\[Zeta]^2-gamma/(2 \[Pi]) arg\[Zeta];
\[Psi] = -(im\[Zeta]+imc\[Zeta]/abs\[Zeta]^2+gamma/(2 \[Pi]) Log[abs\[Zeta]]);


  (*Velocity Components*)
  vex = D[\[Psi], x] /. {x -> x[t], y -> y[t]};
  vey = D[\[Psi], y] /. {x -> x[t], y -> y[t]};
  eq1 = Simplify[Derivative[1][x][t] == -vey];
  eq2 = Simplify[Derivative[1][y][t] == vex];
  sol1[i_] := Flatten[NDSolve[{eq1, eq2, x[T1] == indata[[i,1]], y[T1] == indata[[i,2]]}, {x, y}, {t, T1, T2}, MaxSteps -> steps]];
  pl[i_] := ParametricPlot[Evaluate[{x[t], y[t]} /. sol1[i]], {t, T1, T2}, AspectRatio -> Automatic, PlotRange -> {{a, b}, {c, d}}, PlotPoints -> 60];
  plo = Table[pl[i], {i, 1, Length[indata]}];

  (*Streamlines around the Cylinder*)
  StylePrint["Streamlines around the cylinder.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
 Print[Show[plo]];
Do[
  u1[i]= x[t] /. sol1[i];
  u2[i]= y[t] /. sol1[i],{i,1,Length[indata]}];

  (*Joukowsky's Transformation*)
  $Assumptions = {\[Xi] \[Element] Reals, \[Eta] \[Element] Reals};
  \[Zeta] = \[Xi] + I*\[Eta];
  abs\[Zeta] = Sqrt[\[Xi]^2 + \[Eta]^2];
  co\[Zeta] = \[Xi] - I*\[Eta];
  wJ= \[Zeta] + l^2/abs\[Zeta]^2 co\[Zeta];
  w1 = \[Xi]+l^2/abs\[Zeta]^2 \[Xi];
  w2 = \[Eta]-l^2/abs\[Zeta]^2 \[Eta];
  
  
Do[
  w1x[i]= w1 /. {\[Xi] -> u1[i], \[Eta] -> u2[i]};
  w2y[i]= w2 /. {\[Xi] -> u1[i], \[Eta] -> u2[i]};
  p1[i]= ParametricPlot[Evaluate[{Evaluate[w1x[i]], Evaluate[w2y[i]]}], {t, T1, T2}, PlotRange -> {{a, b}, {c, d}}],{i,1,Length[indata]}];

  (*Streamlines around the Wing profile*)
  Print[""];
  StylePrint["Streamlines around the wing profile.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
  Print[Show[Table[p1[i], {i, 1, Length[indata]}], DisplayFunction -> $DisplayFunction]];

]




(*
 --------------------------------------------------------------------------------------------------------------------------
 -
 -   T H E   P R O G R A M  JoukowskyMap:
 -
 -                 UsageJoukowskyMap[]
 -                 HelpJoukowskyMap[]
 -                 JoukowskyMap[xc_,yc_,curve_,data_,range_,option_]
 -
 --------------------------------------------------------------------------------------------------------------------------
*)

UsageJoukowskyMap[]:=
Module[{},
StylePrint["Aims of the program JoukowskyMap", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["When a (closed, open, or piecewise defined) curve or a point set \[CapitalGamma] are assigned, the program determines the corresponding image under Joukowsky's transformation of a unit circle with a given center C\[Congruent](xc,yc).", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
StylePrint["The command raw is", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["JoukowskyMap[xc, yc, curve, data, range, option],", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["where xc and yc are the coordinates of the center C, curve is the option relative to the curve or a point set \[CapitalGamma], data are the parametric equations of \[CapitalGamma] or the point list, range is the variability range of the parameter t in data or null if \[CapitalGamma] is a point set, option is the option to determine the parametric equations of the image curve or the coordinates of the image points.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["For the input datum curve the following choices are possible:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["\!\(\*StyleBox[\"closed\", FontWeight->\"Bold\"]\) if \[CapitalGamma] is a closed curve defined by two parametric equations x = x(t), y = y(t), where the parameter t varies in the interval [\!\(\[Tau]\_1\),\!\(\[Tau]\_2\)], which in input is given by range = {\!\(\[Tau]\_1\),\!\(\[Tau]\_2\)}", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["\!\(\*StyleBox[\"open\", FontWeight->\"Bold\"]\) if \[CapitalGamma] is an open curve defined by two parametric equations x = x(t), y = y(t), where the parameter t varies in the interval [\!\(\[Tau]\_1\),\!\(\[Tau]\_2\)], which in input is given by range = {\!\(\[Tau]\_1\),\!\(\[Tau]\_2\)};", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["\!\(\*StyleBox[\"piecewise\", FontWeight->\"Bold\"]\) if \[CapitalGamma] is a curve defined by more parametric equations \!\(x\_i\) = \!\(x\_i\)(t),  \!\(y\_i\) = \!\(y\_i\)(t), i = 1,\[CenterEllipsis], n, where the parameter t varies in one or more intervals, which are given by range = {\!\(\[Tau]\_1\),\!\(\[Tau]\_2\)}; or range = {{\!\(\[Tau]\_\(1,1\)\),\!\(\[Tau]\_\(1,2\)\)},...,{\!\(\[Tau]\_\(n,1\)\),\!\(\[Tau]\_\(n,2\)\)}};", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["\!\(\*StyleBox[\"points\", FontWeight->\"Bold\"]\) if \[CapitalGamma] is the finite set of points; in this case, range = null.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["If option = parametric, the program supplies the parametric equations or the coordinates of image of \[CapitalGamma]; for a different choice of option, it shows only the plot.", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]



HelpJoukowskyMap[]:=
Module[{},
StylePrint["How to use the program JoukowskyMap", "Output", FontFamily -> "Times-Bold", FontSize -> 12];
Print[StyleForm["The program determines the corresponding image \[CapitalGamma]' by Joukowsky's transformation of a unit circle with a given center C\[Congruent](xc,yc), when a curve or a point set \[CapitalGamma] are assigned.", "Output", FontFamily -> "Times-Plain", FontSize -> 12]];
StylePrint["To run the program the following input data have to be typed:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["xc = abscissa of the center C of the unit circle (f.e. xc = 0.1;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["yc = ordinate of the center C of the unit circle (f.e. yc = 0.2;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["curve = option relative to the curve \[CapitalGamma] for which the choices closed, open, piecewise, and points are possible (f.e. curve = closed;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["data = parametric equations of \[CapitalGamma] or the point list (f.e. data = {xc+2Cos[t],yc+Sin[t]};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["range = range of the parameter t in data or null, if \[CapitalGamma] is a point set (f.e. range = {0,2\[Pi]};)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["option = if option is equal to parametric, the program supplies the parametric equations or the coordinates of \[CapitalGamma]'; for a different choice, it shows only the plot (f.e. option = null;)", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
StylePrint["JoukowskyMap[xc, yc, curve, data, range, option]", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
]


JoukowskyMap[xc_, yc_, curve_, data_, range_, option_] := 
 Module[{(*l,beta,circle,Tpoint,plN1,plN2,MapPoints,symbol1,symbol2,MLP,pl1T,pl2T,pl1,pl2*)},
  $Assumptions = {Global`t \[Element] Reals, \[Xi] \[Element] Reals, \[Eta] \[Element] Reals};
  l = xc + Sqrt[1 - yc^2];
  beta = ArcSin[yc/1];
  (*Joukowsky's Transformation*)
  \[Zeta] = \[Xi] + I*\[Eta];
  (*Im[\[Xi]]^=0;
  Im[\[Eta]]^=0;*)
  absz = Sqrt[\[Xi]^2 + \[Eta]^2];
  coz = \[Xi] - I*\[Eta];
  w = z + l^2/absz^2 coz;
  w1 = \[Xi] + l^2/absz^2 \[Xi];
  w2 = \[Eta] - l^2/absz^2 \[Eta];
  circle = Graphics[Circle[{xc, yc}, 1], AspectRatio -> Automatic, Axes -> Automatic];
  Which[curve === Global`closed, Goto[closed], curve === Global`open, 
   Goto[open], curve === Global`points, Goto[list], 
   curve === Global`piecewise, Goto[piecewise]];
  
  Label[piecewise];
  Which[{Length[data]} === Dimensions[data] && Length[data] === 2 && {Length[range]} === Dimensions[range] && 
    Length[range] === 2, Goto[open], {Length[data]} =!= Dimensions[data] && Length[Dimensions[data]] > 1 && {Length[range]} === 
     Dimensions[range] && Length[range] === 2, Goto[sp1range], Length[data] === Length[range] === Dimensions[data][[1]], 
   Goto[plusrange]];
  
  Label[sp1range];
  Tpoint = Table[{data[[i, 1]] //. Global`t -> range[[1]], data[[i, 2]] //. Global`t -> range[[1]]}, {i, 1, Length[data]}];
  plN1[i_] := ParametricPlot[{data[[i, 1]], data[[i, 2]]}, {Global`t, range[[1]], range[[2]]}, AspectRatio -> Automatic, PlotRange -> All, 
    PlotPoints -> 40, PlotStyle -> {RGBColor[1, 0, 0]}];
  plN2[i_] := ParametricPlot[{w1 /. {\[Xi] -> data[[i, 1]], \[Eta] -> data[[i, 2]]}, 
     w2 /. {\[Xi] -> data[[i, 1]], \[Eta] -> data[[i, 2]]}}, {Global`t, range[[1]], range[[2]]}, 
    AspectRatio -> Automatic, PlotRange -> All, PlotPoints -> 40, PlotStyle -> {RGBColor[0, 0, 1]}];
  Goto[graphics];
  
  Label[plusrange];
  Tpoint = 
   Table[{data[[i, 1]] //. Global`t -> range[[i, 1]], data[[i, 2]] //. Global`t -> range[[i, 1]]}, {i, 1, 
     Length[data]}];
  plN1[i_] := 
   ParametricPlot[{data[[i, 1]], data[[i, 2]]}, {Global`t, range[[i, 1]], range[[i, 2]]}, AspectRatio -> Automatic, 
    PlotRange -> All, PlotPoints -> 40, PlotStyle -> {RGBColor[1, 0, 0]}];
  plN2[i_] := 
   ParametricPlot[{w1 /. {\[Xi] -> data[[i, 1]], \[Eta] -> data[[i, 2]]}, 
     w2 /. {\[Xi] -> data[[i, 1]], \[Eta] -> data[[i, 2]]}}, {Global`t, range[[i, 1]], range[[i, 2]]}, 
    AspectRatio -> Automatic, PlotRange -> All, PlotPoints -> 40, PlotStyle -> {RGBColor[0, 0, 1]}];
  Goto[graphics];
  Label[graphics];
  MapPoints = 
   Table[{w1 /. {\[Xi] -> Tpoint[[i, 1]], \[Eta] -> Tpoint[[i, 2]]}, w2 /. {\[Xi] -> Tpoint[[i, 1]], \[Eta] -> Tpoint[[i, 2]]}}, {i, 
     1, Length[data]}];
  
  symbol1 = 
   MapThread[Text[#1, #2 + {.1, .1}] &, {Table[Subscript["P", i], {i, 1, Length[data]}], Tpoint}];
  symbol2 = 
   MapThread[Text[#1, #2 + {.1, .1}] &, {Table[Derivative[1][Subscript["P", i]], {i, 1, Length[data]}], MapPoints}];
  MLP = {Graphics[{symbol1, symbol2}], ListPlot[{Tpoint, MapPoints}, PlotStyle -> {Red, Blue}, PlotMarkers -> Automatic]};
  
  pl1T = Table[plN1[i], {i, 1, Length[data]}];
  pl2T = Table[plN2[i], {i, 1, Length[data]}];
  If[option === Global`parametric, StylePrint["Parametric equation of the image curve:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
   Do[Print[Subscript["x", i], "(t) = ", Chop[FullSimplify[ComplexExpand[w1 /. {\[Xi] -> data[[i, 1]], \[Eta] -> data[[i, 2]]}]]]];
    Print[Subscript["y", i], "(t) = ", Chop[FullSimplify[ComplexExpand[w2 /. {\[Xi] -> data[[i, 1]], \[Eta] -> data[[i, 2]]}]]]], {i,
      1, Length[data]}], Goto[graphic]];
  Label[graphic];
  Print[Show[circle, MLP, pl1T, pl2T, AspectRatio -> Automatic, ImageSize -> Automatic]];
  Goto[end];
    
  Label[closed];
  pl1 = ParametricPlot[{data[[1]], data[[2]]}, {Global`t, range[[1]], range[[2]]}, AspectRatio -> Automatic, PlotRange -> All, PlotPoints -> 40, PlotStyle -> {RGBColor[1, 0, 0]}];
  pl2 = ParametricPlot[{Chop[FullSimplify[ComplexExpand[w1 /. {\[Xi] -> data[[1]], \[Eta] -> data[[2]]}]]], 
  Chop[FullSimplify[ComplexExpand[w2 /. {\[Xi] -> data[[1]], \[Eta] -> data[[2]]}]]]}, {Global`t, range[[1]], range[[2]]}, 
    AspectRatio -> Automatic, PlotRange -> All, PlotPoints -> 40, PlotStyle -> {RGBColor[0, 0, 1]}];
  If[option === Global`parametric, StylePrint["Parametric equation of the image curve:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
   Print["x(t) = ", w1 /. {\[Xi] -> data[[1]], \[Eta] -> data[[2]]}];
   Print["y(t) = ", w2 /. {\[Xi] -> data[[1]], \[Eta] -> data[[2]]}], Goto[graphic]];
  Label[graphic];
  Print[Show[circle, pl1, pl2, AspectRatio -> Automatic, ImageSize -> Automatic]];
  Goto[end];
  
  Label[open];
  pl1 = ParametricPlot[{data[[1]], data[[2]]}, {Global`t, range[[1]], range[[2]]}, AspectRatio -> Automatic, PlotRange -> All, 
    PlotPoints -> 40, PlotStyle -> {RGBColor[1, 0, 0]}];
  pl2 = ParametricPlot[{w1 /. {\[Xi] -> data[[1]], \[Eta] -> data[[2]]}, w2 /. {\[Xi] -> data[[1]], \[Eta] -> data[[2]]}}, {Global`t, 
     range[[1]], range[[2]]}, AspectRatio -> Automatic, PlotRange -> All, PlotPoints -> 90, PlotStyle -> {RGBColor[0, 0, 1]}];
  Tpoint = Table[{data[[1]] //. Global`t -> range[[i]], data[[2]] //. Global`t -> range[[i]]}, {i, 1, Length[data]}];
  MapPoints = Table[{Chop[FullSimplify[ComplexExpand[w1 /. {\[Xi] -> Tpoint[[i, 1]], \[Eta] -> Tpoint[[i, 2]]}]]], Chop[FullSimplify[
       ComplexExpand[w2 /. {\[Xi] -> Tpoint[[i, 1]], \[Eta] -> Tpoint[[i, 2]]}]]]}, {i, 1, Length[data]}];
  
  symbol1 = MapThread[Text[#1, #2 + {.1, .1}] &, {Table[Subscript["P", i], {i, 1, Length[data]}], Tpoint}];
  symbol2 = MapThread[Text[#1, #2 + {.1, .1}] &, {Table[Derivative[1][Subscript["P", i]], {i, 1, Length[data]}], MapPoints}];
  MLP = {Graphics[{symbol1, symbol2}], ListPlot[{Tpoint, MapPoints}, PlotStyle -> {Red, Blue}, PlotMarkers -> Automatic]};
  If[option === Global`parametric, StylePrint["Parametric equations of the image curve:", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
   Print["x(t) = ", Chop[FullSimplify[ComplexExpand[w1 /. {\[Xi] -> data[[1]], \[Eta] -> data[[2]]}]]]];
   Print["y(t) = ", Chop[FullSimplify[ComplexExpand[w2 /. {\[Xi] -> data[[1]], \[Eta] -> data[[2]]}]]]], Goto[graphic]];
  Label[graphic];
  Print[Show[MLP, circle, pl1, pl2, AspectRatio -> Automatic, Axes -> Automatic, ImageSize -> Automatic]];
  Goto[end];
  
  Label[list];
  MapPoints = Table[{Chop[FullSimplify[ComplexExpand[w1 /. {\[Xi] -> data[[i, 1]], \[Eta] -> data[[i, 2]]}]]], 
     Chop[FullSimplify[ComplexExpand[w2 /. {\[Xi] -> data[[i, 1]], \[Eta] -> data[[i, 2]]}]]]}, {i, 1, Length[data]}];
  symbol1 = MapThread[Text[#1, #2 + {.15, .15}] &, {Table[Subscript["P", i], {i, 1, Length[data]}], data}];
  symbol2 = MapThread[Text[#1, #2 + {.15, .15}] &, {Table[Derivative[1][Subscript["P", i]], {i, 1, Length[data]}], MapPoints}];
  If[option === Global`parametric, StylePrint["Point list and its image", "Output", FontFamily -> "Times-Plain", FontSize -> 12];
   TabPoints = Table[Subscript["P", i] \[Congruent] data[[i]], {i, 1, Length[data]}];
   TabMapPoints = Table[Derivative[1][Subscript["P", i]] \[Congruent] MapPoints[[i]], {i, 1, Length[data]}]; 
   Print[TableForm[TabPoints, TableDepth -> 2], 
   " \!\(\( \[RightArrow] \&\(\(\\ \\ \\ \\ \\ \\ \\ \\ \
\)\(\(Joukowsky'\) \(s\)\(\\ \)\(Map\)\(\\ \\ \\ \\ \\ \\ \\ \\ \\ \\ \
\\ \)\)\)\)\) ", TableForm[TabMapPoints, TableDepth -> 2]], 
  Goto[graphic]];
  Label[graphic];
  MLP = {Graphics[{symbol1, symbol2}], ListPlot[{data, MapPoints}, PlotStyle -> {Red, Blue}, PlotMarkers -> Automatic]};
  Print[Show[circle, MLP, AspectRatio -> Automatic, PlotRange -> All, ImageSize -> Automatic]];
  
  Label[end];
  Print[Style["———  ", "Text"], Style["Circle of unit radius.", "Text", FontSize -> 12, FontFamily -> "Times-Plain"]]; 
  Print[Style["———  ", "Text", Red], Style["Input curve or list of points.", "Text", FontSize -> 12, FontFamily -> "Times-Plain"]]; 
  Print[Style["———  ", "Text", Blue], Style["Output curve or list of points.", "Text", FontSize -> 12, FontFamily -> "Times-Plain"]];
  
  ]

(* = = = = = = = = = = = = = = =  E N D  P R O G R A M S  = = = = = = = = = = = = = = = = *)


End[]


(* = = = = = = = = = = = = = = =  E N D  P A C K A G E  = = = = = = = = = = = = = = = = = *)


EndPackage[]
