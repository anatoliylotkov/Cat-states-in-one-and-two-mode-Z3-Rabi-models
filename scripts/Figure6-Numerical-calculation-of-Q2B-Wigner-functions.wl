(* Import package *)
AppendTo[$Path, FileNameJoin[{Directory[], "src"}]]
<< WignerFunctions`
<< Z3Rabi`

(* Initialization *)
SetNphotTrunc[50]
InitializeOps[]

(* Parameters *)
maxAbsVal = 1.;
labelStyle = Directive[FontSize -> 18];
energyShift = 100;
NumOfStates = 1;

(* u and v define a 2d section along which we plot the 4d Wigner function *)
u = {1, 1};
v = {I, -I};


(* Functions for computing the Wigner function from a numerical state *)
wignerQutritBosonFromState[z1_, z2_, a_, b_, state_] := 
 ConjugateTranspose[state] . 
  BothBosonDisplacementOperator[2*z1, 
   2*z2] . (QutritDisplacementOperator[2 a, 2 b] . TotalParity ) . 
  state


(* Here we choose the resolution and range of the resulting plots *)
Options[ComputeWignerQutritBoson2DFromState] =
 {"Resolution" -> 
   20,
  "PhaseRange" -> {-2, 2}};

ComputeWignerQutritBoson2DFromState[state_, uSpec_, vSpec_, 
  opts : OptionsPattern[]] := 
 Module[{resolution, phaseRange, qGrid, pGrid, dataArray, omega, 
   coordinateFunc},
	
  resolution = OptionValue["Resolution"];
  phaseRange = OptionValue["PhaseRange"];
  omega = Exp[2*Pi*I/3];
  qGrid = 
   createPhaseGrid[phaseRange[[1]], phaseRange[[2]], resolution];
  pGrid = qGrid;
  coordinateFunc = 
   WignerFunctions`Private`GenerateCoordinateFunction2D[uSpec, 
    vSpec];
  
  dataArray = ParallelTable[Module[{values},
     values = 
      Table[Chop[#, 10^(-6)] &@
        wignerQutritBosonFromState[coordinateFunc[q, p], 
         Quotient[qt, 3], Mod[qt, 3], state], {p, pGrid}, {q, qGrid}];
      WignerData[
      values, <|"Method" -> "QutritBoson from State", 
       "Parameters" -> <|"a" -> Quotient[qt, 3], "b" -> Mod[qt, 3]|>, 
       "PhaseGrid" -> {qGrid, pGrid}, "PhaseRange" -> phaseRange|>]],
    {qt, 0, 8}];
  
  QutritBosonWignerData[Partition[dataArray, 3]]]


(* Parallel compute *)
LaunchKernels[];
DistributeDefinitions[wignerQutritBosonFromState, 
  BosonDisplacementOperator, QutritDisplacementOperator, TotalParity, 
  TwoModeEigenStatesAndEnergies];




(* Actual computation *)

(* 1st Z3 Rabi model parameter set *)
Module[{\[Omega] = 1,
   \[CapitalDelta] = 0.1,
   \[Phi] = 7 Pi/6,
   \[Lambda] = 1,
   spectrum, state},
  spectrum = 
   TwoModeEigenStatesAndEnergies[\[Omega], \[CapitalDelta], \[Phi], \
\[Lambda], energyShift, NumOfStates];
  state = spectrum[[1, 2]];
  data1 = ComputeWignerQutritBoson2DFromState[state, u, v];
  plot1 = PlotQutritBosonRow[data1, 0,
	                     MaxAbsVal -> maxAbsVal,
			     LabelStyle -> labelStyle];
  
  Export[GetOutputPath["Figure6a.png"], plot1];
  Print["Saved the plot of the numerically-computed Q2B Wigner function with lambda = ", \[Lambda] " to ", GetOutputPath["Figure6a.png"], "\n"];
  ];



(* 2nd Z3 Rabi model parameter set *)
Module[{\[Omega] = 1,
   \[CapitalDelta] = 0.1,
   \[Phi] = 7 Pi/6,
   \[Lambda] = 0.5,
	spectrum, state},
  spectrum = 
   TwoModeEigenStatesAndEnergies[\[Omega], \[CapitalDelta], \[Phi], \
\[Lambda], energyShift, NumOfStates];
  state = spectrum[[1, 2]];
  data2 = ComputeWignerQutritBoson2DFromState[state, u, v];
  plot2 = PlotQutritBosonRow[data2, 0,
	                     MaxAbsVal -> maxAbsVal,
			     LabelStyle -> labelStyle];
  
  Export[GetOutputPath["Figure6b.png"], plot2];
  Print["Saved the plot of the numerically-computed Q2B Wigner function with lambda = ", \[Lambda] " to ", GetOutputPath["Figure6b.png"], "\n"];
  ];



(* 3rd Z3 Rabi model parameter set *)
Module[{\[Omega] = 1,
   \[CapitalDelta] = 0.1,
   \[Phi] = 7 Pi/6,
   \[Lambda] = 0.1,
   spectrum, state},
  spectrum = TwoModeEigenStatesAndEnergies[\[Omega], \[CapitalDelta], \[Phi], \
\[Lambda], energyShift, NumOfStates];
  state = spectrum[[1, 2]];
  data3 = ComputeWignerQutritBoson2DFromState[state, u, v];
  plot3 = PlotQutritBosonRow[data3, 0,
	                     MaxAbsVal -> maxAbsVal,
			     LabelStyle -> labelStyle];

  Export[GetOutputPath["Figure6c.png"], plot3];
  Print["Saved the plot of the numerically-computed Q2B Wigner function with lambda = ", \[Lambda] " to ", GetOutputPath["Figure6c.png"], "\n"];
  ];
