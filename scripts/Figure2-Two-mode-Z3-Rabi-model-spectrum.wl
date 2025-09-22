(* Import package *)
AppendTo[$Path, FileNameJoin[{Directory[], "src"}]]
<< WignerFunctions`
<< Z3Rabi`
<< MaTeX`

(* Innitialization *)
SetNphotTrunc[50]
InitializeOps[]

(* Parameters *)
\[Omega]0 = 1;
\[CapitalDelta] = 0.1;
\[Phi] = 7 Pi/6;
\[Lambda]Min = 0;
\[Lambda]Max = 2;
d\[Lambda] = (\[Lambda]Max - \[Lambda]Min)/pointsNum;
GroundStateEnergy[\[Lambda]_] = -2 \[Lambda]^2

pointsNum = 80;
energyShift = 100;
NumOfStates = 3;

(* Parallel compute *)
LaunchKernels[];
DistributeDefinitions[TwoModeEigenStatesAndEnergies];


(* Stylistic options *)
fontSize = 12;
SetOptions[MaTeX, FontSize -> fontSize, "DisplayStyle" -> False];
SetOptions[Plot, 
  BaseStyle -> {FontSize -> fontSize, 
    FontFamily -> "Latin Modern Math"}];


(* Calculation *)
EigenEnergiesAndVectorsList = 
  ParallelTable[
   Prepend[#, \[Lambda]] & /@ 
    TwoModeEigenStatesAndEnergies[\[Omega]0, \[CapitalDelta], \[Phi], \[Lambda], energyShift, NumOfStates],
   {\[Lambda], \[Lambda]Min, \[Lambda]Max, d\[Lambda]}];

shiftEigenEnergies = Transpose@EigenEnergiesAndVectorsList[[All, All, 1 ;; 2]];
For[i = 1, i <= 3, i++,
 For[j = 1, j <= Length[EigenEnergiesAndVectorsList], j++,
  shiftEigenEnergies[[i, j, 2]] = 
  shiftEigenEnergies[[i, j, 2]] - 100 - GroundStateEnergy[d\[Lambda]*(j-1)];
 ]];

numericalPlot = ListPlot[shiftEigenEnergies[[{1, 3, 2}]], 
     PlotRange -> {Automatic, {-0.2, 0.2}}, Frame -> True, 
     FrameStyle -> BlackFrame, 
     FrameLabel -> MaTeX@{"\\lambda", "\\epsilon_k(\\lambda)"}, 
     PlotLegends -> PointLegend[Automatic, 
       MaTeX @ {"\\epsilon_0", "\\epsilon_1", "\\epsilon_2"}, 
       LegendMarkerSize -> fontSize], 
     FrameTicksStyle -> Directive[FontSize -> fontSize - 2], 
     ImageSize -> 246];

analyticalPlot = Plot[Evaluate @
      Table[0.2 Exp[-3 (\[Lambda])^2] Cos[7 Pi/6 + 2 Pi k/3 ], {k, 0, 
        2}],
		{\[Lambda], 0, 2}, PlotRange -> {{0, 2}, {-0.2, 0.2}}];

plot = Show[numericalPlot, analyticalPlot];

Export[GetOutputPath["Figure2.pdf"], plot];
Print["Saved the eigenenergies plot of the two-mode Z3 Rabi model to ", GetOutputPath["Figure2.pdf"], "\n"];
