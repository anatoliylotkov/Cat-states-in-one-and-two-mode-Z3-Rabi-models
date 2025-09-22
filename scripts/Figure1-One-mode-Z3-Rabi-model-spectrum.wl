(* Import package *)
AppendTo[$Path, FileNameJoin[{Directory[], "src"}]];
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

pointsNum = 80;
energyShift = 100;
NumOfStates = 3;
GroundStateEnergy[\[Lambda]_] := \[Lambda]^2


(* Parallel compute *)
LaunchKernels[];
DistributeDefinitions[OneModeEigenStatesAndEnergies];


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
    OneModeEigenStatesAndEnergies[\[Omega]0,
			    \[CapitalDelta],
			    \[Phi],
			    \[Lambda],
			    energyShift,
			    NumOfStates],
   {\[Lambda], \[Lambda]Min, \[Lambda]Max, d\[Lambda]}];



shiftEigenEnergies = Transpose @ EigenEnergiesAndVectorsList[[All, All, 1 ;; 2]];
 
For[i = 1, i <= 3, i++,
 For[j = 1, j <= Length[EigenEnergiesAndVectorsList], j++,
     shiftEigenEnergies[[i, j, 2]] = shiftEigenEnergies[[i, j, 2]] - energyShift + GroundStateEnergy[d\[Lambda]*(j-1)];
     ]];

{shiftEigenEnergies[[2, 30 ;; 81, 2]], 
   shiftEigenEnergies[[3, 30 ;; 81, 2]]} = {shiftEigenEnergies[[3, 
    30 ;; 81, 2]], shiftEigenEnergies[[2, 30 ;; 81, 2]]};

numericalPlot = ListPlot[shiftEigenEnergies[[{1, 3, 2}]], 
    PlotRange -> {Automatic, {-0.2, 0.2}}, Frame -> True, 
    FrameStyle -> BlackFrame, 
    FrameLabel -> MaTeX@{"\\lambda", "\\epsilon_k(\\lambda)"}, 
    PlotLegends -> PointLegend[Automatic, MaTeX@{"\\epsilon_0", "\\epsilon_1", "\\epsilon_2"}, 
    LegendMarkerSize -> fontSize], 
    FrameTicksStyle -> Directive[FontSize -> fontSize - 2], 
    ImageSize -> 246
   ];

analyticalPlot = Plot[
	Evaluate @
          Table[0.2 Exp[-3 \[Lambda]^2/2] Cos[(2 Pi k)/3 + (7 Pi)/6 - 
                Sqrt[3] \[Lambda]^2/2], {k, 0, 2}], {\[Lambda], 0, 2}, 
    PlotStyle -> Thickness[0.007], 
    PlotRange -> {{0, 2}, {-0.2, 0.2}}];


plot = Show[numericalPlot, analyticalPlot];

Export[GetOutputPath["Figure1.pdf"], plot];
Print["Saved the eigenenergies plot of the one-mode Z3 Rabi model to ", GetOutputPath["Figure1.png"], "\n"];
