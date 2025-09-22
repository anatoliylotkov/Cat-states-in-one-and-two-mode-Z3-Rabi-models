(* Import package *)
AppendTo[$Path, FileNameJoin[{Directory[], "src"}]]
<< WignerFunctions`
<< MaTeX`

(* Parameters *)
omega = Exp[2 Pi I/3];
alpha = 3;
catStateNum = 0;

(* u and v define a 2d section along which we plot of the 4d Wigner function *)
uLHSplot = {1, 1}/Sqrt[2];
vLHSplot = {I, -I}/Sqrt[2];
uListRHSplot = Table[{Sin[Pi/4 - 2 b Pi/3]omega^(b), Sin[Pi/4 + 2 b Pi/3]omega^(-b)}, {b,0,2}];
vListRHSplot = Table[{Sin[Pi/4 + 2 b Pi/3] omega^(b), -Sin[Pi/4 - 2 b Pi/3] omega^(-b)}, {b,0,2}];

maxAbsVal = 0.012;
resolution = 20;
labelStyle = Directive[FontSize -> 24];

BlueWhiteRed = Function[x, 
      Which[
        x < -maxAbsVal, Blue,
        x > maxAbsVal, Red,
        True, Blend[{Blue, White, Red}, (x + maxAbsVal)/(2*maxAbsVal)]
      ]
    ];


(* Wigner function evaluation *)

qbDataLHSplots = ComputeWignerQutritBoson2D[alpha, uLHSplot, vLHSplot, "Resolution" -> resolution, "k" -> catStateNum];
qbDataRHSplots = Table[
		ComputeWignerQutritBoson2D[alpha, uListRHSplot[[b]], vListRHSplot[[b]], "Resolution" -> resolution, "k" -> catStateNum],
		{b, 1, 3}];


(* Left Hand Side plots *)
LHSData = (qbDataLHSplots /. (QutritBosonWignerData[data_] :> data));
LHSPlots = Table[
    PlotWigner[LHSData[[a + 1, b + 1]],
      ImageSize -> 300,
      TriangleGuide -> "Triangle",
      PlotLegends -> None,
      FrameLabel -> {"Re[w]", "Im[w]"},
      ColorFunction -> "BlueWhiteRed",
      MaxAbsVal -> maxAbsVal,
      LabelStyle -> labelStyle
    ],
    {a, 0, 2}, {b, 0, 2}];


(* Right Hand Side plots *)
RHSData = Replace[qbDataRHSplots, QutritBosonWignerData[data_] :> data, {1}];
RHSPlots = Table[Module[
  {plotData = RHSData[[i]], plots},
  plots = Table[
    PlotWigner[plotData[[a + 1, b + 1]],
      ImageSize -> 300,
      TriangleGuide -> "Median",
      PlotLegends -> None,
      FrameLabel -> {"Re[w]", "Im[w]"},
      ColorFunction -> "BlueWhiteRed",
      MaxAbsVal -> maxAbsVal,
      LabelStyle -> labelStyle
    ],
    {a, 0, 2}, {b, 0, 2}
  ];

  plots
], {i, 1, 3}];


(* We combine the LHS and RHS plots into the single list plotPairs *)
allPlots = Transpose[Prepend[RHSPlots, LHSPlots], {3, 1, 2}];
plotPairs = Table[allPlots[[i,j,{1, 1 + j}]],{i,1,3},{j,1,3}];

(* The plot legend *)
legend = BarLegend[{BlueWhiteRed[10^-2 #]&,{-1.2, 1.2}}, LegendLabel -> Placed[Style["\!\(\*SuperscriptBox[\(x10\), \(-2\)]\)",GrayLevel[0.4]],Top], LegendMarkerSize->{10,190}, LabelStyle -> labelStyle]


(* For each qutrit variables value (a,b), we combine the LHS and RHS plots into a single object *)
rowedPairs = Table[
	GraphicsRow[plotPairs[[i,j]],
		    Spacings -> Scaled[0.03],
		    BaselinePosition -> Center,
		    ImageSize -> 750,
		    ImagePadding -> {{0, 100}, {0, 30}},
		    Epilog -> {Inset[legend, Scaled[{1.02, .59}], {Left, Center}],
			       Inset[Style[StringForm["(a,b) = (``,``)",i-1,j-1], GrayLevel[0.4], 24], Scaled[{0.55, 1.05}]]}
	            ],
	{i,1,3}, {j,1,3}];

(* We add a frame around each pair of plots *)
framedPairs = Map[
  Framed[#,
      RoundingRadius -> 12,
      FrameMargins  -> 10,
      FrameStyle    -> Directive[GrayLevel[0.4], Thickness[2]],
      Background    -> White
  ] &, rowedPairs, {2}];

(* We combine all the plots into a 3x3 grid *)
gridPairs = Grid[
  framedPairs,
  Spacings   -> {4,4},
  Alignment  -> Center
];

Export[GetOutputPath["Figure4.png"], gridPairs];
Print["Saved the Q2B Wignef function plot to ", GetOutputPath["Figure4.png"], "\n"];
