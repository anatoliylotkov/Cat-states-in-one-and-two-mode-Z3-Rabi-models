(* ::Package:: *)

(* WignerFunctions Package *)
(* Modern implementation for computing and visualizing Wigner functions *)
(* Version 2.0 *)

BeginPackage["WignerFunctions`"]

(* Public symbols *)
WignerData::usage = "WignerData[values, metadata] represents computed Wigner function data.";
QutritBosonWignerData::usage = "QutritBosonWignerData[dataArray] represents qutrit-boson Wigner function data.";

(* A high-dimensional (analytic or numerically callable) Wigner function.  The
   first argument is a pure function that takes the appropriate number of
   real coordinates (e.g. {q1,p1,q2,p2}) and returns the Wigner value.
   metadata is an association analogue to WignerData, storing parameters, the
   system type, etc. *)
WignerFunction::usage = "WignerFunction[f, metadata] represents a callable (typically 4-D) Wigner function together with descriptive metadata.  This wrapper allows visualisation functions such as PlotWigner3D to work with the same data-centric API that PlotWigner uses for 2-D heat maps.";

(* 3-D visualisation *)
PlotWigner3D::usage = "PlotWigner3D[wigner, parameterisation, {u,u1,u2}, {v,v1,v2}, {w,w1,w2}, options] visualises a high-dimensional Wigner function on the 3-D manifold given by the parameterisation. wigner is normally a WignerFunction object (recommended) but a bare Function is also accepted for backward compatibility.";
PlotQutritBoson3D::usage = "PlotQutritBoson3D[\[Alpha], {a,b}, parameterisation, uRange, vRange, wRange, options] is a convenience wrapper around PlotWigner3D for joint qutrit\[Dash]boson cat states.";

(* Type patterns *)
PositionCoordinate::usage = "PositionCoordinate represents a phase space position coordinate (q).";
MomentumCoordinate::usage = "MomentumCoordinate represents a phase space momentum coordinate (p).";
PhasePoint::usage = "PhasePoint represents a phase space point {q, p}.";
QutritIndex::usage = "QutritIndex represents a qutrit index (0, 1, or 2).";
QutritPhasePoint::usage = "QutritPhasePoint represents a qutrit phase point {a, b}.";
ComplexAmplitude::usage = "ComplexAmplitude represents a complex amplitude parameter.";

(* Utility functions *)
SetOutputDirectory::usage = "SetOutputDirectory[dir] sets the default output directory for plots. If dir is None, uses current directory.";
GetOutputPath::usage = "GetOutputPath[filename] returns the full path for saving a file in the output directory.";
$WignerOutputDirectory::usage = "$WignerOutputDirectory is the current output directory for saving plots.";

(* Backend functions *)
ComputeWignerBosonCat::usage = "ComputeWignerBosonCat[alpha, plus, options] computes Wigner function for boson cat states. If plus is True, uses + superposition, otherwise -.";
ComputeWigner2BosonZ3Cat2D::usage = "ComputeWigner2BosonZ3Cat2D[alpha, uSpec, vSpec, options] computes a 2D slice of the Wigner function for 2-boson Z_3 cat states. Option 'k' is responsible for choosing the specific cat state.";
ComputeWignerQutritBoson2D::usage = "ComputeWignerQutritBoson2D[alpha, uSpec, vSpec, options] computes joint qutrit-boson cat state Wigner functions.";

(* Visualization functions *)
PlotWigner::usage = "PlotWigner[wignerData, options] plots a single Wigner function as a heatmap.";
PlotWignerGrid::usage = "PlotWignerGrid[wignerDataList, layout, options] plots multiple Wigner functions in a grid.";
PlotQutritBosonGrid::usage = "PlotQutritBosonGrid[qutritBosonData, options] plots all 9 components of qutrit-boson Wigner function.";
PlotQutritBosonRow::usage = "PlotQutritBosonRow[qutritBosonData, rowNum, options] plots all 3 components of qutrit-boson Wigner function with a fixed a = rowNum.";
PlotQutritBosonFixed::usage = "PlotQutritBosonFixed[qutritBosonData, qutritPoint, options] plots single component for fixed qutrit indices.";

(* Constants *)
DEFAULTGRIDSIZE = 100;

(* Initialize output directory *)
$WignerOutputDirectory = "output_images";

(* Custom patterns for type safety *)
(* Phase space coordinate patterns *)
PhaseCoordinate = _?NumericQ;
PositionCoordinate = _?NumericQ;  (* q coordinate *)
MomentumCoordinate = _?NumericQ;  (* p coordinate *)
PhasePoint = {PositionCoordinate, MomentumCoordinate};

(* Discrete qutrit index patterns *)
QutritIndex = _Integer?((0 <= # <= 2) &);
QutritPhasePoint = {QutritIndex, QutritIndex};

(* Complex amplitude pattern *)
ComplexAmplitude = _?NumericQ;

(* Helper functions *)
createPhaseGrid[minPhasePoint_, maxPhasePoint_, resolution_Integer] := Module[
  {coordMin, coordMax, grid},
  
  coordMin = minPhasePoint;
  coordMax = maxPhasePoint;
  
  grid = Range[coordMin, coordMax, (coordMax - coordMin)/(resolution - 1)];
  grid
]

validatePhasePoint[phasePoint_PhasePoint] := True
validatePhasePoint[___] := False

Begin["`Private`"]

(* Core data structures *)
(* WignerData and QutritBosonWignerData are simple containers - no additional definitions needed *)
(* They will be created using WignerData[values, metadata] in the functions below *)

(* Utility functions for output directory management *)
SetOutputDirectory[dir_] := (
  If[dir === None,
    $WignerOutputDirectory = ".",
    (* Create directory if it doesn't exist *)
    If[!DirectoryQ[dir], CreateDirectory[dir]];
    $WignerOutputDirectory = dir
  ];
  $WignerOutputDirectory
)

GetOutputPath[filename_String] := Module[{path},
  (* Ensure output directory exists *)
  If[!DirectoryQ[$WignerOutputDirectory], CreateDirectory[$WignerOutputDirectory]];
  
  (* Build full path *)
  path = FileNameJoin[{$WignerOutputDirectory, filename}];
  path
]

(* Utility functions*)
RoundToSignificant[x_, n_] := 
  N @ Floor[x, 10^(Floor[Log10[Abs[x]]] - n + 1)]



(* Backend: Boson Cat State Module *)
Options[ComputeWignerBosonCat] = {
  "Resolution" -> DEFAULTGRIDSIZE,
  (* Automatic chooses a range that comfortably contains both coherent
     peaks \[PlusMinus]\[Sqrt]2 \[Alpha] as well as the interference region around the origin. *)
  "PhaseRange" -> Automatic
};

ComputeWignerBosonCat[alpha_?NumericQ, plus_?BooleanQ, opts:OptionsPattern[]] := Module[
  {resolution, phaseRange, norm, qGrid, pGrid, values, metadata},
  
  resolution = OptionValue["Resolution"];
  phaseRange = OptionValue["PhaseRange"];

  (* Derive a reasonable phase-space window automatically when the user
     does not specify one explicitly.  We include both coherent peaks,
     situated at q, p = \[PlusMinus]\[Sqrt]2 Re/Im[\[Alpha]], and leave a margin of 2 units so that
     they are not clipped by the plot frame. *)
  If[phaseRange === Automatic,
    Module[{alphaRe, alphaIm, qRange, pRange, range},
      {alphaRe, alphaIm} = ReIm[alpha];
      qRange = Abs[Sqrt[2]*alphaRe];
      pRange = Abs[Sqrt[2]*alphaIm];
      range = Max[2, qRange, pRange] + 2; (* extra padding *)
      phaseRange = {-range, range};
    ]
  ];


  qGrid = createPhaseGrid[phaseRange[[1]], phaseRange[[2]], resolution];
  pGrid = qGrid;
  
  (* Compute Wigner function values using the provided formula *)
  values = Table[
    wignerBosonCat[{q, p}, alpha, plus],
    {p, pGrid}, {q, qGrid}
  ];
  
  metadata = <|
    "Method" -> "BosonCat",
    "Parameters" -> <|"alpha" -> alpha, "plus" -> plus|>,
    "PhaseGrid" -> {qGrid, pGrid},
    "PhaseRange" -> phaseRange,
    "Normalization" -> norm
  |>;
  
  WignerData[values, metadata]
]

(* Boson cat state Wigner function *)
wignerBosonCat[phasePoint_, alpha_?NumericQ, plus_?BooleanQ] := Module[
  {q, p, term1, term2, term3, sign, norm, alphaPrime, alphaPrimePrime},

  (* Decompose alpha into real and imaginary parts *)
  {alphaPrime, alphaPrimePrime} = ReIm[alpha];
  norm = If[plus, 
      2*(1 + Exp[-2*Abs[alpha]^2]),
      2*(1 - Exp[-2*Abs[alpha]^2])
    ];

  {q, p} = phasePoint;
  sign = If[plus, 1, -1];
  
  term1 = Exp[-(q - Sqrt[2]*alphaPrime)^2 - (p - Sqrt[2]*alphaPrimePrime)^2];
  term2 = Exp[-(q + Sqrt[2]*alphaPrime)^2 - (p + Sqrt[2]*alphaPrimePrime)^2];
  term3 = 2*Exp[-q^2 - p^2]*Cos[2*Sqrt[2]*(alphaPrime*p - q*alphaPrimePrime)];
  
  1/(Pi*norm) * (term1 + term2 + sign*term3)
]



Options[ComputeWigner2BosonZ3Cat2D] = {
  "Resolution" -> DEFAULTGRIDSIZE,
  "PhaseRange" -> {-4, 4},
  "k" -> 0  (* Additional phase parameter *)
};
ComputeWigner2BosonZ3Cat2D[alpha_?NumericQ, uSpec_, vSpec_, opts:OptionsPattern[]] := Module[
  {resolution, phaseRange, k, qGrid, pGrid, dataArray, omega,coordinateFunc},
  
  resolution = OptionValue["Resolution"];
  phaseRange = OptionValue["PhaseRange"];
  If[phaseRange === Automatic,
    Module[{alphaRe, alphaIm, qRange, pRange, range},
      {alphaRe, alphaIm} = ReIm[alpha];
      qRange = Abs[alphaRe];
      pRange = Abs[alphaIm];
      range = Max[2, qRange, pRange] + 2; (* extra padding *)
      phaseRange = {-range, range};
    ]
  ];
  
  k = OptionValue["k"];
  
  qGrid = createPhaseGrid[phaseRange[[1]], phaseRange[[2]], resolution];
  pGrid = qGrid;

  coordinateFunc = GenerateCoordinateFunction2D[uSpec, vSpec];
  
  values = Table[wigner2BosonZ3Cat[coordinateFunc[q, p], alpha, k],
        {p, pGrid}, {q, qGrid}
      ];
      
  (* Create WignerData for this component *)
  WignerData[
     values,
      <|
        "Method" -> "2-boson Z3 cat state",
        "Parameters" -> <|"alpha" -> alpha, "k" -> k|>,
        "PhaseGrid" -> {qGrid, pGrid},
        "PhaseRange" -> phaseRange
      |>
    ]
]


(* Boson cat state Wigner function *)
wigner2BosonZ3Cat[z1_?NumericQ, z2_?NumericQ, alpha_?NumericQ, parity_?IntegerQ] := Module[
  {omega, term1, term2, norm},
  omega = Exp[2 Pi I/3];

  norm = 3 + 2 * Sum[Cos[2 Pi * parity * b / 3] Exp[ - 2 (1 + Cos[2 Pi b / 3]) Abs[alpha]^2], {b, 0, 2}];

  
  term1 = 1/(Pi^2 norm) * Sum[Exp[-2 Abs[alpha - omega^b z1]^2 - 2 Abs[alpha - omega^(-b) z2]^2], {b, 0, 2}];
  term2 = 1/(Pi^2 norm) * Sum[Exp[-Abs[alpha + 2 omega^b z1]^2/2 - Abs[alpha + 2 omega^(-b) z2]^2/2] *
      Cos[2 Sqrt[3] Re[alpha (omega^b Conjugate[z1] - omega^(-b) Conjugate[z2])] + 2 Pi parity/3],
    {b, 0, 2}];
  
  term1 + term2
]



(* Backend: Qutrit-Boson Module *)

(* Utility function that transforms 2D real coordinates into 2D complex coordinates z_1, z_2. uSpec, vSpec are complex 2D vectors *)
GenerateCoordinateFunction2D[uSpec_, vSpec_]:= Function[{q,p}, Sequence @@ (q*uSpec + p*vSpec)]
GenerateCoordinateFunction3D[uSpec_, vSpec_, wSpec_]:= Function[{q,p,r}, Sequence @@ (q*uSpec + p*vSpec + r*wSpec)]

GenerateOpacityFunction[maxAbsVal_]:= Function[x, If[Abs[x] < 0.001, 0., E^(1/(0.001-Abs[x]))/(E^(1/(0.001 - Abs[x])) + E^(-1000/(maxAbsVal - Abs[x])))]]
	
(* Backend: Qutrit-Boson Module *)
Options[ComputeWignerQutritBoson2D] = {
  "Resolution" -> DEFAULTGRIDSIZE,
  "PhaseRange" -> Automatic,
  "k" -> 0  (* Additional phase parameter *)
};

ComputeWignerQutritBoson2D[alpha_, uSpec_, vSpec_, opts:OptionsPattern[]] := Module[
  {resolution, phaseRange, k, qGrid, pGrid, dataArray, omega,coordinateFunc},
  
  resolution = OptionValue["Resolution"];
  phaseRange = OptionValue["PhaseRange"];
  k = OptionValue["k"];
  omega = Exp[2*Pi*I/3];

  If[phaseRange === Automatic,
    Module[{alphaRe, alphaIm, qRange, pRange, range},
      {alphaRe, alphaIm} = ReIm[alpha];
      qRange = Abs[Sqrt[2]*alphaRe];
      pRange = Abs[Sqrt[2]*alphaIm];
      range = Max[2, qRange, pRange] + 2; (* extra padding *)
      phaseRange = {-range, range};
    ]
  ];

  
  qGrid = createPhaseGrid[phaseRange[[1]], phaseRange[[2]], resolution];
  pGrid = qGrid;

  coordinateFunc = GenerateCoordinateFunction2D[uSpec, vSpec];
  
  (* Compute Wigner function for each (a,b) pair *)
  dataArray = Table[
    Module[{values},
      values = Table[
              wignerQutritBoson[coordinateFunc[q, p], {a, b}, alpha, k, omega],
        {p, pGrid}, {q, qGrid}
      ];
      
      (* Create WignerData for this component *)
      WignerData[
        values,
        <|
          "Method" -> "QutritBoson",
          "Parameters" -> <|"alpha" -> alpha, "a" -> a, "b" -> b, "k" -> k|>,
          "PhaseGrid" -> {qGrid, pGrid},
          "PhaseRange" -> phaseRange
        |>
      ]
    ],
    {a, 0, 2}, {b, 0, 2}
  ];
  
  QutritBosonWignerData[dataArray]
]

Options[ComputeWignerQutritBoson3D] = {
  "Resolution" -> DEFAULTGRIDSIZE,
  "PhaseRange" -> {-4, 4},
  "k" -> 0  (* Additional phase parameter *)
};

ComputeWignerQutritBoson3D[alpha_, uSpec_, vSpec_, wSpec_, opts:OptionsPattern[]] := Module[
  {resolution, phaseRange, k, qGrid, pGrid, rGrid, dataArray, omega, coordinateFunc},
  
  resolution = OptionValue["Resolution"];
  phaseRange = OptionValue["PhaseRange"];
  k = OptionValue["k"];
  omega = Exp[2*Pi*I/3];
  
  qGrid = createPhaseGrid[phaseRange[[1]], phaseRange[[2]], resolution];
  pGrid = qGrid;
  rGrid = pGrid;

  coordinateFunc = GenerateCoordinateFunction3D[uSpec, vSpec, wSpec];

  (* Compute Wigner function for each (a,b) pair *)
  dataArray = Table[
    Module[{values},
      values = Table[
              wignerQutritBoson[coordinateFunc[q, p, r], {a, b}, alpha, k, omega],
              {p, pGrid}, {q, qGrid}, {r, rGrid}
      ];
      
      (* Create WignerData for this component *)
      WignerData[
        values,
        <|
          "Method" -> "QutritBoson",
          "Parameters" -> <|"alpha" -> alpha, "a" -> a, "b" -> b, "k" -> k|>,
          "PhaseGrid" -> {qGrid, pGrid},
          "PhaseRange" -> phaseRange
        |>
      ]
    ],
    {a, 0, 2}, {b, 0, 2}
  ];
  
  QutritBosonWignerData[dataArray]
]

(* Joint qutrit-boson Wigner function *)
wignerQutritBoson[z1_?NumericQ, z2_?NumericQ, qutritPoint_, alpha_, k_, omega_] := Module[
  {a, b, term1, term2, arg},
  {a, b} = qutritPoint;
  
  term1 = Exp[-2*Abs[alpha - omega^(-b)*z1]^2 - 2*Abs[alpha - omega^b*z2]^2];
  
  arg = 2*Sqrt[3]*Re[alpha*(Conjugate[omega^b*z1] - Conjugate[omega^(-b)*z2])] + 2*Pi*(a + k)/3;
  term2 = Exp[-Abs[alpha + 2*omega^(-b)*z1]^2/2 - Abs[alpha + 2*omega^b*z2]^2/2] * Cos[arg];
  
  1/(9*Pi^2) * (term1 + term2)
]




(* Visualization Module *)

WignerColor[ColorFunction_, maxAbsVal_] := Which[
    ColorFunction == "BlueWhiteRed",
    (* Balanced blue-white-red centered at zero *)
    Function[x, 
      Which[
        x < -maxAbsVal, Blue,
        x > maxAbsVal, Red,
        True, Blend[{Blue, White, Red}, (x + maxAbsVal)/(2*maxAbsVal)]
      ]
    ],
    ColorFunction == "BlueWhiteOrange",
    (* Alternative balanced scheme *)
    Function[x,
      Which[
        x < -maxAbsVal, RGBColor[0, 0.4, 0.8],
        x > maxAbsVal, RGBColor[1, 0.5, 0],
        True, Blend[{RGBColor[0, 0.4, 0.8], White, RGBColor[1, 0.5, 0]}, (x + maxAbsVal)/(2*maxAbsVal)]
      ]
    ],
    ColorFunction == "PurpleWhiteGreen",
    (* Another balanced scheme *)
    Function[x,
      Which[
        x < -maxAbsVal, Purple,
        x > maxAbsVal, Darker[Green],
        True, Blend[{Purple, White, Darker[Green]}, (x + maxAbsVal)/(2*maxAbsVal)]
      ]
    ],
    True,
    ColorData[ColorFunction]
  ]

(* TriangleShape creates a dashed circles forming a triangle to improve readability of the Z3 cat state Wigner function plots *)
TriangleShape[alpha_, r_:0.75] :=
  Graphics[{Gray, Dashed, 
    Line[{{alpha Sqrt[2], 0}, {-alpha Sqrt[2]/2, 
     alpha Sqrt[6]/2}, {-alpha Sqrt[2]/2, -alpha Sqrt[6]/
       2}, {alpha Sqrt[2], 0}}], FaceForm[Transparent], 
    EdgeForm[Directive[Dashed, Gray]], Disk[{3 Sqrt[2], 0}, r], 
    Disk[{-alpha Sqrt[2]/2, alpha Sqrt[6]/2}, r], 
    Disk[{-alpha Sqrt[2]/2, -alpha Sqrt[6]/2}, r],
    Disk[{-alpha Sqrt[2]/2, 0}, r], 
    Disk[{alpha Sqrt[2]/4, -alpha Sqrt[6]/4}, r], 
    Disk[{alpha Sqrt[2]/4, alpha Sqrt[6]/4}, r]}]

(* MedianShape creates two dashed circles placed on a median of the triangle *)
MedianShape[alpha_, num_, r_:0.75] :=
  Part[{
    Graphics[{FaceForm[Transparent], EdgeForm[Directive[Dashed, Gray]], Disk[{3 Sqrt[2], 0}, r], Disk[{-alpha Sqrt[2]/2, 0}, r]}],
    Graphics[{FaceForm[Transparent], EdgeForm[Directive[Dashed, Gray]], Disk[{-alpha Sqrt[2]/2, alpha Sqrt[6]/2}, r], Disk[{alpha Sqrt[2]/4, -alpha Sqrt[6]/4}, r]}],
    Graphics[{FaceForm[Transparent], EdgeForm[Directive[Dashed, Gray]], Disk[{-alpha Sqrt[2]/2, -alpha Sqrt[6]/2}, r], Disk[{alpha Sqrt[2]/4, alpha Sqrt[6]/4}, r]}]},
  num]




(* PlotWigner *)
Options[PlotWigner] = {
  ColorFunction -> "BlueWhiteRed",
  PlotRange -> All,
  PlotLabel -> None,
  FrameLabel -> {"x", "p"},
  ImageSize -> 400,
  PlotLegends -> Automatic,
  SaveAs -> None,
  MaxAbsVal -> None,
  LabelStyle -> Automatic,
  TriangleGuide -> None
};

PlotWigner[WignerData[values_, metadata_], opts:OptionsPattern[]] := Module[
  {colorFunc, phaseGrid, aspectRatio, minVal, maxVal, maxAbsVal, plot, savePath, alpha},
  
  alpha = metadata["Parameters", "alpha"];
  qutritB = metadata["Parameters", "b"];
  phaseGrid = metadata["PhaseGrid"];
  aspectRatio = If[Length[phaseGrid] == 2,
    (phaseGrid[[2, -1]] - phaseGrid[[2, 1]])/(phaseGrid[[1, -1]] - phaseGrid[[1, 1]]),
    1
  ];
  
  (* Calculate range for balanced color scaling *)
  minVal = Min[values];
  maxVal = Max[values];
  maxAbsVal = If[OptionValue[MaxAbsVal] === None, Max[Abs[minVal], Abs[maxVal]], OptionValue[MaxAbsVal]];
  
  colorFunc = WignerColor[OptionValue[ColorFunction], maxAbsVal];
  plotRange = If[OptionValue[PlotRange] === All,
      {-maxAbsVal, maxAbsVal},
      OptionValue[PlotRange]
    ];
  plotLegends = If[OptionValue[PlotLegends] === Automatic,
		   BarLegend[
			   {colorFunc, {-maxAbsVal, maxAbsVal}},
			   Ticks -> ({#, ToString[#]} & /@ Subdivide[-RoundToSignificant[maxAbsVal, 3], RoundToSignificant[maxAbsVal, 3], 4])],
		   OptionValue[PlotLegends]];


  preliminaryPlot = ListDensityPlot[values,
    ColorFunction -> colorFunc,
    (* Prevent Mathematica from rescaling the data to {0,1} before it is
       passed to our custom ColorFunction.  The colour schemes defined
       above assume the raw Wigner values, where zero is mapped to white
       and \[PlusMinus]maxAbsVal to the fully saturated colours.  With the default
       behaviour (ColorFunctionScaling -> True) everything becomes
       positive and the plot appears uniformly red. *)
    ColorFunctionScaling -> False,
    InterpolationOrder -> 1,
    PlotRange -> plotRange,
    PlotLabel -> OptionValue[PlotLabel],
    FrameLabel -> OptionValue[FrameLabel],
    ImageSize -> OptionValue[ImageSize],
    PlotLegends -> plotLegends,
    AspectRatio -> aspectRatio,
    DataRange -> If[Length[phaseGrid] == 2,
      {{phaseGrid[[1, 1]], phaseGrid[[1, -1]]}, {phaseGrid[[2, 1]], phaseGrid[[2, -1]]}},
      Automatic
    ],
    Frame -> True,
    FrameTicks -> True,
    LabelStyle -> OptionValue[LabelStyle]
  ];

  plot = Switch[OptionValue[TriangleGuide],
		None, preliminaryPlot,
		"Triangle", Show[preliminaryPlot, TriangleShape[alpha]],
		"Median", Show[preliminaryPlot, MedianShape[alpha, qutritB + 1]]];
  
  (* Handle SaveAs option *)
  If[OptionValue[SaveAs] =!= None,
    savePath = GetOutputPath[OptionValue[SaveAs]];
    Export[savePath, plot];
    Print["Plot saved to: ", savePath];
  ];
  
  plot
]



Options[PlotWigner3D] = Join[
  Options[PlotWigner],
  {
    Resolution -> 40,
    EmbeddingFunction -> Automatic,
    AxesLabel -> Automatic,
    BoxRatios -> {1, 1, 1}
  }
];

(* ----------------------------------------------------------------- *)
(*  Internal helper \[Dash] core implementation expecting a pure function   *)
(* ----------------------------------------------------------------- *)
plotWigner3DFromFunction[func_Function, param_,
  uSpec : {u_Symbol, uMin_?NumericQ, uMax_?NumericQ},
  vSpec : {v_Symbol, vMin_?NumericQ, vMax_?NumericQ},
  wSpec : {w_Symbol, wMin_?NumericQ, wMax_?NumericQ},
  opts : OptionsPattern[PlotWigner3D]]:= Module[
  {res, uGrid, vGrid, wGrid, embedFn, pts},
  res = OptionValue[Resolution];
  embedFn = OptionValue[EmbeddingFunction];
  
  If[embedFn === Automatic,
    embedFn = Function[{q1, p1, q2, p2}, {q1, p1, q2}]
  ];

  (* Parameter grids *)
  uGrid = Range[uMin, uMax, (uMax - uMin)/(res - 1)];
  vGrid = Range[vMin, vMax, (vMax - vMin)/(res - 1)];
  wGrid = Range[wMin, wMax, (wMax - wMin)/(res - 1)];

  (* Sample the Wigner function *)
  pts = Reap[
    Do[
      With[{phase = param[uVal, vVal, wVal]},
        With[{x = embedFn @@ phase, value = func @@ phase},
          Sow[Append[x, value]]
        ]
      ],
      {uVal, uGrid}, {vVal, vGrid}, {wVal, wGrid}
    ]
  ][[2, 1]];

  badPlotWigner3D[pts, opts]					 
]

badPlotWigner3D[pts_,
  opts : OptionsPattern[PlotWigner3D]] := Module[
  {vals, minVal, maxVal, maxAbsVal,
   colorFunc, plotRangeLocal, plot, savePath, axesLbl, boxRatios},

  axesLbl = OptionValue[AxesLabel];
  boxRatios = OptionValue[BoxRatios];


  vals = pts[[All, -1]];
  minVal = Min[vals];
  maxVal = Max[vals];
  maxAbsVal = Max[Abs[minVal], Abs[maxVal]];

  colorFunc = WignerColor[OptionValue[ColorFunction], maxAbsVal];
  plotRangeLocal = If[ MemberQ[{"BlueWhiteRed", "BlueWhiteOrange", "PurpleWhiteGreen"}, OptionValue[ColorFunction]] && OptionValue[PlotRange] === All,
    {-maxAbsVal, maxAbsVal},
    OptionValue[PlotRange]
  ];
  plotLegends = If[OptionValue[PlotLegends] === Automatic,
		   BarLegend[
			   {colorFunc, {-maxAbsVal, maxAbsVal}},
			   Ticks -> ({#, ToString[#]} & /@ Subdivide[-RoundToSignificant[maxAbsVal, 3], RoundToSignificant[maxAbsVal, 3], 4])],
		   OptionValue[PlotLegends]];

  plot = ListDensityPlot3D[
    pts,
    PlotRange -> {All, All, All, plotRangeLocal},
    ColorFunction -> colorFunc,
    ColorFunctionScaling -> False,
    PlotLabel -> OptionValue[PlotLabel],
    AxesLabel -> axesLbl,
    BoxRatios -> boxRatios,
    PlotLegends -> plotLegends,
    ImageSize -> OptionValue[ImageSize]
  ];

  If[OptionValue[SaveAs] =!= None,
    savePath = GetOutputPath[OptionValue[SaveAs]];
    Export[savePath, plot];
    Print["3D plot saved to: ", savePath];
  ];

  plot
];


Options[PlotWigner3D] = Join[
  Options[PlotWigner],
  {
    AxesLabel -> Automatic,
    BoxRatios -> {1, 1, 1},
    OpacityFunction -> Automatic
  }
];
(* Rewritten plotWigner3D to take WignerData as an argument*)
plotWigner3D[WignerData[values_, metadata_],
	     opts : OptionsPattern[PlotWigner3D]] := Module[
  {vals, minVal, maxVal, maxAbsVal,
   colorFunc, plotRangeLocal, plot, savePath, axesLbl, boxRatios},

  axesLbl = OptionValue[AxesLabel];
  boxRatios = OptionValue[BoxRatios];


  minVal = Min[values];
  maxVal = Max[values];
  maxAbsVal = Max[Abs[minVal], Abs[maxVal]];

  colorFunc = WignerColor[OptionValue[ColorFunction], maxAbsVal];
  plotRangeLocal = If[ MemberQ[{"BlueWhiteRed", "BlueWhiteOrange", "PurpleWhiteGreen"}, OptionValue[ColorFunction]] && OptionValue[PlotRange] === All,
    {-maxAbsVal, maxAbsVal},
    OptionValue[PlotRange]
  ];
  plotLegends = If[OptionValue[PlotLegends] === Automatic,
		   BarLegend[
			   {colorFunc, {-maxAbsVal, maxAbsVal}},
			   Ticks -> ({#, ToString[#]} & /@ Subdivide[-RoundToSignificant[maxAbsVal, 3], RoundToSignificant[maxAbsVal, 3], 4])],
		   OptionValue[PlotLegends]];
  opacityFunc = GenerateOpacityFunction[maxAbsVal];

  plot = ListDensityPlot3D[
    values,
    PlotRange -> {All, All, All, plotRangeLocal},
    ColorFunction -> colorFunc,
    ColorFunctionScaling -> False,
    PlotLabel -> OptionValue[PlotLabel],
    AxesLabel -> axesLbl,
    BoxRatios -> boxRatios,
    PlotLegends -> plotLegends,
    ImageSize -> OptionValue[ImageSize],
    OpacityFunction -> opacityFunc,
    OpacityFunctionScaling -> False
  ];

  If[OptionValue[SaveAs] =!= None,
    savePath = GetOutputPath[OptionValue[SaveAs]];
    Export[savePath, plot];
    Print["3D plot saved to: ", savePath];
  ];

  plot
];


(* Grid plotting *)
Options[PlotWignerGrid] = Join[
  Options[PlotWigner],
  {GridLayout -> Automatic, Spacings -> 3}
];

PlotWignerGrid[wignerDataList_List, opts:OptionsPattern[]] := Module[
  {layout, n, plots, grid, savePath},
  
  n = Length[wignerDataList];
  layout = If[OptionValue[GridLayout] === Automatic,
    {Ceiling[Sqrt[n]], Ceiling[n/Ceiling[Sqrt[n]]]},
    OptionValue[GridLayout]
  ];
  
  plots = Map[
    PlotWigner[#, 
	       Sequence @@ FilterRules[{opts}, Except[{SaveAs, GridLayout}, Options[PlotWigner]]]
    ] &,
    wignerDataList
  ];
  
  grid = Grid[
    Partition[plots, layout[[2]], layout[[2]], {1, 1}, {}],
    Spacings -> OptionValue[Spacings],
    ItemSize -> Full
  ];
  
  (* Handle SaveAs option *)
  If[OptionValue[SaveAs] =!= None,
    savePath = GetOutputPath[OptionValue[SaveAs]];
    Export[savePath, grid];
    Print["Grid plot saved to: ", savePath];
  ];
  
  grid
]

(* Qutrit-Boson visualization *)
Options[PlotQutritBosonGrid] = Join[
  Options[PlotWigner],
  {Spacings -> 1}
];

PlotQutritBosonGrid[QutritBosonWignerData[dataArray_], opts:OptionsPattern[]] := Module[
  {plots, grid, savePath, r, alpha, triangle},
  r = 0.5;
  alpha = 3.0;
  (*   triangle = Graphics[{Gray, Dashed, 
  Line[{{alpha Sqrt[2], 0}, {-alpha Sqrt[2]/2, 
     alpha Sqrt[6]/2}, {-alpha Sqrt[2]/2, -alpha Sqrt[6]/
       2}, {alpha Sqrt[2], 0}}], FaceForm[Transparent], 
  EdgeForm[Directive[Dashed, Gray]], Disk[{3 Sqrt[2], 0}, r], 
  Disk[{-alpha Sqrt[2]/2, alpha Sqrt[6]/2}, r], 
  Disk[{-alpha Sqrt[2]/2, -alpha Sqrt[6]/2}, r],
  Disk[{-alpha Sqrt[2]/2, 0}, r], 
  Disk[{alpha Sqrt[2]/4, -alpha Sqrt[6]/4}, r], 
  Disk[{alpha Sqrt[2]/4, alpha Sqrt[6]/4}, r]}]; *)
  
  plots = Table[
    PlotWigner[dataArray[[a + 1, b + 1]],
      PlotLabel -> StringForm["(a,b) = (``,``)", a, b],
      ImageSize -> 300,
      PlotLegends -> If[b == 2, Automatic, None],
      Sequence @@ FilterRules[{opts}, Except[{SaveAs}, Options[PlotWigner]]]
    ],
    {a, 0, 2}, {b, 0, 2}
  ];

  (*plots = If[OptionValue[TriangleGuide],
     Map[Show[#, triangle]&, plots, {2}],
	     plots];*)
  
  grid = Grid[
	  plots,
	  Spacings -> OptionValue[Spacings],
	  ItemSize -> Full
	 ];
  
  (* Handle SaveAs option *)
  If[OptionValue[SaveAs] =!= None,
    savePath = GetOutputPath[OptionValue[SaveAs]];
    Export[savePath, grid];
    Print["Qutrit-Boson grid plot saved to: ", savePath];
  ];
  
  grid
]

(* Qutrit-Boson visualization *)
Options[PlotQutritBosonRow] = Join[
  Options[PlotWigner],
  {Spacings -> 1}
];
PlotQutritBosonRow[QutritBosonWignerData[dataArray_], rowNum:QutritIndex, opts:OptionsPattern[]] := Module[
  {plots, row, savePath, r, alpha, triangle},

  
  plots = Table[
    PlotWigner[dataArray[[rowNum + 1, b + 1]],
      PlotLabel -> StringForm["(a,b) = (``,``)", rowNum, b],
      PlotLegends -> If[b == 2, Automatic, None],
      ImageSize -> 300,
	       Sequence @@ FilterRules[{opts}, Except[{SaveAs, RoundingRadius, Frame}, Options[PlotWigner]]]
    ],
    {b, 0, 2}
  ];

  row = GraphicsRow[Rasterize /@ plots];
  
  (* Handle SaveAs option *)
  If[OptionValue[SaveAs] =!= None,
    savePath = GetOutputPath[OptionValue[SaveAs]];
    Export[savePath, row];
    Print["Qutrit-Boson row plot saved to: ", savePath];
  ];
  
  row
]

PlotQutritBosonFixed[QutritBosonWignerData[dataArray_], qutritPoint_, opts:OptionsPattern[]] := Module[
  {a, b},
  {a, b} = qutritPoint;
  PlotWigner[dataArray[[a + 1, b + 1]], opts]
]

(* ================================================================= *)
(*                      3-D Heat-Map Visualisation                    *)
(* ================================================================= *)

(* When we have two bosonic modes the full Wigner function lives in a
   4-dimensional phase space (q\:2081, p\:2081, q\:2082, p\:2082).  A faithful visualisation
   therefore requires us to restrict the function to a 3-manifold
   embedded in that 4-space.  The user specifies this manifold by
   providing a **parameterisation** \[Phi] which maps (u,v,w) \[RightTeeArrow] {q\:2081,p\:2081,q\:2082,p\:2082}.

   PlotWigner3D samples the Wigner function on a uniform (u,v,w) grid,
   evaluates it at the corresponding phase-space points, and produces a
   ListDensityPlot3D heat-map.  An optional *embedding function* \[Chi] can
   be supplied to choose how the 4-dimensional coordinates are
   projected into 3-D plotting space (default: drop the 4th coordinate
   and take {q\:2081,p\:2081,q\:2082}).                                        *)


(* Public wrapper \[Dash] recommended usage with WignerFunction wrapper *)

PlotWigner3D[func_Function, param_, uSpec_, vSpec_, wSpec_, opts : OptionsPattern[]] :=
  plotWigner3DFromFunction[ func, param, uSpec, vSpec, wSpec, opts ];



(* ----------------------------------------------------------------- *)
(* Convenience wrapper for the qutrit-boson case                     *)
(* ----------------------------------------------------------------- *)


Options[PlotQutritBoson3D] = Options[PlotWigner3D];

PlotQutritBoson3D[alpha_?NumericQ, qutritPoint_,
  param_, uSpec_, vSpec_, wSpec_, opts : OptionsPattern[]] := Module[
  {omega, func, wfMetadata, a, b},

  omega = Exp[2 Pi I/3];
  {a, b} = qutritPoint;

  (* Build a numeric Wigner function of four real arguments *)
  func = Function[{q1, p1, q2, p2},
    Module[{z1, z2},
      z1 = (q1 + I p1)/Sqrt[2];
      z2 = (q2 + I p2)/Sqrt[2];
      wignerQutritBoson[z1, z2, qutritPoint, alpha, 0, omega]
    ]
  ];


  PlotWigner3D[func, param, uSpec, vSpec, wSpec, FilterRules[{opts}, Options[PlotWigner3D]]]
]

  PlotQutritBosonFixed3D[QutritBosonWignerData[dataArray_], qutritPoint_, opts:OptionsPattern[]] := Module[{a, b},
  {a, b} = qutritPoint;
  plotWigner3D[dataArray[[a + 1, b + 1]], opts]
]

End[] (* Private *)

EndPackage[]
