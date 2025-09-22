(* Import package *)
AppendTo[$Path, FileNameJoin[{Directory[], "src"}]]
<< WignerFunctions`

(* Parameters *)
alpha = 3.0;
uSpec = {1, 1}/Sqrt[2];
vSpec = {I, -I}/Sqrt[2];

maxAbsVal = 0.012;
labelStyle = Directive[FontSize -> 24];

(* Q2B cat state *)

Print["Computing Q2B cat state with \[Alpha] = ", alpha];

Q2BCatData = ComputeWignerQutritBoson2D[alpha, uSpec, vSpec, "Resolution" -> 50, "k" -> 0];
Q2BCatRow = PlotQutritBosonRow[Q2BCatData, 0,
			       TriangleGuide -> "Triangle",
			       MaxAbsVal -> maxAbsVal,
			       LabelStyle -> labelStyle];
 
Export[GetOutputPath["Figure5a.png"], Q2BCatRow];
Print["Saved a=0 row of the Q2B cat state Wigner function to ", GetOutputPath["Figure5a.png"]];




(* Mixed state *)

Print["Computing mixed state with \[Alpha] = ", alpha];

wignerMixedState[z1_?NumericQ, z2_?NumericQ, qutritPoint_, alpha_, k_, omega_] := Module[
  {a, b, term},
  {a, b} = qutritPoint;
  
  term = Exp[-2*Abs[alpha - omega^(-b)*z1]^2 - 2*Abs[alpha - omega^b*z2]^2];
  
  1/(9*Pi^2) * term
]

MixedStateData = Block[{WignerFunctions`Private`wignerQutritBoson = wignerMixedState},
				ComputeWignerQutritBoson2D[alpha, uSpec, vSpec, "Resolution" -> 50, "k" -> 0]];
MixedStateRow = PlotQutritBosonRow[MixedStateData, 0,
				   TriangleGuide -> "Triangle",
				   MaxAbsVal -> maxAbsVal,
				   LabelStyle -> labelStyle];

Export[GetOutputPath["Figure5b.png"], MixedStateRow];
Print["Saved a=0 row of the mixed state Wigner function to ", GetOutputPath["Figure5b.png"], "\n"];





(* 2B cat state*)

Print["Computing 2B cat state with \[Alpha] = ", alpha];

TwoBCatSinglePlot = ComputeWigner2BosonZ3Cat2D[alpha, uSpec, vSpec, "Resolution" -> 50, "k" -> 0, "PhaseRange" -> {-6,6}];

(* Fill QutritBosonWignerData with zeros *)
metaData = TwoBCatSinglePlot[[2]];
TwoBvalues = TwoBCatSinglePlot[[1]]/3;
values = Table[0, {x, 50}, {y, 50}];
dataArray = Table[Module[{metaData = metaData},
					metaData["Parameters"]["a"] = a;
					metaData["Parameters"]["b"] = b;
					If[a == 0 && b == 0,
					WignerData[TwoBvalues, metaData],
					WignerData[values, metaData]]],
				{a, 0, 2}, {b, 0, 2}];
TwoBCatData = QutritBosonWignerData[dataArray];

TwoBCatRow = PlotQutritBosonRow[TwoBCatData, 0,
				TriangleGuide -> "Triangle",
				MaxAbsVal -> maxAbsVal,
				LabelStyle -> labelStyle];

Export[GetOutputPath["Figure5c.png"], TwoBCatRow];
Print["Saved a=0 row of the 2B cat state Wigner function to ", GetOutputPath["Figure5c.png"], "\n"];

