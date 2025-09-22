(* Import package *)
AppendTo[$Path, FileNameJoin[{Directory[], "src"}]]
<< WignerFunctions`

(* Parameters setup *)
alpha = 3.0;
labelStyle = Directive[FontSize -> 24];

(* Z2 cat state *)

Print["Computing Z2 cat state with \[Alpha] = ", alpha];

catPlus = ComputeWignerBosonCat[alpha, True, "Resolution" -> 50];

(* Plot both states side by side *)
plotPlus = PlotWigner[catPlus, 
		      ColorFunction -> "BlueWhiteRed",
		      PlotLabel -> " ",
		      LabelStyle -> labelStyle
];

Export[GetOutputPath["Figure3a.png"], plotPlus];
Print["Saved comparison plot to ", GetOutputPath["Figure3a.png"], "\n"];


(* Z3 cat state *)

Print["Computing Z3 cat state with \[Alpha] = ", alpha];

uSpec = {1,1}/Sqrt[2];
vSpec = {I, -I}/Sqrt[2];
catParity0 = ComputeWigner2BosonZ3Cat2D[alpha, uSpec, vSpec, "Resolution" -> 50, "k" -> 0, "PhaseRange" -> {-6,6}];

(* Plot both states side by side *)
plot = PlotWigner[catParity0,
		   TriangleGuide -> "Triangle",
		   ColorFunction -> "BlueWhiteRed",
		   PlotLabel -> " ",
   		   LabelStyle -> labelStyle
];


Export[GetOutputPath["Figure3b.png"], plot];
Print["Saved comparison plot to ", GetOutputPath["Figure3b.png"], "\n"];

