(* ::Package:: *)

BeginPackage["Z3Rabi`"]

(* Public symbols - what users can access *)
SetNphotTrunc::usage = "SetNphotTrunc[n] sets the photonic truncation number to n.";
InitializeOps::usage = "InitializeOps[] computes all operators and matrices with the current NphotTrunc value.";
GetNphotTrunc::usage = "GetNphotTrunc[] returns the current photonic truncation number.";
IsInitialized::usage = "IsInitialized[] returns True if the package has been initialized.";

(* Operator symbols *)
a::usage = "Photonic annihilation operator in truncated Hilbert space.";
a†::usage = "Photonic creation operator in truncated Hilbert space.";
idphot::usage = "Identity operator in photonic Hilbert space.";
idel::usage = "Identity operator in electronic/spin Hilbert space.";
Id::usage = "Identity operator in full Hilbert space.";
A::usage = "A[i] photonic annihilation operator in full Hilbert space.";
A†::usage = "A†[i] photonic creation operator in full Hilbert space.";
PhotNum::usage = "PhotNum[i] photon number operator.";
PhotNumTotal::usage = "Total photon number operator.";
J::usage = "J[i,j] spin operators in full Hilbert space.";
TotalParity::usage = "Parity operator for both bosons and qutrit.";

(* Displacement operators *)
BosonDisplacementOperator::usage = "Creates a displacemnt operator with this parameter z."
QutritDisplacementOperator::usage = "Creates a qutrit displacment operator."
BothBosonDisplacementOperator::usage = "Creates a displacement operator acting on both bosons."
				       
(* Hamiltonian functions *)
OneModeHamiltonianOp::usage = "HamiltonianOp[ω0,Δ,φ,λ] returns the Hamiltonian operator of the one-mode Z3 Rabi model.";
OneModeEigenStatesAndEnergies::usage = "EigenStatesAndEnergies[ω0,Δ,φ,λ,shift,n] returns {energies,states} of the one-mode Z3 Rabi model.";
TwoModeHamiltonianOp::usage = "HamiltonianOp[ω0,Δ,φ,λ] returns the Hamiltonian operator of the two-mode Z3 Rabi model.";
TwoModeEigenStatesAndEnergies::usage = "EigenStatesAndEnergies[ω0,Δ,φ,λ,shift,n] returns {energies,states} of the two-mode Z3 Rabi model.";

(* Other functions *)
SineBasis::usage = "SineBasis[i,j] returns sine basis elements.";
NeumannEntropy::usage = "NeumannEntropy[state] computes von Neumann entropy.";

Begin["`Private`"]

(* Private variables *)
$NphotTrunc = Undefined;
$Initialized = False;

(* Error messages *)
SetNphotTrunc::positive = "NphotTrunc must be a positive integer.";
InitializeOps::notset = "NphotTrunc must be set before initialization.";
General::notinit = "Package not initialized. Call InitializeOps[] first.";

(* Helper function to check if initialized *)
checkInitialized[] := If[!$Initialized, Message[General::notinit]; $Failed]

(* Public functions *)
SetNphotTrunc[n_Integer?Positive] := (
  $NphotTrunc = n;
  $Initialized = False; (* Force re-initialization *)
  n
)
SetNphotTrunc[_] := (Message[SetNphotTrunc::positive]; $Failed)

GetNphotTrunc[] := $NphotTrunc

IsInitialized[] := $Initialized

InitializeOps[] := Module[{},
  If[$NphotTrunc === Undefined, 
    Message[InitializeOps::notset]; Return[$Failed]
  ];
  
  Print["Initializing operators with NphotTrunc = ", $NphotTrunc, "..."];
  
  (* Define basic constants *)
  ω = Exp[2Pi I/3];
  X = {{0,0,1},{1,0,0},{0,1,0}};
  Z = {{1,0,0},{0,ω,0},{0,0,ω^2}};
  QParity = {{1,0,0}, {0,0,1}, {0,1,0}};
  SineBasis[i_Integer,j_Integer] := Power[ω,i*j]MatrixPower[X,Mod[i,3]] . MatrixPower[Z,Mod[j,3]];
  
  (* Define operators in photonic space *)
  a = Table[If[j==i+1,Sqrt[i]//N,0],{i,1,$NphotTrunc},{j,1,$NphotTrunc}]//SparseArray;
  a† = Transpose[a]//SparseArray;
  idphot = IdentityMatrix[$NphotTrunc]//SparseArray;
  idel = IdentityMatrix[3]//SparseArray;
  photParity = DiagonalMatrix[Table[(-1)^i, {i, 0, $NphotTrunc-1}]]//SparseArray;
  
  (* Define operators in full Hilbert space *)
  Id = KroneckerProduct[idel,idphot,idphot]//SparseArray;
  A[1] = KroneckerProduct[idel,a,idphot]//SparseArray;
  A[2] = KroneckerProduct[idel,idphot,a]//SparseArray;
  A†[1] = KroneckerProduct[idel,a†,idphot]//SparseArray;
  A†[2] = KroneckerProduct[idel,idphot,a†]//SparseArray;
  PhotNum[1] = KroneckerProduct[idel,a† . a,idphot]//SparseArray;
  PhotNum[2] = KroneckerProduct[idel,idphot,a† . a]//SparseArray;
  PhotNumTotal = PhotNum[1]+PhotNum[2];

  (* Parity operators *)
  Parity[1] = KroneckerProduct[idel, photParity, idphot]//SparseArray;
  Parity[2] = KroneckerProduct[idel, idphot, photParity]//SparseArray;
  TotalParity = KroneckerProduct[QParity, photParity, photParity]//SparseArray;

  (* Define J operators *)
  J[i_,j_] := KroneckerProduct[SineBasis[i,j],idphot,idphot]//SparseArray;
  
  $Initialized = True;
  Print["Initialization complete."];
  True
]


(* Displacement operators *)
BosonDisplacementOperator[z_?NumericQ, bosonNum_] := Module[{},
  If[checkInitialized[] === $Failed, Return[$Failed]];
  MatrixExp[z * A†[bosonNum] - Conjugate[z] * A[bosonNum]]]
QutritDisplacementOperator[a_, b_] := KroneckerProduct[Power[ω, - a * b/2] * MatrixPower[Z, b] . MatrixPower[X, a], idphot, idphot] // SparseArray;
BothBosonDisplacementOperator[z1_?NumericQ, z2_?NumericQ] := Module[{},
  If[checkInitialized[] === $Failed, Return[$Failed]];
  KroneckerProduct[idel, MatrixExp[z1 * a† - Conjugate[z1] * a], MatrixExp[z2 * a† - Conjugate[z2] * a]]]

					  
(* Hamiltonian definitions - these check if initialized *)
(*HamiltonianOp[ω0_,Δ_,φ_,λ_] := Module[{},
  If[checkInitialized[] === $Failed, Return[$Failed]];
  ω0 PhotNumTotal + Δ (Exp[I φ]J[0,1]+Exp[-I φ]J[0,2]) - 
  λ/Sqrt[2] ((A[1]+A†[1])+I(A[2]+A†[2])) . J[1,0] - 
  λ/Sqrt[2] ((A[1]+A†[1])-I(A[2]+A†[2])) . J[2,0]]*)

OneModeHamiltonianOp[\[Omega]0_, \[CapitalDelta]_, \[Phi]_, \[Lambda]_] := Module[{},
  If[checkInitialized[] === $Failed, Return[$Failed]];
  \[Omega]0 KroneckerProduct[a\[Dagger] . a, idel] +
  \[CapitalDelta] KroneckerProduct[idphot,
    (Exp[I \[Phi]] SineBasis[0, 1] + 
     Exp[-I \[Phi]] SineBasis[0, 2])] -
  \[Lambda] KroneckerProduct[a, SineBasis[1, 0]] -
  \[Lambda] KroneckerProduct[a\[Dagger], SineBasis[2, 0]]
]									   

TwoModeHamiltonianOp[ω0_,Δ_,φ_,λ_] := Module[{},
  If[checkInitialized[] === $Failed, Return[$Failed]];
  ω0 PhotNumTotal + Δ (Exp[I φ]J[0,1]+Exp[-I φ]J[0,2]) - 
  λ(A[1] + A†[2]) . J[1,0] - 
  λ(A[2] + A†[1]) . J[2,0]
]


(* Wrapper functions for backwards compatibility *)
OneModeEigenStatesAndEnergies[ω0_,Δ_,φ_,λ_,energyShift_,numOfStates_] := Module[{},
  If[checkInitialized[] === $Failed, Return[$Failed]];
  Reverse@Transpose@Eigensystem[
	  OneModeHamiltonianOp[ω0,Δ,φ,λ] + energyShift KroneckerProduct[idphot, idel], -numOfStates
  ]
]

TwoModeEigenStatesAndEnergies[ω0_,Δ_,φ_,λ_,energyShift_,numOfStates_] := Module[{},
  If[checkInitialized[] === $Failed, Return[$Failed]];
  Reverse@Transpose@Eigensystem[
    TwoModeHamiltonianOp[ω0,Δ,φ,λ] + energyShift Id, -numOfStates
  ]
]

(* Entropy functions *)
ReshapeState[state_] := Module[{},
  If[checkInitialized[] === $Failed, Return[$Failed]];
  ArrayReshape[state,{3,$NphotTrunc,$NphotTrunc}]
]

NeumannEntropy[state_] := Module[{reshapedState, densityMatrix},
  If[checkInitialized[] === $Failed, Return[$Failed]];
  reshapedState = ReshapeState[state];
  densityMatrix = Table[
    Conjugate[Flatten[reshapedState[[All,i]]]] . Flatten[reshapedState[[All,j]]],
    {i,1,$NphotTrunc},{j,1,$NphotTrunc}
  ];
  -Tr[densityMatrix . MatrixLog[densityMatrix]]//Chop
]

End[] (* End Private context *)

EndPackage[]
