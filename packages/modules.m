(* ::Package:: *)

ClearAll[tauMin, tauMax, HI, solutionOfPhi, z, f, epsilon, eta, xi, powerSpectrumAtHorizonCrossingEpoch, kOfY, lnPs];


ClearAll[v, tau, phi];


(* ::Text:: *)
(** The Background Part:*)


(* ::Text:: *)
(*The horizon crossing time of the boundary values of the observable window,*)


tauMin = -1/kMin bunchDavisVacuumEpochFactor;
tauMax = -1/kMax;
(* or, a little wilder: *)
tauMin = -1/kMin bunchDavisVacuumEpochFactor (1+0.1);
tauMax = -1/kMax (1-0.1);


(* ::Text:: *)
(*By slow-rolling conditions, the Hubble parameter in inflation epoch is*)


HI[tau_] := Sqrt[V[phi[tau]]/3];


solutionOfPhi = NDSolve[{phi'[tau] == 1/tau V'[phi[tau]]/V[phi[tau]], phi[-1/kN] == phiN}, phi, {tau, tauMin, tauMax}];


z[tau_] := Evaluate[1/HI[tau] 1/tau V'[phi[tau]]/V[phi[tau]]/.solutionOfPhi];


(* ::Text:: *)
(** The Perturbation Part:*)


(* ::Text:: *)
(*The z''/z =: 1/\[Tau]^2 f(\[Phi]) in Mukharnov equation:*)


(* ::Text:: *)
(*Generally, it has the form:*)


f[phi_] := Evaluate[Simplify[2+10 epsilon[phi]^2+epsilon[phi] (5-9 eta[phi])-3 eta[phi]+eta[phi]^2+xi[phi]]];


(* ::Text:: *)
(*where the slow-rolling parameters are defined as:*)


epsilon[phi_] := 1/2 (V'[phi]/V[phi]);
eta[phi_] := V''[phi]/V[phi];
xi[phi_] := V'[phi]/V[phi] V'''[phi]/V[phi];


u[tau_] := v[tau]/z[tau];


(* ::Text:: *)
(*Thus, the power specturm:*)


powerSpectrumAtHorizonCrossingEpoch[k_] :=
Module[{tauK, tauI, solutionOfMukharnovEqu, powerSpectrumAtHorizonCrossingEpochOfKMode},
	(* the time of horizon-crossing epoch of this k-mode: *)
	tauK = -1/k;
	(* the time of epoch of establishing the Bunch-Davis vacuum: *)
	tauI = bunchDavisVacuumEpochFactor*tauK;
	(* numerically solve the Mukharnov equation: *)
	solutionOfMukharnovEqu = NDSolve[{v''[tau] + (k^2 - 1/tau^2 f[Evaluate[phi[tau]/.solutionOfPhi]]) v[tau] == 0, v[tauI] == 1/Sqrt[2k] Exp[-I k tauI], v'[tauI] == (-I k)/Sqrt[2k] Exp[-I k tauI]}, v, {tau, tauI, tauK}];
	(* thus the power spectrum at horizon-crossing epoch of this k-mode: *)
	powerSpectrumAtHorizonCrossingEpochOfKMode = (Evaluate[k^3 Abs[u[tau]/.solutionOfMukharnovEqu]^2 ]/.{tau -> -1/k})[[1,1]];
	Return@powerSpectrumAtHorizonCrossingEpochOfKMode]


(* ::Text:: *)
(*Define y := ln (k/kN) (thus k = kN Exp(y)),*)
(*g(y) := ln Ps (y):*)


kOfY[y_] := kN Exp[y];
lnPs[y_] := Log[powerSpectrumAtHorizonCrossingEpoch[kOfY[y]]];
