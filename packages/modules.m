(* ::Package:: *)

ClearAll[tauMin, tauMax, HI, solutionOfPhi, z, u, powerSpectrumAtHorizonCrossingEpoch, kOfY, lnPs];


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


z[tau_] := Evaluate[(phi'[tau]/HI[tau]^2 1/tau)/.solutionOfPhi];


(* ::Text:: *)
(** The Perturbation Part:*)


u[tau_] := v[tau]/z[tau];


(* ::Text:: *)
(*Thus, the power specturm:*)


(* ::Text:: *)
(*The z''/z in Mukharnov equation:*)


gravitationalEffect[tau_] := Evaluate[(z''[tau]/z[tau])/.solutionOfPhi];


powerSpectrumAtHorizonCrossingEpoch[k_] :=
Module[{tauK, tauI, solutionOfMukharnovEqu, powerSpectrumAtHorizonCrossingEpochOfKMode},
	(* the time of horizon-crossing epoch of this k-mode: *)
	tauK = -1/k;
	(* the time of epoch of establishing the Bunch-Davis vacuum: *)
	tauI = bunchDavisVacuumEpochFactor*tauK;
	(* numerically solve the Mukharnov equation: *)
	solutionOfMukharnovEqu = NDSolve[{v''[tau] + (k^2 - gravitationalEffect[tau]) v[tau] == 0, v[tauI] == 1/Sqrt[2k] Exp[-I k tauI], v'[tauI] == (-I k)/Sqrt[2k] Exp[-I k tauI]}, v, {tau, tauI, tauK}];
	(* thus the power spectrum at horizon-crossing epoch of this k-mode: *)
	powerSpectrumAtHorizonCrossingEpochOfKMode = (Evaluate[k^3 Abs[u[tau]/.solutionOfMukharnovEqu]^2 ]/.{tau -> -1/k})[[1,1]];
	Return@powerSpectrumAtHorizonCrossingEpochOfKMode]


(* ::Text:: *)
(*Define y := ln (k/kN) (thus k = kN Exp(y)),*)
(*g(y) := ln Ps (y):*)


kOfY[y_] := kN Exp[y];
lnPs[y_] := Log[powerSpectrumAtHorizonCrossingEpoch[kOfY[y]]];
