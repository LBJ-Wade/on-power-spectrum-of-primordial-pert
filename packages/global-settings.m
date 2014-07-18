(* ::Package:: *)

(* ::Text:: *)
(*The pivot scale (wave number):*)


kN = 0.05;


(* ::Text:: *)
(*On potential (do NOT change the name (head) of this function):*)


(* ::Text:: *)
(*Such as,*)


V[phi_] := phi^4;
phiN = 20;


(* ::Text:: *)
(*On solving the background e.o.m.:*)


(* ::Text:: *)
(*The observable window:*)


kMax = 5;
kMin = 0.00005;
