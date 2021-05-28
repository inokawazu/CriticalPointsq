#!/usr/bin/env wolframscript
(* ::Package:: *)

wp = 2*MachinePrecision;
ct = 10^-8;
lowerGridPoint = 40;
higherGridPoint = 45;
kRange = Subdivide[-5,5,250];

If[
	$ScriptCommandLine==={},
	SetDirectory@NotebookDirectory[],
	SetDirectory@DirectoryName@ExpandFileName@First[$ScriptCommandLine]
];

dataDir = FileNameJoin[{Directory[],"data"}];
Get[FileNameJoin[{Directory[],"source","PseudoSpectralUtils.wl"}]]
Get[FileNameJoin[{Directory[],"source","QNMEigenSolver.wl"}]]
Get[FileNameJoin[{Directory[],"source","BulkSolverUtils.wl"}]]
Get[FileNameJoin[{Directory[],"eoms","schwarzschild-black-brane-eoms.wls"}]]

If[Not@DirectoryQ[dataDir],CreateDirectory[dataDir]];

eoms = {tensorEom,vectorEom,scalarEom};
toTitle=Function[{sector},sector<>"-ct"<>ToString@Round@Abs@Log10[ct]<>"-wp"<>ToString@Round@Abs[wp]<>"-gp"<>ToString[lowerGridPoint]<>"-"<>ToString[higherGridPoint]<>".csv"];
dataFileNames = Map[FileNameJoin[{dataDir,#}]&]@Map[toTitle]@{"tensor","vector","scalar"};


Do[
	abMats = funcTable[alphabetamat,
		{{eoms[[eomi]]}, gp, {F}, u, w, gridps, derivativemat, WorkingPrecision->wp},
		{gp,lowerGridPoint,higherGridPoint}
	];
	
	qnm = DeleteDuplicates[
		qnmTable[modespseudo, {abMats/.{k -> kR+I kI},"convergenceTolerance"->ct}, 
		{kR,kRange},{kI,kRange}],
		LessThan[ct*10]@*Norm@*Subtract]//N;
		
	Export[dataFileNames[[eomi]],qnm]
,{eomi,Length@eoms}]
