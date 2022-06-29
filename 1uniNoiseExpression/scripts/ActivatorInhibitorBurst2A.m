(* ::Package:: *)

(* ::Text:: *)
(*In this notebook is simulated the Activator-Inhibitor system regulation with GA and the three-stage model (burst expression). The notebook has two parts:*)
(*1. To evaluate the "sensitivity" to the whole set the parameters. With an arbitrary start point to calculate the estimates (FF, CV2,Mean) (30% of the total number of iterations), and different number of iterations to get a final t of 1500. *)
(*2. To evaluate 4 specific regions or set of parameters*)


(* ::Section:: *)
(*Parameters for the two sessions*)


(* ::Text:: *)
(*These are the ranges of values for parameters variation. The range of each parameter is divided into smaller ranges (subRanges). In this way the complete range is evaluated in parallel files each one with a subRange of parameters values.*)


(* ::Input::Initialization:: *)
kg\[Gamma]sListBySample=Import[(*"/media/lina/DATA/noiseGRNdevelopmentLinaMRuizG_git/BiologicalNoise_DevelopmentalGeneRegulatoryNetworks/1uniNoiseExpression/scripts/*)"kg\[Gamma]gSamples.csv"]//ToExpression;
sample=1;
paraV=kg\[Gamma]sListBySample[[sample]];


(* ::Text:: *)
(*These are the default values of parameters*)


(* ::Input::Initialization:: *)
{kgA,kgI}={0.01158,0.01158};(* Gene activation rate (t^-1)*)
{\[Gamma]gA,\[Gamma]gI}= {2.082,2.082};(* Gene desactivation rate (t^-1) *)
{kmA,kmI}={0.696,0.696};(* Production rate of mRNA (molecules*t^-1) *)
{\[Gamma]mA,\[Gamma]mI}={0.02082  ,0.02082 };(* mRNA Degradation rate (molecules*t^-1) *)
{kpA,kpI}={1.386,1.386};(* protein Production (molecules*t^-1) *)
{\[Gamma]pA,\[Gamma]pI}={0.02082  ,0.02082 }(* protein Degradation rate (molecules*t^-1) *)
h=3 (* Hill's constant *);
kIA=kAA=kAI=100;(* Michaeles-Menten Constant *)

{bmA,bmI}={kmA/(\[Gamma]gA+kgA),kmI/(\[Gamma]gI+kgI)};(* Burst size mean (molecules-mRNA) *)
{bpA,bpI}={kpA/\[Gamma]mA,kpI/\[Gamma]mI};(* Burst size mean (molecules-proteins) *)
kaa=\[Gamma]g/kg(* Michaeles-Menten Constant of regulation of a by a *);
kia=\[Gamma]g/kg(* Michaeles-Menten Constant of regulation of a by i *);
kai=\[Gamma]g/kg(* Michaeles-Menten Constant of regulation of i by a *);
kat=\[Gamma]g/kg(* Michaeles-Menten Constant of regulation of t by a *);




(* ::Section:: *)
(*Functions for the two sessions*)


(* ::Text:: *)
(*These are the reactions (changesByRXN) and propensity functions (pF) of the GA for the specific regulation system: A: Activator -> I : Inhibitor gene*)


(* ::Input::Initialization:: *)
changesByRXN=Function/@<|1->{#1,#2,#3,#4,#5-1,#6+1,#7,#8},2->{#1,#2,#3,#4,#5+1,#6-1,#7,#8},3->{#1+1,#2,#3,#4,#5,#6,#7,#8},4->{#1,#2+1,#3,#4,#5,#6,#7,#8},5->{#1-1,#2,#3,#4,#5,#6,#7,#8},6->{#1,#2-1,#3,#4,#5,#6,#7,#8},7->{#1,#2,#3,#4,#5,#6,#7-1,#8+1},8->{#1,#2,#3,#4,#5,#6,#7+1,#8-1},9->{#1,#2,#3+1,#4,#5,#6,#7,#8},10->{#1,#2,#3,#4+1,#5,#6,#7,#8},11->{#1,#2,#3-1,#4,#5,#6,#7,#8},12->{#1,#2,#3,#4-1,#5,#6,#7,#8}|>;(* Each hash represents a specie: mRNAA-proteinA-mRNAB-proteinB-da0A:promoterInactivatedA-da1A:promoterActivatedA-db0B-db1B
   The numbers are the next reactions: 1-7\[Rule] Activation, 2-8\[Rule] Desactivation, 3-4-9-10\[Rule] Synthesis, 5-6-11-12\[Rule] Degradation*)

pF=Compile[{{h,_Real},{kIA,_Real},{kAA,_Real},{kAI,_Real},{kgA,_Real},{\[Gamma]gA,_Real},{kmA,_Real},{\[Gamma]mA,_Real},{kpA,_Real},{\[Gamma]pA,_Real},{kgI,_Real},{\[Gamma]gI,_Real},{kmI,_Real},{\[Gamma]mI,_Real},{kpI,_Real},{\[Gamma]pI,_Real},{mRNAA,_Integer},{proteinA,_Integer},{mRNAI,_Integer},{proteinI,_Integer},{d0A,_Integer},{d1A,_Integer},{d0I,_Integer},{d1I,_Integer}},{d0A*kgA*(1/(1+(proteinI/kIA)^h))(proteinA^h/(kAA+proteinA ^h))(*La inhibici\[OAcute]n de I hacia A esta aqui*),d1A*\[Gamma]gA,d1A*kmA+0.01,kpA*mRNAA,mRNAA*\[Gamma]mA,\[Gamma]pA*proteinA,d0I*kgI*(proteinA^h/(kAI+proteinA ^h)),d1I*\[Gamma]gI,d1I*kmI+0.01,kpI*mRNAI,mRNAI*\[Gamma]mI,\[Gamma]pI*proteinI},RuntimeAttributes->{Listable},Parallelization->True];


(* ::Text:: *)
(*This is equal for all  sub-circuits*)


(* ::Input::Initialization:: *)
ff[dato_]:=N[Variance[dato]/Mean[dato]]
cv[data_]:=N[StandardDeviation[data]/Mean[data]]
cv2[data_]:=cv[data]^2


(* ::Input::Initialization:: *)
(* This is the algorithm to estimate the Steady State start from Kelly and Hedengren 2013 *)
steadyStateStart[window_,listData_]:=Block[{n,dataByWindow,mPorVentana,\[Mu],sd,probabilidades},
n=Floor[Length[#]/window]&[listData];
dataByWindow=Partition[#,n]&[listData[[All,2]]];mPorVentana=Map[Mean[#[[2;;All]]-#[[1;;-2]]]&,dataByWindow];\[Mu]=MapThread[(1/n)(Total[#1]-#2*Total[Range[1,n](*tPorVentana*)])&,{dataByWindow,mPorVentana}];sd=Table[Sqrt[(1/(n-2))*Sum[(dataByWindow[[i,j]]-mPorVentana[[i]]*j-\[Mu][[i]])^2,{j,1,Length[dataByWindow[[i]]]}]],{i,window}];probabilidades=Map[Total[#]/n&,Table[If[Abs[dataByWindow[[i,j]]-\[Mu][[i]]]<=sd[[i]],1,0],{i,window},{j,1,n}]]//N;
((Position[probabilidades,x_/;x>0.1][[1,1]]-1)*n)+1]


(* ::Input:: *)
(**)


(* ::Section:: *)
(*S1: All range parameters values*)


(* ::Subsection:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
iteraciones1=400000; (* Number of algorithm iterations *)
replicas1=100; (* Number of times the algorithm is run *)


(* ::Subsection:: *)
(*Simulation*)


(* ::Input::Initialization:: *)
(* Launch the number of kernels in which you will parallelize the simulation *)
LaunchKernels[44] 


(* ::Text:: *)
(*This function runs the GA algorithm and processes the output according to the previous functions.*)


(* ::Input::Initialization:: *)
DistributeDefinitions[changesByRXN,pF,ff,expressionNoise,iteraciones,h,kIA,kAA,kAI,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kgI,\[Gamma]gI,kmI,\[Gamma]mI,kpI,\[Gamma]pI,cv,cv2,iteraciones1,replicas1];
(*List by parameters*)
mRNAffcv2ListA={};
proteinffcv2ListA={};
mRNAffcv2ListI={};
proteinffcv2ListI={};

Do[
(*List by replicates*)
 listEstimatesReplicatesA = {};
 listEstimatesReplicatesI = {};

SetSharedVariable[ listEstimatesReplicatesA,listEstimatesReplicatesI,listEstimatemRNAA,listEstimateproteinA,listEstimatemRNAI,listEstimateproteinI, listDataAndTimeT];
(*******************************)
{kgA,\[Gamma]gA}={para[[1]],para[[2]]};
(*******************************)

ParallelDo[

listDataAndTime = CreateDataStructure["OrderedHashSet"];
abundancia={1,1,1,1,1,0,1,0}(*mRNAA,proteinA,mRNAI,proteinI,da0A,da1A,db0I,db1I*);
t=0;

Do[
\[Tau]is=Map[If[#>0,RandomVariate[ExponentialDistribution[#]],\[Infinity]]&,pF[h,kIA,kAA,kAI,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kgI,\[Gamma]gI,kmI,\[Gamma]mI,kpI,\[Gamma]pI,Sequence@@abundancia]];
reactionAnd\[Tau]\[Mu]={FirstPosition[#,Min[#]],Min[#]}&[\[Tau]is];abundancia=Through[{reactionAnd\[Tau]\[Mu][[1]]/.changesByRXN}[[1,1]]@@abundancia][[1]];t=t+reactionAnd\[Tau]\[Mu][[2]];         
 listDataAndTime["Union", {{t, abundancia}}];
            
If[Length[#]>3&&#[[2]]==5&[IntegerDigits[Round[t]]](* para que Break a los 1500*), Break[]],iteraciones1 ];

  listData = Normal[listDataAndTime][[All, 2]];
        ini =Round[(Length[listData]*30)/100.] (*el estado estable se considera depu\[EAcute]s del 30% de iteraciones *);
(*ini = steadyStateStart[window,listData];*)
{mRNAffcv2ItA,proteinffcv2ItA}=Table[{ff[i], cv[i],cv2[i], Mean[i],StandardDeviation[i],Variance[i],Kurtosis[i],Skewness[i],Table[Moment[i,j],{j,10}],Table[CentralMoment[i,j],{j,10}]},{i,{listData[[ini;;All,1]],
listData[[ini;;All,2]]}}];

{mRNAffcv2ItI,proteinffcv2ItI}=Table[{ff[i], cv[i],cv2[i], Mean[i],StandardDeviation[i],Variance[i],Kurtosis[i],Skewness[i],Table[Moment[i,j],{j,10}],Table[CentralMoment[i,j],{j,10}]},{i,{listData[[ini;;All,3]],
listData[[ini;;All,4]]}}];

 AppendTo[listEstimatesReplicatesA, {mRNAffcv2ItA,proteinffcv2ItA}];
AppendTo[listEstimatesReplicatesI, {mRNAffcv2ItI,proteinffcv2ItI}];

,
replicas1
];

  Export["controlAIGA_Asample\[CapitalDelta]A_"<>ToString[sample]<>","<>ToString[para]<>".csv",listEstimatesReplicatesA];
  Export["controlAIGA_Isample\[CapitalDelta]A_"<>ToString[sample]<>","<>ToString[para]<>".csv",listEstimatesReplicatesI];
,
    
 {para,paraV}]; (*CAMBIAR AQUI EN EL .m*)
