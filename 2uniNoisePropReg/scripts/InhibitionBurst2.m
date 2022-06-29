(* ::Package:: *)

(* ::Text:: *)
(*This scripts is almost the same than InhibitionBurst . m but here we only run the simulation for a set of parameters that were chosen randomly for grids*)


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
{kgA,kgB}={0.01158,0.01158};(* Gene activation rate (t^-1)*)
{\[Gamma]gA,\[Gamma]gB}= {2.082,2.082};(* Gene desactivation rate (t^-1) *)
{kmA,kmB}={0.696,0.696};(* mRNA production rate(molecules*t^-1) *)
{\[Gamma]mA,\[Gamma]mB}={0.02082  ,0.02082 } ;(* mRNA degradation rate(molecules*t^-1) *)
{kpA,kpB}={1.386,1.386};(* Protein production rate (molecules*t^-1) *)
{\[Gamma]pA,\[Gamma]pB}={0.02082  ,0.02082 };(* Protein degradation rate (molecules*t^-1) *)
kAB=100;(* Michaelis-Menten constant *)
h=3 (* Hill's constant *);
{bmA,bmB}={kmA/(\[Gamma]gA+kgA),kmB/(\[Gamma]gB+kgB)};(* Burst size mean (molecules-mRNA) *)
{bpA,bpB}={kpA/\[Gamma]mA,kpB/\[Gamma]mB};(* Burst size mean (molecules-proteins) *)



(* ::Section:: *)
(*Functions for the two sessions*)


(* ::Subsection:: *)
(*Functions*)


(* ::Text:: *)
(*These are the reactions (changesByRXN) and propensity functions (pF) of the GA for the specific regulation system: A: Inhibitor -> B : Inhibited gene*)


(* ::Input::Initialization::Plain:: *)
changesByRXN=Function/@<|1->{#1,#2,#3,#4,#5-1,#6+1,#7,#8},2->{#1,#2,#3,#4,#5+1,#6-1,#7,#8},3->{#1+1,#2,#3,#4,#5,#6,#7,#8},4->{#1,#2+1,#3,#4,#5,#6,#7,#8},5->{#1-1,#2,#3,#4,#5,#6,#7,#8},6->{#1,#2-1,#3,#4,#5,#6,#7,#8},7->{#1,#2,#3,#4,#5,#6,#7-1,#8+1},8->{#1,#2,#3,#4,#5,#6,#7+1,#8-1},9->{#1,#2,#3+1,#4,#5,#6,#7,#8},10->{#1,#2,#3,#4+1,#5,#6,#7,#8},11->{#1,#2,#3-1,#4,#5,#6,#7,#8},12->{#1,#2,#3,#4-1,#5,#6,#7,#8}|>;(* Each hash represents a specie: mRNAA-proteinA-mRNAB-proteinB-da0A:promoterActivatedA-da1A:promoterInactivatedA-db0B-db1B
   The numbers are the next reactions:
 1-7\[Rule] Activation,
 2-8\[Rule] Desactivation,
 3-4-9-10\[Rule] Synthesis,
5-6-11-12\[Rule] Degradation*)
pF=Compile[{{h,_Real},{kAB,_Real},{kgA,_Real},{\[Gamma]gA,_Real},{kmA,_Real},{\[Gamma]mA,_Real},{kpA,_Real},{\[Gamma]pA,_Real},{kgB,_Real},{\[Gamma]gB,_Real},{kmB,_Real},{\[Gamma]mB,_Real},{kpB,_Real},{\[Gamma]pB,_Real},{mRNAA,_Integer},{proteinA,_Integer},{mRNAB,_Integer},{proteinB,_Integer},{d0A,_Integer},{d1A,_Integer},{d0B,_Integer},{d1B,_Integer}},{d0A*kgA,d1A*\[Gamma]gA,d1A*kmA+0.01,kpA*mRNAA,mRNAA*\[Gamma]mA,\[Gamma]pA*proteinA,d0B*kgB(1/(1+(proteinA/kAB)^h)),d1B*\[Gamma]gB,d1B*kmB+0.01,kpB*mRNAB,mRNAB*\[Gamma]mB,\[Gamma]pB*proteinB},RuntimeAttributes->{Listable},Parallelization->True];


(* ::Text:: *)
(*This is equal for all  sub-circuits*)


(* ::Input::Initialization:: *)
ff[dato_]:=N[Variance[dato]/Mean[dato]]
cv[data_]:=N[StandardDeviation[data]/Mean[data]]
cv2[data_]:=cv[data]^2


(* ::Input::Initialization:: *)
(* This is the algorithm to estimate the Steady State start from Kelly and Hedengren 2013 *)steadyStateStart[window_,listData_]:=Block[{n,dataByWindow,mPorVentana,\[Mu],sd,probabilidades},
n=Floor[Length[#]/window]&[listData];
dataByWindow=Partition[#,n]&[listData[[All,2]]];mPorVentana=Map[Mean[#[[2;;All]]-#[[1;;-2]]]&,dataByWindow];\[Mu]=MapThread[(1/n)(Total[#1]-#2*Total[Range[1,n](*tPorVentana*)])&,{dataByWindow,mPorVentana}];sd=Table[Sqrt[(1/(n-2))*Sum[(dataByWindow[[i,j]]-mPorVentana[[i]]*j-\[Mu][[i]])^2,{j,1,Length[dataByWindow[[i]]]}]],{i,window}];probabilidades=Map[Total[#]/n&,Table[If[Abs[dataByWindow[[i,j]]-\[Mu][[i]]]<=sd[[i]],1,0],{i,window},{j,1,n}]]//N;
((Position[probabilidades,x_/;x>0.1][[1,1]]-1)*n)+1]


(* ::Section:: *)
(*S1:  All range parameters values*)


(* ::Subsection:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
iteraciones2=400000; (* Number of algorithm iterations *)
replicas2=100; (* Number of times the algorithm is run *)


(* ::Subsection:: *)
(*Simulation*)


(* ::Input::Initialization:: *)
(* Launch the number of kernels in which you will parallelize the simulation *)
LaunchKernels[44] 


(* ::Text:: *)
(*This function runs the GA algorithm and processes the output according to the previous functions.*)


(* ::Input::Initialization:: *)
DistributeDefinitions[changesByRXN,pF,ff,expressionNoise,iteraciones2,replicas2,h,kAB,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kgB,\[Gamma]gB,kmB,\[Gamma]mB,kpB,\[Gamma]pB];
(*List by parameters*)
mRNAffcv2ListA={};
proteinffcv2ListA={};
mRNAffcv2ListI={};
proteinffcv2ListI={};

Do[
(*List by replicates*)
  listEstimatesReplicatesA = {};
  listEstimatesReplicatesB = {};

SetSharedVariable[listEstimatesReplicatesA,listEstimatesReplicatesB];

(*******************************)
{kgB,\[Gamma]gB}={para[[1]],para[[2]]};
(*******************************)

ParallelDo[

listDataAndTime = CreateDataStructure["OrderedHashSet"];
abundancia={1,1,1,1,1,0,1,0}(*mRNAA,proteinA,mRNAI,proteinI,da0A,da1A,db0I,db1I*);
t=0;


Do[
\[Tau]is=Map[If[#>0,RandomVariate[ExponentialDistribution[#]],\[Infinity]]&,pF[h,kAB,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kgB,\[Gamma]gB,kmB,\[Gamma]mB,kpB,\[Gamma]pB,Sequence@@abundancia]];
reactionAnd\[Tau]\[Mu]={FirstPosition[#,Min[#]],Min[#]}&[\[Tau]is];abundancia=Through[{reactionAnd\[Tau]\[Mu][[1]]/.changesByRXN}[[1,1]]@@abundancia][[1]];t=t+reactionAnd\[Tau]\[Mu][[2]];         
 listDataAndTime["Union", {{t, abundancia}}];
            If[Length[#]>3&&#[[2]]==5&[IntegerDigits[Round[t]]](* para que Break a los 1500*), Break[]],iteraciones2];

  listData = Normal[listDataAndTime][[All, 2]];
        ini =Round[(Length[listData]*30)/100.] (*el estado estable se considera depu\[EAcute]s del 30% de iteraciones *);
(*ini = steadyStateStart[window,listData];*)
{mRNAffcv2ItA,proteinffcv2ItA}=Table[{ff[i], cv[i],cv2[i], Mean[i],StandardDeviation[i],Variance[i],Kurtosis[i],Skewness[i],Table[Moment[i,j],{j,10}],Table[CentralMoment[i,j],{j,10}]},{i,{listData[[ini;;All,1]],
listData[[ini;;All,2]]}}];
{mRNAffcv2ItB,proteinffcv2ItB}=Table[{ff[i], cv[i],cv2[i], Mean[i],StandardDeviation[i],Variance[i],Kurtosis[i],Skewness[i],Table[Moment[i,j],{j,10}],Table[CentralMoment[i,j],{j,10}]},{i,{listData[[ini;;All,3]],
listData[[ini;;All,4]]}}];

  AppendTo[listEstimatesReplicatesA, {mRNAffcv2ItA,proteinffcv2ItA}];    
  AppendTo[listEstimatesReplicatesB, {mRNAffcv2ItB,proteinffcv2ItB}];    
,
replicas2];

    Export["controlInhGA_Asample\[CapitalDelta]B_"<>ToString[sample]<>","<>ToString[para]<>".csv",listEstimatesReplicatesA];
    Export["controlInhGA_Bsample\[CapitalDelta]B_"<>ToString[sample]<>","<>ToString[para]<>".csv",listEstimatesReplicatesB];
    ,
 {para,paraV}];

