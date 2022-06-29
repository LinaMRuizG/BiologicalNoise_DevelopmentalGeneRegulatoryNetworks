(* ::Package:: *)

(* ::Text:: *)
(*In this notebook is simulated a One self-Regulated gene with GA and the three-stage model (burst expression). The notebook has two parts:*)
(*1. To evaluate the "sensitivity" to the whole set the parameters. With an arbitrary start point to calculate the estimates (FF, CV2,Mean) (30% of the total number of iterations), and different number of iterations to get a final t of 1500. *)
(*3. To evaluate 4 specific regions*)


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
Clear[kg,\[Gamma]g,km,\[Gamma]m,kp,\[Gamma]p,\[Tau]];
kg=0.01158;(* Gene activation rate (t^-1)*)
\[Gamma]g=2.082;(* Gene desactivation rate (t^-1) *)
km=0.696;(* Production rate of mRNA (molecules*t^-1) *)
\[Gamma]m=0.02082(*degradation mRNA*);
kp=1.386(*syntesis proteins*);
\[Gamma]p=0.02082(*degradation proteins*);
h=3.0; (* Michaelis-Menten constant *)
kAA=100(* Hill's constant *);


(* ::Section:: *)
(*Functions for the two sessions*)


(* ::Text:: *)
(*These are the reactions (changesByRXN) and propensity functions (pF) of the GA for the specific regulation system: self-regulated gene*)


(* ::Input::Initialization:: *)
changesByRXN=Function/@<|1->{#1,#2,#3-1,#4+1},2->{#1,#2,#3+1,#4-1},3->{#1+1,#2,#3,#4},4->{#1-1,#2,#3,#4},5->{#1,#2+1,#3,#4},6->{#1,#2-1,#3,#4}|>;(* Each hash represents a specie: mRNA,protein,da0:InactivatedPromoter,da1:ActivatedPromoter. The numbers are the next reactions: 1-activation, 2-deactivation, 3-5 synthesis, 4-6 degradation *)pF=Compile[{{kg,_Real},{\[Gamma]g,_Real},{km,_Real},{\[Gamma]m,_Real},{kp,_Real},{\[Gamma]p,_Real},{h,_Real},{kAA,_Real},{a,_Integer},{p,_Integer},{da0,_Integer},{da1,_Integer}},{da0*kg*(p ^h/(kAA+p ^h)),da1*\[Gamma]g,da1*km+0.01(*Esto para mantener una expresi\[OAcute]n basal como en el de Langevin (que no tiene ceros)y que el kon no se haga cero*),a*\[Gamma]m,kp*a,\[Gamma]p*p},RuntimeAttributes->{Listable},Parallelization->True];


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


(* ::Input::Initialization:: *)
plotsLogExpression[rep_,ini_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,parametros_]:=Module[{datesListAm,datesListAp},
SetOptions[ListLogPlot,Frame-> True,ImageSize->Medium,Joined->True];

{datesListAm,datesListAp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,2}];{ListLogPlot[datesListAm[[ini;;]],PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA A molecules"},PlotLabel->labelA],
ListLogPlot[datesListAp[[ini;;]],PlotStyle->
RGBColor[0.2,0.8,1],FrameLabel->{"time","protein A molecules"},PlotLabel->labelAp]}]


(* ::Section:: *)
(*S1:  All range parameters values*)


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
(*Here there is a detail that changes by sub-circuit : the number and type of Kmm (put it in the arguments of the function simu, in the DistributeDefinitions, and pF function)*)


(* ::Text:: *)
(*This function runs the GA algorithm and processes the output according to the previous functions.*)


(* ::Input::Initialization:: *)
DistributeDefinitions[pF,h,kAA,kg,\[Gamma]g,km,\[Gamma]m,kp,\[Gamma]p,changesByRXN,ff,cv,cv2,iteraciones1,replicas1];
(*List by parameters*)
mRNAffcv2List={};
proteinffcv2List={};

Do[
(*List by replicates*)
(*listEstimatemRNA={};
     listEstimateprotein={};
 listDataAndTimeT={};*)
 listEstimatesReplicates = {};

SetSharedVariable[ listEstimatesReplicates ,listEstimatemRNA,listEstimateprotein,listDataAndTimeT];

(*******************************)
{kg,\[Gamma]g}={para[[1]],para[[2]]};
(*******************************)

ParallelDo[

listDataAndTime = CreateDataStructure["OrderedHashSet"];
abundancia={1,1,1,0}(*mRNA,protein,da0,da1*);
t=0;

Do[
\[Tau]is=Map[If[#>0,RandomVariate[ExponentialDistribution[#]],\[Infinity]]&,pF[kg,\[Gamma]g,km,\[Gamma]m,kp,\[Gamma]p,h,kAA,Sequence@@abundancia]];reactionAnd\[Tau]\[Mu]={FirstPosition[#,Min[#]],Min[#]}&[\[Tau]is];abundancia=Through[{reactionAnd\[Tau]\[Mu][[1]]/.changesByRXN}[[1,1]]@@abundancia][[1]];t=t+reactionAnd\[Tau]\[Mu][[2]];         
 listDataAndTime["Union", {{t, abundancia}}];
            If[Length[#]>3&&#[[2]]==5&[IntegerDigits[Round[t]]](* para que Break a los 1500*), Break[]],iteraciones1 ];

  listData = Normal[listDataAndTime][[All, 2]];
        ini =Round[(Length[listData]*30)/100.] (*el estado estable se considera depu\[EAcute]s del 30% de iteraciones *);
(*ini = steadyStateStart[window,listData];*)
{mRNAffcv2It,proteinffcv2It}=Table[{ff[i], cv[i],cv2[i], Mean[i],StandardDeviation[i],Variance[i],Kurtosis[i],Skewness[i],Table[Moment[i,j],{j,10}],Table[CentralMoment[i,j],{j,10}]},{i,{listData[[ini;;All,1]],
listData[[ini;;All,2]]}}];

(*AppendTo[listEstimatemRNA,mRNAffcv2It];
AppendTo[listEstimateprotein,proteinffcv2It];
  AppendTo[listDataAndTimeT, listDataAndTime];*)
 AppendTo[listEstimatesReplicates, {mRNAffcv2It,proteinffcv2It}];
,
replicas1
];

(*AppendTo[mRNAffcv2List,Mean[listEstimatemRNA]];
AppendTo[proteinffcv2List,Mean[listEstimateprotein]]*)
  Export["controlRGA_sample_"<>ToString[sample]<>","<>ToString[para]<>".csv",listEstimatesReplicates];
,
    
 {para,paraV}]; (*CAMBIAR AQUI EN EL .m*)
