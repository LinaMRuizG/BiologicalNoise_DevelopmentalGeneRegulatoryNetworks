(* ::Package:: *)

(* ::Text:: *)
(*In this notebook is simulated a One UnR gene with GA and the three-stage model (burst expression). The notebook has two parts:*)
(*1. To evaluate the "sensitivity" to the whole set the parameters. With an arbitrary start point to calculate the estimates (FF, CV2,Mean) (30% of the total number of iterations are removed to estimate these values in the steady-state), and different number of iterations to get a final t of 1500. *)
(*3. To evaluate 4 specific regions*)


(* ::Section:: *)
(*Parameters for the two sessions*)


(* ::Text:: *)
(*These are the ranges of values for parameters variation. The range of each parameter is divided into smaller ranges (subRanges). In this way the complete range is evaluated in parallel files each one with a subRange of parameters values.*)


(* ::Input::Initialization:: *)
name=0.5; (* This is the last value of the range for the parameter being evluated *)
kgs=Range[0.01,1,0.01];
\[Gamma]gs=Range[0.02,name,0.02]; (* subRanges: 0.02-0.5, 0.52-1.0, 1.02-1.5, 1.52-2 *)
\[Gamma]ms=Range[0.001,0.1(*name*),0.001];(* subRanges: 0.001-0.025, 0.026-0.05, 0.051-0.075, 0.076-0.1 *)
\[Gamma]ps=Range[0.001,0.1(*name*),0.001];(* subRanges: 0.001-0.025, 0.026-0.05, 0.051-0.075, 0.076-0.1 *)
parametro="\[Gamma]g";(* With the name of the second parameter is named the output file *)



(* ::Text:: *)
(*These are the default values of parameters*)


(* ::Input::Initialization:: *)
kgA=0.01158;(* Gene activation rate (t^-1)*)
\[Gamma]gA= 2.082;(* Gene desactivation rate (t^-1) *)
kmA=0.696;(* Production rate of mRNA (molecules*t^-1) *)
\[Gamma]mA=0.02082 ;(* mRNA Degradation rate (molecules*t^-1) *)
kpA=1.386;(* Production rate of proteins (molecules*t^-1) *)
\[Gamma]pA=0.02082;(* Degradation rate (molecules*t^-1) *)
(*THESE PARAMETERS ARE A GOOD CHOICE BECAUSE THE NOISE PATTERNS CAN BE SEEN AND NOT SO MANY ITERATIONS ARE REQUIRED TO REACH THE STABLE STATE DUE TO THE DEGRADATION/PRODUCTION RATIO*)


(* ::Section:: *)
(*Functions for the two sessions*)


(* ::Text:: *)
(*These are the reactions (changesByRXN) and propensity functions (pF) of the GA for an unRegulated gene*)


(* ::Input::Initialization:: *)
changesByRXN=Function/@<|1->{#1,#2,#3-1,#4+1},2->{#1,#2,#3+1,#4-1},3->{#1+1,#2,#3,#4},4->{#1-1,#2,#3,#4},5->{#1,#2+1,#3,#4},6->{#1,#2-1,#3,#4}|>;(* Each hash represents a specie: mRNA,protein,da0:promoterInactivated,da1:promoterActivated. The numbers are the next reactions: 1-activation, 2-deactivation, 3-5 synthesis, 4-6 degradation *)

pF=Compile[{{kg,_Real},{\[Gamma]g,_Real},{km,_Real},{\[Gamma]m,_Real},{kp,_Real},{\[Gamma]p,_Real},{a,_Integer},{p,_Integer},{da0,_Integer},{da1,_Integer}},{da0*kg,da1*\[Gamma]g,da1*km+0.01(*This allows a basal expression*),a*\[Gamma]m,kp*a,\[Gamma]p*p},RuntimeAttributes->{Listable},Parallelization->True];


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


(* ::Input::Initialization:: *)
(* This plots the temporal dynamic of genes mRNA and proteins with log scale *)plotsLogExpression[rep_,ini_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,parametros_]:=Module[{datesListAm,datesListAp},
SetOptions[ListLogPlot,Frame-> True,ImageSize->Medium,Joined->True];

{datesListAm,datesListAp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,2}];{ListLogPlot[datesListAm[[ini;;]],PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA A molecules"},PlotLabel->labelA],
ListLogPlot[datesListAp[[ini;;]],PlotStyle->
RGBColor[0.2,0.8,1],FrameLabel->{"time","protein A molecules"},PlotLabel->labelAp]}]


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
LaunchKernels[10] 


(* ::Text:: *)
(*This function runs the GA algorithm and processes the output according to the previous functions.*)


(* ::Input::Initialization:: *)
DistributeDefinitions[changesByRXN, pF, ff, cv,cv2, steadyStateStart,
     iteraciones1, kgA, \[Gamma]gA, kmA, \[Gamma]mA, kpA, \[Gamma]pA]; 
(*List by parameters*)
(* listEstimatesPara= {};*)
   listEstimatesReplicates = {};

Do[
(*List by replicates*)
    (*listEstimatesReplicates = {};*)

SetSharedVariable[   listEstimatesReplicates,listEstimatemRNAA,listEstimateproteinA,listDataAndTimeT];

ParallelDo[
       
 listDataAndTime = CreateDataStructure["OrderedHashSet"];
        t = 0;
        abundancia = {1, 1, 1, 0}(*mRNAA,proteinA,da0A,da1A*);
        
Do[\[Tau]is=Map[If[# > 0,RandomVariate[ExponentialDistribution[#]],\[Infinity]]&  ,pF[kgA, \[Gamma]gA, kmA, \[Gamma]mA, kpA, \[Gamma]pA, Sequence@@ abundancia]];
            reactionAnd\[Tau]\[Mu] = {FirstPosition[#, Min[#]], Min[#]}&[\[Tau]is];
            abundancia = Through[{reactionAnd\[Tau]\[Mu][[1]] /. changesByRXN}[[1, 1]] @@ abundancia][[1]];
            t = t + reactionAnd\[Tau]\[Mu][[2]];
            listDataAndTime["Union", {{t, abundancia}}];
            If[Length[#]>3&&#[[2]]==5&[IntegerDigits[Round[t]]](* para que Break a los 1500*), Break[]],iteraciones1 ];
      
  listData = Normal[listDataAndTime][[All, 2]];
        ini =Round[(Length[listData]*30)/100.] (*el estado estable se considera depu\[EAcute]s del 30% de iteraciones *);
(*ini = steadyStateStart[window,listData];*)
        {mRNAffcv2ItA, proteinffcv2ItA} = Table[{ff[i], cv[i],cv2[i], Mean[i],StandardDeviation[i],Variance[i],Table[Moment[i,j],{j,10}],Table[CentralMoment[i,j],{j,10}]}, {i, {listData[[ini ;; All, 1]], listData[[ini ;; All, 2]]}}];
       
   AppendTo[listEstimatesReplicates, {mRNAffcv2ItA,proteinffcv2ItA}];
       
      ,
        replicas1

    ];
(*AppendTo[listEstimatesPara,listEstimatesReplicates];*)
    ,
 {kgA,kgs},{\[Gamma]gA,\[Gamma]gs}]; (*CAMBIAR AQUI EN EL .m*)


     Export["controlUnRGA_"<>parametro<>ToString[name]<>".csv",listEstimatesReplicates];
