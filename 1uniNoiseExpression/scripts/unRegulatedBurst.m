(* ::Package:: *)

(* ::Text:: *)
(*In this notebook is simulated a One UnR gene with GA and the three-stage model (burst expression). The notebook has two parts:*)
(*1. To evaluate the "sensitivity" to the whole set the parameters. With an arbitrary start point to calculate the estimates (FF, CV2,Mean) (30% of the total number of iterations), and different number of iterations to get a final t of 1500. *)
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

pF=Compile[{{kg,_Real},{\[Gamma]g,_Real},{km,_Real},{\[Gamma]m,_Real},{kp,_Real},{\[Gamma]p,_Real},{a,_Integer},{p,_Integer},{da0,_Integer},{da1,_Integer}},{da0*kg,da1*\[Gamma]g,da1*km+0.01(*Esto para mantener una expresi\[OAcute]n basal como en el de Langevin (que no tiene ceros)y que el kon no se haga cero*),a*\[Gamma]m,kp*a,\[Gamma]p*p},RuntimeAttributes->{Listable},Parallelization->True];


(* ::Text:: *)
(*This is equal for all  sub-circuits*)


(* ::Input::Initialization:: *)
ff[dato_]:=N[Variance[dato]/Mean[dato]];(* Fano Factor for noise measure *)
expressionNoise[data_]:=N[Variance[data]/Mean[data]^2](* Squared Coefficient of Variation for noise measure *)


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
iteraciones2=400000; (* Number of algorithm iterations *)
replicas2=5; (* Number of times the algorithm is run *)


(* ::Subsection:: *)
(*Simulation*)


(* ::Input::Initialization:: *)
(* Launch the number of kernels in which you will parallelize the simulation *)
LaunchKernels[44] 


(* ::Text:: *)
(*This function runs the GA algorithm and processes the output according to the previous functions.*)


(* ::Input::Initialization:: *)
DistributeDefinitions[changesByRXN, pF, ff, expressionNoise, steadyStateStart,
     iteraciones1, kgA, \[Gamma]gA, kmA, \[Gamma]mA, kpA, \[Gamma]pA]; 
(*List by parameters*)
    mRNAffcv2ListA = {}; 
    proteinffcv2ListA = {}; 

Do[
(*List by replicates*)
    listEstimatemRNAA = {};
    listEstimateproteinA = {};
    listDataAndTimeT={};

SetSharedVariable[listEstimatemRNAA,listEstimateproteinA,listDataAndTimeT];

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
        {mRNAffcv2ItA, proteinffcv2ItA} = Table[{ff[Flatten @ i], expressionNoise[Flatten @ i], Mean[Flatten @ i]}, {i, {listData[[ini ;; All, 1]], listData[[ini ;; All, 2]]}}];
       
   AppendTo[listEstimatemRNAA, mRNAffcv2ItA];
        AppendTo[listEstimateproteinA, proteinffcv2ItA];
        AppendTo[listDataAndTimeT, listDataAndTime];
        
        ,
        replicas1
    ];
  
  AppendTo[mRNAffcv2ListA, Mean[listEstimatemRNAA]];
    AppendTo[proteinffcv2ListA, Mean[listEstimateproteinA]];
    
    ,
 {kgA,{1}},{\[Gamma]gA,{0.01}}]; (*CAMBIAR AQUI EN EL .m*)


(* ::Input::Initialization:: *)
Export["controlUnRGA_"<>parametro<>ToString[name]<>".csv",{mRNAffcv2ListA,proteinffcv2ListA}];


(* ::Input:: *)
(*listDataAndTimeTN=Normal/@listDataAndTimeT;*)
(*MapThread[plotsLogExpression[#1,#2,listDataAndTimeTN,"","",""]&,{Range[1,replicas1],{39639,36610,40648,43881}}]*)


(* ::Section:: *)
(*S2: Evaluation by regions*)


(* ::Subsection:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
iteraciones2=400000; (* Number of algorithm iterations *)
replicas2=5; (* Number of times the algorithm is run *)
rep=1; (* Some plots are made with the output of this replicate*)
lagStep=10; (* Lag time for the Autocorrelation plot *)
species=2;(*mRNA genA protein genA*) 


(* ::Subsection:: *)
(*Functions*)


(* ::Text:: *)
(*These functions make different plots of the temporal genes expression *)


(* ::Input::Initialization:: *)
(* This plots the temporal dynamic of genes mRNA and proteins *)plotsExpression[rep_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,parametros_]:=Module[{datesListAm,datesListAp},
SetOptions[ListLinePlot,Frame-> True,ImageSize->Medium];

{datesListAm,datesListAp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,2}];{ListLinePlot[datesListAm,PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA B molecules"},PlotLabel->labelA],
ListLinePlot[datesListAp,PlotStyle->
RGBColor[0.2,0.8,1],FrameLabel->{"time","protein B molecules"},PlotLabel->labelAp]}]


(* ::Input::Initialization:: *)
(* This plots the temporal dynamic of genes mRNA and proteins with log scale *)plotsLogExpression[rep_,ini_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,parametros_]:=Module[{datesListAm,datesListAp},
SetOptions[ListLogPlot,Frame-> True,ImageSize->Medium,Joined->True];

{datesListAm,datesListAp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,2}];{ListLogPlot[datesListAm[[ini;;]],PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA B molecules"},PlotLabel->labelA],
ListLogPlot[datesListAp[[ini;;]],PlotStyle->
RGBColor[0.2,0.8,1],FrameLabel->{"time","protein B molecules"},PlotLabel->labelAp]}]


(* ::Input::Initialization:: *)
(* This plots the auto-correlation of the temporal dynamic of genes mRNA and proteins *)autoCorrPlot[lagStep_,listDataAndTimeTN_,positions_]:=Module[{data,lenght,data2,corr,negatives,cortes},
data=Table[MapThread[Extract[#1,#2]&,{listDataAndTimeTN[[All,All,2]][[All,All,i]],positions}],{i,{1,2}}];
lenght=Min[Length/@#]&/@data;
data2=MapThread[TemporalData[#1[[All,1;;#2]],{Range[1,#2]}]&,{data,lenght}](*This is to get the corr in each lag time as the average between replicates *);
DistributeDefinitions[data,data2];
corr=ParallelMap[CorrelationFunction[#,{1,#["PathLengths"][[1]]-1,lagStep}]&,data2];
negatives=FirstPosition[Normal[#],{_,_?Negative}]*lagStep&/@corr;
cortes=Count[Partition[#//Normal,2,1],{{_,_?Positive},{_,_?Negative}}|{{_,_?Negative},{_,_?Positive}}]&/@corr;
ListPlot[#[[1]],Filling->Axis ,Frame->True,ImageSize->Medium,FrameLabel->{"Lag-"<>ToString[lagStep],#[[2]]},PlotStyle->#[[4]],PlotLabel->"First negative: "<>ToString[#[[3]]]<>"\n cortes: "<>ToString[#[[5]]]]&/@{{corr[[1]],"mRNA A",negatives[[1]],RGBColor[0.2,0.6,1],cortes[[1]]},{corr[[2]],"protein A",negatives[[2]],RGBColor[0.2,0.8,1],cortes[[2]]}}]


(* ::Input::Initialization:: *)
(* This plots the histogram of the steady state gene expression for mRNA and proteins *)distriPlot[ini_,listDataAndTimeTN_]:=Module[{distribucionAm,distribucionAp},
{distribucionAm,distribucionAp}=Table[Flatten[listDataAndTimeTN[[All,All,2]][[All,ini;;All,i]]],{i,2}];Histogram[#[[1]],Automatic(*{Min[#[[1]]]//Round,Max[#[[1]]]//Round,5}*),"Probability",FrameLabel->{#[[2]],"Probability"},PlotLabel->#[[4]]<>" \n FF: "<>ToString[ff[#[[1]]]]<>", \!\(\*SuperscriptBox[\(CV\), \(2\)]\): "<>ToString[expressionNoise[#[[1]]]]<>", Mean: "<>ToString[Mean[#[[1]]]//N],ChartStyle->#[[3]],Frame->True,ImageSize->Medium]&/@{{distribucionAm,"mRNA B molecules",RGBColor[0.2,0.6,1],labelA},{distribucionAp,"protein B molecules",RGBColor[0.2,0.8,1],labelAp}}]


(* ::Input::Initialization:: *)
(* This plots the velocity of the temporal dynamic of  gene expression *)veloPlot[rep_,listDataAndTimeTN_,positionsVelocity_]:=Module[{listS,listT,velocidades},
listS=Table[Extract[#[[All,i]],positionsVelocity],{i,2}]&[listDataAndTimeTN[[rep,All,2]]];
listT=Extract[#,positionsVelocity]&[listDataAndTimeTN[[rep,All,1]]];
velocidades=Table[(listS[[i,2;;]]-listS[[i,;;-2]])/(listT[[2;;]]-listT[[;;-2]]),{i,2}];
Table[ListLinePlot[Transpose[{listT[[2;;]],velocidades[[i[[1]]]]}],PlotStyle->i[[2]],PlotLabel->"Velocities",FrameLabel->{"time",i[[3]]}],{i,{{1,RGBColor[0.2,0.6,1],"mRNA A molecules"},{2,RGBColor[0.2,0.8,1],"protein A molecules"}}}]]


(* ::Text:: *)
(*Here there is a detail that changes by sub-circuit : the number and type of Kmm (put it in the arguments of the function simu, in the DistributeDefinitions, and pF function)*)


(* ::Text:: *)
(*This function runs the GA algorithm and processes the output according to the previous functions.*)


(* ::Input::Initialization:: *)
simu[iteraciones_,replicas_,kgA_,\[Gamma]gA_,kmA_,\[Gamma]mA_,kpA_,\[Gamma]pA_,rep_,parametros_,lagStep_,species_]:=

Module[{ini,listEstimatemRNAA,listEstimateproteinA,listEstimatemRNAI,listEstimateproteinI,listDataAndTimeT,listDataAndTime,t,abundancia,\[Tau]is,reactionAnd\[Tau]\[Mu],listData,mRNAffcv2ItA,proteinffcv2ItA,mRNAffcv2ItI,proteinffcv2ItI,mRNAffcv2MeanA,proteinffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanI,listDataAndTimeTN,positions},

DistributeDefinitions[changesByRXN,pF,ff,expressionNoise,steadyStateStart,iteraciones,ini,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA];

listEstimatemRNAA={};
listEstimateproteinA={};
listDataAndTimeT={};

SetSharedVariable[listEstimatemRNAA,listEstimateproteinA,listDataAndTimeT];

ParallelDo[

listDataAndTime=CreateDataStructure["OrderedHashSet"];
t=0;
abundancia={1,1,1,0}(*mRNAA,proteinA,da0A,da1A*);


Do[
          \[Tau]is=Map[If[#>0,RandomVariate[ExponentialDistribution[#]],\[Infinity]]&,pF[kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,Sequence@@abundancia]];reactionAnd\[Tau]\[Mu]={FirstPosition[#,Min[#]],Min[#]}&[\[Tau]is];abundancia=Through[{reactionAnd\[Tau]\[Mu][[1]]/.changesByRXN}[[1,1]]@@abundancia][[1]];t=t+reactionAnd\[Tau]\[Mu][[2]];listDataAndTime["Union",{{t,abundancia}}];
If[Length[#]>3&&#[[2]]==5&[IntegerDigits[Round[t]]](* para que Break a los 1500*), Break[]],iteraciones];

listData=Normal[listDataAndTime][[All,2]];
 ini =Round[(Length[listData]*30)/100.];

{mRNAffcv2ItA,proteinffcv2ItA}=Table[{ff[Flatten@i],expressionNoise[Flatten@i],Mean[Flatten@i]//N},{i,{listData[[ini;;All,1]],
listData[[ini;;All,2]]}}];

AppendTo[listEstimatemRNAA,mRNAffcv2ItA];
AppendTo[listEstimateproteinA,proteinffcv2ItA];
AppendTo[listDataAndTimeT,listDataAndTime];,

replicas];

mRNAffcv2MeanA=Mean[listEstimatemRNAA];
proteinffcv2MeanA=Mean[listEstimateproteinA];

Export["controlUnRGA_"<>parametros<>".csv",{listEstimatemRNAA,listEstimateproteinA}];

listDataAndTimeTN=Normal/@listDataAndTimeT;

{labelA,labelAp}=Table["FF: "<>ToString[i[[1]]]<>" CV2: "<>ToString[i[[2]]]<>" Mean: "<>ToString[i[[3]]]<>parametros,{i,{mRNAffcv2MeanA,proteinffcv2MeanA}}];

 ini =Round[(Length[listDataAndTimeTN[[rep]]]*30)/100.];

positions=Table[DeleteCases[FirstPosition[listDataAndTimeTN[[i,All,1]]//Round,#]&/@Range[1,1500,1],Missing["NotFound"]],{i,replicas}];(* These are the positions the ~1500 points to get the velocities and the auto-correlation plot*)

{plotsExpression[rep,listDataAndTimeTN,mRNAffcv2MeanA,proteinffcv2MeanA,parametros],
plotsLogExpression[rep,ini,listDataAndTimeTN,mRNAffcv2MeanA,proteinffcv2MeanA,parametros],
autoCorrPlot[lagStep,listDataAndTimeTN,positions],
distriPlot[ini,listDataAndTimeTN],
veloPlot[rep,listDataAndTimeTN,positions[[rep]]],listDataAndTimeT}]


(* ::Subsection:: *)
(*Simulation*)


(* ::Input:: *)
(*s=Flatten[Table[simu[iteraciones2,replicas2,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,rep,"\n kgA: " <>ToString[kgA]<>" \[Gamma]gA: "<>ToString[\[Gamma]gA](*parametros*),lagStep,species],{kgA,{0.05,1}},{\[Gamma]gA,{0.2,2}}],{3}];*)


(* ::Input:: *)
(*s*)
