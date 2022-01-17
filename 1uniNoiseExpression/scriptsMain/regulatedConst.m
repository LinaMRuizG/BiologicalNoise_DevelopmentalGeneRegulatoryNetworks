(* ::Package:: *)

(* ::Text:: *)
(*In this nb is simulated One self-regulated gene A with  GA and constitutive expression. The idea is to evaluate 4 regions of parameters  *)


(* ::Section:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
kmA=0.696;(* mRNA Production rate (molecules*t^-1) *)
\[Gamma]mA=0.02082 ;(* mRNA Degradation rate (molecules*t^-1) *)
kpA=1.386;(* Proteins production rate (molecules*t^-1) *)
\[Gamma]pA=0.02082;(* Proteins degradation rate (molecules*t^-1) *)
h=3.0;(* Hill's constant *)
kAA=100; (* Michaeles-Menten Constant *)


(* ::Input::Initialization:: *)
iteraciones2=400000; (* Number of algorithm iterations *)
replicas2=4;  (* Number of times the algorithm is run *)
rep=1;  (* Some plots are made with the output of this replicate*)
lagStep=10;  (* Lag time for the Autocorrelation plot *)
species=4;(*mRNAA, proteinA*)


(* ::Section:: *)
(*Functions*)


(* ::Text:: *)
(*These are the reactions (changesByRXN) and propensity functions (pF) of the GA for the specific regulation system: one self-regulated gene*)


(* ::Input::Initialization:: *)
changesByRXN=Function/@<|1->{#1+1,#2},2->{#1-1,#2},3->{#1,#2+1},4->{#1,#2-1}|>;(* Each hash represents a specie: mRNA,protein,da0:promoterInactivated,da1: promoterActivated. The numbers are the next reactions: 1-activation, 2-deactivation, 3-5 synthesis, 4-6 degradation *)pF=Compile[{{h,_Real},{kAA,_Real},{km,_Real},{\[Gamma]m,_Real},{kp,_Real},{\[Gamma]p,_Real},{a,_Integer},{p,_Integer}},{km*(p ^h/(kAA+p ^h))+0.01,a*\[Gamma]m,kp*a,\[Gamma]p*p},RuntimeAttributes->{Listable},Parallelization->True];


(* ::Text:: *)
(*These functions are equal for all regulation systems and they are used to estimate the FF and CV2 in the steady state*)


(* ::Input::Initialization:: *)
ff[dato_]:=N[Variance[dato]/Mean[dato]];(* Fano Factor for noise measure *)
expressionNoise[data_]:=N[Variance[data]/Mean[data]^2](* Squared Coefficient of Variation for noise measure *)


(* ::Text:: *)
(*These functions make different plots of the temporal genes expression *)


(* ::Input::Initialization:: *)
(* This plots the temporal dynamic of genes mRNA and proteins *)plotsExpression[rep_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,parametros_]:=Module[{datesListAm,datesListAp},
SetOptions[ListLinePlot,Frame-> True,ImageSize->Medium];

{datesListAm,datesListAp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,2}];{ListLinePlot[datesListAm,PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA A molecules"},PlotLabel->labelA],
ListLinePlot[datesListAp,PlotStyle->
RGBColor[0.2,.8,1],FrameLabel->{"time","protein A molecules"},PlotLabel->labelAp]}]


(* ::Input::Initialization:: *)
(* This plots the temporal dynamic of genes mRNA and proteins with log scale *)plotsLogExpression[rep_,ini_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,parametros_]:=Module[{datesListAm,datesListAp},
SetOptions[ListLogPlot,Frame-> True,ImageSize->Medium,Joined->True];

{datesListAm,datesListAp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,2}];{ListLogPlot[datesListAm[[ini;;]],PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA A molecules"},PlotLabel->labelA],
ListLogPlot[datesListAp[[ini;;]],PlotStyle->
RGBColor[0.2,.8,1],FrameLabel->{"time","protein A molecules"},PlotLabel->labelAp]}]


(* ::Input::Initialization:: *)
(* This plots the auto-correlation of the temporal dynamic of genes mRNA and proteins *)autoCorrPlot[lagStep_,listDataAndTimeTN_,positions_]:=Module[{data,lenght,data2,corr,negatives,cortes},
data=Table[MapThread[Extract[#1,#2]&,{listDataAndTimeTN[[All,All,2]][[All,All,i]],positions}],{i,{1,2}}];
lenght=Min[Length/@#]&/@data;
data2=MapThread[TemporalData[#1[[All,1;;#2]],{Range[1,#2]}]&,{data,lenght}](*This is to get the corr in each lag time as the average between replicates *);
DistributeDefinitions[data,data2];
corr=ParallelMap[CorrelationFunction[#,{1,#["PathLengths"][[1]]-1,lagStep}]&,data2];
negatives=FirstPosition[Normal[#],{_,_?Negative}]*lagStep&/@corr;
cortes=Count[Partition[#//Normal,2,1],{{_,_?Positive},{_,_?Negative}}|{{_,_?Negative},{_,_?Positive}}]&/@corr;
ListPlot[#[[1]],Filling->Axis ,Frame->True,ImageSize->Medium,FrameLabel->{"Lag-"<>ToString[lagStep],#[[2]]},PlotStyle->#[[4]],PlotLabel->"First negative: "<>ToString[#[[3]]]<>"\n cortes: "<>ToString[#[[5]]]]&/@{{corr[[1]],"mRNA A",negatives[[1]],RGBColor[0.2,0.6,1],cortes[[1]]},{corr[[2]],"protein A",negatives[[2]],RGBColor[0.5,0.99,1],cortes[[2]]}}]


(* ::Input::Initialization:: *)
(* This plots the histogram of the steady state gene expression for mRNA and proteins *)
distriPlot[ini_,listDataAndTimeTN_]:=Module[{distribucionAm,distribucionAp},
{distribucionAm,distribucionAp}=Table[Flatten[listDataAndTimeTN[[All,All,2]][[All,ini;;All,i]]],{i,2}];Histogram[#[[1]],Automatic(*{Min[#[[1]]]//Round,Max[#[[1]]]//Round,5}*),"Probability",FrameLabel->{#[[2]],"Probability"},PlotLabel->#[[4]]<>" \n FF: "<>ToString[ff[#[[1]]]]<>", \!\(\*SuperscriptBox[\(CV\), \(2\)]\): "<>ToString[expressionNoise[#[[1]]]]<>", Mean: "<>ToString[Mean[#[[1]]]//N],ChartStyle->#[[3]],Frame->True,ImageSize->Medium]&/@{{distribucionAm,"mRNA A molecules",RGBColor[0.2,0.6,1],labelA},{distribucionAp,"protein A molecules",RGBColor[0.2,.8,1],labelAp}}]


(* ::Input::Initialization:: *)
(* This plots the velocity of the temporal dynamic of  gene expression *)veloPlot[rep_,listDataAndTimeTN_,positionsVelocity_]:=Module[{listS,listT,velocidades},
listS=Table[Extract[#[[All,i]],positionsVelocity],{i,2}]&[listDataAndTimeTN[[rep,All,2]]];
listT=Extract[#,positionsVelocity]&[listDataAndTimeTN[[rep,All,1]]];
velocidades=Table[(listS[[i,2;;]]-listS[[i,;;-2]])/(listT[[2;;]]-listT[[;;-2]]),{i,2}];
Table[ListLinePlot[Transpose[{listT[[2;;]],velocidades[[i[[1]]]]}],PlotStyle->i[[2]],PlotLabel->"Velocities",FrameLabel->{"time",i[[3]]}],{i,{{1,RGBColor[0.2,0.6,1],"mRNA A molecules"},{2,RGBColor[0.2,.8,1],"protein A molecules"}}}]]


(* ::Text:: *)
(*Here there is a detail that changes by sub-circuit : the number and type of Kmm (put it in the arguments of the function simu, in the DistributeDefinitions, and pF function)*)


(* ::Text:: *)
(*This function runs the GA algorithm and processes the output according to the previous functions.*)


(* ::Input::Initialization:: *)
simu[iteraciones_,replicas_,h_,kAA_,kmA_,\[Gamma]mA_,kpA_,\[Gamma]pA_,rep_,parametros_,lagStep_,species_]:=

Module[{listEstimatemRNAA,listEstimateproteinA,listEstimatemRNAI,listEstimateproteinI,listDataAndTimeT,listDataAndTime,t,abundancia,\[Tau]is,reactionAnd\[Tau]\[Mu],listData,mRNAffcv2ItA,proteinffcv2ItA,mRNAffcv2ItI,proteinffcv2ItI,mRNAffcv2MeanA,proteinffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanI,listDataAndTimeTN,ini,positions},

DistributeDefinitions[changesByRXN,pF,ff,expressionNoise,iteraciones,h,kAA,kmA,\[Gamma]mA,kpA,\[Gamma]pA];

listEstimatemRNAA={};
listEstimateproteinA={};
listDataAndTimeT={};

SetSharedVariable[listEstimatemRNAA,listEstimateproteinA,listDataAndTimeT];

ParallelDo[

listDataAndTime=CreateDataStructure["OrderedHashSet"];
t=0;
abundancia={1,1}(*mRNAA,proteinA,da0A,da1A*);


Do[
          \[Tau]is=Map[If[#>0,RandomVariate[ExponentialDistribution[#]],\[Infinity]]&,pF[h,kAA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,Sequence@@abundancia]];reactionAnd\[Tau]\[Mu]={FirstPosition[#,Min[#]],Min[#]}&[\[Tau]is];abundancia=Through[{reactionAnd\[Tau]\[Mu][[1]]/.changesByRXN}[[1,1]]@@abundancia][[1]];t=t+reactionAnd\[Tau]\[Mu][[2]];listDataAndTime["Union",{{t,abundancia}}];
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

Export["controlRGAConst_"<>parametros<>".csv",{listEstimatemRNAA,listEstimateproteinA}];

listDataAndTimeTN=Normal/@listDataAndTimeT;

{labelA,labelAp}=Table["FF: "<>ToString[i[[1]]]<>" CV2: "<>ToString[i[[2]]]<>" Mean: "<>ToString[i[[3]]]<>parametros,{i,{mRNAffcv2MeanA,proteinffcv2MeanA}}];

ini =Round[(Length[listDataAndTimeTN[[rep]]]*30)/100.];

positions=Table[DeleteCases[FirstPosition[listDataAndTimeTN[[i,All,1]]//Round,#]&/@Range[1,1500,1],Missing["NotFound"]],{i,replicas}];(* These are the positions the ~1500 points to get the velocities and the auto-correlation plot*)

{plotsExpression[rep,listDataAndTimeTN,mRNAffcv2MeanA,proteinffcv2MeanA,parametros],
plotsLogExpression[rep,ini,listDataAndTimeTN,mRNAffcv2MeanA,proteinffcv2MeanA,parametros],
autoCorrPlot[lagStep,listDataAndTimeTN,positions],
distriPlot[ini,listDataAndTimeTN],
veloPlot[rep,listDataAndTimeTN,positions[[rep]]]}]


(* ::Subsection:: *)
(*Simulation*)


(* ::Input::Initialization:: *)
s=simu[iteraciones2,replicas2,h,kAA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,rep,"\n Default parameters"(*parametros*),lagStep,species];


(* ::Input::Initialization:: *)
Export["controlRGAConstReg.pdf",s];
