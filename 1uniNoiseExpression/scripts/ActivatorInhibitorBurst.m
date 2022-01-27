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
name=0.5; (* This is the last value of the range for the parameter being evluated *)
kgs=Range[0.01,1,0.01];
\[Gamma]gs=Range[0.02,name,0.02]; (* subRanges: 0.02-0.5, 0.52-1.0, 1.02-1.5, 1.52-2 *)
\[Gamma]ms=Range[0.001,0.1(*name*),0.001];(* subRanges: 0.001-0.025, 0.026-0.05, 0.051-0.075, 0.076-0.1 *)
\[Gamma]ps=Range[0.001,0.1(*name*),0.001];(* subRanges: 0.001-0.025, 0.026-0.05, 0.051-0.075, 0.076-0.1 *)
parametro="\[Gamma]g";(* With the name of the second parameter is named the output file *)



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
ff[dato_]:=N[Variance[dato]/Mean[dato]];(* Fano Factor for noise measure *)
expressionNoise[data_]:=N[Variance[data]/Mean[data]^2](* Squared Coefficient of Variation for noise measure *)


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
DistributeDefinitions[changesByRXN,pF,ff,expressionNoise,iteraciones,h,kIA,kAA,kAI,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kgI,\[Gamma]gI,kmI,\[Gamma]mI,kpI,\[Gamma]pI];
(*List by parameters*)
mRNAffcv2ListA={};
proteinffcv2ListA={};
mRNAffcv2ListI={};
proteinffcv2ListI={};

Do[
(*List by replicates*)
listEstimatemRNAA={};
listEstimateproteinA={};
listEstimatemRNAI={};
listEstimateproteinI={};
 listDataAndTimeT={};

SetSharedVariable[listEstimatemRNAA,listEstimateproteinA,listEstimatemRNAI,listEstimateproteinI, listDataAndTimeT];

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
        ini =Echo@Round[(Length[listData]*30)/100.] (*el estado estable se considera depu\[EAcute]s del 30% de iteraciones *);
(*ini = steadyStateStart[window,listData];*)
{mRNAffcv2ItA,proteinffcv2ItA}=Table[{ff[Flatten@i],expressionNoise[Flatten@i],Mean[Flatten@i]},{i,{listData[[ini;;All,1]],
listData[[ini;;All,2]]}}];
{mRNAffcv2ItI,proteinffcv2ItI}=Table[{ff[Flatten@i],expressionNoise[Flatten@i],Mean[Flatten@i]},{i,{listData[[ini;;All,3]],
listData[[ini;;All,4]]}}];

AppendTo[listEstimatemRNAA,mRNAffcv2ItA];
AppendTo[listEstimateproteinA,proteinffcv2ItA];
AppendTo[listEstimatemRNAI,mRNAffcv2ItI];
AppendTo[listEstimateproteinI,proteinffcv2ItI];
AppendTo[listDataAndTimeT,listDataAndTime];,


replicas1];

AppendTo[mRNAffcv2ListA,Mean[listEstimatemRNAA]];
AppendTo[proteinffcv2ListA,Mean[listEstimateproteinA]];
AppendTo[mRNAffcv2ListI,Mean[listEstimatemRNAI]];
AppendTo[proteinffcv2ListI,Mean[listEstimateproteinI]];,

    
 {kg,{1}},{\[Gamma]g,{0.01}}]; (*CAMBIAR AQUI EN EL .m*)


(* ::Input:: *)
(*(* CAMBIAR AQUI EN EL .m _Cuando cambiamos A *****************)*)


(* ::Input::Initialization:: *)
Export["controlAIGAA_A"<>parametro<>ToString[name],{mRNAffcv2ListA,proteinffcv2ListA},"CSV"];
Export["controlAIGAA_I"<>parametro<>ToString[name],{mRNAffcv2ListI,proteinffcv2ListI},"CSV"];


(* ::Input:: *)
(*(* CAMBIAR AQUI EN EL .m _Cuando cambiamos I *****************)*)


(* ::Input::Initialization:: *)
Export["controlAIGAI_A"<>parametro<>ToString[name],{mRNAffcv2ListA,proteinffcv2ListA},"CSV"];
Export["controlAIGAI_I"<>parametro<>ToString[name],{mRNAffcv2ListI,proteinffcv2ListI},"CSV"];


(* ::Input:: *)
(*listDataAndTimeTN=Normal/@listDataAndTimeT;*)
(*MapThread[plotsLogExpression[#1,#2,listDataAndTimeTN,"","","","",""]&,{Range[1,replicas1],{3864,5708,5625,6087}}]*)


(* ::Section:: *)
(*S2: Evaluation by regions*)


(* ::Subsection:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
iteraciones2=400000; (* Number of algorithm iterations *)
replicas2=500; (* Number of times the algorithm is run *)
rep=1; (* Some plots are made with the output of this replicate*)
lagStep=10; (* Lag time for the Autocorrelation plot *)
species=4;(*mRNA genA/B protein genA/B*) 


(* ::Subsection:: *)
(*Functions*)


(* ::Text:: *)
(*These functions make different plots of the temporal genes expression *)


(* ::Input::Initialization:: *)
(* This plots the temporal dynamic of genes mRNA and proteins *)plotsExpression[rep_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,mRNAffcv2MeanI_,proteinffcv2MeanI_,parametros_]:=Module[{datesListAm,datesListAp,datesListBm,datesListBp},
SetOptions[ListLinePlot,Frame-> True,ImageSize->Medium];

{datesListAm,datesListAp,datesListBm,datesListBp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,4}];{ListLinePlot[datesListAm,PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA A molecules"},PlotLabel->labelA],
ListLinePlot[datesListBm,PlotStyle->
RGBColor[0.5,0.99,1],FrameLabel->{"time","mRNA I molecules"},PlotLabel->labelB],
ListLinePlot[datesListAp,PlotStyle->
RGBColor[0.2,.8,1],FrameLabel->{"time","protein A molecules"},PlotLabel->labelAp],
ListLinePlot[datesListBp,PlotStyle->
RGBColor[0.3,.09,0.7],FrameLabel->{"time","protein I molecules"},PlotLabel->labelBp]}]


(* ::Input::Initialization:: *)
(* This plots the temporal dynamic of genes mRNA and proteins with log scale *)plotsLogExpression[rep_,ini_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,mRNAffcv2MeanI_,proteinffcv2MeanI_,parametros_]:=Module[{datesListAm,datesListAp,datesListBm,datesListBp},
SetOptions[ListLogPlot,Frame-> True,ImageSize->Medium,Joined->True];

{datesListAm,datesListAp,datesListBm,datesListBp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,4}];{ListLogPlot[datesListAm[[ini;;]],PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA A molecules"},PlotLabel->labelA],
ListLogPlot[datesListBm[[ini;;]],PlotStyle->
RGBColor[0.5,0.99,1],FrameLabel->{"time","mRNA I molecules"},PlotLabel->labelB],
ListLogPlot[datesListAp[[ini;;]],PlotStyle->
RGBColor[0.2,.8,1],FrameLabel->{"time","protein A molecules"},PlotLabel->labelAp],
ListLogPlot[datesListBp[[ini;;]],PlotStyle->
RGBColor[0.3,.09,0.7],FrameLabel->{"time","protein I molecules"},PlotLabel->labelBp]}]


(* ::Input::Initialization:: *)
(* This plots the auto-correlation of the temporal dynamic of genes mRNA and proteins *)
autoCorrPlot[lagStep_,listDataAndTimeTN_,positions_]:=Module[{data,lenght,data2,corr,negatives,cortes},
data=Table[MapThread[Extract[#1,#2]&,{listDataAndTimeTN[[All,All,2]][[All,All,i]],positions}],{i,{1,3,2,4}}];
lenght=Min[Length/@#]&/@data;
data2=MapThread[TemporalData[#1[[All,1;;#2]],{Range[1,#2]}]&,{data,lenght}](*This is to get the corr in each lag time as the average between replicates *);
DistributeDefinitions[data,data2];
corr=ParallelMap[CorrelationFunction[#,{1,#["PathLengths"][[1]]-1,lagStep}]&,data2];
negatives=FirstPosition[Normal[#],{_,_?Negative}]*lagStep&/@corr;
cortes=Count[Partition[#//Normal,2,1],{{_,_?Positive},{_,_?Negative}}|{{_,_?Negative},{_,_?Positive}}]&/@corr;
ListPlot[#[[1]],Filling->Axis ,Frame->True,ImageSize->Medium,FrameLabel->{"Lag-"<>ToString[lagStep],#[[2]]},PlotStyle->#[[4]],PlotLabel->"First negative: "<>ToString[#[[3]]]<>"\n cortes: "<>ToString[#[[5]]]]&/@{{corr[[1]],"mRNA A",negatives[[1]],RGBColor[0.2,0.6,1],cortes[[1]]},{corr[[2]],"mRNA I",negatives[[2]],RGBColor[0.5,0.99,1],cortes[[2]]},{corr[[3]],"protein A",negatives[[3]],RGBColor[0.2,.8,1],cortes[[3]]},{corr[[4]],"protein I",negatives[[4]],RGBColor[0.3,.09,0.7],cortes[[4]]}}]


(* ::Input::Initialization:: *)
(* This plots the histogram of the steady state gene expression for mRNA and proteins *)
distriPlot[ini_,listDataAndTimeTN_]:=Module[{distribucionAm,distribucionAp,distribucionBm,distribucionBp},
{distribucionAm,distribucionAp,distribucionBm,distribucionBp}=Table[Flatten[listDataAndTimeTN[[All,All,2]][[All,ini;;All,i]]],{i,4}];Histogram[#[[1]],Automatic(*{Min[#[[1]]]//Round,Max[#[[1]]]//Round,5}*),"Probability",FrameLabel->{#[[2]],"Probability"},PlotLabel->#[[4]]<>" \n FF: "<>ToString[ff[#[[1]]]]<>", \!\(\*SuperscriptBox[\(CV\), \(2\)]\): "<>ToString[expressionNoise[#[[1]]]]<>", Mean: "<>ToString[Mean[#[[1]]]//N],ChartStyle->#[[3]],Frame->True,ImageSize->Medium]&/@{{distribucionAm,"mRNA A molecules",RGBColor[0.2,0.6,1],labelA},{distribucionBm,"mRNA I molecules",RGBColor[0.5,0.99,1],labelB},{distribucionAp,"protein A molecules",RGBColor[0.2,.8,1],labelAp},{distribucionBp,"protein I molecules ",RGBColor[0.3,.09,0.7],labelBp}}]


(* ::Input::Initialization:: *)
(* This plots the temporal dynamics of the regulator (proteins) and the regulated gene (proteins and mRNA) and indicates the Pearson Correlation between temporal dynamics of regulator and regulated gene *)
corrPlotAI[replicas_,rep_,listDataAndTimeTN_]:=Module[{dataBmRNA,dataAp,dataBp,corremRNA,correprotein},
{dataBmRNA,dataAp,dataBp}=Table[#/Max[#]&[listDataAndTimeTN[[All,All,2]][[i,All,j]]],{j,{3,2,4}},{i,replicas}];{corremRNA,correprotein}=Mean[Table[Correlation[#[[1]][[i]],#[[2]][[i]]]//N,{i,replicas}]]&/@{{dataAp,dataBmRNA},{dataAp,dataBp}};{ListLinePlot[{Transpose[{listDataAndTimeTN[[rep,All,1]],dataAp[[rep]]}],Transpose[{listDataAndTimeTN[[rep,All,1]],dataBmRNA[[rep]]}]},PlotStyle->{RGBColor[0.2,.8,1],RGBColor[0.5,0.99,1]},Frame->True,ImageSize->Medium,PlotLabel->"Pearson Correlation: "<>ToString[corremRNA],PlotLegends->{"protein A","mRNA I"}],ListLinePlot[{Transpose[{listDataAndTimeTN[[rep,All,1]],dataAp[[rep]]}],Transpose[{listDataAndTimeTN[[rep,All,1]],dataBp[[rep]]}]},PlotStyle->{RGBColor[0.2,.8,1],RGBColor[0.3,.09,0.7]},Frame->True,ImageSize->Medium,PlotLabel->"Pearson Correlation: "<>ToString[correprotein],PlotLegends->{"protein A","protein I"}]}]


(* ::Input::Initialization:: *)
(* This plots the temporal dynamics of the regulator (proteins) and the regulated gene (proteins and mRNA) and indicates the Pearson Correlation between temporal dynamics of regulator and regulated gene *)corrPlotIA[replicas_,rep_,listDataAndTimeTN_]:=Module[{dataAmRNA,dataAp,dataBp,corremRNA,correprotein},
{dataBp,dataAmRNA,dataAp}=Table[#/Max[#]&[listDataAndTimeTN[[All,All,2]][[i,All,j]]],{j,{4,1,2}},{i,replicas}];{corremRNA,correprotein}=Mean[Table[Correlation[#[[1]][[i]],#[[2]][[i]]]//N,{i,replicas}]]&/@{{dataBp,dataAmRNA},{dataBp,dataAp}};{ListLinePlot[{Transpose[{listDataAndTimeTN[[rep,All,1]],dataBp[[rep]]}],Transpose[{listDataAndTimeTN[[rep,All,1]],dataAmRNA[[rep]]}]},PlotStyle->{RGBColor[0.3,.09,0.7],RGBColor[0.2,0.6,1]},Frame->True,ImageSize->Medium,PlotLabel->"Pearson Correlation: "<>ToString[corremRNA],PlotLegends->{"protein I","mRNA A"}],ListLinePlot[{Transpose[{listDataAndTimeTN[[rep,All,1]],dataBp[[rep]]}],Transpose[{listDataAndTimeTN[[rep,All,1]],dataAp[[rep]]}]},PlotStyle->{RGBColor[0.3,.09,0.7],RGBColor[0.2,.8,1]},Frame->True,ImageSize->Medium,PlotLabel->"Pearson Correlation: "<>ToString[correprotein],PlotLegends->{"protein I","protein A"}]}]


(* ::Input::Initialization:: *)
(* This plots the velocity of the temporal dynamic of  gene expression *)veloPlot[rep_,listDataAndTimeTN_,positionsVelocity_]:=Module[{listS,listT,velocidades},
listS=Table[Extract[#[[All,i]],positionsVelocity],{i,4}]&[listDataAndTimeTN[[rep,All,2]]];
listT=Extract[#,positionsVelocity]&[listDataAndTimeTN[[rep,All,1]]];
velocidades=Table[(listS[[i,2;;]]-listS[[i,;;-2]])/(listT[[2;;]]-listT[[;;-2]]),{i,4}];
Table[ListLinePlot[Transpose[{listT[[2;;]],velocidades[[i[[1]]]]}],PlotStyle->i[[2]],PlotLabel->"Velocities",FrameLabel->{"time",i[[3]]}],{i,{{1,RGBColor[0.2,0.6,1],"mRNA A molecules"},{2,RGBColor[0.2,.8,1],"protein A molecules"},{3,RGBColor[0.5,0.99,1],"mRNA I molecules"},{4,RGBColor[0.3,.09,0.7],"protein I molecules"}}}]]


(* ::Text:: *)
(*This function runs the GA algorithm and processes the output according to the previous functions.*)


(* ::Input::Initialization:: *)
simu[iteraciones_,replicas_,h_,kIA_,kAA_,kAI_,kgA_,\[Gamma]gA_,kmA_,\[Gamma]mA_,kpA_,\[Gamma]pA_,kgI_,\[Gamma]gI_,kmI_,\[Gamma]mI_,kpI_,\[Gamma]pI_,rep_,parametros_,lagStep_,species_]:=

Module[{ini,listEstimatemRNAA,listEstimateproteinA,listEstimatemRNAI,listEstimateproteinI,listDataAndTimeT,listDataAndTime,t,abundancia,\[Tau]is,reactionAnd\[Tau]\[Mu],listData,mRNAffcv2ItA,proteinffcv2ItA,mRNAffcv2ItI,proteinffcv2ItI,mRNAffcv2MeanA,proteinffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanI,listDataAndTimeTN,positions},

DistributeDefinitions[changesByRXN,pF,ff,expressionNoise,iteraciones,h,kIA,kAA,kAI,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kgI,\[Gamma]gI,kmI,\[Gamma]mI,kpI,\[Gamma]pI];

listEstimatemRNAA={};
listEstimateproteinA={};
listEstimatemRNAI={};
listEstimateproteinI={};
listDataAndTimeT={};

SetSharedVariable[listEstimatemRNAA,listEstimateproteinA,listEstimatemRNAI,listEstimateproteinI,listDataAndTimeT];

ParallelDo[

listDataAndTime=CreateDataStructure["OrderedHashSet"];
t=0;
abundancia={1,1,1,1,1,0,1,0}(*mRNAA,proteinA,mRNAI,proteinI,da0A,da1A,db0I,db1I*);


Do[
          \[Tau]is=Map[If[#>0,RandomVariate[ExponentialDistribution[#]],\[Infinity]]&,pF[h,kIA,kAA,kAI,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kgI,\[Gamma]gI,kmI,\[Gamma]mI,kpI,\[Gamma]pI,Sequence@@abundancia]];reactionAnd\[Tau]\[Mu]={FirstPosition[#,Min[#]],Min[#]}&[\[Tau]is];abundancia=Through[{reactionAnd\[Tau]\[Mu][[1]]/.changesByRXN}[[1,1]]@@abundancia][[1]];t=t+reactionAnd\[Tau]\[Mu][[2]];listDataAndTime["Union",{{t,abundancia}}];
If[Length[#]>3&&#[[2]]==5&[IntegerDigits[Round[t]]](* para que Break a los 1500*), Break[]],iteraciones];

listData=Normal[listDataAndTime][[All,2]];
ini =Round[(Length[listData]*30)/100.] (*el estado estable se considera depu\[EAcute]s del 30% de iteraciones *);
(*ini = steadyStateStart[window,listData];*)

{mRNAffcv2ItA,proteinffcv2ItA}=Table[{ff[Flatten@i],expressionNoise[Flatten@i],Mean[Flatten@i]//N},{i,{listData[[ini;;All,1]],
listData[[ini;;All,2]]}}];
{mRNAffcv2ItI,proteinffcv2ItI}=Table[{ff[Flatten@i],expressionNoise[Flatten@i],Mean[Flatten@i]//N},{i,{listData[[ini;;All,3]],
listData[[ini;;All,4]]}}];

AppendTo[listEstimatemRNAA,mRNAffcv2ItA];
AppendTo[listEstimateproteinA,proteinffcv2ItA];
AppendTo[listEstimatemRNAI,mRNAffcv2ItI];
AppendTo[listEstimateproteinI,proteinffcv2ItI];
AppendTo[listDataAndTimeT,listDataAndTime];,

replicas];


mRNAffcv2MeanA=Mean[listEstimatemRNAA];
proteinffcv2MeanA=Mean[listEstimateproteinA];
mRNAffcv2MeanI=Mean[listEstimatemRNAI];
proteinffcv2MeanI=Mean[listEstimateproteinI];

Export["controlAIGA_A"<>parametros<>".csv",{listEstimatemRNAA,listEstimateproteinA}];
Export["controlAIGA_I"<>parametros<>".csv",{listEstimatemRNAA,listEstimateproteinA}];


listDataAndTimeTN=Normal/@listDataAndTimeT;

{labelA,labelB,labelAp,labelBp}=Table["FF: "<>ToString[i[[1]]]<>" CV2: "<>ToString[i[[2]]]<>" Mean: "<>ToString[i[[3]]]<>parametros,{i,{mRNAffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanA,proteinffcv2MeanI}}];

ini =Round[(Length[listDataAndTimeTN[[rep]]]*30)/100.] ;

positions=Table[DeleteCases[FirstPosition[listDataAndTimeTN[[i,All,1]]//Round,#]&/@Range[1,1500,1],Missing["NotFound"]],{i,replicas}];(* These are the positions the ~1500 points to get the velocities and the auto-correlation plot*)

{plotsExpression[rep,listDataAndTimeTN,mRNAffcv2MeanA,proteinffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanI,parametros],plotsLogExpression[rep,ini,listDataAndTimeTN,mRNAffcv2MeanA,proteinffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanI,parametros],
autoCorrPlot[lagStep,listDataAndTimeTN,positions],
distriPlot[ini,listDataAndTimeTN],
corrPlotAI[replicas,rep,listDataAndTimeTN],
corrPlotIA[replicas,rep,listDataAndTimeTN],
veloPlot[rep,listDataAndTimeTN,positions[[rep]]],
listDataAndTimeT}]


(* ::Subsection:: *)
(*Simulation*)


(* ::Input:: *)
(*s=Flatten[Table[simu[iteraciones2,replicas2,h,kIA,kAA,kAI,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kgI,\[Gamma]gI,kmI,\[Gamma]mI,kpI,\[Gamma]pI,rep,"\n kgA: " <>ToString[kgA]<>" \[Gamma]gA: "<>ToString[\[Gamma]gA](*parametros*),lagStep,species],{kgA,{0.05,1}},{\[Gamma]gA,{0.2,2}}],{3}];*)


(* ::Input:: *)
(*s*)


(* ::Input:: *)
(*s=Flatten[Table[simu[iteraciones2,replicas2,h,kIA,kAA,kAI,kgA,\[Gamma]gA,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kgI,\[Gamma]gI,kmI,\[Gamma]mI,kpI,\[Gamma]pI,rep,"\n kgA: " <>ToString[kgA]<>" \[Gamma]gA: "<>ToString[\[Gamma]gA](*parametros*),lagStep,species],{kgI,{0.05,1}},{\[Gamma]gI,{0.2,2}}],{3}];*)
