(* ::Package:: *)

(* ::Text:: *)
(*This notebook is to understand how the noise is propagated from one gene A to  gene B in a direct link of activation. Using Gillespie's Algorithm. Both genes have CONSTITUTIVE expression.*)


(* ::Section:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
{kmA,kmB}={0.696,0.696};(* mRNA production rate  (molecules*t^-1) *)
{\[Gamma]mA,\[Gamma]mB}={0.02082,0.02082} ;(* mRNA Degradation rate (molecules*t^-1) *)
{kpA,kpB}={1.386,1.386};(* Proteins production rate (molecules*t^-1) *) 
{\[Gamma]pA,\[Gamma]pB}={0.02082 ,0.02082; (* Proteins degradation rate (molecules*t^-1) *)
kAB=100;(* Michaelis-Menten constant *)
h=3 (* Hill's constant *);


(* ::Input::Initialization:: *)
iteraciones2=400000;(* Number of algorithm iterations *)
replicas2=4; (* Number of times the algorithm is run *)
rep=1;  (* Some plots are made with the output of this replicate*)
lagStep=10;  (* Lag time for the Autocorrelation plot *)
species=4;(*mRNAA, proteinA*)  



(* ::Section:: *)
(*Functions*)


(* ::Text:: *)
(*These are the reactions (changesByRXN) and propensity functions (pF) of the GA for the specific regulation system: A: Inhibitor -> B : Inhibited gene*)


(* ::Input::Initialization:: *)
changesByRXN=Function/@<|1->{#1+1,#2,#3,#4},2->{#1,#2+1,#3,#4},3->{#1-1,#2,#3,#4},4->{#1,#2-1,#3,#4},5->{#1,#2,#3+1,#4},6->{#1,#2,#3,#4+1},7->{#1,#2,#3-1,#4},8->{#1,#2,#3,#4-1}
|>;(* Each hash represents a specie: mRNAA-proteinA-mRNAB-proteinB-da0A:promoterActivatedA-da1A:promoterInactivatedA-db0B-db1B
   The numbers are the next reactions:
1-2-5-6\[Rule] Synthesis,  3-4-7-8\[Rule] Degradation*)
pF=Compile[{{h,_Real},{kAB,_Real},{kmA,_Real},{\[Gamma]mA,_Real},{kpA,_Real},{\[Gamma]pA,_Real},{kmB,_Real},{\[Gamma]mB,_Real},{kpB,_Real},{\[Gamma]pB,_Real},{mRNAA,_Integer},{proteinA,_Integer},{mRNAB,_Integer},{proteinB,_Integer}},{kmA+0.01,kpA*mRNAA,mRNAA*\[Gamma]mA,\[Gamma]pA*proteinA,kmB*1/(1+(proteinA/kAB)^h)+0.01,kpB*mRNAB,mRNAB*\[Gamma]mB,\[Gamma]pB*proteinB},RuntimeAttributes->{Listable},Parallelization->True];



(* ::Text:: *)
(*These functions are equal for all regulation systems and they are used to estimate the FF and CV2 in the steady state*)


(* ::Input::Initialization:: *)
ff[dato_]:=N[Variance[dato]/Mean[dato]];(* Fano Factor for noise measure *)
expressionNoise[data_]:=N[Variance[data]/Mean[data]^2](* Squared Coefficient of Variation for noise measure *)


(* ::Text:: *)
(*These functions make different plots of the temporal genes expression *)


(* ::Input::Initialization:: *)
(* This plots the temporal dynamic of genes mRNA and proteins *)
plotsExpression[rep_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,mRNAffcv2MeanI_,proteinffcv2MeanI_,parametros_]:=Module[{datesListAm,datesListAp,datesListBm,datesListBp},
SetOptions[ListLinePlot,Frame-> True,ImageSize->Medium];

{datesListAm,datesListAp,datesListBm,datesListBp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,4}];{ListLinePlot[datesListAm,PlotStyle->
RGBColor[0.5,0.99,1],FrameLabel->{"time","mRNA A molecules"},PlotLabel->labelA],
ListLinePlot[datesListBm,PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA B molecules"},PlotLabel->labelB],
ListLinePlot[datesListAp,PlotStyle->
RGBColor[0.3,.09,0.7],FrameLabel->{"time","protein A molecules"},PlotLabel->labelAp],
ListLinePlot[datesListBp,PlotStyle->
RGBColor[0.2,.8,1],FrameLabel->{"time","protein B molecules"},PlotLabel->labelBp]}]


(* ::Input::Initialization:: *)
(* This plots the temporal dynamic of genes mRNA and proteins with log scale *)
plotsLogExpression[rep_,ini_,listDataAndTimeTN_,mRNAffcv2MeanA_,proteinffcv2MeanA_,mRNAffcv2MeanI_,proteinffcv2MeanI_,parametros_]:=Module[{datesListAm,datesListAp,datesListBm,datesListBp},
SetOptions[ListLogPlot,Frame-> True,ImageSize->Medium,Joined->True];

{datesListAm,datesListAp,datesListBm,datesListBp}=Table[Transpose[{listDataAndTimeTN[[All,All,1]][[rep]],listDataAndTimeTN[[All,All,2]][[rep,All,i]]}],{i,4}];{ListLogPlot[datesListAm[[ini;;]],PlotStyle->
RGBColor[0.5,0.99,1],FrameLabel->{"time","mRNA A molecules"},PlotLabel->labelA],
ListLogPlot[datesListBm[[ini;;]],PlotStyle->
RGBColor[0.2,0.6,1],FrameLabel->{"time","mRNA B molecules"},PlotLabel->labelB],
ListLogPlot[datesListAp[[ini;;]],PlotStyle->
RGBColor[0.3,.09,0.7],FrameLabel->{"time","protein A molecules"},PlotLabel->labelAp],
ListLogPlot[datesListBp[[ini;;]],PlotStyle->
RGBColor[0.2,.8,1],FrameLabel->{"time","protein B molecules"},PlotLabel->labelBp]}]


(* ::Input::Initialization:: *)
(* This plots the auto-correlation of the temporal dynamic of genes mRNA and proteins *)
autoCorrPlot[lagStep_,listDataAndTimeTN_,positions_]:=Module[{data,lenght,data2,corr,negatives,cortes},
data=Table[MapThread[Extract[#1,#2]&,{listDataAndTimeTN[[All,All,2]][[All,All,i]],positions}],{i,{1,3,2,4}}];
lenght=Min[Length/@#]&/@data;
data2=MapThread[TemporalData[#1[[All,1;;#2]],{Range[1,#2]}]&,{data,lenght}];(*This is to get the corr in each lag time as the average between replicates *)
DistributeDefinitions[data,data2];
corr=ParallelMap[CorrelationFunction[#,{1,#["PathLengths"][[1]]-1,lagStep}]&,data2];
negatives=FirstPosition[Normal[#],{_,_?Negative}]*lagStep&/@corr;
cortes=Count[Partition[#//Normal,2,1],{{_,_?Positive},{_,_?Negative}}|{{_,_?Negative},{_,_?Positive}}]&/@corr;
ListPlot[#[[1]],Filling->Axis ,Frame->True,ImageSize->Medium,FrameLabel->{"Lag-"<>ToString[lagStep],#[[2]]},PlotStyle->#[[4]],PlotLabel->"First negative: "<>ToString[#[[3]]]<>"\n cortes: "<>ToString[#[[5]]]]&/@{{corr[[1]],"mRNA A",negatives[[1]],RGBColor[0.5,0.99,1],cortes[[1]]},{corr[[2]],"mRNA B",negatives[[2]],RGBColor[0.2,0.6,1],cortes[[2]]},{corr[[3]],"protein A",negatives[[3]],RGBColor[0.3,.09,0.7],cortes[[3]]},{corr[[4]],"protein B",negatives[[4]],RGBColor[0.2,.8,1],cortes[[4]]}}]


(* ::Input::Initialization:: *)
(* This plots the histogram of the steady state gene expression for mRNA and proteins *)
distriPlot[ini_,listDataAndTimeTN_]:=Module[{distribucionAm,distribucionAp,distribucionBm,distribucionBp},
{distribucionAm,distribucionAp,distribucionBm,distribucionBp}=Table[Flatten[listDataAndTimeTN[[All,All,2]][[All,ini;;All,i]]],{i,4}];Histogram[#[[1]],Automatic(*{Min[#[[1]]]//Round,Max[#[[1]]]//Round,5}*),"Probability",FrameLabel->{#[[2]],"Probability"},PlotLabel->#[[4]]<>" \n FF: "<>ToString[ff[#[[1]]]]<>", \!\(\*SuperscriptBox[\(CV\), \(2\)]\): "<>ToString[expressionNoise[#[[1]]]]<>", Mean: "<>ToString[Mean[#[[1]]]//N],ChartStyle->#[[3]],Frame->True,ImageSize->Medium]&/@{{distribucionAm,"mRNA A molecules",RGBColor[0.5,0.99,1],labelA},{distribucionBm,"mRNA B molecules",RGBColor[0.2,0.6,1],labelB},{distribucionAp,"protein A molecules",RGBColor[0.3,.09,0.7],labelAp},{distribucionBp,"protein B molecules ",RGBColor[0.2,.8,1],labelBp}}]


(* ::Input::Initialization:: *)
(* This plots the temporal dynamics of the regulator (proteins) and the regulated gene (proteins and mRNA) and indicates the Pearson Correlation between temporal dynamics of regulator and regulated gene *)corrPlot[replicas_,rep_,listDataAndTimeTN_]:=Module[{dataBmRNA,dataAp,dataBp,corremRNA,correprotein},
{dataBmRNA,dataAp,dataBp}=Table[#/Max[#]&[listDataAndTimeTN[[All,All,2]][[i,All,j]]],{j,{3,2,4}},{i,replicas}];{corremRNA,correprotein}=Mean[Table[Correlation[#[[1]][[i]],#[[2]][[i]]]//N,{i,replicas}]]&/@{{dataAp,dataBmRNA},{dataAp,dataBp}};{ListLinePlot[{Transpose[{listDataAndTimeTN[[rep,All,1]],dataAp[[rep]]}],Transpose[{listDataAndTimeTN[[rep,All,1]],dataBmRNA[[rep]]}]},PlotStyle->{RGBColor[0.3,.09,0.7],RGBColor[0.2,0.6,1]},Frame->True,ImageSize->Medium,PlotLabel->"Pearson Correlation: "<>ToString[corremRNA],PlotLegends->{"protein A","mRNA B"}],ListLinePlot[{Transpose[{listDataAndTimeTN[[rep,All,1]],dataAp[[rep]]}],Transpose[{listDataAndTimeTN[[rep,All,1]],dataBp[[rep]]}]},PlotStyle->{RGBColor[0.3,.09,0.7],RGBColor[0.2,.8,1]},Frame->True,ImageSize->Medium,PlotLabel->"Pearson Correlation: "<>ToString[correprotein],PlotLegends->{"protein A","protein B"}]}]


(* ::Input::Initialization:: *)
(* This plots the velocity of the temporal dynamic of  gene expression *)veloPlot[rep_,listDataAndTimeTN_,positionsVelocity_]:=Module[{listS,listT,velocidades},
listS=Table[Extract[#[[All,i]],positionsVelocity],{i,4}]&[listDataAndTimeTN[[rep,All,2]]];
listT=Extract[#,positionsVelocity]&[listDataAndTimeTN[[rep,All,1]]];
velocidades=Table[(listS[[i,2;;]]-listS[[i,;;-2]])/(listT[[2;;]]-listT[[;;-2]]),{i,4}];
Table[ListLinePlot[Transpose[{listT[[2;;]],velocidades[[i[[1]]]]}],PlotStyle->i[[2]],PlotLabel->"Velocities",FrameLabel->{"time",i[[3]]}],{i,{{1,RGBColor[0.5,0.99,1],"mRNA A molecules"},{2,RGBColor[0.3,.09,0.7],"protein A molecules"},{3,RGBColor[0.2,0.6,1],"mRNA B molecules"},{4,RGBColor[0.2,.8,1],"protein B molecules"}}}]]


(* ::Text:: *)
(*This function runs the GA algorithm and processes the output according to the previous functions.*)


(* ::Input::Initialization:: *)
simu[iteraciones_,replicas_,h_,kAB_,kmA_,\[Gamma]mA_,kpA_,\[Gamma]pA_,kmB_,\[Gamma]mB_,kpB_,\[Gamma]pB_,rep_,parametros_,lagStep_,species_]:=

Module[{listEstimatemRNAA,listEstimateproteinA,listEstimatemRNAI,listEstimateproteinI,listDataAndTimeT,listDataAndTime,t,abundancia,\[Tau]is,reactionAnd\[Tau]\[Mu],listData,mRNAffcv2ItA,proteinffcv2ItA,mRNAffcv2ItI,proteinffcv2ItI,mRNAffcv2MeanA,proteinffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanI,listDataAndTimeTN,positions,ini},

DistributeDefinitions[changesByRXN,pF,ff,expressionNoise,iteraciones,ini,h,kAB,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kmB,\[Gamma]mB,kpB,\[Gamma]pB];

listEstimatemRNAA={};
listEstimateproteinA={};
listEstimatemRNAI={};
listEstimateproteinI={};
listDataAndTimeT={};

SetSharedVariable[listEstimatemRNAA,listEstimateproteinA,listEstimatemRNAI,listEstimateproteinI,listDataAndTimeT];

ParallelDo[

listDataAndTime=CreateDataStructure["OrderedHashSet"];
t=0;
abundancia={1,1,1,1}(*mRNAA,proteinA,mRNAI,proteinI*);


Do[
          \[Tau]is=Map[If[#>0,RandomVariate[ExponentialDistribution[#]],\[Infinity]]&,pF[h,kAB,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kmB,\[Gamma]mB,kpB,\[Gamma]pB,Sequence@@abundancia]];reactionAnd\[Tau]\[Mu]={FirstPosition[#,Min[#]],Min[#]}&[\[Tau]is];abundancia=Through[{reactionAnd\[Tau]\[Mu][[1]]/.changesByRXN}[[1,1]]@@abundancia][[1]];t=t+reactionAnd\[Tau]\[Mu][[2]];listDataAndTime["Union",{{t,abundancia}}];
 If[Length[#]>3&&#[[2]]==5&[IntegerDigits[Round[t]]](* para que Break a los 1500*), Break[]],iteraciones];

listData=Normal[listDataAndTime][[All,2]];
ini =Echo@Round[(Length[listData]*30)/100.] ;

{mRNAffcv2ItA,proteinffcv2ItA}=Table[{ff[Flatten@i],expressionNoise[Flatten@i],Mean[Flatten@i]},{i,{listData[[ini;;All,1]],
listData[[ini;;All,2]]}}];
{mRNAffcv2ItI,proteinffcv2ItI}=Table[{ff[Flatten@i],expressionNoise[Flatten@i],Mean[Flatten@i]},{i,{listData[[ini;;All,3]],
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

Export["controlABGAInhibConst_"<>parametros<>".csv",{listEstimatemRNAA,listEstimateproteinA,listEstimatemRNAI,listEstimateproteinI}];

listDataAndTimeTN=Normal/@listDataAndTimeT;

{labelA,labelB,labelAp,labelBp}=Table["FF: "<>ToString[i[[1]]]<>" CV2: "<>ToString[i[[2]]]<>" Mean: "<>ToString[i[[3]]]<>parametros,{i,{mRNAffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanA,proteinffcv2MeanI}}];

ini =Round[(Length[listDataAndTimeTN[[rep]]]*30)/100.] ;

positions=Table[DeleteCases[FirstPosition[listDataAndTimeTN[[i,All,1]]//Round,#]&/@Range[1,1500,1],Missing["NotFound"]],{i,replicas}];(* These are the positions the ~1500 points to get the velocities and the auto-correlation plot*)

{plotsExpression[rep,listDataAndTimeTN,mRNAffcv2MeanA,proteinffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanI,parametros],
plotsLogExpression[rep,ini,listDataAndTimeTN,mRNAffcv2MeanA,proteinffcv2MeanA,mRNAffcv2MeanI,proteinffcv2MeanI,parametros],
autoCorrPlot[lagStep,listDataAndTimeTN,positions],
distriPlot[ini,listDataAndTimeTN],
corrPlot[replicas,rep,listDataAndTimeTN],
veloPlot[rep,listDataAndTimeTN,positions[[rep]]]}]


(* ::Subsection:: *)
(*Simulation*)


(* ::Input::Initialization:: *)
s=simu[iteraciones2,replicas2,h,kAB,kmA,\[Gamma]mA,kpA,\[Gamma]pA,kmB,\[Gamma]mB,kpB,\[Gamma]pB,rep,"\n Default para"(*parametros*),lagStep,species];


(* ::Input::Initialization:: *)
Export["graficasABGAInhibConstReg.pdf",s];
