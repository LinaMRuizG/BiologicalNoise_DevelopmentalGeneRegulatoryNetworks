(* ::Package:: *)

(* ::Text:: *)
(*This notebook simulates a self-regulated gene with the Chemical Langevin Equation (CLE) implemented by Yan et al 2017. Particularly, the CLE for burst in mRNA but no in protein (CLE1).*)


(* ::Section:: *)
(*Parameters*)


(* ::Text:: *)
(*These are the ranges of values for parameters variation . The range of each parameter is divided into smaller ranges (subRanges) . In this way the complete range is evaluated in parallel files each one with a subRange of parameters values .*)


(* ::Input::Initialization:: *)
iteraciones=5000; (* Number of algorithm iterations *)
replicas=10;  (* Number of times the algorithm is run *)
name=0.5;(* This is the last value of the range for the parameter being evluated *)
kgs=Range[0.01,1,0.01];
\[Gamma]gs=Range[0.02,name,0.02]; (* subRanges: 0.02-0.5, 0.52-1.0, 1.02-1.5, 1.52-2 *)
\[Gamma]ms=Range[0.001,0.1(*name*),0.001];(* subRanges: 0.001-0.025, 0.026-0.05, 0.051-0.075, 0.076-0.1 *)
\[Gamma]ps=Range[0.001,0.1(*name*),0.001];(* subRanges: 0.001-0.025, 0.026-0.05, 0.051-0.075, 0.076-0.1 *)
parametro="\[Gamma]g";


(* ::Input::Initialization:: *)
Clear[kg,\[Gamma]g,bm,km,\[Gamma]m,kp,\[Gamma]p,\[Tau],\[Epsilon],g,nc];
kg=0.01158;(* Gene activation rate (t^-1)*)
\[Gamma]g= 2.082;(* Gene desactivation rate (t^-1) *)
km=0.696;(* Production rate of mRNA (molecules*t^-1) *)
bm=km/(\[Gamma]g+kg);(* Burst size mean (molecules-mRNA) *)
\[Gamma]m=0.02082;(* Degradation rate (molecules*t^-1) *)
kp=1.386;(* Production rate of proteins (molecules*t^-1) *)
\[Gamma]p=0.02082 ;(* Degradation rate (molecules*t^-1) *)
bp=kp/\[Gamma]m;(* Burst size mean (molecules-protein) *)
kaa=\[Gamma]g/kg; (* Michaelis-Mente constant *)
h=3;  (* Hill's constant *)
\[Epsilon]=0.03; (* First control parameter_ the change in the propensity function or in the number of molecules in a step \[Tau] is bounded by this number *)
g={1};(* Highest order of reaction in which specie i appears as a reactant, for these cases production and degradation are the Fisrt Order*)
nc=10; (* second control parameters_ indicates the maximum number of fires for a reaction be considered as a critical one *)
species:={mRNA,protein}; (* The values of these are changing in each step \[Tau]*)


(* ::Input:: *)
(*(* NOTA IMPORTANTE: kp y \[Gamma]p no pueden ser exactamente iguales *)*)


(* ::Section:: *)
(*Burst mRNA and non-protein burst*)


(* ::Subsection:: *)
(*Functions*)


(* ::Input::Initialization:: *)
numero=0.01; (* Basal transcription *)
autoActivation[activator_,\[Gamma]g_,kg_,h_]:=(activator^h/((\[Gamma]g/kg)+activator ^h)); (* self-regulation Hill function *)
ff1[dato_]:=N[Variance[dato]/Mean[dato]];  (* FF for noise meadure *)
expressionNoise[data_]:=N[Variance[data]/Mean[data]^2](* CV2 for noise meadure *)


(* ::Input::Initialization:: *)
(* CLE for mRNA expression *)equmRNA[mRNA_,bm_,kg_,\[Gamma]m_,regulation_,\[Tau]_]:=Module[{rnd1=RandomVariate[NormalDistribution[0,1]],rnd2=RandomVariate[NormalDistribution[0,1]],mRNAs},mRNAs=mRNA+(kg \[Tau] bm regulation+(kg \[Tau] bm (2 bm + 1)regulation)^(1/2) rnd1 )-(\[Gamma]m mRNA \[Tau]+(\[Gamma]m mRNA \[Tau])^(1/2)   rnd2); If[mRNAs>0,mRNAs,numero]];

(* CLE for protein expression*)
equProtein[protein_,mRNA_,kp_,\[Gamma]p_,\[Tau]_]:=Module[{rnd1=RandomVariate[NormalDistribution[0,1]],rnd2=RandomVariate[NormalDistribution[0,1]],proteins},proteins=protein+(kp \[Tau] mRNA+(kp \[Tau] mRNA)^(1/2) rnd1 )-(\[Gamma]p protein \[Tau]+(\[Gamma]p protein \[Tau])^(1/2)  rnd2); If[proteins>0,proteins,numero]];

(* Algorithm to select \[Tau] or \[CapitalDelta]t *)
\[Tau]Choicee[mRNA_,protein_,regulation_,kg_,bm_,\[Gamma]m_,kp_,\[Gamma]p_]:=Module[{all\[Mu]varNC,\[Tau]1,mRNAi,proteini},
mRNAi=If[mRNA==0,numero,mRNA];
proteini=If[protein==0,numero,protein];
(* reactionsVandPF={{{1,kg*bm},{-1,\[Gamma]m*mRNAi}},{{1,kp*proteini},{-1,\[Gamma]p*proteini}}}For each specie introduce the state change value (v) and the propensity (PF) of synthesis and degradation reactions: {{v,PF}-synthesis,{v,PF}-degradation}-specie *);
all\[Mu]varNC={{mRNAi,kg*bm*regulation,kg*bm*regulation},{mRNAi,-\[Gamma]m*mRNAi,\[Gamma]m*mRNAi},{proteini,kp*mRNAi,kp*mRNAi},{proteini,-\[Gamma]p*proteini,\[Gamma]p*proteini}};
\[Tau]1={Sequence@@Table[Min[Max[\[Epsilon] i[[1]]/g,1]/Abs[i[[2]] ],Max[\[Epsilon] i[[1]]/g,1]^2/i[[3]] ],{i,all\[Mu]varNC[[1;;2]]}],Min[Max[\[Epsilon] proteini/g,1]/Abs[Total[all\[Mu]varNC[[3;;4,2]]] ],Max[\[Epsilon] proteini/g,1]^2/Total[all\[Mu]varNC[[3;;4,3]]]]};(*ESTO SALE DEL ART\[CapitalIAcute]CULO DE Cao 2016 y CHING-YAN 2017. Para cada especie se calcula el tao. Cuando ambas reacciones son no-crticas este se calcula con la suma de la media y la varianza de ambas reacciones. Cuando la degradaci\[OAcute]n se torna critica estas se separan. Cuando hay burst se calcula una tipo de tao para el burst y otro para la degradaci\[OAcute]n (uno cuando es no-critica y otro cuando es critica) *)
Min[\[Tau]1[[1]],
If[mRNAi>=10,\[Tau]1[[2]],1/(\[Gamma]m*mRNAi) Log[1/(1-RandomReal[{0,1}])]] (*ESTO SALE DEL ART\[CapitalIAcute]CULO DE Cao 2016*),
If[proteini>=10,\[Tau]1[[3]],Sequence@@{Min[Max[\[Epsilon] proteini/g,1]/Abs[all\[Mu]varNC[[3,2]]],Max[\[Epsilon] proteini/g,1]^2/all\[Mu]varNC[[3,3]]],1/(\[Gamma]p*proteini) Log[1/(1-RandomReal[{0,1}])]}]]]


(* ::Subsection:: *)
(*Solution-iteration*)


(* ::Input::Initialization:: *)
mRNAffcv2List={};
proteinffcv2List={};
Table[
listEstimatemRNA={};
listEstimateprotein={};
listDataTmRNA={};
listDataTprotein={};
timesT={};
Table[
mRNA=1;
mRNAlist={1};
protein=1;
proteinlist={1};
t=0;
tlist={0};
Do[
\[Tau]=\[Tau]Choicee[mRNA,protein,autoActivation[protein,\[Gamma]g,kg,h],kg,km/(\[Gamma]g+kg),\[Gamma]m,kp,\[Gamma]p];
t=t+\[Tau];
AppendTo[tlist,t];
mRNA=equmRNA[mRNA,km/(\[Gamma]g+kg),kg,\[Gamma]m,autoActivation[protein,\[Gamma]g,kg,h],\[Tau]];
AppendTo[mRNAlist,mRNA];
protein=equProtein[protein,mRNA,kp,\[Gamma]p,\[Tau]];
AppendTo[proteinlist,protein],iteraciones];
{mRNAffcv2It,proteinffcv2It}=Table[{ff1[Flatten@i],expressionNoise[Flatten@i]},{i,{mRNAlist[[10000;;All]],
proteinlist[[10000;;All]]}}];
AppendTo[listEstimatemRNA,mRNAffcv2It];
AppendTo[listEstimateprotein,proteinffcv2It],replicas];
AppendTo[mRNAffcv2List,Mean[listEstimatemRNA]];
AppendTo[proteinffcv2List,Mean[listEstimateprotein]],{kg,kgs},{\[Gamma]g,\[Gamma]gs}(*,{\[Gamma]m,\[Gamma]ms}*)];


(* ::Input::Initialization:: *)
Export["burstmRCLE_"<>parametro<>ToString[name],{mRNAffcv2List,proteinffcv2List},"CSV"];
