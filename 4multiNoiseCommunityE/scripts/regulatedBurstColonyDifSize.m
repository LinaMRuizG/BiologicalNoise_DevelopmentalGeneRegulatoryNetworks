(* ::Package:: *)

(* ::Text:: *)
(*Here it is evaluated the noise in a cell colony where each cell expresses a self-regulated gene. The cells are connected by diffusion of the genetic product (i.e. paracrine signal). In each cell, the regulatory system is represented two Chemical Langeving Equations (CLE).  One for mRNA expression and the other one for protein expression according to Yan et al 2017.*)


SetDirectory[]


(* ::Section:: *)
(*Parameters*)


(* ::Text:: *)
(*These are the ranges of values for parameters variation. The range of each parameter is divided into smaller ranges (subRanges). In this way the complete range is evaluated in parallel files each one with a subRange of parameters values.*)


(* ::Input::Initialization:: *)
radius=5; (* Coloby cell radius *)
iteraciones=1000; (* Number of algorithm iterations. An error could be get if there is less than 500 iterations *)
replicas=2;  (* Number of times the algorithm is run *)

name=0.009;   (* This is the last value of the range for the parameter being evluated *)
kgs=Range[0.01,1,0.02];
\[Gamma]gs=Range[0.02,2,0.04]; (* subRanges: 0.02-0.18, 0.22-0.38, 0.42-0.58, 0.62-0.78, 0.82-0.98, 1.02-1.18, 1.22-1.38, 1.42-1.58, 1.62-1.78, 1.82-1.98 *)
\[Gamma]ms=Range[0.001,0.1,0.002];(* subRanges: 0.001-0.009, 0.011-0.019, 0.021-0.029, 0.031-0.039, 0.041-0.049, 0.051-0.059, 0.061-0.069, 0.071-0.079, 0.081-0.089, 0.091-0.099 *)
\[Gamma]ps=Range[0.001,name,0.002];(* subRanges: 0.001-0.009, 0.011-0.019, 0.021-0.029, 0.031-0.039,0.041-0.049, 0.051-0.059, 0.061-0.069, 0.071-0.079,0.081-0.089, 0.091-0.099 *)
dfs=Range[0.01,0.1,0.01]; (* range values for diffussion coeficient*)

(* variables to name the output data*)
size="default1_"<>ToString[radius];
para1="kon: ";
para1v=kg;
para2=" koff: ";
para2v=\[Gamma]g;
parametro1="_\[Gamma]m";
parametro2="_\[Gamma]p";(* Poner aqui el segundo parametro. Con este se nombra el archivo de salida*)
(*NOTA: Para hacer las graficas se debe ir a cambiar los parametros en abajo*)


(* ::Input::Initialization:: *)
LaunchKernels[8]
dtGifs=10; 
ventanas=4;


Clear[kg,\[Gamma]g,bm,km,\[Gamma]m,kp,\[Gamma]p,\[Tau],\[Epsilon],g,nc];
(*************************TO VARIATE THE REGIONS OF EXPRESSION*)
kg=0.01;(* Gene activation rate (t^-1)*)
\[Gamma]g=2;(* Gene desactivation rate (t^-1) *)
(*************************PARAMETERS DEFAULT*)
km=0.696;(* Production rate of mRNA (molecules*t^-1) *)
\[Gamma]m=0.02082;(* Degradation rate (molecules*t^-1) *)
kp=1.386;(* Production rate of proteins (molecules*t^-1) *)
\[Gamma]p=0.02082 ;(* Degradation rate (molecules*t^-1) *)
bp=kp/\[Gamma]m;(* Burst size mean (molecules-protein) *)
bm=km/(\[Gamma]g+kg);(* Burst size mean (molecules-mRNA) *)
kaa=\[Gamma]g/kg; (* Michaelis-Mente constant *)
df=0.01;(* pixeles/minuto , como en Muller 0.7\[Mu]m2/s*)
h=3; (* Hill's constant *)
dfs=0.01;(* diffusion coeficient: pixeles/minuto , as in Muller 0.7\[Mu]m2/s*)
\[Epsilon]=0.03; (* First control parameter_ the change in the propensity function or in the number of molecules in a step \[Tau] is bounded by this number *)
g={1};(* Highest order of reaction in which specie i appears as a reactant, for these cases production and degradation are the Fisrt Order*)
nc=10; (* second control parameters_ indicates the maximum number of fires for a reaction be considered as a critical one *)
species:={mRNA,protein}; (* The values of these are changing in each step \[Tau]*) 


(* ::Input:: *)
(*(* NOTA IMPORTANTE: kp y \[Gamma]p no pueden ser exactamente iguales *)*)


(* ::Section::Initialization:: *)
(*Functions*)


(* ::Input::Initialization:: *)
numero=0.01;  (* Basal transcription *)
autoActivation[activator_,\[Gamma]g_,kg_,h_]:=(activator^h/((\[Gamma]g/kg)+activator ^h));(* self-regulation function *)
ff1[dato_]:=N[Variance[dato]/Mean[dato]];  (* FF for noise meadure *)
expressionNoise[data_]:=N[Variance[data]/Mean[data]^2](* CV2 for noise meadure *)


(* ::Input::Initialization:: *)
(* CLE for mRNA expression*)equmRNA[mRNA_,bm_,kg_,\[Gamma]m_,regulation_,\[Tau]_]:=Module[{rnd1=RandomVariate[NormalDistribution[0,1]],rnd2=RandomVariate[NormalDistribution[0,1]],mRNAs},mRNAs=mRNA+(kg \[Tau] bm regulation+(kg \[Tau] bm (2 bm + 1)regulation)^(1/2) rnd1 )-(\[Gamma]m mRNA \[Tau]+(\[Gamma]m mRNA \[Tau])^(1/2)   rnd2); If[mRNAs>0,mRNAs,numero]];

(* CLE for protein expression*)equProtein[entradaDif_,protein_,mRNA_,kp_,\[Gamma]p_,df_,\[Tau]_]:=Module[{rnd1=RandomVariate[NormalDistribution[0,1]],rnd2=RandomVariate[NormalDistribution[0,1]],proteins},proteins=entradaDif*\[Tau]+protein+(kp \[Tau] mRNA+(kp \[Tau] mRNA)^(1/2) rnd1 )-(\[Gamma]p protein \[Tau]+(\[Gamma]p protein \[Tau])^(1/2)  rnd2)-df protein *\[Tau]; If[proteins>0,proteins,numero]];

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


(* ::Input::Initialization:: *)
(* This function gives  the new abundance of cells expressing the self-regulated system and the \[Tau] *)newAbundanceUnos[kg_,\[Gamma]g_,bm_,\[Gamma]m_,kp_,\[Gamma]p_,h_,abundancia_,entradaDiff_,df_]:=Module[{\[Tau],mRNA,protein},
\[Tau]=\[Tau]Choicee[abundancia[[1]](*mRNA*),abundancia[[2]](*protein*),autoActivation[abundancia[[2]](*protein*),\[Gamma]g,kg,h],kg,bm,\[Gamma]m,kp,\[Gamma]p];mRNA=equmRNA[abundancia[[1]](*mRNA*),bm,kg,\[Gamma]m,autoActivation[abundancia[[2]](*protein*),\[Gamma]g,kg,h],\[Tau]];
protein=equProtein[entradaDiff,abundancia[[2]](*protein*),abundancia[[1]](*mRNA*),kp,\[Gamma]p,df,\[Tau]];{{mRNA,protein},\[Tau]}]


(* ::Input::Initialization:: *)
(* This function gives  the new abundance of cells non-expressing the self-regulated system (border cells) *)newAbundanceCeros[abundancia_,entradaDiff_,df_,kp_]:=Module[{abundancia0=abundancia},
abundancia0[[2]]=abundancia0[[2]]+entradaDiff-df*abundancia0[[2]]-kp*abundancia0[[2]];
abundancia0]
newAbundanceCeros2[abundancia_,entradaDiff_,df_,kp_]:=Module[{abundancia0=abundancia},
abundancia0[[2]]=abundancia0[[2]]+entradaDiff-df*abundancia0[[2]]-kp*abundancia0[[2]];
{abundancia0,0}]


(* ::Input::Initialization:: *)
(* This is the algorithm to estimate the Steady State start from Kelly and Hedengren 2013. It gives the FF since the position estimated as the start of the steady-sate *)ffs[data_]:=Module[{first,ff,n,ffPorVentana,mPorVentana,\[Mu],sd,probabilidades,posicion},
first=Position[data,{__,x_,__}/;x>0][[1,1]];
ff=Table[N[Variance[Flatten[data[[i]]]]/Mean[Flatten[data[[i]]]]],{i,first,data//Length}];
n=(ff//Length)/ventanas;
ffPorVentana=Partition[ff,n];
mPorVentana=Map[Mean,Map[#[[2;;All]]-#[[1;;-2]]&,ffPorVentana](*diferencias por ventana*)];
\[Mu]=Table[(1/n)(Total[ffPorVentana[[i]]]-mPorVentana[[i]]*Total[Range[1,n](*tPorVentana*)]),{i,ventanas}];
sd=Table[Sqrt[(1/(n-2))*Sum[(ffPorVentana[[i,j]]-mPorVentana[[i]]*j-\[Mu][[i]])^2,{j,1,Length[ffPorVentana[[i]]]}]],{i,ventanas}];
probabilidades=Map[Total[#]/n&,Table[If[Abs[ffPorVentana[[i,j]]-\[Mu][[i]]]<=sd[[i]],1,0],{i,ventanas},{j,1,n}]]//N;
posicion=((Position[probabilidades,x_/;x>0.1][[1,1]]-1)*n)+1;{ff[[posicion(*con este optengo la posicion a paritir de la cual considerar el ff como en el estado estable*);;All]],posicion}];


(* ::Input::Initialization:: *)
(* This is the algorithm to estimate the Steady State start from Kelly and Hedengren 2013 *)posicionF[data_]:=Module[{first,ff,n,ffPorVentana,mPorVentana,\[Mu],sd,probabilidades,posicion},
first=Position[data,{__,x_,__}/;x>0][[1,1]];
ff=Table[N[Variance[Flatten[data[[i]]]]/Mean[Flatten[data[[i]]]]],{i,first,data//Length}];
n=(ff//Length)/ventanas;
ffPorVentana=Partition[ff,n];
mPorVentana=Map[Mean,Map[#[[2;;All]]-#[[1;;-2]]&,ffPorVentana](*diferencias por ventana*)];
\[Mu]=Table[(1/n)(Total[ffPorVentana[[i]]]-mPorVentana[[i]]*Total[Range[1,n](*tPorVentana*)]),{i,ventanas}];
sd=Table[Sqrt[(1/(n-2))*Sum[(ffPorVentana[[i,j]]-mPorVentana[[i]]*j-\[Mu][[i]])^2,{j,1,Length[ffPorVentana[[i]]]}]],{i,ventanas}];
probabilidades=Map[Total[#]/n&,Table[If[Abs[ffPorVentana[[i,j]]-\[Mu][[i]]]<=sd[[i]],1,0],{i,ventanas},{j,1,n}]]//N;
posicion=((Position[probabilidades,x_/;x>0.1][[1,1]]-1)*n)+1];


(* ::Input::Initialization:: *)
(* With this function the mean and sd of the FFs are calculated for each of the cells of the colony and in the stable state *)msdFFbyCell[data_,posicion_]:=Module[{ffPorCelulaTodoT,disffPorCelulaTodoT},
ffPorCelulaTodoT=Map[ff1[#]&,Transpose[data][[All,posicion;;All]]];
disffPorCelulaTodoT=FindDistribution[ffPorCelulaTodoT,1,"AIC"][[1,1]];
{Mean[disffPorCelulaTodoT],StandardDeviation[disffPorCelulaTodoT]}]


(* ::Input::Initialization:: *)
(* With this function we calculate the mean of the mean and the sd of the number of molecules for all cells in each unit of time within the stable state *)
MeanMolecules[data_,posicion_]:=Mean[Map[Mean[#]&,data[[posicion;;All]]]]


(* ::Input::Initialization:: *)
(*With this function the entropy of the FFs of each circle is calculated in time within the steady state*)FFbyCirculoModuleF[resultados_,matrixCirculos_,posicion_]:=Module[{ffCirculoTiempo,meanFFCirculo},
ffCirculoTiempo=Table[Map[ff1[Pick[#,matrixCirculos,i]]&,resultados[[All,1]][[All,All,2]]],{i,radius}];
Round@Map[Mean,ffCirculoTiempo[[All,posicion;;All]]]
(*Esto es exactamente la Entropia de Shannon en base 10*)]


(* ::Input::Initialization:: *)
(* This gives the mean \[PlusMinus]  SD of some data *)meanSd[datos_]:=Transpose@Map[Through[{Around[Mean@#,StandardDeviation@#]&}@#]&]@Map[Flatten,Transpose[datos]];


(* ::Input::Initialization:: *)
solucion2D[step_,protein_,comunidadEstadolistTiempos_]:=Module[{celulasValue},
celulasValue=MapThread[#1-> #2&,{Flatten[Pick[matrixPos,matrix,1],1],comunidadEstadolistTiempos[[step,All,protein]]}];
ReplacePart[matrix,celulasValue]];(* To get the graphics *)


(* ::Input::Initialization:: *)
plot2DF[resultados_]:=Module[{colonyStateTime1,maxprotein},
colonyStateTime1=resultados[[All,1]][[All,posicionUnos]];
maxprotein=Max[Flatten[colonyStateTime1[[All,All,2]]]];
Table[ArrayPlot[solucion2D[i,2,colonyStateTime1],ColorFunctionScaling->False,PlotLegends->Automatic,PlotLabel->"protein molecules, step: "<>ToString[i],PlotRange->{0,maxprotein}],{i,Range[1,iteraciones,dtGifs]}]]
(* To get the graphics-gifs *)


(* ::Input::Initialization:: *)
ffCircleAcuTimeF[resultados_,posicion_,radius_]:=Module[{ffCirculoTiempo,moleculesCirculoAcuTime},
ffCirculoTiempo=Table[Map[Pick[#,matrixCirculos,i]&,resultados[[All,1]][[posicion;;All,All,2]]],{i,Reverse@Range[1,radius]}];
moleculesCirculoAcuTime=Table[Transpose[ffCirculoTiempo[[1;;i]]],{i,1,radius}];Table[Map[ff1[Flatten[#]]&,moleculesCirculoAcuTime[[j]]],{j,radius}]]


(* ::Input::Initialization:: *)
(*Plots the temporal FFs by circles *)plot1F[ffCirculoAcuTime_,par1_,par2_,resultados_]:=ListLinePlot[Transpose[{Accumulate[resultados[[All,2]]],
ffCirculoAcuTime}],Frame-> True,FrameLabel-> {"time (min)","FF steady state"},PlotLabel->para1<>ToString[par1]<>para2<>ToString[par2] ,PlotLegends->Range[1,radius],PlotRange-> All,ImageSize->Medium,PlotStyle->Flatten@Table[RGBColor[r,gr,b],{r,{0.5,1}},{gr,{0.1,0.9,.5}},{b,{0.1,0.8,0.2,0.6,0.4}}]]


(* ::Input::Initialization:: *)
(*Plots the temporal FFs by circles *)plot1Fotro[ffCirculoAcuTime_,par1_,par2_,times_]:=ListLinePlot[Transpose[{Accumulate[times],#}]&/@ffCirculoAcuTime,Frame-> True,FrameLabel-> {"time (min)","FF steady state"},PlotLabel->para1<>ToString[par1]<>para2<>ToString[par2] ,PlotLegends->Range[1,radius],PlotRange-> All,ImageSize->Medium,PlotStyle->Flatten@Table[RGBColor[r,gr,b],{r,{0.5,1}},{gr,{0.1,0.9,.5}},{b,{0.1,0.8,0.2,0.6,0.4}}]]


(* ::Input::Initialization:: *)
(*Plots the mean of the temporal FFs by circles *)plot2F[data1_,par1_,par2_]:=ListLinePlot[data1,PlotLabel-> para1<>ToString[par1]<>para2<>ToString[par2] ,Frame->True,FrameLabel->{"Disk number","Mean FF steady state"},PlotRange->All,ImageSize->Medium,PlotStyle->RGBColor[0,0.1,1]]


(* ::Input::Initialization:: *)
(*Plots the sd of the temporal FFs by circles *)plot3F[data2_,par1_,par2_]:=ListLinePlot[data2,PlotLabel-> "sd FF by acumulate circles from inside to outside \n "<>parametro1<>ToString[par1]<>parametro2<>ToString[par2] ,Frame->True,FrameLabel->{"N\[Degree] acumulate circle","SD FF steady state"},PlotRange->All,ImageSize->Medium,PlotStyle->RGBColor[0,0.1,.5]]


(* ::Input::Initialization:: *)
ffCircleAcuQF[ffCirculoAcuTime_]:=Module[{uno},
uno=Table[Select[Flatten[i],Quantile[Flatten[i],{0.1}][[1]]<#<Quantile[Flatten[i],{0.95}][[1]]&],{i,ffCirculoAcuTime}];
{Map[Mean,uno],Map[StandardDeviation,uno]}]


meanTempBox[resultados_]:=Module[{times},
times=Range[1,iteraciones,100];
BoxWhiskerChart[resultados[[times,1,All,2]],"Diamond",ChartStyle->RGBColor[0.2,.8,1],
ChartLabels->Placed[Accumulate[resultados[[times,2]]], {{0.5, 0}, {1,1.1}},Rotate[#,Pi/4]&]]]
meanTemp[resultados_]:=ListLinePlot[Transpose[{Accumulate[resultados[[All,2]]],
Map[Mean,resultados[[All,1,All,2]]]}],PlotLabel->para1<>ToString[para1v]<>para2<>ToString[para2v],
PlotRange->{{0,1500},All},Frame->True,FrameLabel->{"Time(min)","Mean molecules"},PlotStyle->RGBColor[0.2,.8,1]]


(* ::Section:: *)
(*Colonia*)


(* ::Text:: *)
(*Here it is created the colony and the neighbors of each cell are identified*)


(* ::Input::Initialization:: *)
capaCeros=3;
totalFila=radius*2+1+capaCeros*2;
matrix=Join[Table[Table[0,radius*2+1+capaCeros*2],3],Table[{0,0,0,Sequence@@i,0,0,0},{i,DiskMatrix[radius]}],Table[Table[0,radius*2+1+capaCeros*2],3]];
noSeconsideran=Position[Flatten[matrix],0];
neighbors=Flatten[MapThread[List,Map[RotateRight[matrix,#]&,{{0,0},{1,0},{0,-1},{-1,0},{0,1},{1,-1},{-1,-1},{-1,1},{1,1}}],2],1];
matrixPos=Partition[Tuples[Range[1,totalFila],2],totalFila];
neighborsPos=Flatten[MapThread[List,Map[RotateRight[matrixPos,#]&,{{0,0},{1,0},{0,-1},{-1,0},{0,1},{1,-1},{-1,-1},{-1,1},{1,1}}],2],1][[All,2;;All]];
filas=Flatten[matrixPos,1];
asociasiones1=Map[Association,Replace[neighborsPos,x_->( x-> 1),{2}]](* Each association is the neighboors for each cell *);
asociasiones0=Replace[filas,x_->( x-> 0),{1}]//Association (* Is the default, to put 0 in the matrix where there is not conection-neighboors *);
matrizCompleta=Map[#/.asociasiones0&,Map[filas/.#&,asociasiones1]];
totalEqu=Length[matrizCompleta];
vecinosT=Map[Count[#,1]&,matrizCompleta];(* Todos tiene 8 vecinos porque es *)
prVectorEntradaA=Table["df* colonyStateTime[[-1,"<>ToString[i]<>",2]]/"<>ToString[vecinosT[[i]]],{i,1,totalEqu}];
terminosEntradaA=Table[prVectorEntradaA . matrizCompleta[[i]],{i,Length[prVectorEntradaA]}](*dot product*);


(* ::Section::Initialization:: *)
(*Solution-iteration*)


(* ::Input::Initialization:: *)
(* Esto es la condicion inicial para cada c\[EAcute]lula de cada una de las especies {mRNA,protein,da0,da1} *)
colonyStateTime={Table[{10,10},totalEqu]};
posicionCeros=Position[Flatten[matrix],0]//Flatten;
posicionUnos=Position[Flatten[matrix],1]//Flatten;posicionesUnoMatrix=DeleteCases[Pick[matrixPos,matrix,1],{}];
posicionesCirculos=Union[Flatten@Table[Replace[Union[Flatten[posicionesUnoMatrix[[{i,-i}]],1],Flatten[posicionesUnoMatrix[[i;;-i,{i,-i}]],1]],x_->(x-> i),{1}] ,{i,radius}],{{radius+capaCeros+1,radius+capaCeros+1}-> 0}];
matrixCirculos=Flatten[ReplacePart[matrix,posicionesCirculos]];


(* ::Input::Initialization:: *)
Table[
allCell=Range[1,totalEqu];
funcionesCeros=Map[Hold[newAbundanceCeros2[colonyStateTime[[-1,#]],Map[ToExpression,terminosEntradaA[[#]]],df,kp]]&,posicionCeros];
funcionesUnos=Map[Hold[newAbundanceUnos[kg,\[Gamma]g,km/(\[Gamma]g+kg)(*bm*),\[Gamma]m,kp,\[Gamma]p,h,colonyStateTime[[-1,#]],Map[ToExpression,terminosEntradaA[[#]]],df]]&,posicionUnos];
asociacionCeros=MapThread[#1-> #2&,{posicionCeros,funcionesCeros}];
asociacionUnos=MapThread[#1-> #2&,{posicionUnos,funcionesUnos}]//Association;
asociacionCeros=MapThread[#1-> #2&,{posicionCeros,funcionesCeros}]//Association;
funcionesAllCell=(allCell/.asociacionCeros)/.asociacionUnos;
(*ojo*)
DistributeDefinitions[posicionF,msdFFbyCell,MeanMolecules,plot2DF,solucion2D ,ffCircleAcuTimeF,ffs,ff1,radius,posicionUnos,iteraciones,dtGifs,matrixPos,matrix,matrixCirculos,FFbyCirculoModuleF,newAbundanceUnos,newAbundanceCeros,newAbundanceCeros2,kg,\[Gamma]g,km,bm,\[Gamma]m,kp,\[Gamma]p,h,df,colonyStateTime,terminosEntradaA,plot1F,plot2F,plot3F,ffCircleAcuQF];
(*********************************to plots 6 and 7*)

(*ojo*)
ffPorCellReplicas={};
meanMoleculesReplicas={};
ffbyCirculoReplicas={};
meanffbyCirculoQReplicas={};
sdffbyCirculoQReplicas={};
SetSharedVariable[resultadosR,ffPorCellReplicas,meanMoleculesReplicas,ffbyCirculoReplicas,meanffbyCirculoQReplicas,sdffbyCirculoQReplicas];

ParallelDo[
Block[{colonyStateTime,resultados,colonyAndSpeciesState,data,posicion,mFFbyCell,sFFbyCell,meanMolecules,ffCirculoAcuTime,meanffCirculoAcuTimeQ,sdffCirculoAcuTimeQ},
colonyStateTime={Table[{10,10},totalEqu]};
resultados=Table[
colonyAndSpeciesState=Map[ReleaseHold,funcionesAllCell];
colonyStateTime={colonyAndSpeciesState[[All,1]]};
{colonyStateTime[[1]],Max[colonyAndSpeciesState[[All,2]]]},iteraciones];
data=resultados[[All,1]][[All,posicionUnos]][[All,All,2]];
posicion=posicionF[data];
{mFFbyCell,sFFbyCell}=msdFFbyCell[data,posicion];
meanMolecules=MeanMolecules[data,posicion];

(*Export["burstmRCLEColony"<>size<>"_plot1_df"<>ToString[df],plot2DF[resultados],"GIF"];*)


ffCirculoAcuTime=ffCircleAcuTimeF[resultados,posicion,radius];
{meanffCirculoAcuTimeQ,sdffCirculoAcuTimeQ}=ffCircleAcuQF[ffCirculoAcuTime];
AppendTo[ffPorCellReplicas,{mFFbyCell,sFFbyCell}];
AppendTo[meanMoleculesReplicas,meanMolecules];
AppendTo[ffbyCirculoReplicas,ffCirculoAcuTime];
AppendTo[meanffbyCirculoQReplicas,meanffCirculoAcuTimeQ];
AppendTo[sdffbyCirculoQReplicas,sdffCirculoAcuTimeQ];
resultadosR=resultados;],

replicas];

{meanFFbyCellrep,sdFFbyCellrep}=Map[Mean,Transpose[ffPorCellReplicas]];
meanMoleculesrep=Mean[Flatten[meanMoleculesReplicas]];
meanEntoripasRep=Mean[Map[N@Entropy[Round[#]]&,meanffbyCirculoQReplicas]];
meanRepFFByCircle=Map[Mean,Transpose[meanffbyCirculoQReplicas]];




,{df,{0.01,0.1}}];





(* ::Section:: *)
(*References*)


(* ::Text:: *)
(*Yan C - CS, Chepyala SR, Yen C - M, Hsu C - P . Efficient and flexible implementation of Langevin simulation for gene burst production . Sci Rep . 2017; 7 : 16851. doi : 10.1038/s41598 - 017 - 16835 - y*)


(*Export["bRCLEcolony"<>size<>"_df:"<>ToString[df]<>"_todos",{meanFFbyCellrep,sdFFbyCellrep,meanMoleculesrep,meanEntoripasRep,meanRepFFByCircle},"CSV"];*)
Export["burstmRCLEColony_plot2"<>size<>"_df"<>ToString[df]<>".png",plot1F[ffbyCirculoReplicas[[1]],para1v,para2v,resultadosR],"PNG"];

Export["burstmRCLEColony_plot3"<>size<>"_df"<>ToString[df]<>".png",plot2F[meanSd[meanffbyCirculoQReplicas][[1]],para1v,para2v],"PNG"];

(*Export["burstmRCLEColony_plot4"<>size<>"_df"<>ToString[df],plot3F[meanSd[sdffbyCirculoQReplicas][[1]],\[Gamma]m,\[Gamma]p],"PNG"];*)

Export["burstmRCLEColony_plot5.png"<>size<>"_df"<>ToString[df]<>".png",Histogram[Round@meanffbyCirculoQReplicas[[1]],
"Probability",PlotLabel->para1<>ToString[para1v]<>para2<>ToString[para2v]<>" \n Entropy: "<>ToString[meanEntoripasRep],Frame->True,FrameLabel->{"FF Mean by disk","Probability"}],"PNG"];
(***************************************************)

(*Export["burstmRCLEColony_plot6"<>size<>"_df"<>ToString[df]<>".png",meanTempBox[resultadosR],"PNG"];*)
(*Export["burstmRCLEColony_plot7"<>size<>"_df"<>ToString[df]<>".png",meanTemp[resultadosR],"PNG"];*)
(*************************************************)



