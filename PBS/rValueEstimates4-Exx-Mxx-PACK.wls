#!/usr/bin/env wolframscript 
(* ::Package:: *)

(* ::Section:: *)
(**)


Print["(*Preliminaries*)"];

$LOG=$DBG=$RVL=False;


LL="31";

MMlist={"10","12","14", "16", "18", "20", "22", "24"};
MMlist=ReadList["!cat",Number];
(*MMlist = $ScriptInputString <> Import["!cat","String"]*)

Print[MMlist];

(*Print[size];
If[IntegerQ[MMsize], MMlist={ToString[size]}];
Print[MMlist];*)

(*maxsamples=Infinity;maxdir=Infinity;*)
maxsamples=100;maxdir=2;


maindir=Directory[]
Print[maindir];


(* ::Section:: *)
(**)


Print["(*Analysis*)"];


SetDirectory[maindir];


Do[
MM=ToString[MMlist[[iMM]]];
Print["--- working on M=",MM];

SetDirectory[maindir];
alldirs=FileNames["L"<>LL<>"-"<>MM<>"-*/L"<>LL<>"_M"<>MM<>"*"];
Print[alldirs];
(*alldirs=FileNames["QPr-ABPN2*"]*)
allavglist={};

starttimeD=AbsoluteTime[];
lendirs=Min[maxdir,Length[alldirs]];
Do[
dirname=alldirs[[idir]];
Print["--- "<>dirname];
SetDirectory[maindir<>"/"<>dirname];
spos=Union[Flatten[StringPosition[dirname,"/"]]][[1]];
topdirname=StringTake[dirname,{1,spos-1}];
bottomdirname=StringTake[dirname,{spos+1,StringLength[dirname]}];

HubDis=ToExpression[StringTake[bottomdirname,{11,16}]]/100.;
RimDis=ToExpression[StringTake[bottomdirname,{20,25}]]/100.;
TarEng=StringTake[bottomdirname,-7];
If[StringContainsQ[TarEng,"-"],
TarEng=ToExpression[TarEng]/100.,
TarEng=ToExpression[StringDrop[TarEng,1]]/100.
];

Print[{Directory[](*,TarEng,HubDis,RimDis*)}];

allfiles=FileNames["EVal*.raw"];
allenglist=allrlist=allstatrlist=allzlist=allstatzlist=allseeds={};
ndone=0;

starttime=AbsoluteTime[];
lensamples=Min[maxsamples,Length[allfiles]];
Do[
filename=allfiles[[ifile]];
If[Length[allfiles]<=0, Continue[]];

If[$LOG,PrintTemporary["--- --- "<>filename]];
rawdata=ReadList[filename,Number];
If[Length[rawdata]<=3, Continue[], ndone++];

AppendTo[allseeds,ToExpression[StringTake[StringDrop[filename,-4],-5]]];

(* drop the length item at the start of the list *)
len=rawdata[[1]];

eng=Sort[DeleteCases[rawdata, _Integer]];
AppendTo[allenglist,eng];

(* raw spacings, removing all negative/zero energies *)
tiny=10^-10;
engpos0=Sort[Select[eng,#>tiny&]];

(* remove flat band energies beyond E=0 *)
Switch[LL,
"31",engpos=engpos0,
"32",engpos=Sort[Select[engpos0,Abs[Abs[#]-1.0]>tiny&]],
"33",engpos=Sort[Select[engpos0,Abs[Abs[#]-Sqrt[2.0]]>tiny&]],
"34",engpos1=Sort[Select[engpos0,Abs[Abs[#]-(1.0+Sqrt[5.0])/2.0]>tiny&]];
engpos=Sort[Select[engpos1,Abs[Abs[#]-(-1.0+Sqrt[5.0])/2.0]>tiny&]]
];

(* compute spacings *)
eng1=Drop[engpos,1];eng0=Drop[engpos,-1];

delta=eng1-eng0;
If[$DBG,Print[delta]];
delta=DeleteCases[delta,x_/;x==0.];
If[$DBG,Print[delta]];

(* Oganesyan r-values *)
rlist={};
Do[
rn=Min[delta[[ind]],delta[[ind-1]]]/Max[delta[[ind]],delta[[ind-1]]];
If[$DBG,Print[{ind,delta[[ind]],delta[[ind-1]],Min[delta[[ind]],delta[[ind-1]]],Max[delta[[ind]],delta[[ind-1]]],rn}]];
AppendTo[rlist,rn],
{ind,2,Length[delta]}
];
AppendTo[allrlist,rlist];
AppendTo[allstatrlist,{Mean[rlist],StandardDeviation[rlist]/Sqrt[Length[rlist]]}];

(* Shindou z-values *)
zlist=zcollec={};

(* lower boundary *)
AppendTo[zcollec,{delta[[1]],delta[[1]]+delta[[2]]}]; (*first delta*)
AppendTo[zcollec,{delta[[1]],delta[[2]],delta[[2]]+delta[[3]]}];(*second delta*)

(* bulk of data *)
Do[           
AppendTo[zcollec,{delta[[ind-2]]+delta[[ind-1]],delta[[ind-1]],delta[[ind]],delta[[ind]]+delta[[ind+1]]}];
,{ind,3,Length[delta]-1}];

(* upper boundary *)

AppendTo[zcollec,{delta[[Length[delta]-2]]+delta[[Length[delta]-1]],delta[[Length[delta]-1]],delta[[Length[delta]]]}]; (*2nd to last delta*)
AppendTo[zcollec,{delta[[Length[delta] - 1]] + delta[[Length[delta]]], delta[[Length[delta]]]}]; (*last delta*)

If[$DBG,Print["ncolec=",zcollec]];

Do[
denn=Min[zcollec[[k]]];
dennn=Min[DeleteCases[zcollec[[k]],denn]];
sk=denn/dennn;
If[$DBG,Print["sk=",sk]];
AppendTo[zlist,sk],
{k,1,Length[zcollec]}];
AppendTo[allzlist,zlist];
AppendTo[allstatzlist,{Mean[zlist],StandardDeviation[zlist]/Sqrt[Length[zlist]]}];

If[
Mod[ndone*100/lensamples,10]==0,
Print[{ifile,N[ndone/lensamples],N[ndone/lensamples],(AbsoluteTime[]-starttime),(AbsoluteTime[]-starttime)/N[ndone/lensamples]}]]
,{ifile,1,lensamples}
];

maxbins=200;

MeanEng=Mean[Flatten[allenglist]];
MedianEng=Median[Flatten[allenglist]];

(* r value *)
If[Length[Flatten[allrlist]]<=2,Continue[]];
 
(*binwidth=N[Min[0.1,Max[1/maxbins,maxbins/Length[Flatten[allrlist]]]]];
rawhistolist=HistogramList[Flatten[allrlist],{binwidth},"PDF"];bins=Drop[rawhistolist[[1]]-rawhistolist[[1,2]]/2,1];
histo=rawhistolist[[2]];
histolist=Transpose[{bins,histo}];

If[$LOG,Print[{dirname,
Length[Flatten[allenglist]],
TarEng,MeanEng,MedianEng,HubDis,RimDis,
Length[Flatten[allrlist]],Mean[Flatten[allrlist]]
}]];*)

AppendTo[allavglist,{dirname,
Length[Flatten[allenglist]],Length[allenglist],
Length[Flatten[allrlist]],Length[allrlist],
Chop[TarEng],MeanEng,MedianEng,Chop[HubDis],Chop[RimDis],
Mean[Flatten[allrlist]],
StandardDeviation[Flatten[allrlist]]/Sqrt[Length[Flatten[allrlist]]],
Mean[Transpose[allstatrlist][[1]]],
StandardDeviation[Transpose[allstatrlist][[1]]]/Sqrt[Length[Transpose[allstatrlist][[1]]]],
Mean[Flatten[allzlist]],
StandardDeviation[Flatten[allzlist]]/Sqrt[Length[Flatten[allzlist]]],
Mean[Transpose[allstatzlist][[1]]],
StandardDeviation[Transpose[allstatzlist][[1]]]/Sqrt[Length[Transpose[allstatzlist][[1]]]]
}];

If[$RVL,
rvlname=bottomdirname<>"_mE"<>StringDrop[ToString[PaddedForm[MeanEng,{6,3},NumberPadding->{"0","0"},NumberPoint->""]],1]<>"_"<>ToString[Min[allseeds]]<>"-"<>ToString[Max[allseeds]]<>"_"<>ToString[Length[allseeds]]<>"_rvl.txt";
PrintTemporary["     writing "<>rvlname];
Export[rvlname,Flatten[allrlist],"Table"]
];

If[
Mod[idir*100/lendirs,10]==0,
Print[{MM,N[idir/lendirs],N[idir/lendirs],(AbsoluteTime[]-starttimeD),(AbsoluteTime[]-starttimeD)/N[idir/lendirs]}]]
,{idir,1,lendirs}
];

(*Print[TableForm[allrmsdlist]];*)
SetDirectory[maindir];

(*Print[allavglist];*)

phasedata=Sort[Transpose[{
Transpose[allavglist][[7]],
Chop[Transpose[allavglist][[9]]],
Transpose[allavglist][[11]],
Transpose[allavglist][[12]],
Transpose[allavglist][[13]],
Transpose[allavglist][[14]],
Transpose[allavglist][[15]],
Transpose[allavglist][[16]],
Transpose[allavglist][[17]],
Transpose[allavglist][[18]]
}]];

phasedata=Sort[phasedata,#1[[2]]<#2[[2]] &];

SetDirectory[maindir];
Export["RZstat_E"<>ToString[Floor[TarEng*10]]<>"_M"<>MM<>"_rvl.txt",allavglist,"Table"];
Export["RZstat_E"<>ToString[Floor[TarEng*10]]<>"_M"<>MM<>"_rvl-phase.txt",phasedata,"Table"];

(*Export["Rstat_E"<>ToString[Floor[TarEng*10]]<>"_M"<>MM<>"_rvl-phase.csv",phasedata,"CSV"];*)

,
{iMM,Length[MMlist]}
];



