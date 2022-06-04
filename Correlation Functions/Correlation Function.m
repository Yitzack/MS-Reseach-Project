(*Import the data exported by Data Collector.nb or Data Collector.m*)
{parameters["24.1"],ImG1["24.1"],ImGC1["24.1"],ImGL1["24.1"],ImGQ1["24.1"],ImGV1["24.1"],ImG2["24.1"],ImGC2["24.1"],ImGL2["24.1"],ImGQ2["24.1"],ImGV2["24.1"],ImG3["24.1"],ImGC3["24.1"],ImGL3["24.1"],ImGQ3["24.1"],ImGV3["24.1"],ReGC1["24.1"],ReGL1["24.1"],ReGQ1["24.1"],ReGV1["24.1"],ReGC2["24.1"],ReGL2["24.1"],ReGQ2["24.1"],ReGV2["24.1"],ReGC3["24.1"],ReGL3["24.1"],ReGQ3["24.1"],ReGV3["24.1"],C0["24.1"],Fraction["24.1"],SpatialConstInter["24.1"],SpatialTotalInter["24.1"],SpatialNon["24.1"],EuclideanConstInter["24.1"],EuclideanTotalInter["24.1"],EuclideanNon["24.1"],EuclideanExt["24.1"]}=Import["Spectralcc24.97.1.m"];

(*Conversion between i and j to P and s or back*)
i1[s_,P_]:=10Sqrt[s+P^2]
j1[s_,P_]:=10(P-Sqrt[s+P^2])
i2[s_,P_]:=936/5+Sqrt[s+P^2]
j2[s_,P_]:=10(P-Sqrt[s+P^2])
P[i_,j_]:=If[j<=150,If[i<=208,(i+j)/10,i+j/10-187.2],(4i)/5]
s[i_,j_]:=If[j<=150,If[i<=208,(-i j)/50-j^2/100,-(j/5)(i-187.2)-j^2/100],If[j<=181,((j-151)/10)^2,If[j<=381,((j-181)/100+3)^2,((j-381)/10+5)^2]]]

(*Stiching the interpolation functions together*)
ImG[s_,P_,version_]:=If[s<-P^2,0,If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImG1[version][i1[s,P],j1[s,P]],ImG2[version][i2[s,P],j2[s,P]]],If[s>=0,ImG3[version][P,s],0]]];
ImGC[s_,P_,version_]:=If[s<-P^2,0,If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImGC1[version][i1[s,P],j1[s,P]],ImGC2[version][i2[s,P],j2[s,P]]],If[s>=0,ImGC3[version][P,s],0]]];
ImGL[s_,P_,version_]:=If[s<-P^2,0,If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImGL1[version][i1[s,P],j1[s,P]],ImGL2[version][i2[s,P],j2[s,P]]],If[s>=0,ImGL3[version][P,s],0]]];
ImGQ[s_,P_,version_]:=If[s<-P^2,0,If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImGQ1[version][i1[s,P],j1[s,P]],ImGQ2[version][i2[s,P],j2[s,P]]],If[s>=0,ImGQ3[version][P,s],0]]];
ImGV[s_,P_,version_]:=If[s<-P^2,0,If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImGV1[version][i1[s,P],j1[s,P]],ImGV2[version][i2[s,P],j2[s,P]]],If[s>=0,ImGV3[version][P,s],0]]];
ReGC[s_,P_,version_]:=If[s<-P^2,0,If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ReGC1[version][i1[s,P],j1[s,P]],ReGC2[version][i2[s,P],j2[s,P]]],If[s>=0,ReGC3[version][P,s],0]]];
ReGL[s_,P_,version_]:=If[s<-P^2,0,If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ReGL1[version][i1[s,P],j1[s,P]],ReGL2[version][i2[s,P],j2[s,P]]],If[s>=0,ReGL3[version][P,s],0]]];
ReGQ[s_,P_,version_]:=If[s<-P^2,0,If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ReGQ1[version][i1[s,P],j1[s,P]],ReGQ2[version][i2[s,P],j2[s,P]]],If[s>=0,ReGQ3[version][P,s],0]]] ;
ReGV[s_,P_,version_]:=If[s<-P^2,0,If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ReGV1[version][i1[s,P],j1[s,P]],ReGV2[version][i2[s,P],j2[s,P]]],If[s>=0,ReGV3[version][P,s],0]]];

(*Assemble the contributions to the spectral function into a spetcral function*)
\[Sigma]Non[s_,P_,version_]:=-(3/\[Pi])ImG[s,P,version]
\[Sigma]C[s_,P_,version_]:=3/\[Pi] C0[version]parameters[version][[1]] Im[(1/4 (ReGC[s,P,version]+I ImGC[s,P,version])^2)/(1-(ReGV[s,P,version]+I ImGV[s,P,version]))]
\[Sigma]L[s_,P_,version_]:=3/\[Pi] C0[version] parameters[version][[1]]Im[(ReGL[s,P,version]+I ImGL[s,P,version])^2/(1-(ReGV[s,P,version]+I ImGV[s,P,version]))]
\[Sigma]Q[s_,P_,version_]:=3/\[Pi] C0[version] parameters[version][[1]]Im[(ReGQ[s,P,version]+I ImGQ[s,P,version])^2/(1-(ReGV[s,P,version]+I ImGV[s,P,version]))]
\[Sigma]ConstInter[s_,P_,version_]:=8\[Sigma]C[s,P,version]
\[Sigma]TotalInter[s_,P_,version_]:=\[Sigma]C[s,P,version]+\[Sigma]L[s,P,version]+\[Sigma]Q[s,P,version]
\[Sigma]Const[s_,P_,version_]:=\[Sigma]Non[s,P,version]+\[Sigma]ConstInter[s,P,version]
\[Sigma]Total[s_,P_,version_]:=\[Sigma]Non[s,P,version]+\[Sigma]TotalInter[s,P,version]

(*Correlation kernels*)
Kernel["Cutoff"][s_,z_,P0_]:=-((z (-2 P0 z (12852-5868 s z^2+1737 s^2 z^4-78 s^3 z^6+P0^8 z^8+s^4 z^8+3 P0^4 z^4 (579+26 s z^2+2 s^2 z^4)+P0^6 (78 z^6+4 s z^8)+P0^2 (5868 z^2+2418 s z^4-78 s^2 z^6+4 s^3 z^8)) Cos[P0 z]+(-5400+33948 s z^2-16074 s^2 z^4+2301 s^3 z^6-88 s^4 z^8+P0^10 z^10+s^5 z^10+P0^8 z^8 (84+5 s z^2)+P0^6 z^6 (2037+164 s z^2+10 s^2 z^4)+P0^4 z^4 (10782+3975 s z^2-12 s^2 z^4+10 s^3 z^6)+P0^2 z^2 (27756-4716 s z^2+4239 s^2 z^4-180 s^3 z^6+5 s^4 z^8)) Sin[P0 z]))/(P0^12 z^12+6 P0^10 z^10 (15+s z^2)+3 P0^8 z^8 (819+90 s z^2+5 s^2 z^4)+(-36+216 s z^2-45 s^2 z^4+s^3 z^6)^2+4 P0^6 z^6 (4878+1593 s z^2+45 s^2 z^4+5 s^3 z^6)+3 P0^4 z^4 (16632+6120 s z^2+2610 s^2 z^4-60 s^3 z^6+5 s^4 z^8)+6 P0^2 z^2 (2592+12312 s z^2-3060 s^2 z^4+1062 s^3 z^6-45 s^4 z^8+s^5 z^10)))
Kernel["Lorentz"][s_,z_]:=\[Pi]/(2Sqrt[s]) Exp[-Sqrt[s]z]
Kernel["General"][s_,z_,P_]:=Cos[P z]/(s+P^2)
Kernel["General Neg"][E_,z_,P_]:=2 Cos[P z]/E
Kernel["Euclidean"][s_,P_,\[Tau]_,T_]:=Cosh[Sqrt[s+P^2](\[Tau]-1/(2T))]/(Sinh[Sqrt[s+P^2]/(2T)]2Sqrt[s+P^2])

(*How to divide up the spectral function for doing the Euclidean correlation function. Space-like mesons are not considered as they are a small contribution.*)
Dimensions[sTab=Join[Range[0,3,.25],{3.6,7,23.5}]^2]

(*Launch parallel kernels. The default is half of your logical CPU cores. I've directed it to launch 22, 2 short of my 24 logical cores*)
LaunchKernels[22]

(*Euclidean correlation function for the interacting spectral function that uses the constant trace everywhere. The backup was to save work so that it isn't lost if my computer locks up.*)
EuclideanConstInter["24.1"]=Table[If[\[Tau]==0,Print[{Now,\[Tau],P}]];{\[Tau],P,ParallelTable[{sTab[[Mod[i,15]+1]],NIntegrate[Kernel["Euclidean"][s,P,{.194,.258,.32,.4}^-1 \[Tau],{.194,.258,.32,.4}]\[Sigma]ConstInter[s,P,"24.1"],{s,sTab[[Mod[i,15]+1]],sTab[[Mod[i,15]+2]]},Method->"GaussKronrodRule"]},{i,0,15},Method->"FinestGrained"]},{P,0,3,.5},{\[Tau],0,1/2,1/40}];
Export["Correlation Function Backup.m",{EuclideanConstInter["24.1"]}]
NotebookSave[];

(*Euclidean correlation function for the interacting spectral function that uses the complete trace everywhere.*)
EuclideanTotalInter["24.1"]=Table[If[\[Tau]==0,Print[{Now,\[Tau],P}]];{\[Tau],P,ParallelTable[{sTab[[Mod[i,15]+1]],NIntegrate[Kernel["Euclidean"][s,P,{.194,.258,.32,.4}^-1 \[Tau],{.194,.258,.32,.4}]\[Sigma]TotalInter[s,P,"24.1"],{s,sTab[[Mod[i,15]+1]],sTab[[Mod[i,15]+2]]},Method->"GaussKronrodRule"]},{i,0,15},Method->"FinestGrained"]},{P,0,3,.5},{\[Tau],0,1/2,1/40}];
Export["Correlation Function Backup.m",{EuclideanConstInter["24.1"],EuclideanTotalInter["24.1"]}]
NotebookSave[];

(*Euclidean correlation function for the non-interacting spectral function.*)
EuclideanNon["24.1"]=Table[If[\[Tau]==0,Print[{Now,\[Tau],P}]];{\[Tau],P,ParallelTable[{sTab[[Mod[i,15]+1]],NIntegrate[Kernel["Euclidean"][s,P,{.194,.258,.32,.4}^-1 \[Tau],{.194,.258,.32,.4}]\[Sigma]Non[s,P,"24.1"],{s,sTab[[Mod[i,15]+1]],sTab[[Mod[i,15]+2]]},Method->"GaussKronrodRule"]},{i,0,15},Method->"FinestGrained"]},{P,0,3,.5},{\[Tau],0,1/2,1/40}];
Export["Correlation Function Backup.m",{EuclideanConstInter["24.1"],EuclideanTotalInter["24.1"],EuclideanNon["24.1"]}]
NotebookSave[];

(*Euclidean correlation function for the non-interacting spectral function using an extension to s=400^2 GeV^2.*)
EuclideanExt["24.1"]=ParallelTable[If[\[Tau]==0,Print[{Now,\[Tau],P}]];{\[Tau],P,NIntegrate[Kernel["Euclidean"][s,P,{.194,.258,.32,.4}^-1 \[Tau],{.194,.258,.32,.4}] (3s)/(8\[Pi]^2) Sqrt[1-(4parameters["24.1"][[3]]^2)/s],{s,552.25,400^2},Method->"GaussKronrodRule"]},{P,0,3,.5},{\[Tau],0,1/2,1/40},Method->"FinestGrained"];
Export["Correlation Function Backup.m",{EuclideanConstInter["24.1"],EuclideanTotalInter["24.1"],EuclideanNon["24.1"],EuclideanExt["24.1"]}]
NotebookSave[];

(*Store the results*)
Export["Spectralcc24.97.1.m",{parameters["24.1"],ImG1["24.1"],ImGC1["24.1"],ImGL1["24.1"],ImGQ1["24.1"],ImGV1["24.1"],ImG2["24.1"],ImGC2["24.1"],ImGL2["24.1"],ImGQ2["24.1"],ImGV2["24.1"],ImG3["24.1"],ImGC3["24.1"],ImGL3["24.1"],ImGQ3["24.1"],ImGV3["24.1"],ReGC1["24.1"],ReGL1["24.1"],ReGQ1["24.1"],ReGV1["24.1"],ReGC2["24.1"],ReGL2["24.1"],ReGQ2["24.1"],ReGV2["24.1"],ReGC3["24.1"],ReGL3["24.1"],ReGQ3["24.1"],ReGV3["24.1"],C0["24.1"],Fraction["24.1"],SpatialConstInter["24.1"],SpatialTotalInter["24.1"],SpatialNon["24.1"],EuclideanConstInter["24.1"],EuclideanTotalInter["24.1"],EuclideanNon["24.1"],EuclideanExt["24.1"]}]

(*How to divide up the time-like mesons*)
Dimensions[sTab=Join[Range[0,7,.25]^2,{23.5^2}]]

(*Spatial correlation function for the interacting spectral function that uses the constant trace*)
SpatialConstInter["24.1"]=Table[Print[{Now,z}];(*check progress*)
{z,
ParallelTable[{sTab[[Mod[i,29]+1]],Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]],NIntegrate[Kernel["General"][s,z,P]{\[Sigma]TotalConst[s,P,"24.1"],\[Sigma]TotalConst[s,0,"24.1"]},{s,sTab[[Mod[i,29]+1]],sTab[[Mod[i,29]+2]]},{P,Join[{0},Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&]][[Floor[i/29]+1]],Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]]},Method->"GaussKronrodRule"]},{i,0,29(Dimensions[Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&]][[1]])-1}],(*time-like contributions*)
ParallelTable[{sTab[[Mod[i,29]+1]],Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]],NIntegrate[Kernel["Cutoff"][s,z,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]]]{\[Sigma]TotalConst[s,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]],"24.1"],\[Sigma]TotalConst[s,0,"24.1"]},{s,sTab[[Mod[i,29]+1]],sTab[[Mod[i,29]+2]]},Method->"GaussKronrodRule"]},{i,0,29(Dimensions[Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&]][[1]])-1}],(*time-like cutoff contributions*)
NIntegrate[Kernel["General Neg"][\[Omega],z,P]\[Sigma]ConstInter[\[Omega]^2-P^2,P,"24.1"],{P,0,15},{\[Omega],0,P},Method->"GaussKronrodRule"]+NIntegrate[Kernel["General Neg"][\[Omega],z,P]\[Sigma]ConstInter[\[Omega]^2-P^2,P,"24.1"],{P,15,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]},{\[Omega],P-15,P},Method->"GaussKronrodRule"]+NIntegrate[2\[Omega] Kernel["Cutoff"][\[Omega]^2-Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]^2,z,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]]\[Sigma]ConstInter[\[Omega]^2-Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]^2,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]],"24.1"],{\[Omega],Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]-15,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]},Method->"GaussKronrodRule"],(*space-like contributions w/ cutoff*)
ParallelTable[{sTab[[i]],NIntegrate[Kernel["Lorentz"][s,z]\[Sigma]TotalConst[s,0,"24.1"],{s,sTab[[i]],sTab[[i+1]]},Method->"GaussKronrodRule"]},{i,29}](*P=0 contributions for comparision with the momentum dependent result and check to on the numerical result*)},{z,.25,6.25,.25}];
Export["Correlation Function Backup.m",{EuclideanConstInter["24.1"],EuclideanTotalInter["24.1"],EuclideanNon["24.1"],EuclideanExt["24.1"],SpatialConstInter["24.1"]}]

(*Spatial correlation function for the interacting spectral function that uses the complete trace*)
SpatialTotalInter["24.1"]=Table[Print[{Now,z}];
{z,
ParallelTable[{sTab[[Mod[i,29]+1]],Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]],NIntegrate[Kernel["General"][s,z,P]{\[Sigma]TotalInter[s,P,"24.1"],\[Sigma]TotalInter[s,0,"24.1"]},{s,sTab[[Mod[i,29]+1]],sTab[[Mod[i,29]+2]]},{P,Join[{0},Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&]][[Floor[i/29]+1]],Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]]},Method->"GaussKronrodRule"]},{i,0,29(Dimensions[Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&]][[1]])-1}],
ParallelTable[{sTab[[Mod[i,29]+1]],Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]],NIntegrate[Kernel["Cutoff"][s,z,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]]]{\[Sigma]TotalInter[s,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[Floor[i/29]+1]],"24.1"],\[Sigma]TotalInter[s,0,"24.1"]},{s,sTab[[Mod[i,29]+1]],sTab[[Mod[i,29]+2]]},Method->"GaussKronrodRule"]},{i,0,29(Dimensions[Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&]][[1]])-1}],
NIntegrate[Kernel["General Neg"][\[Omega],z,P]\[Sigma]TotalInterRev[\[Omega]^2-P^2,P,"24.1"],{P,0,15},{\[Omega],0,P},Method->"GaussKronrodRule"]+NIntegrate[Kernel["General Neg"][\[Omega],z,P]\[Sigma]TotalInterRev[\[Omega]^2-P^2,P,"24.1"],{P,15,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]},{\[Omega],P-15,P},Method->"GaussKronrodRule"]+NIntegrate[2\[Omega] Kernel["Cutoff"][\[Omega]^2-Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]^2,z,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]]\[Sigma]TotalInterRev[\[Omega]^2-Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]^2,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]],"24.1"],{\[Omega],Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]-15,Select[Range[\[Pi]/(2z),(2 32\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=32&][[-1]]},Method->"GaussKronrodRule"],
ParallelTable[{sTab[[i]],NIntegrate[Kernel["Lorentz"][s,z]\[Sigma]TotalInter[s,0,"24.1"],{s,sTab[[i]],sTab[[i+1]]},Method->"GaussKronrodRule"]},{i,29}]},{z,.25,6.25,.25}];
Export["Correlation Function Backup.m",{EuclideanConstInter["24.1"],EuclideanTotalInter["24.1"],EuclideanNon["24.1"],EuclideanExt["24.1"],SpatialConstInter["24.1"],SpatialTotalInter["24.1"]}]

(*Spatial correlation function for the non-interacting spectral function*)
SpatialNon["24.1"]=Table[Print[{Now,z}];
{z,
ParallelTable[{sTab[[Mod[i,29]+1]],Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[Floor[i/29]+1]],NIntegrate[Kernel["General"][s,z,P]{\[Sigma]Non[s,P,"24.1"],\[Sigma]Non[s,0,"24.1"]},{s,sTab[[Mod[i,29]+1]],sTab[[Mod[i,29]+2]]},{P,Join[{0},Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&]][[Floor[i/29]+1]],Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[Floor[i/29]+1]]},Method->"GaussKronrodRule"]},{i,0,29(Dimensions[Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&]][[1]])-1}],
ParallelTable[{sTab[[Mod[i,29]+1]],Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[Floor[i/29]+1]],NIntegrate[Kernel["Cutoff"][s,z,Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[Floor[i/29]+1]]]{\[Sigma]Non[s,Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[Floor[i/29]+1]],"24.1"],\[Sigma]Non[s,0,"24.1"]},{s,sTab[[Mod[i,29]+1]],sTab[[Mod[i,29]+2]]},Method->"GaussKronrodRule"]},{i,0,29(Dimensions[Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&]][[1]])-1}],
NIntegrate[Kernel["General Neg"][\[Omega],z,P]\[Sigma]Non[\[Omega]^2-P^2,P,"24.1"],{P,0,15},{\[Omega],0,P},Method->"GaussKronrodRule"]+NIntegrate[Kernel["General Neg"][\[Omega],z,P]\[Sigma]Non[\[Omega]^2-P^2,P,"24.1"],{P,15,Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[-1]]},{\[Omega],P-15,P},Method->"GaussKronrodRule"]+NIntegrate[2\[Omega] Kernel["Cutoff"][\[Omega]^2-Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[-1]]^2,z,Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[-1]]]\[Sigma]Non[\[Omega]^2-Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[-1]]^2,Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[-1]],"24.1"],{\[Omega],Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[-1]]-15,Select[Range[\[Pi]/(2z),(2 70.4\[Pi])/z+\[Pi]/(2z),\[Pi]/z],#<=70.4&][[-1]]},Method->"GaussKronrodRule"],
ParallelTable[{sTab[[i]],NIntegrate[Kernel["Lorentz"][s,z]\[Sigma]Non[s,0,"24.1"],{s,sTab[[i]],sTab[[i+1]]},Method->"GaussKronrodRule"]},{i,29}]},{z,.25,6.25,.25}];
Export["Correlation Function Backup.m",{EuclideanConstInter["24.1"],EuclideanTotalInter["24.1"],EuclideanNon["24.1"],EuclideanExt["24.1"],SpatialConstInter["24.1"],SpatialTotalInter["24.1"],SpatialNon["24.1"]}]

(*Export results to mathematica file*)
Export["Spectralcc24.97.1.m",{parameters["24.1"],ImG1["24.1"],ImGC1["24.1"],ImGL1["24.1"],ImGQ1["24.1"],ImGV1["24.1"],ImG2["24.1"],ImGC2["24.1"],ImGL2["24.1"],ImGQ2["24.1"],ImGV2["24.1"],ImG3["24.1"],ImGC3["24.1"],ImGL3["24.1"],ImGQ3["24.1"],ImGV3["24.1"],ReGC1["24.1"],ReGL1["24.1"],ReGQ1["24.1"],ReGV1["24.1"],ReGC2["24.1"],ReGL2["24.1"],ReGQ2["24.1"],ReGV2["24.1"],ReGC3["24.1"],ReGL3["24.1"],ReGQ3["24.1"],ReGV3["24.1"],C0["24.1"],Fraction["24.1"],SpatialConstInter["24.1"],SpatialTotalInter["24.1"],SpatialNon["24.1"],EuclideanConstInter["24.1"],EuclideanTotalInter["24.1"],EuclideanNon["24.1"],EuclideanExt["24.1"]}]

CloseKernels[]
NotebookSave[]
Exit[]
