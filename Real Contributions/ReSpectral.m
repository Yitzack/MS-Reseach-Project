assignment = ToExpression[$CommandLine[[4]]];
Temp = $CommandLine[[5]];
version = $CommandLine[[6]];
Master = Import["~/data/Tables/Spectralcc24.97.P_Dep4."<>Temp<>"_"<>version<>".csv"];
If[Master==$Failed,Exit[]];
dim=Dimensions[Master][[1]];
parameters = Select[Master, Dimensions[#][[1]] == 5 &];
Master = Select[Master, Dimensions[#][[1]] != 5 &];
upper =Master[[-1,1]]; (*upper edge of the s>=0 table*)

i1[s_,P_]:=10Sqrt[s+P^2]
j1[s_,P_]:=10(P-Sqrt[s+P^2])
i2[s_,P_]:=936/5+Sqrt[s+P^2]
j2[s_,P_]:=10(P-Sqrt[s+P^2])
P[i_,j_]:=If[j<=150,If[0<=i&&i<=208,(i+j)/10,i+j/10-187.2],If[i>=0,(4i)/5,(Mod[i-1,7]/8-7/8-Floor[i/7])*.8]]
s[i_,j_]:=If[j<=150,If[0<=i&&i<=208,(-i j)/50-j^2/100,-(j/5)(i-187.2)-j^2/100],If[j<=181,((j-151)/10)^2,If[j<=381,((j-181)/100+3)^2,((j-381)/10+5)^2]]]
ReG[f_,i_,j_]:=(NIntegrate[-f[y,P[i,j]]/(s[i,j]-y),{y,Piecewise[{{-P[i,j]^2,P[i,j]<15}},-30P[i,j]+225],Max[Piecewise[{{-P[i,j]^2,P[i,j]<15}},-30P[i,j]+225],s[i,j]-.001]}]+NIntegrate[(f[s[i,j],P[i,j]]-f[y,P[i,j]])/(s[i,j]-y),{y,s[i,j]-.001,s[i,j]}]+NIntegrate[(f[s[i,j],P[i,j]]-f[y,P[i,j]])/(s[i,j]-y),{y,s[i,j],s[i,j]+.001}]-Re[Integrate[(f[s[i,j],P[i,j]])/(s[i,j]-y),{y,s[i,j]-.001,s[i,j]+.001},PrincipalValue->True]]+NIntegrate[-f[y,P[i,j]]/(s[i,j]-y),{y,s[i,j]+.001,729.096309407015895}])/\[Pi]
(*ReG[f_,i_,j_]:=(NIntegrate[-f[y,P[i,j]]/(s[i,j]-y),{y,0,Max[0,s[i,j]-.001]}]+NIntegrate[(f[s[i,j],P[i,j]]-f[y,P[i,j]])/(s[i,j]-y),{y,Max[0,s[i,j]-.001],s[i,j]}]+NIntegrate[(f[s[i,j],P[i,j]]-f[y,P[i,j]])/(s[i,j]-y),{y,s[i,j],s[i,j]+.001}]-Re[Integrate[(f[s[i,j],P[i,j]])/(s[i,j]-y),{y,Max[0,s[i,j]-.001],s[i,j]+.001},PrincipalValue->True]]+NIntegrate[-f[y,P[i,j]]/(s[i,j]-y),{y,s[i,j]+.001,729.096309407015895}])/\[Pi]*)

Dimensions[Sub1=Select[Master[[;;-2,{1,2,5,6,7,8,9}]],#[[2]]<=150&&#[[1]]<=208&]]
Dimensions[Sub2=Select[Master[[;;-2,{1,2,5,6,7,8,9}]],#[[2]]<=150&&#[[1]]>=208&]]
Dimensions[Sub3=Select[Master,#[[2]]>150&&#[[1]]<=upper&][[All,3;;]]]

If[Temp != "0",
Sub1=Table[If[MachineNumberQ[N[Sub1[[i,j]]]],N[Sub1[[i,j]]],0],{i,31559},{j,7}];
dim=Dimensions[Sub2][[1]];
Sub2=Table[If[MachineNumberQ[N[Sub2[[i,j]]]]&&Abs[Sub2[[i,j]]]>$MinMachineNumber,N[Sub2[[i,j]]],0],{i,dim},{j,7}],
Sub1=Flatten[Table[{i,j,0,0,0,0,0},{i,0,208},{j,0,150}],1];
Sub2=Flatten[Table[{i,j,0,0,0,0,0},{i,208,258},{j,0,150}],1]];
dim=Dimensions[Sub3][[1]];
Sub3=Table[If[MachineNumberQ[N[Sub3[[i,j]]]],N[Sub3[[i,j]]],0],{i,dim},{j,7}];

Off[Interpolation::mspl]
If[Sub1=={},ImG1=0&,ImG1=Interpolation[Transpose[{Transpose[{Sub1[[All,1]],Sub1[[All,2]]}],Sub1[[All,3]]}],Method->"Spline",InterpolationOrder->1]];
If[Sub1=={},ImGvC1=0&,ImGvC1=Interpolation[Transpose[{Transpose[{Sub1[[All,1]],Sub1[[All,2]]}],Sub1[[All,4]]}],Method->"Spline",InterpolationOrder->1]];
If[Sub1=={},ImGvL1=0&,ImGvL1=Interpolation[Transpose[{Transpose[{Sub1[[All,1]],Sub1[[All,2]]}],Sub1[[All,5]]}],Method->"Spline",InterpolationOrder->1]];
If[Sub1=={},ImGvQ1=0&,ImGvQ1=Interpolation[Transpose[{Transpose[{Sub1[[All,1]],Sub1[[All,2]]}],Sub1[[All,6]]}],Method->"Spline",InterpolationOrder->1]];
If[Sub1=={},ImGV1=0&,ImGV1=Interpolation[Transpose[{Transpose[{Sub1[[All,1]],Sub1[[All,2]]}],Sub1[[All,7]]}],Method->"Spline",InterpolationOrder->1]];
If[Sub2=={},ImG2=0&,ImG2=Interpolation[Transpose[{Transpose[{Sub2[[All,1]],Sub2[[All,2]]}],Sub2[[All,3]]}],Method->"Spline",InterpolationOrder->1]];
If[Sub2=={},ImGvC2=0&,ImGvC2=Interpolation[Transpose[{Transpose[{Sub2[[All,1]],Sub2[[All,2]]}],Sub2[[All,4]]}],Method->"Spline",InterpolationOrder->1]];
If[Sub2=={},ImGvL2=0&,ImGvL2=Interpolation[Transpose[{Transpose[{Sub2[[All,1]],Sub2[[All,2]]}],Sub2[[All,5]]}],Method->"Spline",InterpolationOrder->1]];
If[Sub2=={},ImGvQ2=0&,ImGvQ2=Interpolation[Transpose[{Transpose[{Sub2[[All,1]],Sub2[[All,2]]}],Sub2[[All,6]]}],Method->"Spline",InterpolationOrder->1]];
If[Sub2=={},ImGV2=0&,ImGV2=Interpolation[Transpose[{Transpose[{Sub2[[All,1]],Sub2[[All,2]]}],Sub2[[All,7]]}],Method->"Spline",InterpolationOrder->1]];
ImG3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,3]]}],Method->"Spline"];
ImGvC3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,4]]}],Method->"Spline"];
ImGvL3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,5]]}],Method->"Spline"];
ImGvQ3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,6]]}],Method->"Spline"];
ImGV3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,7]]}],Method->"Spline"];
ImG[s_,P_]:=If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImG1[i1[s,P],j1[s,P]],ImG2[i2[s,P],j2[s,P]]],If[s>=0,ImG3[P,s],0]];
ImGvC[s_,P_]:=If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImGvC1[i1[s,P],j1[s,P]],ImGvC2[i2[s,P],j2[s,P]]],If[s>=0,ImGvC3[P,s],0]];
ImGvL[s_,P_]:=If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImGvL1[i1[s,P],j1[s,P]],ImGvL2[i2[s,P],j2[s,P]]],If[s>=0,ImGvL3[P,s],0]];
ImGvQ[s_,P_]:=If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImGvQ1[i1[s,P],j1[s,P]],ImGvQ2[i2[s,P],j2[s,P]]],If[s>=0,ImGvQ3[P,s],0]];
ImGV[s_,P_]:=If[s<0&&P<Sqrt[s+P^2]+15,If[s+P^2<=20.8^2,ImGV1[i1[s,P],j1[s,P]],ImGV2[i2[s,P],j2[s,P]]],If[s>=0,ImGV3[P,s],0]];

Off[NIntegrate::izero,NIntegrate::ncvb,NIntegrate::slwcon,NIntegrate::precw,NIntegrate::eincr,NIntegrate::zeroregion,Integrate::idiv];

If[assignment>=0&&Temp!="0",
	ReGTab1=Table[{assignment,j,ReG[ImGvC,assignment,j],ReG[ImGvL,assignment,j],ReG[ImGvQ,assignment,j],ReG[ImGV,assignment,j]},{j,0,150}];,
	ReGTab1={};
]

ReGTab2=If[assignment<=upper,
	Sub3=Select[Sub3,Abs[#[[1]]-P[assignment,200]]<4&];
	ImG3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,3]]}],Method->"Spline"];
	ImGvC3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,4]]}],Method->"Spline"];
	ImGvL3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,5]]}],Method->"Spline"];
	ImGvQ3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,6]]}],Method->"Spline"];
	ImGV3=Interpolation[Transpose[{Transpose[{Sub3[[All,1]],Sub3[[All,2]]}],Sub3[[All,7]]}],Method->"Spline"];
	Table[{assignment,j,ReG[ImGvC,assignment,j],ReG[ImGvL,assignment,j],ReG[ImGvQ,assignment,j],ReG[ImGV,assignment,j]},{j,151,566}]
	,{}
]

ReGTab=Join[ReGTab1,ReGTab2]

Export["~/data/ReSpectralcc24.97.P_Dep4."<>Temp<>"_"<>version<>"."<>ToString[assignment]<>".csv",ReGTab]

Exit[]
