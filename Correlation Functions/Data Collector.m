(*Import contributions to the spectral function*)
Dimensions[ImMaster=Import["~/Physics/Research/MS/data/Spectralcc22.97.1.csv"]]
Dimensions[ReMaster=Import["~/Physics/Research/MS/data/ReSpectralcc22.97.1.csv"]]

(*Seperate the input parameters from the imaginary contributions*)
parameters["curr"]=Select[ImMaster, Dimensions[#]!={9}&][[1]]
Dimensions[ImMaster=Select[ImMaster, Dimensions[#]=={9}&]]

(*Seperate out the inputs into the respective tables for interpolation*)
Dimensions[Sub1=Select[ImMaster[[All, {1, 2, 5, 6, 7, 8, 9}]], #[[2]]<=150&&#[[1]]<=208&]]
Dimensions[Sub2=Select[ImMaster[[All, {1, 2, 5, 6, 7, 8, 9}]], #[[2]]<=150&&#[[1]]>=208&]]
Dimensions[Sub3=Select[ImMaster, #[[2]]>150&&#[[1]]<=40&][[All, 3;;]]]
Dimensions[Sub4=Select[ReMaster, #[[2]]<=150&&#[[1]]<=208&]]
Dimensions[Sub5=Select[ReMaster, #[[2]]<=150&&#[[1]]>=208&]]
Dimensions[Sub6=Select[ReMaster, #[[2]]>150&&#[[1]]<=126&]]

(*Force all numbers in the imaginary contributions to be machine real numbers*)
dim=Dimensions[Sub1][[1]];
Sub1=Table[If[MachineNumberQ[N[Sub1[[i, j]]]], N[Sub1[[i, j]]], 0], {i, dim}, {j, 7}];
dim=Dimensions[Sub2][[1]];
Sub2=Table[If[MachineNumberQ[N[Sub2[[i, j]]]]&&Abs[Sub2[[i, j]]]>$MinMachineNumber, N[Sub2[[i, j]]], 0], {i, dim}, {j, 7}];
dim=Dimensions[Sub3][[1]];
Sub3=Table[If[MachineNumberQ[N[Sub3[[i, j]]]], N[Sub3[[i, j]]], 0], {i, dim}, {j, 7}];

(*Create the interpolation functions*)
ImG1["curr"]=Interpolation[Transpose[{Transpose[{Sub1[[All, 1]], Sub1[[All, 2]]}], Sub1[[All, 3]]}], Method->"Spline"];
ImGC1["curr"]=Interpolation[Transpose[{Transpose[{Sub1[[All, 1]], Sub1[[All, 2]]}], Sub1[[All, 4]]}], Method->"Spline"];
ImGL1["curr"]=Interpolation[Transpose[{Transpose[{Sub1[[All, 1]], Sub1[[All, 2]]}], Sub1[[All, 5]]}], Method->"Spline"];
ImGQ1["curr"]=Interpolation[Transpose[{Transpose[{Sub1[[All, 1]], Sub1[[All, 2]]}], Sub1[[All, 6]]}], Method->"Spline"];
ImGV1["curr"]=Interpolation[Transpose[{Transpose[{Sub1[[All, 1]], Sub1[[All, 2]]}], Sub1[[All, 7]]}], Method->"Spline"];
ImG2["curr"]=Interpolation[Transpose[{Transpose[{Sub2[[All, 1]], Sub2[[All, 2]]}], Sub2[[All, 3]]}], Method->"Spline"];
ImGC2["curr"]=Interpolation[Transpose[{Transpose[{Sub2[[All, 1]], Sub2[[All, 2]]}], Sub2[[All, 4]]}], Method->"Spline"];
ImGL2["curr"]=Interpolation[Transpose[{Transpose[{Sub2[[All, 1]], Sub2[[All, 2]]}], Sub2[[All, 5]]}], Method->"Spline"];
ImGQ2["curr"]=Interpolation[Transpose[{Transpose[{Sub2[[All, 1]], Sub2[[All, 2]]}], Sub2[[All, 6]]}], Method->"Spline"];
ImGV2["curr"]=Interpolation[Transpose[{Transpose[{Sub2[[All, 1]], Sub2[[All, 2]]}], Sub2[[All, 7]]}], Method->"Spline"];
ImG3["curr"]=Interpolation[Transpose[{Transpose[{Sub3[[All, 1]], Sub3[[All, 2]]}], Sub3[[All, 3]]}], Method->"Spline"];
ImGC3["curr"]=Interpolation[Transpose[{Transpose[{Sub3[[All, 1]], Sub3[[All, 2]]}], Sub3[[All, 4]]}], Method->"Spline"];
ImGL3["curr"]=Interpolation[Transpose[{Transpose[{Sub3[[All, 1]], Sub3[[All, 2]]}], Sub3[[All, 5]]}], Method->"Spline"];
ImGQ3["curr"]=Interpolation[Transpose[{Transpose[{Sub3[[All, 1]], Sub3[[All, 2]]}], Sub3[[All, 6]]}], Method->"Spline"];
ImGV3["curr"]=Interpolation[Transpose[{Transpose[{Sub3[[All, 1]], Sub3[[All, 2]]}], Sub3[[All, 7]]}], Method->"Spline"];
ReGC1["curr"]=Interpolation[Transpose[{Transpose[{Sub4[[All, 1]], Sub4[[All, 2]]}], Sub4[[All, 3]]}], Method->"Spline"];
ReGL1["curr"]=Interpolation[Transpose[{Transpose[{Sub4[[All, 1]], Sub4[[All, 2]]}], Sub4[[All, 4]]}], Method->"Spline"];
ReGQ1["curr"]=Interpolation[Transpose[{Transpose[{Sub4[[All, 1]], Sub4[[All, 2]]}], Sub4[[All, 5]]}], Method->"Spline"];
ReGV1["curr"]=Interpolation[Transpose[{Transpose[{Sub4[[All, 1]], Sub4[[All, 2]]}], Sub4[[All, 6]]}], Method->"Spline"];
ReGC2["curr"]=Interpolation[Transpose[{Transpose[{Sub5[[All, 1]], Sub5[[All, 2]]}], Sub5[[All, 3]]}], Method->"Spline"];
ReGL2["curr"]=Interpolation[Transpose[{Transpose[{Sub5[[All, 1]], Sub5[[All, 2]]}], Sub5[[All, 4]]}], Method->"Spline"];
ReGQ2["curr"]=Interpolation[Transpose[{Transpose[{Sub5[[All, 1]], Sub5[[All, 2]]}], Sub5[[All, 5]]}], Method->"Spline"];
ReGV2["curr"]=Interpolation[Transpose[{Transpose[{Sub5[[All, 1]], Sub5[[All, 2]]}], Sub5[[All, 6]]}], Method->"Spline"];
dim=Dimensions[Sub6][[1]];ReGC3["curr"]=Interpolation[Table[{{P[Sub6[[i, 1]], Sub6[[i, 2]]], s[Sub6[[i, 1]], Sub6[[i, 2]]]}, Sub6[[i, 3]]}, {i, dim}], Method->"Spline"];
ReGL3["curr"]=Interpolation[Table[{{P[Sub6[[i, 1]], Sub6[[i, 2]]], s[Sub6[[i, 1]], Sub6[[i, 2]]]}, Sub6[[i, 4]]}, {i, dim}], Method->"Spline"];
ReGQ3["curr"]=Interpolation[Table[{{P[Sub6[[i, 1]], Sub6[[i, 2]]], s[Sub6[[i, 1]], Sub6[[i, 2]]]}, Sub6[[i, 5]]}, {i, dim}], Method->"Spline"];
ReGV3["curr"]=Interpolation[Table[{{P[Sub6[[i, 1]], Sub6[[i, 2]]], s[Sub6[[i, 1]], Sub6[[i, 2]]]}, Sub6[[i, 6]]}, {i, dim}], Method->"Spline"];

(*Record the vacuum coupling constant and the fraction of the coupling constant in-medium*)
C0["curr"]={116.253`, 65.8549`, 38.4541`, 50.3627}[[2]]
Fraction["curr"]=parameters["curr"][[1]]

(*A previous iteration may have correlation functions that are close or otherwise still good,  Import those*)
{parameters["22.1"], ImG1["22.1"], ImGC1["22.1"], ImGL1["22.1"], ImGQ1["22.1"], ImGV1["22.1"], ImG2["22.1"], ImGC2["22.1"], ImGL2["22.1"], ImGQ2["22.1"], ImGV2["22.1"], ImG3["22.1"], ImGC3["22.1"], ImGL3["22.1"], ImGQ3["22.1"], ImGV3["22.1"], ReGC1["22.1"], ReGL1["22.1"], ReGQ1["22.1"], ReGV1["22.1"], ReGC2["22.1"], ReGL2["22.1"], ReGQ2["22.1"], ReGV2["22.1"], ReGC3["22.1"], ReGL3["22.1"], ReGQ3["22.1"], ReGV3["22.1"], C0["22.1"], Fraction["22.1"], SpatialConstInter["22.1"], SpatialTotalInter["22.1"], SpatialNon["22.1"], EuclideanConstInter["22.1"], EuclideanTotalInter["22.1"], EuclideanNon["22.1"], EuclideanExt["22.1"]}=Import["~/Physics/Research/MS/data/Spectralcc24.97.1.Rev.m"];

(*Export everything to a Mathematica file*)
Export["~/Physics/Research/MS/data/Spectralcc24.97.1.Rev.m", {parameters["curr"], ImG1["curr"], ImGC1["curr"], ImGL1["curr"], ImGQ1["curr"], ImGV1["curr"], ImG2["curr"], ImGC2["curr"], ImGL2["curr"], ImGQ2["curr"], ImGV2["curr"], ImG3["curr"], ImGC3["curr"], ImGL3["curr"], ImGQ3["curr"], ImGV3["curr"], ReGC1["curr"], ReGL1["curr"], ReGQ1["curr"], ReGV1["curr"], ReGC2["curr"], ReGL2["curr"], ReGQ2["curr"], ReGV2["curr"], ReGC3["curr"], ReGL3["curr"], ReGQ3["curr"], ReGV3["curr"], C0["curr"], Fraction["curr"], SpatialConstInter["22.1"], SpatialTotalInter["22.1"], SpatialNon["22.1"], EuclideanConstInter["22.1"], EuclideanTotalInter["22.1"], EuclideanNon["22.1"], EuclideanExt["22.1"]}]
