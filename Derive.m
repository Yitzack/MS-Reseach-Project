Spectral=Import["~/Physics/Research/MS/Spectral.1.xml"] (*Import the table for interpolation extrapolation*)
STable1 = STable[[1]]; (*Extract the subtables*)
STable2 = STable[[2]];
STable3 = STable[[3]];

s[j_] := If[j <= 425, N[(j - 400)/10], If[j <= 626, 2.5040308 + (j - 426)/200, If[j <= 667, 3.55 + (11 (j - 627))/800, 4.1 + (j - 667)/10]]]  (*Remove the enumeration of positive s values*)
For[i = 1, i <= 348176, i++, STable3[[i]][[1]][[2]] = s[STable3[[i]][[1]][[2]]]^2; STable3[[i]][[1]][[1]] = 4/5 STable3[[i]][[1]][[1]]]

e1=Spectral[[1;;462,1]][[All,2]] (*Extract the (invariant mass)^2 abscissa*)
Clear[i] (*Prepare i for use in the derivitive calculation*)

f = Interpolation[STable1, Method -> "Spline"] (*Turn the table into an interpolation*)
g = Interpolation[STable2, Method -> "Spline"]
h = Interpolation[STable3, Method -> "Spline"]
List1 = D[f[i, j], i] /. i -> Range[0, 13] /. j -> Range[0, 400]; (*Calculate the derivitives*)
List2 = D[f[i, j], j] /. i -> Range[0, 13] /. j -> Range[0, 400];
List3 = D[f[i, j], i, j] /. i -> Range[0, 13] /. j -> Range[0, 400];
List4 = D[g[i, j], i] /. i -> Range[13, 751] /. j -> Range[0, 400];
List5 = D[g[i, j], j] /. i -> Range[13, 751] /. j -> Range[0, 400];
List6 = D[g[i, j], i, j] /. i -> Range[13, 751] /. j -> Range[0, 400];
List7 = D[h[P, e], P] /. P -> Range[0, 600.8, .8] /. e -> e1;
List8 = D[h[P, e], e] /. P -> Range[0, 600.8, .8] /. e -> e1;
List9 = D[h[P, e], P, e] /. P -> Range[0, 600.8, .8] /. e -> e1;

Export["~/Physics/Research/MS/SpectralLarge.1.xml", {List1, List2, List3, List4, List5, List6, List7, List8, List9}] (*Export the table*)

Exit[]
