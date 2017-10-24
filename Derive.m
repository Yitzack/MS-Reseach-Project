input = $CommandLine

Master = Import[StringJoin[Directory[],"/",input[[Dimensions[input]]],".xml"]]; (*Import the table for interpolation*)
Sub1 = N[Transpose[{Master[[1,All,1]],Master[[1,All,2,2]]}]];(*Extract the subtables*)
Sub2 = N[Transpose[{Master[[1,All,1]],Master[[1,All,2,3]]}]];
Sub3 = N[Transpose[{Master[[2,All,1]],Master[[2,All,2,2]]}]];
Sub4 = N[Transpose[{Master[[2,All,1]],Master[[2,All,2,3]]}]];
Sub5 = N[Transpose[{Master[[3,All,1]],Master[[3,All,2,2]]}]];
Sub6 = N[Transpose[{Master[[3,All,1]],Master[[3,All,2,3]]}]];
e1 = Sub1[[;;465,1,2]]; (*Pull the choice of s from the table*)

ImGV1a = Interpolation[Sub1, Method -> "Spline", InterpolationOrder -> 2];
ImGV2a = Interpolation[Sub2, Method -> "Spline", InterpolationOrder -> 2];
ImGV1b = Interpolation[Sub3, Method -> "Spline", InterpolationOrder -> 2];
ImGV2b = Interpolation[Sub4, Method -> "Spline", InterpolationOrder -> 2];
ImGV1c = Interpolation[Sub5, Method -> "Spline", InterpolationOrder -> 2];
ImGV2c = Interpolation[Sub6, Method -> "Spline", InterpolationOrder -> 2];

List1 = Flatten[D[ImGV1a[i, j], i] /. i -> Range[0, 208] /. j -> Range[0, 150];
List2 = Flatten[D[ImGV1a[i, j], j] /. i -> Range[0, 208] /. j -> Range[0, 150];
List3 = Flatten[D[ImGV1a[i, j], i, j] /. i -> Range[0, 208] /. j -> Range[0, 150];
List4 = Flatten[D[ImGV2a[i, j], i] /. i -> Range[0, 208] /. j -> Range[0, 150];
List5 = Flatten[D[ImGV2a[i, j], j] /. i -> Range[0, 208] /. j -> Range[0, 150];
List6 = Flatten[D[ImGV2a[i, j], i, j] /. i -> Range[0, 208] /. j -> Range[0, 150];
List7 = Flatten[D[ImGV1b[i, j], i] /. i -> Range[208, 788] /. j -> Range[0, 150];
List8 = Flatten[D[ImGV1b[i, j], j] /. i -> Range[208, 788] /. j -> Range[0, 150];
List9 = Flatten[D[ImGV1b[i, j], i, j] /. i -> Range[208, 788] /. j -> Range[0, 150];
List10 = Flatten[D[ImGV2b[i, j], i] /. i -> Range[208, 788] /. j -> Range[0, 150];
List11 = Flatten[D[ImGV2b[i, j], j] /. i -> Range[208, 788] /. j -> Range[0, 150];
List12 = Flatten[D[ImGV2b[i, j], i, j] /. i -> Range[208, 788] /. j -> Range[0, 150];
List13 = Flatten[D[ImGV1b[P, e], P] /. P -> Range[0, 600.8, .8] /. e -> e1];
List14 = Flatten[2*Sqrt[e]*D[ImGV1c[P, e], e] /. P -> Range[0, 600.8, .8] /. e -> e1]; (*2*Sqrt[e]* is to change the derivative from d/ds to d/d(sqrt(s))*)
List15 = Flatten[2*Sqrt[e]*D[ImGV1c[P, e], P, e] /. P -> Range[0, 600.8, .8] /. e -> e1];
List16 = Flatten[D[ImGV2c[P, e], P] /. P -> Range[0, 600.8, .8] /. e -> e1];
List17 = Flatten[2*Sqrt[e]*D[ImGV2c[P, e], e] /. P -> Range[0, 600.8, .8] /. e -> e1];
List18 = Flatten[2*Sqrt[e]*D[ImGV2c[P, e], P, e] /. P -> Range[0, 600.8, .8] /. e -> e1];

Export[StringJoin[Directory[],"/",input[[Dimensions[input]]],"ImGV1.csv"], Transpose[{List1,List2,List3,List7,List8,List9,List13,List14,List15}]] (*Export the table*)
Export[StringJoin[Directory[],"/",input[[Dimensions[input]]],"ImGV2.csv"], Transpose[{List4,List5,List6,List10,List11,List12,List16,List17,List18}]]

Exit[]
