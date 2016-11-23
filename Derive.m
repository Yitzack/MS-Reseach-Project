input = $CommandLine

Master = Import[StringJoin[Directory[],"/",input[[Dimensions[input]]],".xml"]]; (*Import the table for interpolation*)
Sub1 = N[Transpose[{Master[[2,All,1]],Master[[2,All,2,2]]}]];(*Extract the subtables*)
Sub2 = N[Transpose[{Master[[2,All,1]],Master[[2,All,2,3]]}]];
e1 = Sub1[[;;465,1,2]]; (*Pull the choice of s from the table*)

ImGV1 = Interpolation[Sub1, Method -> "Spline"];
ImGV2 = Interpolation[Sub2, Method -> "Spline"];

List1 = Flatten[D[ImGV1[P, e], P] /. P -> Range[0, 600.8, .8] /. e -> e1];
List2 = Flatten[2*Sqrt[e]*D[ImGV1[P, e], e] /. P -> Range[0, 600.8, .8] /. e -> e1]; (*2*Sqrt[e]* is to change the derivative from d/ds to d/d(sqrt(s))*)
List3 = Flatten[2*Sqrt[e]*D[ImGV1[P, e], P, e] /. P -> Range[0, 600.8, .8] /. e -> e1];
List4 = Flatten[D[ImGV2[P, e], P] /. P -> Range[0, 600.8, .8] /. e -> e1];
List5 = Flatten[2*Sqrt[e]*D[ImGV2[P, e], e] /. P -> Range[0, 600.8, .8] /. e -> e1];
List6 = Flatten[2*Sqrt[e]*D[ImGV2[P, e], P, e] /. P -> Range[0, 600.8, .8] /. e -> e1];

Export[StringJoin[Directory[],"/",input[[Dimensions[input]]],"ImGV1.csv"], Transpose[{List1,List2,List3}]] (*Export the table*)
Export[StringJoin[Directory[],"/",input[[Dimensions[input]]],"ImGV2.csv"], Transpose[{List4,List5,List6}]]

Exit[]
