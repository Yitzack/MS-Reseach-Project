input = $CommandLine

Master = Import[StringJoin[Directory[],"/",input[[Dimensions[input]]],".xml"]]; (*Import the table for interpolation*)
Sub1 = N[Transpose[{Master[[1,All,1]],Master[[1,All,2,1]]}]];(*Extract the subtables*)
Sub2 = N[Transpose[{Master[[1,All,1]],Master[[1,All,2,2]]}]];
Sub3 = N[Transpose[{Master[[1,All,1]],Master[[1,All,2,3]]}]];
Sub4 = N[Transpose[{Master[[2,All,1]],Master[[2,All,2,1]]}]];
Sub5 = N[Transpose[{Master[[2,All,1]],Master[[2,All,2,2]]}]];
Sub6 = N[Transpose[{Master[[2,All,1]],Master[[2,All,2,3]]}]];
Sub7 = N[Transpose[{Master[[3,All,1]],Master[[3,All,2,1]]}]];
Sub8 = N[Transpose[{Master[[3,All,1]],Master[[3,All,2,2]]}]];
Sub9 = N[Transpose[{Master[[3,All,1]],Master[[3,All,2,3]]}]];
e1 = Sub7[[;;465,1,2]]; (*Pull the choice of s from the table*)

ImGa = Interpolation[Sub1, Method -> "Spline"];
ImGV1a = Interpolation[Sub2, Method -> "Spline"];
ImGV2a = Interpolation[Sub3, Method -> "Spline"];
ImGb = Interpolation[Sub4, Method -> "Spline"];
ImGV1b = Interpolation[Sub5, Method -> "Spline"];
ImGV2b = Interpolation[Sub6, Method -> "Spline"];
ImGc = Interpolation[Sub7, Method -> "Spline"];
ImGV1c = Interpolation[Sub8, Method -> "Spline"];
ImGV2c = Interpolation[Sub9, Method -> "Spline"];

List1=Flatten[ImGa[i,j]/.i->Range[0,208]/.j->Range[0,150]];
List2=Flatten[D[ImGa[i,j],i]/.i->Range[0,208]/.j->Range[0,150]];
List3=Flatten[D[ImGa[i,j],j]/.i->Range[0,208]/.j->Range[0,150]];
List4=Flatten[D[ImGa[i,j],i,j]/.i->Range[0,208]/.j->Range[0,150]];
List5=Flatten[ImGV1a[i,j]/.i->Range[0,208]/.j->Range[0,150]];
List6=Flatten[D[ImGV1a[i,j],i]/.i->Range[0,208]/.j->Range[0,150]];
List7=Flatten[D[ImGV1a[i,j],j]/.i->Range[0,208]/.j->Range[0,150]];
List8=Flatten[D[ImGV1a[i,j],i,j]/.i->Range[0,208]/.j->Range[0,150]];
List9=Flatten[ImGV2a[i,j]/.i->Range[0,208]/.j->Range[0,150]];
List10=Flatten[D[ImGV2a[i,j],i]/.i->Range[0,208]/.j->Range[0,150]];
List11=Flatten[D[ImGV2a[i,j],j]/.i->Range[0,208]/.j->Range[0,150]];
List12=Flatten[D[ImGV2a[i,j],i,j]/.i->Range[0,208]/.j->Range[0,150]];
List13=Flatten[ImGb[i,j]/.i->Range[208,788]/.j->Range[0,150]];
List14=Flatten[D[ImGb[i,j],i]/.i->Range[208,788]/.j->Range[0,150]];
List15=Flatten[D[ImGb[i,j],j]/.i->Range[208,788]/.j->Range[0,150]];
List16=Flatten[D[ImGb[i,j],i,j]/.i->Range[208,788]/.j->Range[0,150]];
List17=Flatten[ImGV1b[i,j]/.i->Range[208,788]/.j->Range[0,150]];
List18=Flatten[D[ImGV1b[i,j],i]/.i->Range[208,788]/.j->Range[0,150]];
List19=Flatten[D[ImGV1b[i,j],j]/.i->Range[208,788]/.j->Range[0,150]];
List20=Flatten[D[ImGV1b[i,j],i,j]/.i->Range[208,788]/.j->Range[0,150]];
List21=Flatten[ImGV2b[i,j]/.i->Range[208,788]/.j->Range[0,150]];
List22=Flatten[D[ImGV2b[i,j],i]/.i->Range[208,788]/.j->Range[0,150]];
List23=Flatten[D[ImGV2b[i,j],j]/.i->Range[208,788]/.j->Range[0,150]];
List24=Flatten[D[ImGV2b[i,j],i,j]/.i->Range[208,788]/.j->Range[0,150]];
List25=Flatten[ImGc[P,e]/.P->Range[0,600.8,.8]/.e->e1];
List26=Flatten[D[ImGc[P,e],P]/.P->Range[0,600.8,.8]/.e->e1];
List27=Flatten[2*Sqrt[e]*D[ImGc[P,e],e]/.P->Range[0,600.8,.8]/.e->e1];
List28=Flatten[2*Sqrt[e]*D[ImGc[P,e],P,e]/.P->Range[0,600.8,.8]/.e->e1];
List29=Flatten[ImGV1c[P,e]/.P->Range[0,600.8,.8]/.e->e1];
List30=Flatten[D[ImGV1c[P,e],P]/.P->Range[0,600.8,.8]/.e->e1];
List31=Flatten[2*Sqrt[e]*D[ImGV1c[P,e],e]/.P->Range[0,600.8,.8]/.e->e1];
List32=Flatten[2*Sqrt[e]*D[ImGV1c[P,e],P,e]/.P->Range[0,600.8,.8]/.e->e1];
List33=Flatten[ImGV2c[P,e]/.P->Range[0,600.8,.8]/.e->e1];
List34=Flatten[D[ImGV2c[P,e],P]/.P->Range[0,600.8,.8]/.e->e1];
List35=Flatten[2*Sqrt[e]*D[ImGV2c[P,e],e]/.P->Range[0,600.8,.8]/.e->e1];
List36=Flatten[2*Sqrt[e]*D[ImGV2c[P,e],P,e]/.P->Range[0,600.8,.8]/.e->e1];

Export[StringJoin[Directory[],"/",input[[Dimensions[input]]],"ImG.csv"], Join[Transpose[{List1,List2,List3,List4}],Transpose[{List13,List14,List15,List16}],Transpose[{List25,List26,List27,List28}]]] (*Export the table*)
Export[StringJoin[Directory[],"/",input[[Dimensions[input]]],"ImGV1.csv"], Join[Transpose[{List5,List6,List7,List8}],Transpose[{List17,List18,List19,List20}],Transpose[{List29,List30,List31,List32}]]]
Export[StringJoin[Directory[],"/",input[[Dimensions[input]]],"ImGV2.csv"], Join[Transpose[{List9,List10,List11,List12}],Transpose[{List21,List22,List23,List24}],Transpose[{List33,List34,List35,List36}]]]

Exit[]
