Spectral=Import["~/T.0/Spectral.xml"] (*Import the table for interpolation extrapolation*)
e1=Spectral[[1;;462,1]][[All,2]] (*Extract the energy abscissa*)

f=Interpolation[Spectral,Method->"Spline"] (*Turn the table into an interpolation*)
List1=D[f[P,e],P]/.P->Range[0,600.8,.8]/.e->e1 (*Calculate the derivitives*)
List2=D[f[P,e],e]/.P->Range[0,600.8,.8]/.e->e1
List3=D[f[P,e],P,e]/.P->Range[0,600.8,.8]/.e->e1

Export["~/T.0/Spectral.xml",{List1,List2,List3}] (*Export the table*)

Parameters=Table[err,{i,462},{j,4}] (*Declare a table for parameterization of the extrapolation*)

For[i = 1, i <= 462, i++, 
  Parameters[[i]] = {err, a, b, n} /. 
  FindFit[Transpose[{Spectral[[i ;; ;; 462, 1]][[All, 1]], Spectral[[i ;; ;; 462, 2]]}][[251 ;;]],
    b P^n, {{a, 1}, b, n}, P, MaxIterations -> 1000]] (*Calculate the parameterization of the extrapolation*)
For[i = 1, i <= 462, i++, Check[
  Parameters[[i]] = {1, a, b, n} /. 
    FindFit[Transpose[{Spectral[[i ;; ;; 462, 1]][[All, 1]], Spectral[[i ;; ;; 462, 2]]}][[251 ;;]], 
    a + b P^n, {a, {b, Parameters[[i, 3]]}, {n, Parameters[[i, 4]]}}, P, MaxIterations -> 1000], 
  Parameters[[i, 1]] = err, FindFit::cvmit]]

For[i = 1, i <= 462, i++, 
  If[Or[Parameters[[i, 1]] == err, Parameters[[i,4]] > 1], 
    Parameters[[i]] = {err, a, b, n} /. 
    FindFit[Transpose[{Spectral[[i ;; ;; 462, 1]][[All, 1]], Spectral[[i ;; ;; 462, 2]]}][[251 ;;]], 
    Spectral[[462 ;; ;; 462, 2]][[752]] + b P^n, {{a, Spectral[[462 ;; ;; 462, 2]][[752]]}, b, n}, P, MaxIterations -> 1000]]] (*Calculate the parameterization of the extrapolation*)
For[i = 1, i <= 462, i++, 
  If[Parameters[[i, 1]] == err, 
    Check[Parameters[[i]] = {1, a, b, n} /. 
    FindFit[Transpose[{Spectral[[i ;; ;; 462, 1]][[All, 1]], Spectral[[i ;; ;; 462, 2]]}][[251 ;;]], 
    a + b P^n, {{a, Parameters[[i, 2]]}, {b, Parameters[[i, 3]]}, {n, Parameters[[i, 4]]}}, P, MaxIterations -> 1000], 
  Parameters[[i, 1]] = err, FindFit::cvmit]]]

For[i = 1, i <= 462, i++,
 If[Or[Parameters[[i, 1]] == err, Parameters[[i,4]] > 1], 
  Holder[[1]] = {err, a, b, n} /. 
    FindFit[Transpose[{Spectral[[i ;; ;; 462, 1]][[All, 1]], Spectral[[i ;; ;; 462, 2]]}][[251 ;;]], 
    b P^n, {{a, 1}, b, n}, P, MaxIterations -> 1000]; 
  Holder[[1]] = {3, a, b, n} /. 
    FindFit[Transpose[{Spectral[[i ;; ;; 462, 1]][[All, 1]], Spectral[[i ;; ;; 462, 2]]}][[251 ;;]], 
    a + b P^n, {a, {b, Holder[[1, 3]]}, {n, Holder[[1, 4]]}}, P, MaxIterations -> 1000];
  Holder[[2]] = {err, a, b, n} /. 
    FindFit[Transpose[{Spectral[[i ;; ;; 462, 1]][[All, 1]], Spectral[[i ;; ;; 462, 2]]}][[251 ;;]], 
    Spectral[[462 ;; ;; 462, 2]][[752]] + b P^n, {{a, Spectral[[462 ;; ;; 462, 2]][[752]]}, b, n}, P, MaxIterations -> 1000];
  Holder[[2]] = {3, a, b, n} /. 
    FindFit[Transpose[{Spectral[[i ;; ;; 462, 1]][[All, 1]], Spectral[[i ;; ;; 462, 2]]}][[251 ;;]], 
    a + b P^n, {{a, Holder[[1, 2]]}, {b, Holder[[1, 3]]}, {n, Holder[[1, 4]]}}, P, MaxIterations -> 1000];
  If[Holder[[1, 4]] > 1,
   If[Holder[[2, 4]] > 1,
    Parameters[[i, 1]] = 4; 
    Parameters[[i, 2]] = (Holder[[1, 2]] + Holder[[2, 2]])/2.; 
    Parameters[[i, 3]] = (Holder[[1, 3]] + Holder[[2, 3]])/2.; 
    Parameters[[i, 4]] = (Holder[[1, 4]] + Holder[[2, 4]])/2.,
    Parameters[[i]] = Holder[[2]]],
   If[Holder[[2, 4]] > 1,
    Parameters[[i]] = Holder[[1]],
    Parameters[[i, 1]] = 3; 
    Parameters[[i, 2]] = (Holder[[1, 2]] + Holder[[2, 2]])/2.; 
    Parameters[[i, 3]] = (Holder[[1, 3]] + Holder[[2, 3]])/2.; 
    Parameters[[i, 4]] = (Holder[[1, 4]] + Holder[[2, 4]])/2.]]]]
    
a=Interpolation[Transpose[{e1,Parameters[[All,2]]}],Method->"Spline"](*Interpolation of parameters over energy*)
b=Interpolation[Transpose[{e1,Parameters[[All,3]]}],Method->"Spline"]
n=Interpolation[Transpose[{e1,Parameters[[All,4]]}],Method->"Spline"]

List1 = a[e] /. e -> e1 (*Place the parameterization in the list*)
List2 = D[a[e], e] /. e -> e1 (*Place the derivative of the parameterization in the list*)
List3 = b[e] /. e -> e1
List4 = D[b[e], e] /. e -> e1
List5 = n[e] /. e -> e1
List6 = D[n[e], e] /. e -> e1

Export["~/T.0/SpectralExtrapolate.xml",{List1,List2,List3,List4,List5,List6}] (*Export the table to a second file. Can't
go in the first file as it isn't of the same size.*)

Exit[]
