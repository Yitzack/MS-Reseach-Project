T=Import["~/T.5/Peak.xml"] (*Import the table*)
Gammas=Range[752] (*Declare an array, not sure about necessity*)
For[i=1,i<=752,i++,Gammas[[i]]=({.8*i-.8, a, gamma, x0}/.If[Dimensions[T[[i]]][[1]] <= 20,{a->1.6,gamma->0,x0->3.040308},FindFit[T[[i]],{a/Pi*(gamma/((x-x0)^2+gamma^2))},{{a,1.6},{x0,3.040308},gamma},x]])] (*Find a fit if there there enough valid data to do so and store the parameters to the array*)
Gammas[[All,2]]=Abs[Gammas[[All,2]]] (*Change the signs, sometimes Mathematica will put negatives on both*)
Gammas[[All,3]]=Abs[Gammas[[All,3]]]
Gammas[[All,3]]*=2 (*Change the gamma, which is half-width-half-max, to full-width-half-max*)
Picture=ListPlot[Gammas[[All,3]],PlotRange->Full]
Export["~/T.5/Epsilon Plot.png",Picture]
Export["~/T.5/Epsilon List.xml",Gammas]
Exit[]
