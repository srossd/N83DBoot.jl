$ZcacheLoaded = ! FileExistsQ["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"];
If[FileExistsQ["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m"], <<"/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m"];

Clear[Hinvseries, Iseries, Rseries];
Hinvseries[k_, l_, n_, order_] := 
  Hinvseries[k, l, n, order] = 
   Series[1/(2  Cosh[((m1 - m2)  n)/2] - 
       Exp[(2  \[Pi]   I   l)/k]  (-1)^
         n  2  Cosh[((m1 + m2)  n)/2]), {m1, 0, order}, {m2, 0, 
     order}];
Iseries[k_, l_, n_, order_] := 
  Iseries[k, l, n, order] = 
   Series[Exp[((m1 - m2)  n)/2] - 
     Exp[(2  \[Pi]   I   l)/k]  (-1)^n  Exp[((m1 + m2)  n)/2], {m1, 0,
      order}, {m2, 0, order}];
Rseries[order_] := 
  Series[I   Exp[(I   m1   m2)/(2  \[Pi])], {m1, 0, order}, {m2, 0, 
    order}];

prefactorseries[k_, M_, NN_, order_] := 
  With[{series = 
     Normal@Series[
       Exp[I   \[Pi]  (M^3/(3  k) + 
            M^2/2 + (-11/12  k + 1/2 - 
               1/(6  k))  M)]  (-Exp[-((m1 + m2)/2)])^(M   NN), {m1, 
        0, order}, {m2, 0, order}]}, 
   Series[ComplexExpand[Normal@series], {m1, 0, order}, {m2, 0, 
     order}]];
invprefactorseries[k_, M_, NN_, order_] := 
  With[{series = 
     Normal@Series[(Exp[
           I   \[Pi]  (M^3/(3  k) + 
              M^2/2 + (-11/12  k + 1/2 - 
                 1/(6  k))  M)]  (-Exp[-((m1 + m2)/
                 2)])^(M   NN))^-1, {m1, 0, order}, {m2, 0, order}]}, 
   Series[ComplexExpand[Normal@series], {m1, 0, order}, {m2, 0, 
     order}]];

Clear[Zseries];
Zseries[k_, M_, 0, order_] := 
  Zseries[k, M, 0, order] = 
   With[{series = 
      Series[I^(M^2/2 - M)  Exp[-I   Pi   M  (M^2 - 1)/(6  k)]  k^(-M/
           2)  Product[
         2  Sin[(\[Pi]  (r - s))/k], {s, 1, M}, {r, s + 1, M}], {m1, 
        0, order}, {m2, 0, order}]}, 
    Series[ComplexExpand[Normal@series], {m1, 0, order}, {m2, 0, 
      order}]];
Zseries[k_, 0, 1, order_] := 
  Zseries[k, 0, 1] = 
   Series[1/(4  k   Cosh[m1/2]  Cosh[m2/2]), {m1, 0, order}, {m2, 0, 
     order}];

Clear[recur];
recur[1, 0, NN_, order_, recOrder_] := 
  recur[1, 0, NN, order, recOrder] = Module[{ans},
    ans = 
     Hinvseries[1, 0, NN, order]   Simplify[
       Sum[Rseries[recOrder]^(2   n - NN + 1) Zseries[1, 0, n, 
           recOrder]   Zseries[1, 0, NN - 1 - n, recOrder], {n, 0, 
          NN - 1}] - 
        Sum[Iseries[1, 0, 2   n - NN, recOrder]   Zseries[1, 0, n, 
           recOrder]   Zseries[1, 0, NN - n, recOrder], {n, NN - 1}]];
    ans
    ];

Zseries[k_, M_, NN_, order_] /; OddQ[order] := (Zseries[k, M, NN, order + 1]; Zseries[k, M, NN, order]);
Zseries[1, 0, NN_, order_] := 
  Module[{recOrder = If[EvenQ[NN], order + 2, order], rec, poly},
  (* If[!TrueQ[$ZcacheLoaded], <<"/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"; $ZcacheLoaded = True;]; *)
   Block[{$saveBlock = True},
    rec = recur[1, 0, NN, order, recOrder];
    ];
   poly = 
    Simplify@
      CoefficientList[rec, {m1, m2}][[;; order + 1, ;; order + 1]] . 
     m1^Range[0, order] . m2^Range[0, order];
   Do[Zseries[1, 0, NN, o] = Series[poly, {m1, 0, o}, {m2, 0, o}], {o,
      0, order}];
   (* If[!TrueQ[$saveBlock],
    If[FileExistsQ["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"],
      DeleteFile["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"]
    ];
    Save["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m", {Zseries, Zpseries}];
   ]; *)
   Zseries[1, 0, NN, order]
   ];

Zpseries[k_, M_, NN_, order_] := Zpseries[k, M, NN, order] =
  Simplify[
   prefactorseries[k, M, NN, order]   Zseries[k, M, NN, order]];

recur[k_, 0, NN_, order_, recOrder_] /; k >= 2 := 
  recur[k, 0, NN, order, recOrder] = Module[{ans},
    ans = Hinvseries[k, 0, NN, order]   Simplify[
       Sum[Rseries[recOrder]^(k   n)  Zpseries[k, 1, NN - 1 - n, 
           recOrder]   Zpseries[k, k - 1, n, recOrder], {n, 0, 
          NN - 1}] - 
        Sum[Iseries[k, 0, 2   n - NN, recOrder]   Zseries[k, 0, n, 
           recOrder]   Zseries[k, 0, NN - n, recOrder], {n, 
          NN - 1}]];
    ans
    ];
Zseries[k_, 0, NN_, order_] /; k >= 2 := 
 Module[{recOrder = If[EvenQ[NN], order + 2, order], rec, poly},
  (* If[!TrueQ[$ZcacheLoaded], <<"/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"; $ZcacheLoaded = True;]; *)
  Block[{$saveBlock = True},
   rec = Simplify@recur[k, 0, NN, order, recOrder]
   ];
  poly = (Simplify@CoefficientList[rec, {m1, m2}])[[;; order + 1, ;; 
       order + 1]] . m1^Range[0, order] . m2^Range[0, order];
  Do[Zseries[k, 0, NN, o] = Series[poly, {m1, 0, o}, {m2, 0, o}], {o, 
    0, order}];
   (* If[!TrueQ[$saveBlock],
    If[FileExistsQ["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"],
      DeleteFile["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"]
    ];
    Save["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m", {Zseries, Zpseries}];
   ]; *)
  Zseries[k, 0, NN, order]]

recur[2, 1, NN_, order_, recOrder_] := 
  recur[2, 1, NN, order, recOrder] = Module[{ans},
    ans = Hinvseries[2, 1, NN, order]   invprefactorseries[2, 1, NN, 
       recOrder]   (1/Zpseries[2, 1, 0, order])   Simplify[
       Sum[Rseries[recOrder]^(-2   n)  Zpseries[2, 0, n, 
           recOrder]   Zseries[2, 0, NN - n, recOrder], {n, 0, NN}] - 
        Sum[Iseries[2, 1, 2   n - NN, recOrder]   Zpseries[2, 1, n, 
           recOrder]   Zpseries[2, 1, NN - n, recOrder], {n, NN - 1}]];
    ans
    ];

Zseries[2, 1, NN_, order_] := 
 Module[{recOrder = If[OddQ[NN], order + 2, order], rec, poly},
  (* If[!TrueQ[$ZcacheLoaded], <<"/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"; $ZcacheLoaded = True;]; *)
  Block[{$saveBlock = True},
   rec = Simplify@recur[2, 1, NN, order, recOrder];
   ];
  poly = (Simplify@CoefficientList[rec, {m1, m2}])[[;; order + 1, ;; 
       order + 1]] . m1^Range[0, order] . m2^Range[0, order];
  Do[Zseries[2, 1, NN, o] = Series[poly, {m1, 0, o}, {m2, 0, o}], {o, 
    0, order}];
   (* If[!TrueQ[$saveBlock],
    If[FileExistsQ["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"],
      DeleteFile["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m"]
    ];
    Save["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache.m", {Zseries, Zpseries}];
   ]; *)
  Zseries[2, 1, NN, order]]

Zderiv[k_, M_, NN_, order1_, order2_] := 
  If[order1 == order2 == 0, Normal[Zseries[k, M, NN, 0]], 
   Pi^(order1 + order2)   (order1!)   (order2!)   Coefficient[
      Zseries[k, M, NN, Max[order1, order2]], 
      m1^order1   m2^order2] /. {m1 -> 0, m2 -> 0}];

cTfull[k_, M_, N_] := If[Head[cTfullcache[k, M, N]] =!= cTfullcache,
  cTfullcache[k, M, N],
  cTfullcache[k, M, N] = 
    -(64/\[Pi]^2)   (Zderiv[k, M, N, 2, 0]/
          Zderiv[k, M, N, 0, 
          0] - (Zderiv[k, M, N, 1, 0]/Zderiv[k, M, N, 0, 0])^2) // 
      Simplify;
  (* If[FileExistsQ[
  "/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m"],DeleteFile[
  "/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m"]];
  Save["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m",
  {cTfullcache, \[Lambda]B2fullcache, rhscache}]; *)
  
  cTfullcache[k, M, N]
];

\[Lambda]B2full[k_, M_, N_] := If[Head[\[Lambda]B2fullcache[k, M, N]] =!= \[Lambda]B2fullcache,
  \[Lambda]B2fullcache[k, M, N],
  \[Lambda]B2fullcache[k, M, N] = 
    16  (-1 - 256/(9 cTfull[k, M, N]) + 
        5/9  ((Zderiv[k, M, N, 4, 0]/Zderiv[k, M, N, 0, 0] - 
            4 (Zderiv[k, M, N, 1, 0] Zderiv[k, M, N, 3, 0])/
              Zderiv[k, M, N, 0, 0]^2 + 
            6 (Zderiv[k, M, N, 1, 0]^2 Zderiv[k, M, N, 2, 0])/
              Zderiv[k, M, N, 0, 0]^3 - 
            3 Zderiv[k, M, N, 1, 0]^4/
              Zderiv[k, M, N, 0, 0]^4)/((-\[Pi]^2/64) cTfull[k, M, 
              N])^2)) // Simplify;
  (* If[FileExistsQ[
  "/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m"],DeleteFile[
  "/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m"]];
  Save["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m",
  {cTfullcache, \[Lambda]B2fullcache, rhscache}]; *)
  
  \[Lambda]B2fullcache[k, M, N]
];
  

cT[k_, M_, N_] := cTfull[k, M, N] - 16   Boole[k == 1] // Simplify;
\[Lambda]B2[k_, M_, N_] := 
  If[k > 1, \[Lambda]B2full[k, M, 
    N], (\[Lambda]B2full[k, M, N]   cTfull[k, M, N]^2 - 
       2   (32/3)   (16)   cT[k, M, N])/cT[k, M, N]^2 // Simplify];
\[Lambda]Bp[k_, M_, 
   N_] := (4   (256/cT[k, M, N]) + \[Lambda]B2[k, M, N] + 16)/5 // 
   Simplify;

rhs[k_, M_, N_] := If[Head[rhscache[k, M, N]] =!= rhscache,
  rhscache[k, M, N],
  rhscache[k, M, N] = 
    64/(\[Pi]^2  cT[k, M, 
        N]^2)   (-((-((6 Zderiv[k, M, N, 0, 1]^2)/
            Zderiv[k, M, N, 0, 0]^4) + (2 Zderiv[k, M, N, 0, 2])/
          Zderiv[k, M, N, 0, 0]^3)  Zderiv[k, M, N, 1, 0]^2) - (
      8 Zderiv[k, M, N, 0, 1] Zderiv[k, M, N, 1, 0] Zderiv[k, M, N, 1, 
        1])/Zderiv[k, M, N, 0, 0]^3 + (
      2 Zderiv[k, M, N, 1, 1]^2 + 
      2 Zderiv[k, M, N, 1, 0] Zderiv[k, M, N, 1, 2])/
      Zderiv[k, M, N, 0, 
      0]^2 - ((2 Zderiv[k, M, N, 0, 1]^2)/Zderiv[k, M, N, 0, 0]^3 - 
        Zderiv[k, M, N, 0, 2]/Zderiv[k, M, N, 0, 0]^2)  Zderiv[k, M, N,
        2, 0] + (2 Zderiv[k, M, N, 0, 1] Zderiv[k, M, N, 2, 1])/
      Zderiv[k, M, N, 0, 0]^2 - Zderiv[k, M, N, 2, 2]/
      Zderiv[k, M, N, 0, 0]) // Simplify;
  (* If[FileExistsQ[
  "/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m"],DeleteFile[
  "/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m"]];
  Save["/home/sdempsey/bootstrap/abjm/inputs/nosaka_cache_results.m",
  {cTfullcache, \[Lambda]B2fullcache, rhscache}]; *)
  
  rhscache[k, M, N]
];