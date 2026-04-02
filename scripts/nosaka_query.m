(* nosaka_query.m
   Command-line wrapper for localization data in 3DN8Boot.jl.

   Usage:
       math -script nosaka_query.m N k M

   Outputs a single space-separated line:
       cT  lambda2_Bp0020  lambda2_Bp0040  lambda2_Btwo0200  rhs

   where:
       lambda2_Bp0020 = 256/cT(k,M,N)                  [stress-tensor OPE coeff]
       lambda2_Bp0040 = (4*lambda2_Bp0020 + lambda2_Btwo0200 + 16)/5  [Ward identity]
       lambda2_Btwo0200 = \[Lambda]B2[k,M,N]            [twist-2 OPE coeff]
       cT  = 256 / lambda2_Bp0020                       [Weyl anomaly coefficient]
       rhs = localization integral RHS from nosaka.m

   This script is called by N83DBoot._compute_localization_data when N is
   outside the tabulated range of the pre-computed CSV files.
*)

{abjmN, abjmK} = ToExpression /@ $CommandLine[[4;;5]];
abjmM = If[Length[$CommandLine] >= 6, ToExpression[$CommandLine[[6]]], 0];

(* Load the matrix model routines *)
scriptDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{scriptDir, "nosaka.m"}]];

Off[N::meprec];
Off[Series::ztest1];

(* Compute the three independent OPE coefficients *)
lam2_Bp0020   = N[256/cT[abjmK, abjmM, abjmN], 150];
lam2_Btwo0200 = N[\[Lambda]B2[abjmK, abjmM, abjmN], 150];
lam2_Bp0040   = N[\[Lambda]Bp[abjmK, abjmM, abjmN], 150];
cT_val        = N[cT[abjmK, abjmM, abjmN], 150];
rhs_val       = N[rhs[abjmK, abjmM, abjmN], 150];

Print[StringRiffle[
    ToString /@ {cT_val, lam2_Bp0020, lam2_Bp0040, lam2_Btwo0200, rhs_val},
    " "]]
