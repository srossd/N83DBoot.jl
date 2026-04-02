Get[FileNameJoin[{scriptsDir, "superblocks.m"}]];
Get[FileNameJoin[{scriptsDir, "scalar_blocks.m"}]];

If[TrueQ[$loadIntegrals],
  Get[FileNameJoin[{$n83dDataDir, "block_expansions", "block-expansion-Lmax70-20.m"}]];
  rOrderInt = 20;
];
If[TrueQ[$loadIntegralCache],
  Check[Get[#], Print[#]]& /@ FileNames[FileNameJoin[{$n83dDataDir, "gint_cache", "*"}]];
  Set @@@ (DownValues[GintCache] /. GintCache -> Gint);
];

If[!ValueQ[rOrder], rOrder = 30];
maxDerivs = 20;
scalarBlocks[spin_] := scalarBlocks[spin] = LoadScalarBlocks[spin, "Dimension" -> 3, "Order" -> rOrder, "KeptPoles" -> rOrder, "MaxDerivs" -> maxDerivs, "Directory" -> FileNameJoin[{$n83dDataDir, "blocks"}]];

cross[1] = Subscript[F, 1, 0] + Subscript[F, 1, 1] + (5/3) Subscript[F, 2, 0] - (2/5) Subscript[F, 2, 1] - (14/3) Subscript[F, 2, 2];
cross[2] = Subscript[F, 0, 0] - (1/4) Subscript[F, 1, 1] - (20/21) Subscript[F, 2, 0] + Subscript[F, 2, 1] - (14/15) Subscript[F, 2, 2];

Frule = Subscript[F, n_, m_] :> uvmultiply[Subscript[A, n, m], 0, 1];

expansionCoeff[a_ + b_, zOrder_, zbOrder_, prec_ : 300] := expansionCoeff[a, zOrder, zbOrder, prec] + expansionCoeff[b, zOrder, zbOrder, prec];
expansionCoeff[a_, zOrder_, zbOrder_, prec_ : 300] /; FreeQ[a, G | uvmultiply] := a Boole[zOrder == zbOrder == 0];
expansionCoeff[a_?NumericQ b_, zOrder_, zbOrder_, prec_ : 300] /; FreeQ[a, G] := a expansionCoeff[b, zOrder, zbOrder, prec];

uvCoeffs[uOrder_, vOrder_] := uvCoeffs[uOrder, vOrder] = CoefficientList[Series[((zp + 1/2) (zbp + 1/2))^uOrder ((1/2 - zp) (1/2 - zbp))^vOrder, {zp, 0, 100}, {zbp, 0, 100}], {zp, zbp}];
uvmultiply[0, __] = 0;
expansionCoeff[uvmultiply[expr_, uOrder_, vOrder_], zOrder_, zbOrder_, prec_ : 300] := With[{coeffs = uvCoeffs[uOrder, vOrder]},
  Sum[coeffs[[ii + 1, jj + 1]] expansionCoeff[expr, zOrder - ii, zbOrder - jj, prec], 
    {ii, 0, Min[zOrder, Length[coeffs] - 1]}, 
    {jj, 0, Min[zbOrder, Length[coeffs] - 1]}
  ]
]

expansionCoeff[Subscript[G, \[CapitalDelta]_, j_], zOrder_, zbOrder_, prec_ : 300] := Gcoeff[\[CapitalDelta], j, zOrder, zbOrder, prec];

Clear[Gcoeff];
Gcoeff[0, 0, zOrder_, zbOrder_, prec_ : 300] := Boole[zOrder == zbOrder == 0];
Gcoeff[\[CapitalDelta]_, j_, zOrder_, zbOrder_, prec_ : 300] /; zOrder >= zbOrder := Gcoeff[\[CapitalDelta], j, zOrder, zbOrder, prec] = 
  If[zOrder < 0 || zbOrder < 0, 0, 
  With[{num = SetPrecision[((3 - 2 Sqrt[2]))^(\[CapitalDelta]) (zzbDeriv[zOrder, zbOrder] /. scalarBlocks[Round@j]), prec],
        denom = SetPrecision[Times @@ Flatten[Table[(x - poles[i])^i, {i, Length[shiftedPoles /. scalarBlocks[Round@j]]}] /. poles[i_] :> (shiftedPoles /. scalarBlocks[Round@j])[[i]]], prec]},
  1/(zOrder! zbOrder!)  If[(\[CapitalDelta] == j + 1 && j > 0) || (\[CapitalDelta] == 0 && j == 0), D[num, x]/D[denom, x], num/denom] /. x -> SetPrecision[\[CapitalDelta] - (j + 1), prec]
  ]
  ];
Gcoeff[\[CapitalDelta]_, j_, zOrder_, zbOrder_, prec_ : 300] /; zOrder < zbOrder := 
  Gcoeff[\[CapitalDelta], j, zbOrder, zOrder, prec];


measure = (1/96) (64 r^2 (r^2 - 1))/((1 + r^2 + 2 r \[Eta])^2 (1 + r^4 + 2 r^2 (1 - 2 \[Eta]^2)));
ufunc = 16 r^2/(r^2 + 2 r \[Eta] + 1)^2;
vfunc = (r^2 - 2 r \[Eta] + 1)^2 / (r^2 + 2 r \[Eta] + 1)^2;
rmax = 2 + Abs[\[Eta]] - Sqrt[\[Eta]^2 + 4 Abs[\[Eta]] + 3];
rSolve[t_] := 1/(4 (t + t^3)) ((1 - 10 t^2 + t^4) \[Eta] + Sqrt[64 (t + t^3)^2 + (-1 + t^2)^2 (1 - 34 t^2 + t^4) \[Eta]^2] - Sqrt[2] \[Sqrt](24 (t + t^3)^2 + (1 - 28 t^2 + 86 t^4 - 28 t^6 + t^8) \[Eta]^2 + \[Eta] Sqrt[64 (t + t^3)^2 + (-1 + t^2)^2 (1 - 34 t^2 + t^4) \[Eta]^2] - 10 t^2 \[Eta] Sqrt[64 (t + t^3)^2 + (-1 + t^2)^2 (1 - 34 t^2 + t^4) \[Eta]^2] + t^4 \[Eta] Sqrt[64 (t + t^3)^2 + (-1 + t^2)^2 (1 - 34 t^2 + t^4) \[Eta]^2]));
rSolve2[t_] := ((-1 + t (12 + t (10 - (-12 + t) t))) \[Eta])/((-1 + t)^2 (1 + t (6 + t))) + (2 Sqrt[(-1 + t)^4 (1 + t (6 + t))^2 - 8 t (1 + t)^2 (1 + (-6 + t) t) (1 + t^2) \[Eta]^2])/((-1 + t)^2 (1 + t (6 + t))) - 1/2 \[Sqrt](12 + 1/((-1 + t)^4 (1 + t (6 + t))^2) 4 \[Eta] (\[Eta] + t (-56 + t (252 + t (504 + t (646 + t (504 + t (252 + (-56 + t) t)))))) \[Eta] - 4 (1 + t (-12 + t (-10 + (-12 + t) t))) \[Sqrt]((-1 + t)^4 (1 + t (6 + t))^2 - 8 t (1 + t)^2 (1 + (-6 + t) t) (1 + t^2) \[Eta]^2)));

masterIntegral[prefactor_] := masterIntegral[prefactor] = Integrate[measure r^power (prefactor /. {u -> ufunc, v -> vfunc}), r];
masterIntegral[1/u] = -1/48*(r^(1 + power)*((\[Eta] + Sqrt[-1 + \[Eta]^2])*
     Hypergeometric2F1[1, (1 + power)/2, (3 + power)/2, 
      -(r^2/(1 - 2*\[Eta]^2 + 2*\[Eta]*Sqrt[-1 + \[Eta]^2]))] + 
    (\[Eta] - Sqrt[-1 + \[Eta]^2])*Hypergeometric2F1[1, (1 + power)/2, 
      (3 + power)/2, r^2/(-1 + 2*\[Eta]^2 + 2*\[Eta]*Sqrt[-1 + \[Eta]^2])]))/
  ((1 + power)*\[Eta]);
masterIntegral[1/u^2] = -1/384*(r^(-1 + power)*(power + power^2 + power*r^2 - 
    power^2*r^2 + 4*r*\[Eta] - 4*power^2*r*\[Eta] - 8*power*\[Eta]^2 - 
    8*power^2*\[Eta]^2 - 4*power*(1 + power)*\[Eta]*
     (-\[Eta] + Sqrt[-1 + \[Eta]^2])*Hypergeometric2F1[1, -1 + power, 
      power, r/(\[Eta] - Sqrt[-1 + \[Eta]^2])] + 4*power*(1 + power)*\[Eta]*
     (\[Eta] + Sqrt[-1 + \[Eta]^2])*Hypergeometric2F1[1, -1 + power, 
      power, r/(\[Eta] + Sqrt[-1 + \[Eta]^2])]))/(power*(-1 + power^2));
masterIntegral[v/u^2] = -1/384*(r^(-1 + power)*(power + power^2 + power*r^2 - 
    power^2*r^2 - 4*r*\[Eta] + 4*power^2*r*\[Eta] - 8*power*\[Eta]^2 - 
    8*power^2*\[Eta]^2 - 4*power*(1 + power)*\[Eta]*
     (-\[Eta] + Sqrt[-1 + \[Eta]^2])*Hypergeometric2F1[1, -1 + power, 
      power, r/(-\[Eta] + Sqrt[-1 + \[Eta]^2])] + 4*power*(1 + power)*\[Eta]*
     (\[Eta] + Sqrt[-1 + \[Eta]^2])*Hypergeometric2F1[1, -1 + power, 
      power, -(r/(\[Eta] + Sqrt[-1 + \[Eta]^2]))]))/(power*(-1 + power^2));
      
rint[pow_, prefactor_, b_] := If[b === 0, (# /. r -> rmax), (# /. r -> rSolve[b]) + (# /. r -> rSolve2[b]) - (# /. r -> rmax)]& @ (masterIntegral[prefactor] /. power -> pow);

GintIntegrand[\[CapitalDelta]p_, j_, prefactor_, b_, prec_] := Sum[
    (blockCoeff[j, order] /. \[CapitalDelta] -> SetPrecision[\[CapitalDelta]p, 2 prec]) rint[SetPrecision[\[CapitalDelta]p + order, 2 prec], prefactor, SetPrecision[b, 2 prec]],
  {order, 0, rOrderInt, 2}]; 

(* Gint[\[CapitalDelta]p_, j_, prefactor_, b_, prec_] := Gint[\[CapitalDelta]p, j, prefactor, b, prec] = 6 Re@NIntegrate[
  Sum[
    (blockCoeff[j, order] /. \[CapitalDelta] -> SetPrecision[\[CapitalDelta]p, 2 prec]) rint[SetPrecision[\[CapitalDelta]p + order, 2 prec], prefactor, SetPrecision[b, 2 prec]],
  {order, 0, rOrderInt, 2}], 
  {\[Eta], 0, 1}, WorkingPrecision -> prec]; *)

integral[superblock_, b_, prec_] := 6 Re@NIntegrate[
  SetPrecision[With[{rule = superblockRule @@ superblock},
    (Simplify[(Subscript[A, 0, 0] - Subscript[A, 1, 1]/4 +Subscript[A, 2, 0]/21 + Subscript[A, 2, 2]/15) /. rule] /. Subscript[G, \[CapitalDelta]_, j_] :> GintIntegrand[\[CapitalDelta], j, 1/u, b, prec])
  + (Simplify[(Subscript[A, 2, 0] + Subscript[A, 2, 1] + Subscript[A, 2, 2]) /. rule] /. Subscript[G, \[CapitalDelta]_, j_] :> GintIntegrand[\[CapitalDelta], j, 1/u^2, b, prec])
  + (Simplify[(Subscript[A, 2, 0] - Subscript[A, 2, 1] + Subscript[A, 2, 2]) /. rule] /. Subscript[G, \[CapitalDelta]_, j_] :> GintIntegrand[\[CapitalDelta], j, v/u^2, b, prec])
  - 3 If[superblock === {Id}, GintIntegrand[0, 0, 1/u, b, prec], 0]
], 2 prec], {\[Eta], 0, 1}, WorkingPrecision -> prec];
