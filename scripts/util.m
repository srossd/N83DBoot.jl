(* Utility functions for N83DBoot Mathematica scripts *)

TruncateZeros[s_String] := StringReplace[s, RegularExpression["(\\d+\\.\\d*?)(0+)$"] -> "$1"];
SetAttributes[canonicalTwist, Listable];
canonicalTwist[tau_Integer] := tau;
canonicalTwist[tau_] := ToExpression@ToString[Quiet@DecimalForm[SetPrecision[tau, 30], {30, 15}]];

ProcessArgs[args_, defaults_, exprVars_, canonVars_] := AssociationMap[
  (#[[1]] -> (If[MemberQ[canonVars, #[[1]]], canonicalTwist, Identity] @ If[MemberQ[exprVars, #[[1]]], ToExpression, Identity] @ If[StringQ[#[[2]]], #[[2]], ToString[#[[2]], InputForm]])) &,
  Fold[If[!KeyExistsQ[#1, #2[[1]]],Append[#1, #2], #1] &, args, Normal[defaults]]
];

Options[SaveCSV] = {"CanonicalColumns" -> {1}};
SaveCSV[table_, filename_, prec_, opt : OptionsPattern[]] := Module[{out},
  Quiet@CreateDirectory[DirectoryName[filename], CreateIntermediateDirectories -> True];
  out = OpenWrite[filename];
  Do[
      If[!AllTrue[row, NumericQ],
        Print["Skipping non-numeric row ", row],
        WriteString[out, StringRiffle[TruncateZeros@*ToString /@ (Table[If[MemberQ[OptionValue["CanonicalColumns"], i], Quiet@DecimalForm[canonicalTwist[row[[i]]], {30, 15}], Quiet@DecimalForm[SetPrecision[row[[i]], prec], prec]], {i, Length[row]}]), ","] <> "\n"]
      ],
      {row, If[AssociationQ[table], Flatten[List @@ #] & /@ Normal[table], table]}
  ];
  Close[out];
];

LoadCSV[filename_] := If[FileExistsQ[filename],
  Module[{table},
    table = Map[ToExpression, If[StringQ[#[[1]]], StringSplit[#[[1]], ","], #] & /@ Import[filename, "Table"], {2}];
    Which[
      table == {}, <||>,
      Length[Dimensions[table]] == 2 && Dimensions[table][[2]] == 2, Association[Rule @@@ table],
      Dimensions[table] == {1,1}, First[First[table]],
      True, Association[Table[Most[row] -> Last[row], {row, table}]]
    ]
  ],
  <||>
];

Options[MonitoredTable] = {"UpdateFrequency" -> 1, "Title" -> "Progress"};
SetAttributes[MonitoredTable, HoldAllComplete];
MonitoredTable[expr_, {var_, n_}, OptionsPattern[]] := Module[{startTime, uniqueVar, colStrings, val, ctr, len, print, table},
  uniqueVar = Unique[ReleaseHold[ToString@var]];
  startTime = AbsoluteTime[];
  ctr = 0;
  len = If[IntegerQ[n], n, Length[n]];
  print[x___] := Print[x];
  table = Table;
  If[Length[Kernels[]] > 0,
    SetSharedVariable[ctr];
    SetSharedFunction[print];
    table = ParallelTable;
  ];
  table[
    val = ReleaseHold[Hold[expr] /. var -> uniqueVar];
    
    ctr += 1;
    If[Mod[ctr, OptionValue["UpdateFrequency"]] == 0,
    colStrings = {
      OptionValue["Title"]<>": "<>ToString[ctr]<>"/"<>ToString[len]<>" ("<>ToString@NumberForm[100. ctr/len, {3, 1}]<>"%)",
      "Elapsed time: "<>ToString@NumberForm[AbsoluteTime[] - startTime, {5, 1}]<>"s", 
      "ETA: "<>ToString@NumberForm[(AbsoluteTime[] - startTime) (len - ctr)/ctr, {5, 1}]<>"s"
    };
    print @@ (StringPadRight[#, Max[StringLength /@ colStrings] + 5]& /@ colStrings);
    ];
    
    val,
    {uniqueVar, n}
  ]
];
