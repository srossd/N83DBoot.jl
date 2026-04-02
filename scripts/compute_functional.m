(* Unified entry point for computing N83DBoot functional values.

   Invocation:
     math -script compute_functional.m '{type -> "Crossing", cross -> 1, m -> 4, n -> 0,
       precision -> 1024, rOrder -> 100,
       operators -> {"Id", "Bp0020"},
       semishorts -> {{"Ap0020", {0, 2, 4}}, {"Atwo0100", {1, 3}}},
       long -> {{0, {1.2, 1.5}}, {2, {1.0}}},   (* twists τ = Δ−J *)
       filename -> "/path/to/cache",
       dataDir -> "/path/to/abjm/inputs"}'

   type "Crossing": crossing equation derivatives (CrossingFunctional)
   type "Integral":  localization measure integrals (IntegralFunctional)
*)

(* === Parse args and load utilities === *)
args = KeyMap[ToString, Association[ToExpression[$CommandLine[[4]]]]];
scriptsDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{scriptsDir, "util.m"}]];

If[!KeyExistsQ[args, "type"],     Print["Error: must specify type"]; Abort[]];
If[!KeyExistsQ[args, "filename"], Print["Error: must specify filename"]; Abort[]];
If[!KeyExistsQ[args, "dataDir"],  Print["Error: must specify dataDir"]; Abort[]];

(* Set data dir before loading abjm_defs.m *)
$n83dDataDir = args["dataDir"];

(* Map from operator name strings to fixed-multiplet symbols *)
opNameToSymbol = <|
  "Id"       -> Id,
  "Bp0020"   -> Bp0020,
  "Bp0040"   -> Bp0040,
  "Btwo0200" -> Btwo0200
|>;

opNameToFile = <|
  "Id"       -> "identity.csv",
  "Bp0020"   -> "Bp0020.csv",
  "Bp0040"   -> "Bp0040.csv",
  "Btwo0200" -> "Btwo0200.csv"
|>;

(* Launch multiple kernels if SLURM_CPUS_PER_TASK is defined *)
If[Environment["SLURM_CPUS_PER_TASK"] =!= $Failed,
    numKernels = ToExpression[Environment["SLURM_CPUS_PER_TASK"]];
    If[numKernels > 1,
        Print["SLURM_CPUS_PER_TASK detected: launching ", numKernels, " kernels"];
        LaunchKernels[numKernels];
    ];
];

(* ============================================================
   CROSSING type
   ============================================================ *)
If[args["type"] == "Crossing",

  args = ProcessArgs[
    args,
    <|"cross" -> "1", "m" -> "0", "n" -> "0", "precision" -> "1024",
      "rOrder" -> "100", "operators" -> "{}", "semishorts" -> "{}", "long" -> "{}"|>,
    {"cross", "m", "n", "precision", "rOrder", "operators", "semishorts", "long"},
    {}
  ];
  crossIdx = args["cross"];
  m        = args["m"];
  n        = args["n"];
  prec     = Ceiling[args["precision"] * Log[10, 2]]; (* convert bits to decimal digits *)
  rOrder   = args["rOrder"];
  cacheDir = args["filename"];

  Get[FileNameJoin[{scriptsDir, "abjm_defs.m"}]];

  (* Helper: compute one crossing derivative for a given superblock substitution rule *)
  computeCrossing[rule_] :=
    SetPrecision[expansionCoeff[cross[crossIdx] /. Frule /. rule, m, n, prec], prec];

  (* --- Fixed operators --- *)
  Do[
    outfile = FileNameJoin[{cacheDir, opNameToFile[opName]}];
    If[FileExistsQ[outfile] && LoadCSV[outfile] =!= <||>,
      Print["Skipping cached: ", opName],
      val = computeCrossing[superblockRule[opNameToSymbol[opName]]];
      Print["Writing: ", opName, " -> ", N[val, 10]];
      Quiet@CreateDirectory[cacheDir, CreateIntermediateDirectories -> True];
      SaveCSV[{{val}}, outfile, prec, "CanonicalColumns" -> {}]
    ],
    {opName, args["operators"]}
  ];

  (* --- Semishort operators: {{"Ap0020", {0,2,4}}, {"Atwo0100", {1,3}}} --- *)
  Do[
    {familyName, spins} = familyEntry;
    multipletSym = ToExpression[familyName];
    Do[
      outfile = FileNameJoin[{cacheDir, "semishort_" <> familyName <> "_J_" <> ToString[J] <> ".csv"}];
      If[FileExistsQ[outfile] && LoadCSV[outfile] =!= <||>,
        Print["Skipping cached: ", familyName, " J=", J],
        val = (-1)^J computeCrossing[superblockRule[multipletSym, J]];
        Print["Writing: ", familyName, " J=", J, " -> ", N[val, 10]];
        Quiet@CreateDirectory[cacheDir, CreateIntermediateDirectories -> True];
        SaveCSV[{{val}}, outfile, prec, "CanonicalColumns" -> {}]
      ],
      {J, spins}
    ],
    {familyEntry, args["semishorts"]}
  ];

  (* --- Long operators: {{J, {tau1, tau2, ...}}, ...} where tau = delta - J --- *)
  Do[
    {J, taus} = longEntry;
    outfile = FileNameJoin[{cacheDir, "spin_" <> ToString[J] <> ".csv"}];
    existingData = LoadCSV[outfile];
    remainingTaus = Complement[canonicalTwist /@ taus, Keys[existingData]];
    If[remainingTaus == {},
      Print["All twists cached for spin=", J, ". Skipping."],
      newData = Association @ MonitoredTable[
        canonicalTwist[tau] ->
          SetPrecision[expansionCoeff[
            cross[crossIdx] /. Frule /. superblockRule[Azero0000, SetPrecision[tau + J, prec], J],
            m, n, prec
          ], prec],
        {tau, remainingTaus}
      ];
      combined = SortBy[Join[List @@@ Normal[existingData], List @@@ Normal[newData]], First];
      Print["Writing: spin_", J, ".csv (", Length[combined], " rows)"];
      Quiet@CreateDirectory[cacheDir, CreateIntermediateDirectories -> True];
      SaveCSV[combined, outfile, prec, "CanonicalColumns" -> {1}]
    ],
    {longEntry, args["long"]}
  ];
];

(* ============================================================
   INTEGRAL type
   ============================================================ *)
If[args["type"] == "Integral",

  args = ProcessArgs[
    args,
    <|"b" -> "0", "precision" -> "256", "rOrder" -> "20",
      "operators" -> "{}", "semishorts" -> "{}", "long" -> "{}"|>,
    {"b", "precision", "rOrder", "operators", "semishorts", "long"},
    {"b"}
  ];
  b        = args["b"];
  prec     = Ceiling[args["precision"] * Log[10, 2]]; (* convert bits to decimal digits *)
  cacheDir = args["filename"];

  $loadIntegrals = True;
  Get[FileNameJoin[{scriptsDir, "abjm_defs.m"}]];

  (* --- Fixed operators --- *)
  Do[
    outfile = FileNameJoin[{cacheDir, opNameToFile[opName]}];
    If[FileExistsQ[outfile] && LoadCSV[outfile] =!= <||>,
      Print["Skipping cached: ", opName],
      val = integral[{opNameToSymbol[opName]}, b, prec];
      Print["Writing: ", opName, " -> ", N[val, 10]];
      Quiet@CreateDirectory[cacheDir, CreateIntermediateDirectories -> True];
      SaveCSV[{{val}}, outfile, prec, "CanonicalColumns" -> {}]
    ],
    {opName, args["operators"]}
  ];

  (* --- Semishort operators --- *)
  Do[
    {familyName, spins} = familyEntry;
    multipletSym = ToExpression[familyName];
    Do[
      outfile = FileNameJoin[{cacheDir, "semishort_" <> familyName <> "_J_" <> ToString[J] <> ".csv"}];
      If[FileExistsQ[outfile] && LoadCSV[outfile] =!= <||>,
        Print["Skipping cached: ", familyName, " J=", J],
        val = (-1)^J integral[{multipletSym, J}, b, prec];
        Print["Writing: ", familyName, " J=", J, " -> ", N[val, 10]];
        Quiet@CreateDirectory[cacheDir, CreateIntermediateDirectories -> True];
        SaveCSV[{{val}}, outfile, prec, "CanonicalColumns" -> {}]
      ],
      {J, spins}
    ],
    {familyEntry, args["semishorts"]}
  ];

  (* --- Long operators: {{J, {tau1, tau2, ...}}, ...} where tau = delta - J --- *)
  Do[
    {J, taus} = longEntry;
    outfile = FileNameJoin[{cacheDir, "spin_" <> ToString[J] <> ".csv"}];
    existingData = LoadCSV[outfile];
    remainingTaus = Complement[canonicalTwist /@ taus, Keys[existingData]];
    If[remainingTaus == {},
      Print["All twists cached for spin=", J, ". Skipping."],
      newData = Association @ MonitoredTable[
        canonicalTwist[tau] -> integral[{Azero0000, tau + J, J}, b, prec],
        {tau, remainingTaus}
      ];
      combined = SortBy[Join[List @@@ Normal[existingData], List @@@ Normal[newData]], First];
      Print["Writing: spin_", J, ".csv (", Length[combined], " rows)"];
      Quiet@CreateDirectory[cacheDir, CreateIntermediateDirectories -> True];
      SaveCSV[combined, outfile, prec, "CanonicalColumns" -> {1}]
    ],
    {longEntry, args["long"]}
  ];
];
