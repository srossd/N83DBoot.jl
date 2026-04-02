(* ::Package:: *)

Options[RunScalarBlocks] = Options[LoadScalarBlocks] = {"Dimension" -> 4, "Order" -> 20, "MaxDerivs" -> 10, "KeptPoles" -> 20, "Delta12" -> 0, "Delta34" -> 0, "Threads" -> 4, "Precision" -> 1000, "Directory" -> "/tigress/sdempsey/bootstrap/n4sym/blocks"};

RunScalarBlocks[spins_, opt:OptionsPattern[]] := Check[LoadScalarBlocks[#,opt]& /@ spins,Run[TemplateApply[StringTemplate["scalar_blocks --dim `Dimension` --order `Order` --max-derivs `MaxDerivs` --spin-ranges `Spins` --poles `KeptPoles` --delta-12 `Delta12` --delta-34 `Delta34` --output-poles --num-threads `Threads` --precision `Precision` -o `Directory`"], Association@Join[Map[# -> OptionValue[#] &, Options[RunScalarBlocks][[;;,1]]],{"Spins" -> StringRiffle[ToString /@ spins, ","]}]]]];

LoadScalarBlocks[spin_, opt:OptionsPattern[]] := LoadScalarBlocks[spin, opt] = Get[TemplateApply[StringTemplate["`Directory`/zzbDerivTable-d`Dimension`-delta12-`Delta12`-delta34-`Delta34`-L`Spin`-nmax`MaxDerivs`-keptPoleOrder`KeptPoles`-order`Order`.m"], Association@Join[Map[# -> OptionValue[#] &, Options[LoadScalarBlocks][[;;,1]]],{"Spin" -> spin}]]];
