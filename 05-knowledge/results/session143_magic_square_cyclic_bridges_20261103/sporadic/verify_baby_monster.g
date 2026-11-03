# Optional independent data rerun: gap -q verify_baby_monster.g
# NOT EXECUTED in the session: no GAP binary was installed.
if LoadPackage("ctbllib", false) = fail then
  Error("GAP CTblLib package is required");
fi;
t := CharacterTable("B");;
names := ClassNames(t, "ATLAS");;
pos := List(["5A", "5B"], x -> Position(names, x));;
chars := Filtered(Irr(t), chi -> chi[1] = 4371);;
if Length(chars) <> 1 then Error("4371-character uniqueness failed"); fi;
chi := chars[1];;
if chi{pos} <> [21,-4] then Error("unexpected order-five values"); fi;
if not ForAll(pos, i -> ForAll([1..4], j -> PowerMap(t,j)[i] = i)) then
  Error("order-five classes are not individually rational");
fi;
if not ForAll(chi, IsInt) then Error("4371-character is not integer-valued"); fi;
if Minimum(Filtered(List(Irr(t), chi -> chi[1]), d -> d > 1)) <> 4371 then
  Error("unexpected minimum nontrivial ordinary degree");
fi;
Print("GAP ", GAPInfo.Version, "\n");
Print("CTblLib ", PackageInfo("ctbllib")[1].Version, "\n");
Print("class positions = ",pos,"\nvalues = ",chi{pos},"\n");
Print("centralizers = ",SizesCentralizers(t){pos},"\n");
for i in pos do
  mult := List([0..4], r -> Sum([0..4], j ->
    chi[PowerMap(t,j)[i]] * E(5)^(-r*j))/5);
  Print(names[i]," multiplicities = ",mult,"\n");
od;
Print("PASS\n");
QUIT;
