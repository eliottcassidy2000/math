# A three-arc repair reaches the next ready target, with an exact module law

**PROVED with an exact finite run-skeleton certificate:** the arbitrary-module formulas, literal repair family, and stated closure obstructions. **FINITE-EXACT:** complete local star probes and the selected graph counts. The universal indegree-one unit-ready405 ray is REFUTED at k17; Chan-Pak remains OPEN.

**Incoming correction to the research target (commit c85ab9a188):** an all-k indegree-one unit-ready405 ray is **REFUTED**. The [complete intrinsic inverse](../../target_directed_cores_20261112/TARGET_INVERSE.md) and its [independent exhaustive audit](../../target_directed_cores_20261112/ports/COMPLETE_INVERSE_AUDIT.md) exclude every tournament order at53,084,231=405*2^17+71. The ordinary Hamilton count is attained, and Chan-Pak remains OPEN. The finite constructions and module formulas here remain valid; future continuation must change the required interface or skip an obstructed target.

## Inheritance and the operation being changed

The anchor is a physical return from H881 to a unit-ready H1691 source. The niche is the inherited [exact-next continuation fibre](../../h17_exact_repairs_20261112/history/CONTINUATION_FIBRES.md); the wildcard is a selective star switch whose intermediate full count drops sharply. The board is: actual deletion-parent count; unique-start face; a retained module; old-arc repair; final readiness; closure under another application.

The closest positive mechanism is the [A2 core repair](../../changing_core_returns_20261113/preparation/NEUTRAL_881.md#n2-a-genuinely-different-core-repairs-a-with-only-two-arcs). The canonical hostile is F>=H: a unit-ready canonical newborn has a unique-start parent and count H+F-1>=2H-1. Thus no induced growth from H881 can directly end in a unit-ready1691 newborn. The corrected near miss is its B9 sibling, whose birth lands at1689 with large readiness defect. The least-used sidecar is the full count of the graph obtained by deleting the final root: it can decrease under an old-arc repair even when the final returned target increases.

Unlike the earlier same-H preparation, the operation below deliberately changes the old full count. No claim of neutrality or a native-only trace is made for its three individual arc reversals. All old vertices and labels are retained through that switch; one canonical newborn is then added at a specified physical port.

## S1. The literal 881 ->153 ->1691 return

Start with the actual A2 source

```
A=(518,388,456,1011,487,455,387,257,1,502),    H(A)=881.
```

Select old vertex0 as the preparation port. Reverse precisely the three arcs on unordered pairs {0,2},{0,7},{0,9}. The resulting strong parent is

```
P=(130,388,457,1011,487,455,387,256,1,503),
(H,s,q,F)(P,0)=(153,1,27,1539).
```

Add a new vertex beaten only by0, and mark the newborn. Its literal rows and state are

```
T=(1154,388,457,1011,487,455,387,256,1,503,1022),
root10,     (H,s,q,F)=(1691,152,153,1692).
```

The new root has indegree1 and q-s=F-H=1, so this is an actual strong unit-ready source at the requested next target. The intermediate count drop881->153 pays the previously unavoidable parent-count budget. Selecting port0 is part of the declared operation; the original A2 active root9 is not asserted to remain active.

Among all old-arc patches of widths0..2 on A, at every selected physical port, none can give a unique-start parent with H+F=1692. At width3 exactly this patch succeeds. The C++ search separately shows that preserving the old active root9 leaves no start-count1 parent even through width4. These finite minimum statements concern a patch before one canonical birth, not arbitrary helper programs.

## S2. Replace the surviving transitive pair by any tournament

Contract the untouched homogeneous pair {4,5} in A to a module vertex8, placing the other old names in order(0,1,2,3,6,7,8,9). The resulting source quotient and repaired quotient are

```
Q =(134,100,120,499,99,65,1,374,119),
Q'=(34,100,121,499,99,64,1,375,119).
```

The switch in these indices is {0,2},{0,5},{0,7}. Replace quotient vertex8 by an arbitrary nonempty tournament B with uniform exterior incidence. Keep port0. Both quotients are strong, and substitution of a nonempty tournament into a vertex of a strong quotient preserves strongness.

Let K_j(B) count ordered spanning covers of B by j nonempty directed paths. Write u=K_1(B), v=K_2(B), w=K_3(B). The exact states before and after the switch are

```
source:
 H=369u+256v+52w,    s=8u+3v,
 q=15u+6v,          F=766u+630v+156w;

repaired parent:
 H=85u+34v,         s=u,
 q=15u+6v,          F=591u+474v+114w.                 (1)
```

Consequently the canonical newborn has

```
H_new=675u+508v+114w,
s_new=84u+34v,       q_new=85u+34v,
F_new=676u+508v+114w,       q_new-s_new=u.             (2)
```

The switch always lowers H by284u+222v+52w>0. The full switch-plus-birth macro always raises the starting H by306u+252v+62w>0. Thus the same physical operation converts a controlled temporary decrease into a larger returned count. Its final readiness defect is exactly H(B), rather than an inferred scalar property.

### Proof by run skeletons

Every entry into B comes either at the start of a path or after quotient vertex3 or7. Therefore a full path has at most three B-runs. In Q', vertex7 has sole in-neighbour3. This prevents independent internal entries3->B and7->B unless one path starts at7; hence repaired full paths and root-deleted full paths have at most two runs. A source path starting at0 has at most two runs; a repaired path starting at0 has at most one.

For F, the first path ends at0. In Q the only in-neighbours of7 are0 and3; vertex0 cannot precede7 because it ends the first path. Thus7 either begins one of the two paths or is preceded by3. In the first case one of the two path-start opportunities for a B-run is lost; in the second case3 cannot also precede a B-run. The same argument applies to Q'. Therefore every relevant two-cover has at most three B-runs.

Contract each maximal B-run to an identical module symbol. Number the symbols by their occurrence in the ordered path or ordered pair. For any fixed legal skeleton, expansions are exactly K_j(B); contraction and expansion are inverse. This gives a finite coefficient problem, independent of B and its order. The coefficients are

| j | source (H,s,q,F) | repaired (H,s,q,F) |
|---:|---|---|
|1|(369,8,15,766)|(85,1,15,591)|
|2|(256,3,6,630)|(34,0,6,474)|
|3|(52,0,0,156)|(0,0,0,114)|

The checked recurrence in `star_family.py` counts those skeletons directly. Its state is(outside subset,number of module symbols,last symbol); consecutive module symbols within one path are forbidden. A rooted two-cover is obtained by splitting the outside subset and the module-symbol budget between a first path ending at0 and a second arbitrary path, allowing that second path to be empty. This is an exact finite combinatorial certificate, not an extrapolation from sample module widths. A separate literal-graph subset counter verifies the resulting formulas for every labelled B through order3 and the transitive family through width6. QED.

## S3. An infinite family of actual unit-ready returns

A tournament has H(B)=1 exactly when it is transitive. In (2), unit readiness therefore holds exactly for B=TT_m, m>=1. Substituting

```
K_1=1,  K_2=2^m-2,  K_3=3^m-3*2^m+3
```

gives the source and returned counts

```
A_m=52*3^m+100*2^m+13,
B_m=114*3^m+166*2^m+1.                              (3)
```

The repaired parent has

```
(H,s,q,F)=(34*2^m+17,1,6*2^m+3,
           114*3^m+132*2^m-15).
```

The first returned unit-ready targets are675,1691,4407,11891,33015,93731. The same three specified old arcs work at every width, followed by the same single newborn operation. This is an all-parameter physical family, rather than a replacement by a pre-enumerated target graph.

The returned unit-ready family is physically renewable: insert another vertex into its retained transitive module, with the same uniform outside incidence and its prescribed place in the total order. This gives the next returned graph while retaining the active root, its unique in-neighbour, and readiness1. Its deletion count grows as t_m=34*2^m+17; the first step from1691 reaches4407 rather than3311.

This renewal does not reuse the original three-arc repair source type. For every m>=1,

```
B_m-A_m=62*3^m+66*2^m-12>0,
A_(m+1)-B_m=42*3^m+34*2^m+12>0.                    (4)
```

Thus its output falls strictly between consecutive source counts. In fact every tournament B' of order m+1 has K_j(B')>=K_j(TT_(m+1)), since each ordered vertex partition has at least one path in each induced part. All source coefficients in(1) are nonnegative. Hence the returned graph cannot re-enter even the arbitrary-B source family at its new order: the smallest source at that order already exceeds B_m.

The405-ray coincidence is also isolated. A_m-71 has two-adic valuation exactly1, so a source on405*2^k+71 requires k1 and then m2. For outputs, reduction modulo3 first forces m even. For every even m>=2, B_m-71 has valuation exactly2, forcing k2; strict growth then forces m2. Therefore this family intersects the requested consecutive ray targets precisely at881->1691. It does not provide the next3311 return by increasing m.

## S4. Complete local failure at the next target

For any fixed graph T and selected port x, put U=T-x. If the new star makes x beat O, then

```
s_x=sum_(u in O) S_U(u).
```

Consequently s_x=1 exactly when O contains one vertex with S_U=1, any subset of vertices with S_U=0, and no other vertex. This exhausts every star orientation at every width; the coordinator's separate directed-cut kernel gives its exact target response.

For A, only deletion at port0 has any such choices: its deletion count is27, the one-start entries are1,2,6, and zero-start vertices are7,8. The12 faces give exactly the S1 repair above.

For the produced1691 graph, only deletions at ports2 and10 qualify. Each has three one-start entries and no zero-start vertices: six total choices. None gives H_parent+F_parent-1=3311. Allowing one arbitrary old-vertex deletion before the star switch gives656 face choices across all11 deletions, again with no3311 hit. These are complete obstructions for these declared menus, including weak parent/output possibilities. A subsequent exhaustive two-star scan also gives no3311 hit:56,320 orientations of a second old star inside a root-deleted graph, followed by all target-eligible unique-start faces, account for14,816 full parent counts. Weak parents and outputs were allowed. These obstructions do not forbid other selective edge sets, extra helpers, or a different current1691 representative.

## Reproduction and scope

Run `star_family.py` normally and with `-O -B`; the two receipts agree. The literal star-face searches are `star_successor_probe.py` and `delete_star_3311.py`. The two C++ patch probes retain their exact input and width universes in their raw receipts; the consolidated `search_replay.py` compiles and exactly repeats every declared search, then independently recounts every positive parent and newborn.

Source -> target is a two-port source with a uniform module -> a selective old-star repair -> a canonical unit-ready newborn when that module is transitive. The operation preserves all surviving old names and the uniform module boundary, changes the selected star and the intermediate H, and adds one vertex. Its required sidecars are the selected port, three nonempty-cover counts, and final readiness. The arithmetic family is not closed under the returned graph, and no Chan-Pak coverage claim follows.
