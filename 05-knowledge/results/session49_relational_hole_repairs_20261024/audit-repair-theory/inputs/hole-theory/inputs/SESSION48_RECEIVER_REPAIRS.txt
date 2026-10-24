# Adaptive repairs after a fourth frozen cheap record

Status: **PROVED conditional routing construction; FINITE-EXACT native receiver census and full-record authentication.** The underlying native packet is pinned at commit `c09b642a629b7d12b529b88c2697abaee78e73e7`. This is an extension at order eight, not an all-order supply theorem.

## 1. What is retained

There are the same 28 core sources, 54 helper sources, and 164 reserved normalized options as in the inherited packet. Every arrival occupies 82 units, or 328 inverse/converse-closed words; the reserved alternative store contains 656 words. No frozen source is deleted. A frozen record stays occupied and contributes its constant charge and population.

Keep the actual external owner equal to **4**, not 3. The literal shared-owner equations force helper368 to U and helpers776,3139 to V. Exactly two of the other 51 helpers are V. Their 1,275 pairs cover only 972 charges, so neither a complete helper-only phase nor an owner3 excursion is available.

The full paired dyadic charge is

    Q = Z/128 × F2^4,    lambda_s = (X_s mod128, Y1_s,...,Y4_s mod2).

The original exact integer core and helper populations remain separate. All costs below count actual changed per-slot owner-clock occurrences, not a union of distinct global owners. Events are simultaneous full-record replacements, not asserted decompositions into individual metagraph edges.

For clarity, write delta_s for the change of the V bit and f for an original potential. The changed source row is2*delta_s*(f_U-f_V), zero at its appropriate modulus by the authenticated target codes. The V-target row changes by-2*sum_s delta_s*(f_s-f_V), with the opposite U-target change. The exact population equality gives sum_s delta_s=0; preservation of totalQ then kills these target rows, including their diagonal terms. X uses modulus256, the radical lifts modulus4, and the complete paired dyadic theorem justifies passing from these authenticated generators to the whole original dyadic bank.

The twenty sources with actual U/V path distance two are

    176,1394,1425,1457,1477,1480,1549,1638,1691,1921,
    1974,2071,2176,2185,2581,2837,2905,2945,3000,3369.

The eight outside sources, paired labels, and switching costs are:

| Source | Label | Cost |
|---:|---:|---:|
|1532|725|4|
|1841|1592|3|
|1923|176|4|
|2087|1720|3|
|2553|1672|4|
|2835|1347|3|
|3015|487|3|
|3360|1173|3|

`verify_native.py` independently enumerates all Hamiltonian paths in all82 canonical source cells. It checks that every guarded U/V target menu is exactly the declared two full records; computes every lift in both actual resource fibres; authenticates target classes from canonical U/V and converse tournaments rather than solely trusting an owner-map entry; and checks actual inverse, converse, kernel, indexed-clock and moving-block data. It independently reconstructs the652-owner source frames, all shared source/partner edges, the complete eight-row guard truth table, and disjointness of all164 normalized options. Full native source/target black capacities equal twice the independently counted resource units on all164 used source/target edges. Thus a charge label is never substituted for an unauthenticated resource.

## 2. Exact one-slot criterion and the three-freeze audit

For a receiver with grade-r support A_r and holes H_r=Q minus A_r, adding one compatible existing binary slot of label lambda gives

    A'_r = A_r union (lambda + A_(r-1)),
    H'_r = H_r intersect (lambda + H_(r-1)).

Consequently the repair works at a prescribed set of grades if and only if

    lambda not in union_r (H_r - H_(r-1)).

This is an equality of actual record supports, not an estimate. A forbidden difference retains the literal grade, hole and previous-hole witness. It requires the extra slot to be free and both options to glue with every retained receiver record being counted.

The independent fresh support calculation reproduces every one of the1,351 frozen subsets of size at most three, not merely exceptional examples. For all1,140 triples, all eight actual outside sources repair grades8,9,10. The six exceptional central sections and every hole at grades7,8,9,10 agree with the original receipt; their union of forbidden differences has exactly32 labels. Four triple erasures have no complete grade; two have complete grades7 and10 but not8 or9. The old claims about nonmonotone grades and eight universal triple repairs are therefore preserved.

There is a small quantitative improvement independent of the new fourth-pin work. The14 largest actual core switching costs sum to

    3*4 + 5*3 + 6*2 = 39,

not the coarse bound14*4=56. The four largest free-helper costs sum to29. At a boundary grade where at most seven free core records are selected or unselected, direct endpoint replacement therefore costs at most68 clocks and changes at most18 slots. Frozen restrictions only reduce this bound. Combining with the inherited repaired interior gives **20-slot/68-clock global connectivity** under three cheap freezes and any additional outside freezes leaving a repair source available. The inherited20/85 statement remains valid, but is not sharp as an upper estimate.

## 3. The fourth-freeze repair hypergraph

For each four-element cheap frozen set F, let G(F) be the subset of the eight outside sources for which adjoining that source to the remaining sixteen cheap sources supplies **both grades8 and9 completely**. Every one of the4,845 sets F is checked. The exact distribution is:

| Number of sources in G(F) | Number of four-freeze sets |
|---:|---:|
|8|4522|
|7|241|
|6|47|
|5|17|
|4|9|
|3|5|
|2|4|

Thus every F has at least two choices, but no single outside source works for all F. The only two-element sets meeting every G(F) are

    {1923,2087} and {1923,2835}.

These are **protected availability pairs**: after F is known, choose one member that belongs to G(F). The theorem does not assert that both members individually work at every F, nor does it require both to be put into the receiver.

The size-two lower bound is sharp, because every singleton fails for some F. The uniform cardinality guarantee is also sharp: any seven remaining outside sources meet G(F), but six need not. For example,

    F = {1425,1480,2071,2837},    G(F) = {1923,3015}.

Freezing those two outside sources leaves six available but removes every successful single-slot repair. This is failure of the specified17-slot complete-grade repair mechanism, not a proof that the entire original fibre is disconnected.

The other three four-freeze sets with only two successful sources are

    {1425,1638,1974,2071},
    {1425,1638,2071,2837},
    {1425,2071,2837,3369}.

Their successful pair is also{1923,3015}. `certificate.json` retains all4,845 actual G(F), and a literal residual hole for every failed actual source.

An arbitrary two-slot enlargement does not silently fix the failure. Across all4,845*28=135,660 choices of F and two actual outside repair slots, exactly six enlarged18-slot receivers fail completeness at some grade8,9,10. The certificate gives each pair and hole. Any pair containing a successful singleton does work: the17-slot grades8 and9 then imply the18-slot grades8,9,10 by the same adjacent-grade union. The remaining pairs are tested directly.

## 4. Conditional routing theorem after four freezes

Freeze any four cheap full records, and any a of the eight outside full records, with0<=a<=7. Let the available outside set meet G(F), and choose an actual available p in that intersection. Write c for its literal cost, either3 or4, and e=7-a. Reassign p to the receiver; the free core now consists of a17-slot receiver R and e outside slots E. Nothing is newly occupied or duplicated. Let t be the number of V records among these free core slots, after subtracting the frozen population and charge.

**Theorem.** Every nonempty fibre of the fixed external-owner4, fixed-helper-population-four packet connects with at most **18 changed slots and68 clocks per event**. For

    8 <= t <= 9+e,

all residual Q charges occur, and the stronger event bound is **18 slots and45+c clocks**, namely48 if c=3 and49 if c=4.

In particular, preserving either certified availability pair allows any other six outside records to be frozen, in addition to the four cheap records: ten frozen records with a uniform guarantee. If the frozen pattern is chosen after selecting a good repair p, seven outside records may be frozen: eleven in total, with a17-slot free receiver complete at grades8 and9. Unlike the triple-freeze theorem, an arbitrary sole surviving outside source is not enough.

**Proof.** The free-helper skeleton is the Johnson graph on two-element subsets of51 private helpers; each edge changes two helpers and costs at most15. Every core frame is private, and all alternative resource units are disjoint, so fixing additional actual U/V records creates no cross-factor. The census supplies every receiver debt at exact grades8 and9.

Apply the receiver migration construction with(h,l,u,e)=(17,8,9,e). Normalize both literal endpoints to receiver grade r0=max(8,t-e), which is8 or9. If necessary, enter the complete band by moving selections between receiver and outside while refreshing the receiver at the exact new residual charge. Move down one grade when needed, then use outside Johnson exchanges, then helper Johnson exchanges, always selecting an actual receiver witness at the arrival. Finally restore the exact target receiver record and reverse its normalization. Both populations, the shared owner and the physical paired charge hold at every arrival; an optional same-charge receiver reset avoids composition of unrelated pairwise minima.

For receiver grades r,s, the maximum support change is

    min(r+s,34-r-s).

At r=s=8 or9 this is16, with actual clock sum at most30+c. A helper edge gives16+2=18 slots and45+c clocks. A consecutive9-to8 receiver change uses at most17 receiver slots plus one outside slot, costing at most32+c+4=36+c. An outside exchange costs at most30+c+8=38+c, also at most18 slots. For entry from r<8, use the bound W_R(r+8)+4(8-r); for r>9 use W_R(25-r)+4(r-9), where W_R sums the largest clamped number of receiver costs. Since e<=7, these are at most44+c, and their slot support is at most16. All are bounded by45+c. These estimates apply to every chosen section, so the construction composes with exact endpoints.

Outside the interior band, t<=7 or(17+e)-t<=7. Direct replacement of two feasible endpoints then changes at most14 free core records and four free helpers. Their actual costs are at most39+29=68. Both endpoints satisfy every guard and invariant, so this is a legal atomic event. Complementation is used only to bound symmetric difference; it is not itself a charge-preserving move. QED.

A sharper source-dependent global bound is

    max(45+c, W_free_core(14)+29).

When exactly the protected pair remains free (the ten-frozen guarantee), its costs are4 and3, so W_free_core(14)=4+3+12*2=31. The global bound then specializes to **18 slots/60 clocks**. When just one successful source p remains free after eleven freezes, it is **18 slots/(55+c) clocks**, hence58 or59. The displayed uniform68 is convenient and is not claimed optimal. None of these boundary statements asserts charge coverage at every boundary grade.

## 5. A sharp mixed-relation warning and the correct generalization

For an arbitrary fixed boundary b, a binary repair record may constrain the old receiver differently in its U and V branches. Let A^U_(b,r) and A^V_(b,r-1) be the charges of old full receiver records that really glue with those respective branches. The exact support is then

    A^U_(b,r) union (lambda + A^V_(b,r-1)).

It is complete exactly when

    (Q minus A^U_(b,r)) intersect
    (lambda + (Q minus A^V_(b,r-1))) is empty.

Thus the forbidden differences are between the **two distinct conditional hole sets**. They may depend on the actual repair source and its full records. The unary formula is recovered only when both branch restrictions leave the same old receiver domain. One must not reuse a unary charge certificate after adding a relation.

A literal native full-record counterexample uses the cheap freezes

    {176,1425,2071},

and the universally successful triple repair1841. Add the record relation

    target(1841) = target(1457).

This is a same-target relation on actual U/V records, not asserted to be a new indexed-owner equality in the fixed652-owner host. If B is the bank of the remaining sixteen cheap slots, the linked pair has options00 and11, so its exact grade-r support is

    B_r union (lambda_1841 + lambda_1457 + B_(r-2)).

Its holes are exactly{477} at grade8, none at grade9, and{1098} at grade10. Both the group-algebra calculation and a separate literal enumeration of all131,072 full assignments verify this. Merely treating the linked pair as a freely available grade-one repair would manufacture a forbidden population/record state.

## 6. Scope and evidence

Fresh calculations here: all82 native source menus/kernels and host/resource facts; every erasure set through size four; all actual single-slot repairs; the exact availability transversals; two-slot failure witnesses; the relational counterexample; and every finite grade/cost case used by the route bounds. These are finite universal quantifiers because all declared cases are exhausted.

The original full dyadic saturation theorem is not inferred solely from checking LX=0 and LY=0. The earlier independent complete filtration audit is bundled unchanged, with byte-identical original merged edges, marked potentials and label data; `replay_filtration.py` re-executes it. Its independently computed2-primary group is(Z/2)^3358+(Z/4)^4+Z/256, whose paired image isQ. It proves all paired dyadic characters, not odd-primary preservation or the full critical-group exponent. The Y0/Y4 coordinate relation is explicitly verified rather than identifying incompatible four-bit displays.

These results concern the original order-eight aggregate physical rows, including target diagonals and the paired-black factor two. They do not imply preservation of ordinary per-port phases, exact integer-X currents, arbitrary extra relations, all-order native supply, or a Hamiltonian observer for the whole652-owner host. Resource capacity is global at each arrival and allows a released option to be used again later. **Strict temporal nonreturn is a different constraint and is not proved.**
