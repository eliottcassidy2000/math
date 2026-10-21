# Same-parent annulus obstruction and three-pair temporal guards

2026-10-20. **PROVED, scoped:** complete same-parent classification in the
64-scale new annulus, its six charge bands, absence of an exact cancelling
same-parent pair, and an all-horizon guard for three selected interaction
terms. **FINITE-EXACT:** the declared native controls. **OPEN:** the bound
on the remaining interactions and universal Collatz coverage. No priority
claim. Files: [checker](check.py), [receipt](receipt.json).

## 1. Inheritance and the intended transfer

The nearest mechanism is S2 in [source orientation and boundary
contexts](../source_orientation_boundary_contexts_20261018k.md): near
siblings x,2x+1 share a first-arrival depth and retain clocks t,t+1, once
neither is already in the stopping frontier. The complete sibling tensor
in S1 has an indefinite leading term; a trace cancellation is not a
positive contraction. The exact153N source bank and owned interactions
are inherited from [session18](../session18_native_annulus/README.md).
The canonical hostile is N7 at depth8, and the least-used coordinate is
the birth horizon of a particular unordered source-pair interaction.
The inverse formulas and their period-three packaging are inherited from
[THM-4690](../../../01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md)
and the incoming [native phase lifts](../session18_native_phase_lifts/README.md),
sections2 and5. The additional conclusion here is their restriction to
one new annulus, followed by the exact charge/interaction classification.

The board is: original source; one-step parent; charge band; selected pair
address; actual clock/depth; shared old core. Anchor: guarded interaction
payment. Niche: complete inverse-parent classification. Wildcard: the
paper's common-space interpolation, used to ask for small exact tables
rather than to alter native source weights. The cheap probe rejects the
attractive conjecture that same-parent sources give ready-made zero blocks.

## 2. Every same-parent pair in the new annulus

Primitive sources are positive odd s with s !=5 mod8. The canonical map P
uses a=2 at s=1 mod8, otherwise a=1, then removes all B(x)=4x+1 copies
from (3s+1)/2^a until a primitive parent u remains. If j copies are removed,
the actual edge clock is a+2j. Inverting gives

    s_(1,j)=((6u+2)4^j-5)/9,
    s_(2,j)=((12u+4)4^j-7)/9.                         (1)

The allowable j>=0 have these residue classes modulo3:

| u mod3 | a=1 | a=2 |
|---:|---:|---:|
|0|2|1|
|1|1|0|
|2|0|2|

These follow by divisibility modulo9; the formulas also give the required
source residues modulo8. Moving j by3 sends a source to64s+35 for a=1,
or64s+49 for a=2. Consequently each branch has at most one member of
J_N=[4N,256N), for every integer N>=1. A parent has at most two new sources.

If both branches occur, j_(a=1)-j_(a=2)=1 modulo3. Adjacent
members have one of the following forms; skipping to another branch
member multiplies size by at least64 and leaves J_N:

* near: x=1 mod8, y=2x+1, clocks t,t+1;
* far: x=3 or7 mod8, y=32x+17, clocks t,t+5.

Both statements can also be checked directly using P and B. These are a
complete classification, not just two sufficient constructions. The
boundary ROOT self-loop is harmless in (1), but is not silently traversed:
the first-arrival argument below explicitly excludes sources already in V.

## 3. The exact charge table rules out same-parent cancellation

Use L=((1,2),(-2,-1),(1,-1)) and U=((0,1),(-1,-1),(1,0)), with indices
modulo3. For s in J_N put r=s mod3. The enlarged charge A_(64N)(s) is:

| Source interval | Charge |
|---|---|
|[4N,8N)|0|
|[8N,16N)|L_r+L_(r+1)+U_(r+2)|
|[16N,32N)|L_r+L_(r+1)|
|[32N,64N)|L_r+U_(r+1)|
|[64N,128N)|L_r|
|[128N,256N)|U_r|

Proof: the B-chain has at most three points. For j=0,1,2 its lower and
upper band tests are 4^j s+(4^j-1)/3 <128N or256N. Since s is odd and
each relevant threshold divided by4^j is even, these are exactly the
displayed interval tests, with no rounded endpoint exception. The first
row uses L_0+L_1+L_2=0. This proves the table for every N, independently
of any assertion about orbit convergence.

For a near pair, x traverses the first five bands and y the next band.
Their respective norms are (0,1),(1,3),(3,1),(1,3),(3,1). A far pair has
norms (0,1). They are never equal. Thus **no same-parent pair of new
sources can be an exact opposite-charge zero packet**. When neither is
in V their equal first-arrival depths and opposite clock signs turn the
pair vector into a signed A_x-A_y, which is nonzero. This identifies the
failure mechanism: synchronized transport exists, but charge equality does not.

## 4. Three owned interactions can pay together at every horizon

For a near pair outside V the selected within-pair interaction is

    e(x)=-beta(A_x,A_(2x+1)),
    beta(v,w)=Q(v+w)-Q(v)-Q(w),  Q(a,b)=a^2-ab+b^2.

It is born at the common actual arrival depth d(x). The exact table is

| x interval | r=0 | r=1 | r=2 |
|---|---:|---:|---:|
|[4N,8N)|0|0|0|
|[8N,16N)|3|0|-3|
|[16N,32N)|0|3|-3|
|[32N,64N)|0|3|-3|
|[64N,128N)|3|0|-3|

Choose x,x+8,x+16 in a single band, all1 mod8, with all six source
endpoints outside V. They give three disjoint pair addresses and all
three residues modulo3. Except in the zero band, the sum of their
selected interactions at horizon d is exactly

    3*1[d_+<=d] - 3*1[d_-<=d].                         (2)

It is nonpositive at **every finite horizon iff d_-<=d_+**. If the plus
pair never arrives, the guard holds without assuming arrival of the
minus pair. If the plus pair does arrive, the certificate must supply the
minus pair's actual arrival no later. This is the strongest unconditional
form available from these three terms alone. Completed cancellation of
the terms does not suffice: the reverse arrival order violates (2).

The certificate preserves six original sources, their actual terminals,
full integer clocks, arrival depths and compatible mod15/mod35 registers.
It does not identify their different terminals or infer clocks from a
register. It retires only three selected unordered-pair terms in I_d.
Every interaction between different sibling pairs, with other sources,
or with the old core remains in the ledger. In particular this is not a
claim that the six-source vector packet meets its source budget.

**Composition with chronological opposite matching.** Let D_d<=0 be
the sum of disjoint guarded sibling-interaction triples. Let Z_d be the
mature opposite-pair credit from [temporal owned
blocks](../owned_pair_blocks_20261020i/README.md), and define
J_d=I_d+Z_d-D_d. The selected nonzero sibling interactions have source
norms1 and3, whereas an opposite pair has equal source norms. Their
unordered pair addresses are therefore distinct, even if they share a
source endpoint. Hence the exact remaining obligation is

    I_d<=M_d  iff  J_d<=M_d+Z_d-D_d.                    (3)

No source diagonal is copied. All counterbalancing cross terms remain
in J_d; (3) alone gives no new numerical bound on the total energy.

## 5. Exact controls and next obstruction

The checker independently enumerates every primitive new source for
1<=N<=64 and checks its direct B-chain, inverse-parent classification,
and the table. At N4096 it groups disjoint consecutive near pairs within
each band, breaks runs at removed frontier sources, and checks every
actual horizon. Its21,152 triple blocks comprise676 zero blocks,
11,624 nonzero blocks satisfying the arrival guard, and8,852 failing it;
22 ungrouped pairs are retained. At depth16 the certified selected
interaction sum is -14,274. It returns to zero after all selected sources
arrive. This is a temporary identified interaction credit, not an extra
addition to153N or a new final-energy saving.

One passing block has x=32769,32777,32785, selected terms+3,-3,0 and
arrival depths9,6,5. One hostile block has x=32793,32801,32809 and terms
+3,-3,0 arriving at7,22,27: it owes3 at every cut7 through21.
The receipt stores full six-source records for both.

Reproduce with

```text
python 05-knowledge/results/sibling_annulus_guards_20261020i/check.py
python -O 05-knowledge/results/sibling_annulus_guards_20261020i/check.py
```

The next meaningful target is the remaining interaction J_d, especially
the8,852 reversed-order blocks. A bridge to a larger paid block must keep
its pair addresses and arrival interval; treating their final zero as
present credit repeats the temporal hostile. The source norm bank stays
exactly153N throughout.

The incoming phase-lift packet gives a complementary positive mechanism:
three successive same-branch lifts can be a synchronized zero packet when
their sparse band defects vanish. They span two scale annuli. This is
compatible with our one-annulus obstruction, where there is at most one
member per branch. A phase triple which mixes old and new sources still
requires the inherited old-baseline correction; its three norms cannot
all be credited to this annulus's153N bank.

## 6. A synchronized mixed packet pays two old/new terms

The incoming three-phase mechanism also gives a positive use of the mixed
baseline. Suppose an admitted, source-retained phase triple has signed
vectors v,Rv,R^2v, the first two old and the last new, and satisfies the
outside-V and exact phase-defect guards. All three have the same actual
first-arrival depth. Put w=R^2v. Their selected old/new interactions are

    beta(v,w)+beta(Rv,w)=beta(-w,w)=-2Q(w).             (4)

They appear together at the common arrival horizon. The packet's own
energy increment is -Q(w), while its new source norm is only Q(w).
If Z is the rest of the old arrived vector in that terminal, the actual
increment on adjoining w is -Q(w)+beta(Z,w). That last term is retained,
as are interactions with all other new sources. This is the inherited
mixed-packet identity specialized to a useful phased configuration;
it never spends the two old source norms against the new bank.

Two incoming literal triples become such certificates at annulus base
N=16384 (old cutoff65536, enlarged source cutoff4194304):

| Old sources | New source | Common terminal | Depth | Full clocks |
|---|---:|---:|---:|---|
|1009,64625|4136049|47|1|6,12,18|
|1015,64995|4159715|571|2|4,10,16|

Each new source has norm1 and selected interaction -2. Our checker
reauthenticates these actual sources, bands, clocks and first arrivals.
They give two additional all-horizon guarded blocks, not a new scale
census. These old/new pair addresses are automatically distinct from
the new/new addresses in sections4–5 and from chronological opposite
matching. Within any proposed larger family, pair-address reuse must
still be excluded. The untouched term beta(Z,w) is the next obligation.

Normal and optimized saved-receipt runs agree (3,841,679 mathematical
checks, plus the read-only receipt comparison). LF SHA256:
`30611d37a7c7509087c127f8dc27f172c0c6b2b7990fe8e0466eb9a74f331237`.
Independent peer proof/API review passed the complete inverse/band and
pair-address arguments. A separate charge/route reader reproduced both
mixed packets and their old-baseline increments; no repair was required.
