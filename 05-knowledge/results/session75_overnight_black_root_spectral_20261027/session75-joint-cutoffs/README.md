# Joint native-depth and integer-clock cutoffs

**Final finite result:** all245826 breakpoint rectangles at the19 explicitly
registered scales satisfy `D_plus<=S_N/2`. This is **not** a scan of every
integer N between4096 and262144. The full denominator retains1000908 original
source labels across these separate scale-specific families. The largest
observed ratio is `4132/44381` at N18275. Every scan and fresh literal receipt
matches its `python -O` replay; no tested source is unresolved.

Joint cutoffs are genuinely stronger: their maximum exceeds both separate
axis maxima at N4096,4097,65536,131072,262144. The authentic clipping hostile
shows a debt increase56 after excluding norm22; the exact boundary-pair
repair pays82. The repair is proved generally and is sharp on an explicitly
non-native scalar family. The depth-only R3 target and the stronger joint
target remain **OPEN**. See [finite coverage](COVERAGE.md) and `SUMMARY.json`.

## Scope and source binding

This tests a **stronger extension**, not a reformulation, of the open depth-only
R3 target. The input repository pin is
`96e35642110104cc547516a63f37eebcf930f707`; complete source extractions and their
hashes are preserved in `SOURCE_PINS.json`. The 627-state frontier is unchanged.
The starting N4096 family is authenticated against producer revision 2 of the
Session75 dyadic packet, manifest SHA-256
`ad8d60b20884230f7bd2180e22115b16f970d4cd0aa284a1149d45cd51f321e8`.
No newer source law has been silently substituted.

The first registration names N4096; ten structured scales; and optional powers
131072,262144. The separately frozen phase-2 registration names the inherited
full-arrival adverse scale18275 and B-chain control11002 with their immediate
neighbors. It was written before any affected joint scan. It does not claim
these were selected independently of the inherited one-resource data.

For every original odd seed `s<4N`, `s !=5 mod8`, retain
`(s,u,p,t,d,A_N(s),v_s,m)` with `m=(s−1)/2` and `v_s=(-1)^t A_N(s)`.
The source charge is the complete original B-chain sum. The literal map uses
`a=2` at `s=1 mod8`, `a=1` otherwise, followed by inverse braids until primitive;
the exact integer edge clock is `a+2j`. It checks
`2^t(3P(s)+1)=3(3s+1)+2^a`. Terminal and entry-port labels remain in every record.
Python integer arithmetic has no fixed-width overflow.

For a rectangle `(D,T)`, arrived means **both** `d_s<=D` and `t_s<=T`.
For every `q=1,2,4,...` through the first `q>=2N`, form

`E_q(D,T) = sum_(u,r modq) Q(sum_(s arrived,u_s=u,m_s=r modq) v_s)`.

There is no cross-terminal vector summation. `D_plus=sum_q max(E_q-E_2q,0)`
takes its positive part **after the entire level sum**, not per bin, port,
terminal, or pair. The denominator is always the original
`S_N=sum_(all source s) Q(A_N(s))`, including moving sources and zero labels.

Every distinct global depth breakpoint is crossed with every distinct global
integer-clock breakpoint. Empty cutoffs have zero debt. Between breakpoints,
all arrived sets are constant; above the maxima they remain constant. Therefore
this is all integer depth/clock rectangles for each completed scale, not a
sample of clock endpoints for selected depths. Distinct cutoffs may have the
same selected source set and are deliberately retained in the receipt.

## Exact two-resource pair identity

For distinct sources `s,z` at the same terminal, let

- `k(s,z)=v_2(m_s-m_z)`;
- `b(s,z)=Q(v_s+v_z)-Q(v_s)-Q(v_z)`;
- `birth(s,z)=(max(d_s,d_z),max(t_s,t_z))`.

Then, at every rectangle,

`Delta_(2^k)(D,T) = sum_(same-terminal pairs with k(s,z)=k,
                        birth_1<=D,birth_2<=T) b(s,z)`.

**Proof.** Expanding the two quadratic bin energies cancels diagonal source
norms. A named pair occurs in the coarse bin but not the next fine bin exactly
at the first differing low bit. Both labels are selected exactly when the
rectangle contains their coordinatewise maximum resource pair. Summing all
such signed terms gives the identity. No positivity estimate is inserted.

This gives a fresh exact verifier: enter each pair's signed interaction at its
birth corner and take two-dimensional prefix sums per level. This verifier
checks every N4096 rectangle against the incrementally maintained energy
scan; it does not reuse the latter's energy-update recurrence. The identity
itself is algebraic and does not prove the proposed bound: its pair measures
are signed, and positive parts do not commute with rectangle restriction.

## Concrete native obstruction to reducing joint scans to axis scans

For N4096, all5917 registered depth×clock breakpoint rectangles satisfy the
stronger proposed bound. The maximum is

`D_plus=495 at (D,T)=(36,64) or (36,65), S_N=9949`.

The separate depth-only maximum is494; the separate clock-only maximum is469.
Thus joint restriction really can increase the maximum beyond **both**
separate scans. This is not a counterexample to `D_plus<=S_N/2`.

The exact clock-border mechanism is particularly small. At depth36, crossing
from clock64/65 to clock66 admits exactly source13759:

| s | terminal u | entry port p | d | t | A | v | m |
|---:|---:|---:|---:|---:|---|---|---:|
|13759|9|6371|35|66|(-1,-1)|(-1,-1)|6879|

Its only nonzero pair interactions with the previously eligible family are:

| other source | u | p | d | t | v | low-bit q | beta |
|---:|---:|---:|---:|---:|---|---:|---:|
|9|9|0|0|0|(0,1)|1|-1|
|2831|9|6371|3|13|(-1,0)|8|+1|
|10065|9|6371|5|18|(0,1)|1|-1|

Consequently `Delta_1:278→276`, `Delta_8:25→26`, and every other
whole-level delta is unchanged. Both affected deltas remain positive,
so debt drops by1, `495→494`. The first interaction crosses entry ports;
the other two retain the same port6371. This isolates the original-source
mechanism without deleting any source from `S_N`. The four labels are an
explanation of the update, **not** a replacement family or denominator.
All same-terminal zero-beta interactions remain represented in the source
ledger and the pair audit.

## Semantic hostile controls

At the N4096 joint maximizer:

- Correct whole-level debt is495.
- Taking positives at each parent bin instead gives11227, a **false violation**
  against the same S9949.
- Coherently merging different terminals gives56715, another false violation.
- Losing the exact clock sign gives209, a different current.
- Substituting odd address `s` for `m=(s−1)/2` changes the whole energy path,
  even though in this instance it happens to leave the total debt495 unchanged.
  Comparing debt alone would miss that mutation.
- Literal full-ledger comparison rejects clock changes by2 even though these
  preserve parity, entry-port changes even though debt ignores that sidecar,
  terminal/depth/address/charge mutations, deletion of a zero-charge source,
  and deletion of a nonarrived nonzero source.

The saved maximum receipt also decomposes each level's signed interaction
into same-entry-port and cross-entry-port contributions. Neither component
can be replaced by its positive part before the level sum. These are
source-resolved diagnostics, not a new source law or a probability measure.

## Verification and limitations

`scan.py` uses the prescribed O(depths × sources × dyadic levels) sweep.
`literal_verify.py` does **not** import it: it reconstructs every source from
three original real rows and uncached literal routes, checks every stored
source field, and directly re-bins the maximum plus selected controls.
For N4096 it additionally checks every rectangle by the pair-birth identity.
Detailed B-chain/route streams are hashed; complete compact source ledgers
are retained. These two implementations were authored by the same producer:
the second implementation is not claimed to be an independent-agent audit.

Normal and `python -O` receipts are compared byte-for-byte. Every check uses
explicit exception guards, not `assert`. A cycle or unfinished route is a
hard incomplete-run error, never a reason to erase its mass. Completed
receipts report zero unresolved sources. Final numerical coverage belongs
in `SUMMARY.json`; no registered but incomplete N is counted as tested.

All results remain finite. No monotonicity in N or cutoff, universal-entry
claim, asymptotic inference, or proof of depth-only R3 follows from them.
The actual forcing cost is not substituted for coherent source debt. No
external publication or changes to shared repository files are made here.

### Port and resource diagnostics at the larger strict joint maximum

At N65536 the joint maximum is4790 at `(D,T)=(60,103)`, while the
separate depth and clock maxima are4779 and4773. The resource partition
retains97071 labels of norm157380 inside,48 labels of norm52 excluded only
by clock,17 labels of norm33 excluded only by depth, and1168 labels of
norm1692 excluded by both. These norms sum to the unchanged S159157.
Thus both resource constraints genuinely matter at this corner.

The same-port/cross-port level signs are mixed. At q128 the positive delta
3104 splits as2340+764; at q1024 the positive delta1239 splits as1095+144;
at q2048 the tiny positive delta2 splits as97−95. Already at N4096,q8,
the positive delta25 splits as−18+43. Therefore a proof cannot simply
discard all cross-port interactions as negative, or infer the sign of the
whole level from the same-port component. These are exact local diagnostic
obstructions, not a claim that no other port-based proof is possible.

### A native boundary norm alone does not pay for clipping

A post-hoc exact consequence of the registered N4096 surface gives a second
useful, separately scoped hostile. At fixed clock64, the unrestricted-depth
ledger has debt439. Restricting depth to36 excludes exactly16 original labels
of total norm22 but increases the debt to495, a gain56>22. Hence

`D_plus(after restriction) <= D_plus(before restriction) + S(excluded)`

is **false even for this actual native source family**. The complete source
norm is still9949 in both ledgers; the norm22 is diagnostic, not a replacement
denominator. The saved hostile lists all16 excluded source records and both
literal energy paths. This observation does not test, compute, or refute a
paid-forcing identity. It rules out paying an arbitrary depth-clipping debt
increase solely by the removed sources' original norm without interaction
terms. It is not a violation of the joint half-source target or depth-only R3.

## Exact repair: pay the removed boundary interactions

For any one fixed terminal-labelled refinement tree, let `F` be the full
selected family, `R` the removed labels, and `K=F\R` the retained family.
Keep removed labels as zero leaves in the same partition when comparing.
Write `Delta_k(F)` and `Delta_k(K)` for their whole-level signed drops and

`delta_k = Delta_k(F) − Delta_k(K)`.

Then the **general exact repair** is

`D_plus(K) − D_plus(F) <= sum_k max(−delta_k,0)`.

**Proof.** For every real x,delta,
`max(x−delta,0)−max(x,0)<=max(−delta,0)`; apply this with
`x=Delta_k(F)` and sum. By the named-pair/LCA identity, `delta_k` is exactly
the sum of `beta(v_s,v_z)` over same-terminal pairs allocated to level k
with at least one label in R. Pairs with both labels in R appear once,
not twice. Source diagonals cancel from each level difference. Thus the
right-hand side is a **boundary-pair invoice**, not removed-source self-energy.
It retains both terminal and entry-port metadata and charges each pair at
its one actual low-bit LCA.

For the actual N4096 depth-clipping hostile, the exact delta vector is

`[-34,-19,-1,-11,13,-6,-8,2,2,0,-3,0,0]`.

The repaired invoice is82, the actual debt increase is56, and the removed
source norm is22. The receipt explicitly enumerates every boundary pair,
including zero-beta pairs, with terminal, both ports, low-bit q and sign.
The slack26 in the repaired inequality comes from positive-part saturation,
not an unaccounted source or changed denominator.

### Sharp generic family: no uniform removed-self-norm replacement

This paragraph concerns a **generic scalar family**, not the prescribed
native `A_N(s)` family and not another native numerical scale sweep.
Let `n=2^h>=4` leaves occupy a complete balanced low-bit tree, all in one
terminal, with one scalar charge−1 and n−1 scalar charges+1. Remove the−1
by replacing it by zero while retaining its conceptual leaf. At any bin
size a, every bin has a labels; one bin has sum a−2 before removal and
a−1 after, while all other bins have sum a. Therefore

`E_F(a)=(n/a−1)a^2+(a−2)^2=(n−4)a+4`,

`E_K(a)=(n/a−1)a^2+(a−1)^2=(n−2)a+1`.

All coarse-to-fine level drops are nonnegative. Consequently

`D_plus(F)=(n−4)(n−1)`,

`D_plus(K)=(n−2)(n−1)`,

`D_plus(K)−D_plus(F)=2(n−1)`.

The removed norm is1. At the level with bin size a, `delta=-a`, so
`sum max(-delta,0)=n+n/2+...+2=2(n−1)`: the repaired inequality is sharp.
No constant times removed self-energy can replace this boundary invoice
uniformly over arbitrary labelled families. The proof holds for every
h>=2; direct exact checks cover h=2,...,12. The native N4096 example above
is the separate evidence that source-norm-only clipping already fails in
the authentic source family; the generic family supplies unbounded
sharpness only in its explicitly broader class.

## Reproduction

The code uses Python standard library only. Complete rerun:

```sh
sh run_all.sh
```

This deliberately repeats the full registered finite scan under normal and
optimized Python and can take substantial time. For a smaller regression:

```sh
python scan.py 4096
python -O scan.py 4096 --suffix .optimized
python literal_verify.py 4096
python -O literal_verify.py 4096 --suffix .optimized
cmp results/receipt_N4096.json results/receipt_N4096.optimized.json
cmp results/literal_receipt_N4096.json results/literal_receipt_N4096.optimized.json
```

The additionally saved `pair_check.py 65536` receipt independently recomputes
all23160 N65536 rectangle debts from21818296 same-terminal named pairs,
of which17979937 have nonzero beta. Its normal and optimized receipts agree.
This complements the all5917-rectangle pair check at N4096 and does not
silently claim a full pair-enumeration audit at every larger scale.

At the largest declared scale262144, literal reconstruction records an
exterior state28167930305. This exceeds32-bit unsigned range; arbitrary
precision is material here, not only a theoretical precaution. No tested
source is unresolved. Retaining its complete original norm in S_N would
remain mandatory even if a future route were unresolved; the checker would
then stop as incomplete rather than silently delete it.
