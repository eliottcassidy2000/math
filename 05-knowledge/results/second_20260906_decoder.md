# A genuine balanced decoder entry defeats automatic minimum-coordinate closure

**Status: FINITE-EXACT actual-entry witnesses + INDEPENDENTLY AUDITED, with
elementary proofs of entry and failure of the stated sufficient comparisons.** Two primitive
thirteen-speed rows have the actual `6+7` decoder partition, satisfy
`W_(Q,3)=V_dec` inside the physical box, contain a literal unit in the larger
primitive component, and pass every retained joint-shadow gcd profile.
Nevertheless their smaller primitive minimum exceeds the current balanced
cutoff, and the stronger native unit inequality fails too. Both rows are
strictly safe. This is an obstruction to forcing those particular numerical
conditions, not to other safety proofs or to LRC(14), which remains **OPEN**.

The [independent root audit](second_20260906_root_audit.md) accepts the
actual-entry construction, full inherited profiles, and the complete
maximum-endpoint native-failure checks.

## 1. Inheritance, connection, and hostile-first search

The closest mechanism is the independently audited
[overnight15 larger-unit theorem](overnight15_20260906_lrc_larger_unit.md).
For actual components `tV union gU` of sizes `a+b=13`, `a<=b`, with a unit
in `U`, it closes the balanced split when

    min V <= floor(3Q/28)=60,843,134,147,
    Q=91^6=567,869,252,041.

Its exact native comparison is stronger but remains sufficient. The same
report's five-vertex star stress is the recovered construction mechanism.
The least-used coordinate is common content of rooted tree cofactors: it
can leave a primitive positive shape with a large minimum.

Actual entry uses **THM-3818**,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`,
and the complete support distinction emphasized by the audited
[overnight13 entry decoder](overnight13_20260906_lrc_entry_decoder.md).
The inherited relation budget and actual physical box are those of
**THM-2052**,
`01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md`.
The exact support checks use the audited
[overnight12 signed-box theorem](overnight12_20260906_lrc_gcd_semigroup.md).
The [creative walk compiler](creative_20260906_lrc_bridge.md) retains
collective versus pair gcds but does not force a low component minimum.

The canonical hostile is the primitive primorial normalization example of
**THM-4052**, `01-canon/theorems/THM-4052-lrc14-affine-component-width-escape-cones.md`:
primitive gcd one does not control the maximum or supply a unit. Unlike
that older hostile, the rows below pass every currently retained large-subset
gcd profile. The corrected near misses are a merely displayed partition,
a selected relation rank instead of the full span, and scalar gcd caps
without the complementary gcd word. Recent unit-placement, physical-box,
and profile corrections in `01-canon/MISTAKES.md` were read.

Anchor: balanced actual decoder entry. Niche: star cofactor content and
prime separation. Wildcard: force opposite mixed-support orientations to
fail for different exact reasons. The live board has six objects: primitive
minimum, bounded decoder edges, distinguished coefficient, two-label
amplitude, hereditary gcd word, and literal safe phase. The new witness
shows the first five do not force the current balanced minimum/native
inequality; the sixth prevents a failed sufficient test being called unsafe.

Exact-coordinate, cofactor-star, minimum-bound, and current-cutoff searches
on canon and the current overnight results found the five-vertex precursor
but no identical six-vertex completion. No external novelty claim is made.

## 2. The shapes and physical scales

Put

    q=(179,181,183,185,187),  P=product(q)=205,114,343,115,
    V={P} union {(356-q_i)P/q_i : i=1,...,5},
    U={1,49,109,331,331^2,331^3,331^4}.

Explicitly,

    V=(185370716505,189592176609,193905909065,
       198314972625,202822562745,205114343115),
    U=(1,49,109,331,109561,36264691,12003612721).

The five denominators `q_i` are pairwise coprime and each is coprime to
`356-q_i`. Consequently the gcd of the six cofactors is one. Their minimum
is `169P/187=185,370,716,505`, more than three times the balanced cutoff.
Every star edge between `P` and one leaf has primitive coefficient sum
`356=2^2*89`, which belongs to the actual decoder atlas.

The seven-shape has the actual spanning tree

    1--49, 1--109, 1--331, 331--331^2--331^3--331^4,

whose primitive sums are `50=2*5^2`, `110=2*5*11`, and `332=2^2*83`.
All their prime factors are `2 mod3`, with exponents at most two.
The unit in `U` makes this shape primitive.

Take `g=1` and either of the scales

    t_odd=36,883,259,177,
    t_mix=36,883,259,192.

The first is the smallest integer strictly above
`Q*(max U+second-largest U)/min V`; the second has exact 2-adic valuation
three. Both are coprime to `49*109*331`, hence to every label of `U`.
For `n(t)=tV union U`, the physical sums are respectively

    43,342,280,629,195,004,440,991,
    43,342,280,646,821,814,650,951,

strictly below

    Q^2=322,475,487,413,604,782,665,681.

All thirteen labels are distinct and positive. The full row is primitive,
because `U` contains one. The mixed row is the principal hostile; the odd
row supplies a simpler independent positive safety control.

## 3. Proof of actual equality entry, not merely a partition

The two displayed internal spanning trees survive their respective common
scalings. There are no cross-component decoder edges: for `v in V,u in U`,
the larger reduced coordinate of `(tv,u)` is at least `tv/u>Q>355`.
Thus the actual graph has exactly these two connected components. Its
weighted incidence rows span the eleven-dimensional componentwise kernel

    V_dec={z in Q^13 : sum_(i in tV)n_i z_i=0,
                         sum_(j in U)n_j z_j=0}.

Every decoder row is itself a support-two relation of height at most 355.
It remains to exclude **all** mixed bounded relations.

For two physical `V` labels `tv_i,tv_j` and one `U` label `u`, a relation
with nonzero distinguished coefficient `C` forces

    t*gcd(v_i,v_j)/gcd(t*gcd(v_i,v_j),u) divides C.

The minimum of this cleared coefficient over all `15*7=105` selected
supports is

    218,681,470,675,839,009       for t_odd,
    218,681,470,764,774,264       for t_mix.

Both exceed `Q`, excluding this entire mixed orientation. The distinguished
coefficient cannot vanish in a relation whose support actually meets both
components. Zero coefficients on one of the pair labels are allowed, so
support-two crossings are included.

For one physical `V` label and two distinct `U` labels, the selected scales
satisfy

    t min V > Q*(max U+second-largest U).

A nonzero `V` contribution therefore has larger absolute value than any
sum of two allowed `U` contributions. This excludes the other `6*21=126`
mixed supports, again including zero pair coefficients.

Hence every support-at-most-three relation of height at most `Q` lies
entirely in one component and thus belongs to `V_dec`. The opposite
containment was supplied by the decoder edge rows. Therefore

    W_(Q,3)=V_dec,    dim W=11.

The verifier independently reconstructs all graph edges, obtains rank eleven
by rational elimination of the literal edge rows, verifies every internal
pair height, and checks all 231 mixed supports by exact signed-box membership.
This cross-check uses the actual coefficient budget rather than inferring
full entry from the graph rank alone.

## 4. All current hereditary profiles survive, but the selected gates fail

In fact every subset of at least seven physical labels has gcd **one**.
Here is the structural reason. A seven-subset containing just one `U` label
contains all six `V` labels; its gcd is `gcd(t,u)=1`. If it contains two
coprime `U` labels its gcd is already one. The only noncoprime pair among
`U` is a pair of positive powers of 331. A seven-subset whose selected `U`
labels are exclusively such powers contains a `V` label, and all `V` labels
and both scales are coprime to 331. Thus its gcd is one too. Larger subsets
inherit the conclusion from any seven-subset.

The completed, independently audited supplier
[lrc14_joint_shadow_empty_core_next_sep06.md](lrc14_joint_shadow_empty_core_next_sep06.md)
is the source of the current retained profiles; no reserved theorem namespace
is used as a dependency. The source checks all 4,095 subsets of sizes
7 through 12 in each row against the complete pinned profile table,
including complementary gcd words. Every signature is `(1,(1,...,1))`.

Nevertheless the balanced minimum condition fails:

    min V=185370716505 > floor(3Q/28)=60843134147.

The native larger-unit condition has `delta=gcd(g,min V)=1`; its coefficient
gate passes because `g=1`, but its other inequality fails:

    56*(max U)*(min V) > 6*Q*(max U+1).

The concurrently proved
[unit-free maximum-endpoint criterion](second_20260906_entry.md) also fails
for every one of the six partners of `L=max U`, in both physical rows.
Indeed, the three partners `1,49,109` have endpoint gcd `D=1`, while the
other three have `D=331,331^2,331^3`. All six distinguished gcds are
`delta=gcd(D,t min V)=1`; hence their coefficient gates `c=D<=Q` pass.
Nevertheless their exact signed radii satisfy

    R=Q*(u/D+L/D)-(u/D-1)*(L/D-1) <= Q*(L+109),
    6*Q*(L+109) < 56*L*min V.

Thus none satisfies the native phase comparison
`6 delta R>=56 L min V`. For `D=1`, the radius bound follows from
`u<=109`; for `D>1`, the normalized pair sum is at most `1+331^3`, which
is below `L+109`. The frozen output records every `u,D,delta,c,R` and the
source checks all six failed comparisons separately for both rows. This
preserves the distinction between an eligible coefficient gate and a
sufficient physical phase bound. The reversed component orientation has
already been blocked by its cleared coefficients exceeding `Q`.

Both direct universal Lipschitz-arc grid comparisons also fail:

    6t < 56 max U,            7g < 49 max V.

These are failed **sufficient** comparisons. They are not equivalent to
unsafe behavior, and no phase-optimality claim is attached to them.
The main row has the literal safe phase

    x=11/23,          min_(n in n(t_mix)) ||nx||=3/23>1/14.

An exact rational-grid check finds no weakly safe phase with denominator
2 through 22; this is a bounded denominator observation, not a global
optimality theorem. All odd sixteenth-grid phases have clearance `1/16`
and miss the target. The odd positive-control row is safe at `x=1/2`,
with clearance exactly `1/2`.

**The exact refuted implication** is: actual balanced equality entry plus a
unit in the larger shape and every retained hereditary gcd profile forces
the current smaller-minimum inequality (or the stronger native unit
comparison). The first failed step is upgrading connected primitive shape
and hereditary divisibility restrictions to that small numerical minimum.
The strongest survivor is the inherited conditional theorem itself. A
repaired research target must retain another phase/shape coordinate or use
a different sufficient condition. No claim is made that all balanced
unit-containing rows cannot be closed by some other argument.

The first exploratory construction used a `355`-power `U` chain. Although
all scalar caps passed, a nine-body gcd-five word with four coprime outside
labels failed the complete profile test. That rejected version is not a
surviving control. Changing the chain prime to 331 and adjoining 49 and
109 preserves actual connectivity while separating the repeated prime
from the six-shape. This is the mechanism that repairs the profile failure.

## 5. Exact reproduction

### Consequence for the incoming cross-divisor score

The concurrent [cross-divisor criterion](open_frontier_sep06_decoder.md),
commit `8e560f2142`, permits any smaller label v and an endpoint-pair gcd D.
Its second sufficient condition is `lcm(D,v)<=3Q/28` at the balanced split.
Every such choice in either saved row fails, because
`lcm(D,v)>=v>=min V>3Q/28`. Thus these actual, fully filtered entries answer
the [incoming forced-score question](open_frontier_sep06_board.md) negatively.
No new finite scan is needed for this universal-over-choices implication;
the saved minimum and the elementary lcm inequality decide it. Safety
remains witnessed by the displayed phases.

### Frozen controls

[Standalone source](../../04-computation/second_20260906_decoder.py) and
[frozen output](second_20260906_decoder.out) use standard-library exact
integer and rational arithmetic, with no producer imports. The universe
is the two explicit physical rows, all 231 mixed supports and 4,095 complete
large-subset profiles per row, and the main row's complete denominator
2..22 and odd-sixteenth controls. Every mathematical conclusion beyond
these fixed rows is the explicitly stated elementary mechanism above.

    python3 -B 04-computation/second_20260906_decoder.py
    python3 -B -O 04-computation/second_20260906_decoder.py

Both runs pass **17,617 always-active gates** with byte-identical LF output.
Source SHA256:
`31f0c395c1c41eaa944f0efb60ac0c605cec2376e49daa11a82718cb688574cd`.
Output SHA256:
`74d1bc5bbae7dcadc833121a3efd70abb87d97f08bbe1ad2c3345823c67447da`.
Inherited profile JSON SHA256:
`935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`.
The script includes a consistency comparison with the concurrently derived
stronger six-vertex real-cofactor bound; that comparison is not used to prove
entry, safety, or the failed balanced gate, and claims no independent proof
of that incoming theorem.
