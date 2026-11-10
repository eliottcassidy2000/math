# Joint deletion cores and the sharp candidate threshold 95

**PROVED reductions, relative to the cited Moon–Busch minimum; FINITE-EXACT
small marked-class census and inherited audited larger census. The proposed
all-height branching rule is OPEN.** Local names J1–J4 belong only to this file.

Here `H(T)` counts ordered directed Hamiltonian vertex words, with
`H(empty)=1`. A marked class is a tournament with one distinguished vertex;
isomorphisms preserve that vertex. Converse classes are not merged.

## Inheritance and the procedural choice

The closest mechanism is ordered-SCC multiplication. The exact strong floor,
logarithmic unmarked core, and rational **single-count** class series are
inherited from [the spectrum-basins packet](../../spectrum_basins_20261020e/spectrum.md),
[THM-1370, all-order omissions](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
and [THM-4744, interval chains and bounded strong inventories](../../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md).
These are not new here.

The corrected near miss is the finite-only statement that no realiser of13
of order at most10 extends by one vertex to27, in
[THM-5127's source, section8](../../tilted_ensemble_local_limit_20261110/README.md).
The useful neglected coordinate is **which SCC contains the new vertex**.
The new refinement applies the inherited decomposition simultaneously to
`H(T)` and `H(T-v)`, then to their marked isomorphism classes.

The anchor is ternary target coverage; the niche is exact marked-class
enumeration; the wildcard is replacing a large all-realisers search by a
small strong-child deletion test. The live board is: scalar target / marked
SCC / common factor / strong floor / ternary branch / retained attachments.
A successful scalar reduction must still state which attachment information
it forgets. No connection to LRC(14) is asserted.

## J1. The exact common-factor reduction

Let `v` be a vertex of a tournament `T`, let `C` be its strongly connected
component, and set

    N=H(T),  h=H(T-v),  r=product_{D != C, D an SCC of T} H(D).

Then

    N=r H(C),                  h=r H(C-v).                 (J1)

In particular `r | gcd(h,N)`. If `N>h` and `gcd(h,N)=1`, then all SCCs of
`T` outside `C` are singletons and `C` is nontrivial and strong. Therefore

    T = TT_a -> C -> TT_b,
    H(C)=N,   H(C-v)=h.                                   (J2)

Conversely every marked strong core with these two counts, with arbitrary
transitive prefix/suffix, gives such an extension pair. Every nontrivial SCC
of the source `T-v` is inside `C-v`: the single new vertex has to join all
of them into its strong component.

**Proof.** Tournament SCCs are totally ordered, and every Hamiltonian path
concatenates one spanning path from each SCC. Deleting `v` changes only the
factor of its own SCC. Although `C-v` can split into several SCCs, their
product is precisely `H(C-v)`. This proves(J1). A nontrivial strong tournament
has at least three Hamiltonian paths; thus outside product1 forces every
outside SCC to be a singleton. If `C` were a singleton, deleting it would
leave `N=h`; hence it is nontrivial in the strict case. The converse is the
same concatenation argument. QED.

This is an **existence reduction**, not a legal operation that deletes old
required names. It retains the marked core and the two scalar counts, but
forgets the names and positions of outside singleton vertices. To restore
the original graph retain the ordered named prefix and suffix. If future
new vertices may connect nonuniformly to those names, their incidences are
also required. Even the two counts together are not a closed induction
state; neither the full joint17-query observer nor its F49/F343 charts are
identified with this pair.

## J2. Finite cutoff and the marked rational series

The classical exact strong minimum is

    f(3q)=3*5^(q-1), f(3q+1)=5^q, f(3q+2)=9*5^(q-1).      (J3)

Thus a coprime pair with `N>h` exists iff some strong tournament `C` with
`f(|C|)<=N` has a vertex deletion of count `h`. In particular its core has
`|C| <= 1+3 log_5 N`. This follows from the **cited** strong minimum, not
from the finite tests below. The [primary Busch paper, Theorem1 and
Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
was rechecked on2026-11-10. The unmarked logarithmic reduction was already
proved in the inherited spectrum-basins packet.

For ordinary tournament classes let the inherited series be

    F_a(z)=sum_[A]:H(A)=a z^|A|,       F_1(z)=1/(1-z).

Let `Q_(b,a)(z)` count marked **strong** classes `(C,v)` with
`H(C)=a,H(C-v)=b`. Include the marked singleton as `Q_(1,1)=z`; all other
coefficients of `Q_(1,1)` vanish. Let `G_(h,N)(z)` count all marked classes
with those two counts. Then the exact all-order refinement is

    G_(h,N)(z) = sum_{r | gcd(h,N)} Q_(h/r,N/r)(z)
                           sum_{ab=r} F_a(z) F_b(z).       (J4)

Every `Q` is a finite polynomial by(J3), and every `F` is an inherited
rational series, so `G` is rational. For a coprime pair,

    G_(h,N)(z) = Q_(h,N)(z)/(1-z)^2.                       (J5)

In particular the number of marked classes at order `n` is eventually a
linear polynomial, with slope `Q_(h,N)(1)`, or is identically zero when
`Q=0`. The fixed count `N` may have a much higher-degree class polynomial;
the coprime deletion constraint removes every nontrivial outside factor.

**Proof.** The SCC containing the mark is unique. Its predecessor SCCs form
one arbitrary tournament `A`, and its successor SCCs one arbitrary tournament
`B`, each possibly empty. Their ordered triple `(A,(C,v),B)` determines the
marked tournament, and every marked tournament has exactly one such triple.
Automorphisms preserve SCC order and the marked component. The two counts
are `(ab H(C-v),ab H(C))`, proving(J4) without labelled multiplicities or
overcounting. Formula(J5) follows at `gcd(h,N)=1`. QED.

For example the complete core polynomials for targets below45 include

| deletion h | child N | Q_(h,N)(z) |
|---:|---:|---|
|3|11|0|
|3|13|2z^5|
|5|17|2z^6|
|5|19|2z^6|
|9|29|10z^6+2z^7|
|9|31|4z^6+5z^7|
|11|35|2z^7|
|11|37|6z^6+2z^7|
|13|27|0|
|13|41|2z^7|
|13|43|4z^6+2z^7|

These polynomials are **FINITE-EXACT and complete at all orders** because
`f(8)=45`. For instance there are exactly `2(n-4)` marked classes with
pair `(3,13)` at every `n>=5`. The coefficient counts root automorphism
orbits, not raw choices of a vertex or labelled tournaments.

## J3. Nine ternary edge failures are permanent

For odd `h`, both `gcd(h,3h+2)` and `gcd(h,3h+4)` are1. Consequently J1
applies to every proposed one-vertex branch of the desired induction.

The separately audited [THM-5127 census](../../tilted_ensemble_local_limit_20261110/README.md)
examines all base classes of orders at most10 and all one-vertex cuts.
Its frozen `out/pA_le10.bin.gz`, with `h<1024`, has exactly the following
missing branch pairs among realised source counts:

| offset | absent pairs (h,3h+offset) |
|---|---|
|+2|(3,11), (15,47), (35,107), (93,281)|
|+4|(1,7), (17,55), (35,109), (41,127), (87,265)|

**All nine absences hold at every order.** Indeed their largest child
count is281 and `f(12)=375`. If any missing pair had a witness at any order,
J1 would strip it to a strong child of order at most11. Its deletion would
be a base of order at most10, already in the complete inherited census:
a contradiction. This is a new all-order consequence of the **inherited**
census and J1, not a new larger census.

The separate missing edge `13 -> 27` is also impossible at every order:
its strong core has order at most7 and the independently rebuilt small
census above has `Q_(13,27)=0`. This upgrades the explicitly finite-only
statement in the inherited note. It does not say13 or27 is absent from
the global count spectrum; both are attained.

The table is a complete failure list only for the declared source range
`h<1024`. It proves no absence beyond that range. Every odd
`95<=h<=1023` has both branches in the inherited census:465 source values,
930 positive branch pairs. Across all realised `h<1024`, there are1,011
positive and9 missing branch pairs. The row `(93,281)` shows why a proposed
universal one-vertex rule cannot start at any odd threshold below95.

The earlier adaptive native constructions allow seed changes, several
births, and old-arc rewrites. J3 does not refute those broader languages.
It rules out these particular edges when the old source remains an
induced tournament and exactly one new vertex is added.

## J4. A smaller, precise sufficient induction target

Consider the following **OPEN** two-target version of the window route.
It is stronger than scalar spectrum completeness and does not demand an
entire interval at each source:

> For every realised odd `h>=95`, there is a one-vertex extension pair
> `(h,3h+2)` and one `(h,3h+4)`. Their source tournaments may be different.

This rule, together with the already certified odd values through287 and
ordered-join multiplication by3, would prove that the only missing positive
odd counts are7 and21.

**Proof of the implication.** The required base follows from
[THM-1370](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
which already gives all allowed odd counts through609. For an odd
`N>=289`, choose `delta` in `{0,2,4}` such that `h=(N-delta)/3` is an odd
integer. It satisfies `95<=h<N`. Strong induction supplies a realiser of
`h`. If `delta=0`, join it to a directed triangle. Otherwise apply the
open rule to a possibly different realiser of the same count. This gives
`N`. The inherited all-order exclusions handle7 and21. QED.

By J2 the open rule can equivalently demand two strong children, each with
a deletion of count `h`, at the explicit logarithmic order cutoffs. Thus
individual membership tests terminate without an arbitrary-order search.
The unresolved work is an **all-height existence theorem**, not deciding
one finite instance. For a constructive proof retaining actual old names,
ports and next-step responses, one also needs a compatible source-selection
and transport rule; this scalar implication alone does not provide it.

## Verification and dependency boundary

Run from repository root:

    python3 -B 05-knowledge/results/decomposition16_odd_spectrum_20261110/induction/check.py
    python3 -O -B 05-knowledge/results/decomposition16_odd_spectrum_20261110/induction/check.py

`check.json` and `check.optimized.json` are byte-identical. The script:

- rebuilds every class through order6 from1,100 Hamiltonian-path frame
  presentations, canonically comparing all actual Hamiltonian words;
- verifies the inherited456 order-seven representatives are distinct and
  independently regenerates their complete universe by all3,584 extensions
  of the56 six-vertex classes;
- compares direct path-word enumeration with subset DP on all532 nonempty
  classes through7, checks3,613 marked SCC factorizations, and checks(J4)
  through order7 on all3,403 marked classes and614 distinct count pairs;
- reads and hashes the separately audited inherited larger table to extract
  the nine missing pairs and finite95..1023 positive band. **It does not rerun
  the order-ten census.** The all-order promotions additionally use J1 and
  the cited strong floor.

A separate [root-agent audit](../core_audit.py), with
[normal](../core_audit.json) and [optimized](../core_audit.optimized.json)
receipts, also checks all ten negative pairs directly against the older,
independently audited complete strong inventory through count375. It
recounts263 selected strong classes and2,171 vertex deletions by literal
path-word enumeration, with five positive controls. This is independent
of the inherited order-ten extension table; completeness of the older
strong inventory remains a stated inherited audited input.

No numerical value-coverage record is claimed, and no canonical theorem ID
is reserved or promoted here. The proofs are new scoped refinements of
inherited SCC machinery; no external priority claim is made.
