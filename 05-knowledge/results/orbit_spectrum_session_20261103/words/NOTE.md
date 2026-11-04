# Strong source tests and exact two apex continuation limits

**PROVED scoped append rigidity / FINITE-EXACT bases with CITED all-order
cutoffs. The full Hamiltonian spectrum remains OPEN.**

The inherited two-apex language is an actual tournament: named apices x→y,
and a transitive chain with letters0,1,2,3. Bit0 says x beats the chain
vertex; bit1 says y beats it. A suffix adds vertices after the entire old
chain. This is a restrictive physical attachment rule, even when the two
apices can be chosen afresh before growth.

The [Session134 decoder](../../session134_constructive_spectrum_frontier_20261101/two_apex/THEORY.md)
already gives rank17, exact actual-word recovery, and an all-order scalar
membership algorithm. Its inherited finite coverage through99999 is not
extended here. The new points are a strong-source spanning certificate,
an all-length classification of universal affine suffixes, and complete
target menus under every allowed marking of small actual sources.

## 1 Strong sources already test every affine suffix identity

Let P(w) be the seventeen continuation counts with suffixes

    empty,0,1,2,00,10,20,30,01,11,21,31,000,100,200,010,210.

Every H(ws) is a linear functional of P(w), by the inherited closed
observer. The following seventeen words all give strong tournaments:

    2,12,20,21,22,30,32,102,112,120,121,122,130,132,200,201,210.

**Proposition W1.** Their seventeen profile rows have determinant150994944.
The constant function1 also belongs to the span of these seventeen query
functions. Consequently an affine identity

    H(ws)=a H(w)+b

for a fixed suffix s and rational constants a,b holds for every strong
word w iff it holds for every word w.

**Certificate and mechanism.** The checker independently implements the
gap bijection P,a,Q,b,R for the two orders of the apices. It merges the
empty and good-last states of P, yielding24 nonnegative coordinates.
This is an implementation variant of the inherited36-state transfer,
not a new minimal dimension. The reachable physical span has dimension20.
The seventeen query rows close under all four letters on that span;68
exact rational identities certify closure. Summing the first-order gap
states with Q=R=empty counts the unique assignment of all chain vertices
to P, so its value is1. This row is invariant under every letter and is
an exact linear combination of the seventeen queries; the coefficients
are printed in `check.json`. Exact elimination gives the displayed nonzero
determinant. Thus a linear combination of these query rows vanishing on
the seventeen strong rows vanishes on all actual words. This proves W1.

No numerical rank tolerance or finite-depth extrapolation is used. The
physical span is closed by a queue until linear closure; each transition
is checked on every basis vector. The certificate is independently
reproduced under normal Python and optimization.

**Corollary W2, complete fixed-suffix classification.** Such a universal
affine identity on strong words exists exactly when s consists only of
the letter3. In that case a=1 and b=0.

Proof. By W1 the identity holds also on w=3^m, whose count is1. If s
contains a letter c≠3, its vertex beats at least one apex r. The induced
tournament on r, the m old chain vertices, and this new vertex has the
cyclic block pattern r→TT_m→new→r. Its count is2^m+1. Hamiltonian-path
counts cannot decrease on induced extension, so H(3^m s) is unbounded,
contradicting the constant a+b. Conversely each3 is a universal sink,
and therefore preserves H. The strong counts3 and5 force a=1,b=0.

This excludes every fixed suffix of arbitrary length, not just a bounded
macro search. It does not exclude a source-dependent choice of suffix,
reselection of the apices, interior chain insertion, or native editing.
The [moving-mark construction](../coverage/NOTE.md) deliberately leaves
the fixed-source condition.

## 2 Every suffix is finitely testable from a strong source

Starting with a strong word w, append an arbitrary suffix s and remove
its trailing3s. Those removed vertices are universal sinks and preserve H.
If the remaining suffix is nonempty, its last vertex beats an apex and
is reached by every old chain vertex. The old source is strong, so every
intervening appended chain vertex joins its strong component. Therefore
the remaining entire tournament is strong.

Write f(n) for the cited exact Moon–Busch strong minimum:

    f(3k)=3*5^(k−1), f(3k+1)=5^k, f(3k+2)=9*5^(k−1).

For target h, put nmax(h)=max{n:f(n)≤h}. A successful suffix, with final
sinks removed, has length at most nmax(h)−|w|−2. Enumerating this finite
set therefore decides **all-length suffix reachability** from the given
actual source and fixed apices. Monotone-H pruning is valid; rejecting a
source does not reject the target globally.

This is a specialization of the inherited bounded-H decision principle,
with the smaller strong cutoff justified by this suffix language. The
primary source is [Busch, Theorem1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf);
the corrected formula is also recorded in
[Session133 core reduction](../../session133_fractran_spectrum_refinement_20261101/core_reduction/THEORY.md).

## 3 Adaptive marking repairs some sources and fails at others

Two presentations of the same four-vertex strong tournament illustrate
why a scalar H or an unmarked class is insufficient:

| Fixed marked word | H | Target3H+2 | Target3H+4 |
|---|---:|---|---|
|21|5|suffix11 gives17|no suffix gives19|
|32|5|no suffix gives17|suffix32 gives19|

Both exclusions are all-length statements: the strong cutoff allows at
most two relevant births. Other markings of this same class do supply
both branches. Reselecting a port is a real change of the construction
state, even when the underlying tournament and H are unchanged.

The checker then groups **all** strong two-apex words at chain lengths1–4
by ordinary isomorphism, and considers every word in each class. Every
ordered apex pair whose deletion is transitive appears in this list.
Converse classes remain separate.

| Source order | Strong classes admitting two apices | Classes supplying both ternary targets after some marking |
|---:|---:|---:|
|3|1|1|
|4|1|1|
|5|6|3|
|6|34|5|

There are35 strong classes at order6; one has no two-apex presentation and
is outside this table. The table is a complete source census in the stated
language, not a full tournament-growth classification.

**Corollary W3.** No strong tournament with H=15 can reach47 by choosing
two apices with transitive deletion and then appending an arbitrary
transitive-chain suffix while keeping those apices fixed.

Proof. Since f(7)=25, every strong H15 source has order at most6. Independent
Hamilton-frame enumeration at orders5 and6, followed by literal subset-DP
counts and isomorphism reduction, gives exactly four classes. All four
occur in the word census. Every permitted apex pair misses47 in its
complete suffix menu; nmax(47)=8 makes each menu finite. This proves the
all-order claim. Both order-five H15 classes also miss49 under every
such marking; one of the order-six classes reaches49. QED.

Thus an arbitrary sequence of **strong H15-neutral resets**, followed by
one such suffix phase, still cannot supply47. This includes the earlier
five-vertex native reset, which repaired a continuation divisor but did
not promise this attachment rule. Native operations interleaved with
growth, weak intermediate sources, moving marks during growth, and larger
feedback sets are not covered by W3.

## 4 A literal interior-insertion repair and its boundary

The H9 class represented by112 cannot reach31 under any of its two-apex
markings followed only by suffix growth. However

    112 → 10102,       H9 → H31

is an actual two-birth extension retaining both apices and all old arcs:
the old chain occupies positions0,2,4 of the final chain. Insert a0 after
each of the first two old chain vertices. A new0 loses to its chain
prefix and beats both apices and its chain suffix. The old chain cut is
therefore additional attachment data; the suffix observer alone does not
authorize this operation.

The tempting repeated-pattern extension does not give a uniform+4 law.
For m≥3, direct gap transfer gives

    H(1010 1^(m−3) 2)=11*2^(m−2)+9.

It equals3(2^m+1)+4 only at m=3. The scalar success exposes a useful
operation but does not close the induction. This family identity is
checked by a closed constant-letter matrix recurrence in the checker.

Positive controls1202,1210,11122,10212 have counts29,31,47,49 by independent
adjacency DP. Thus every displayed missing arrow has an attained target.
The missing ingredient is compatible attachment and continuation, not a
new forbidden count.

## 5 Reproduction

Run `python3 -B check.py` and `python3 -O -B check.py` in this directory.
The receipts agree byte for byte. They include341 literal word checks,
the strong-language criterion,68 row-closure identities, the exact strong
minor, complete small-source suffix menus, and1088 independent H15
Hamilton frames. The map from word to class deliberately forgets the
apex names; retaining the complete word menu restores the information
needed for the stated existential attachment test.
