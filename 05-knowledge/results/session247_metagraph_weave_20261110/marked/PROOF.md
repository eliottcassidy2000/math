# Marked interfaces activate invisible cubes, but carry a fixed-divisor charge

**Status: new proofs below, conditional only on the explicitly inherited
Session170 seed profile. Small literal checks concern the new interface,
not an audit of that earlier session. Chan–Pak remains OPEN.**

## 1. The historical weave and the precise change of object

Session170's `metagraph_topology/NOTE.md`, section 4, gives strong six-vertex
tournaments A and B with lexicographic pair masks 18 and 22 (bit 1 means
`i -> j` for `i<j`). They differ in one arc, have H=15, and share

    pc(c) = (15,93,138,73,15,1), c=1,...,6.

The topology compiler takes ordinal sums of these modules. Its whole cube
is invisible to H and to the complete unmarked path-cover profile. The
older metagraph-stalk reflection asks what extra incidence makes an
invisible change observable. THM-3121's walk-content kernel says exactly
why an unmarked substitution cannot help: it factors through that profile.

We make the incidence asymmetric with one actual vertex. This is a physical
operation on a literal arc-flip cube, not a complement-line move and not a
claim about cube interiors. Separately, the coefficient-tiling and
prime-irreducible-address reflections suggest retaining the fixed divisor
of the integer-valued slot polynomial. Here that suggestion becomes an
exact theorem about how many marked vertices are needed to remove a prime
factor from a modular construction.

## 2. General fixed-divisor charge and bounded-interface theorem

For a nonempty tournament M define

    f_M(c)=c! pc_M(c),
    delta(M)=gcd(f_M(1),...,f_M(|M|)).

Equivalently, if F_M(s) counts allocations into s labelled path slots,
allowing empty slots, then

    F_M(s)=sum_c binom(s,c) f_M(c),
    delta(M)=gcd{F_M(s):s=0,1,2,...}.

The equality follows in both directions from Newton finite differences:
values are integer combinations of the f_M(c), and f_M(c)=Delta^c F_M(0)
is an integer combination of values. Thus polynomial coefficient content
in the monomial basis is not the right charge.

Let M_1,...,M_k be disjoint nonempty modules, with every intermodule pair
oriented uniformly according to any quotient tournament Q. Add a set X of
q further vertices. Arcs inside X and between X and individual module
vertices are arbitrary. Those marked arcs may destroy the module property
in the final tournament T.

**Theorem 1 (interface valuation survival).** For every prime p, put
`e_i=v_p(delta(M_i))`, sorted `e_(1)<=...<=e_(k)`. Then

    v_p(H(T)) >= sum_(i=1)^(max(0,k-2q)) e_(i).             (1)

In particular, if p divides every delta(M_i), then

    p^max(0,k-2q) divides H(T).                            (2)

**Proof.** Partition Hamilton paths by their skeleton: replace every
maximal consecutive run within a module by its module label, and retain
each individual vertex of X. Consecutive equal module labels are forbidden.
For a fixed skeleton, call a module touched if one of its occurrences is
adjacent to an X vertex. There are at most 2q touched modules, because a
Hamilton path has at most two incidences at each X vertex.

For each untouched M_i occurring c_i times, its internal runs form an
unordered c_i-path cover. Bijection of its distinct path components with
its ordered skeleton occurrences gives exactly f_(M_i)(c_i) choices. These
choices are independent: all neighboring module-to-module arcs are uniform,
and no endpoint of these runs meets X. For touched modules, all endpoint
constraints can be counted in a nonnegative integer factor; it need not
be divisible by delta(M_i). Hence the skeleton contribution is divisible
by the product of delta(M_i) over untouched modules. Its p-valuation is at
least the sum of the smallest k-2q of the e_i. Every skeleton contribution
has this same lower bound, so their sum does too. This is an ordinary
integer divisibility proof, with no cancellation assumption. QED.

For an ordinal quotient, each segment between X vertices moves forward
through the modules. Thus c_i<=q+1 and the theorem strengthens on replacing
`delta(M_i)` by

    delta_(q+1)(M_i)=gcd{c! pc_(M_i)(c):1<=c<=min(q+1,|M_i|)}.

A more detailed skeleton-dependent boundary-incidence budget can strengthen
(1); no sharpness for every q,k is claimed.

**Prime-output consequence.** If p divides all module charges and H(T) is a
prime other than p, necessarily `k<=2q`. More generally (1) bounds how many
p-charged modules a low-p-adic-valuation target can accommodate. This applies
even when Q is strongly connected; it is not merely SCC multiplication.

For the inherited A/B pair,

    f=(15,186,828,1752,1800,720), delta=3,
    F(s)=s^6+8s^4+6s^2.

Although F is a primitive monomial polynomial, F(s)=0 mod 3 for every integer
s: mod 3 it is s^4(s^2-1). Also gcd(F(1),F(2))=gcd(15,216)=3.
Consequently every q-marked realization of k such modules has
`3^max(0,k-2q) | H`. The full profile, not only H=15, carries the surviving
charge. To build arbitrary prime counts one must let the interface grow,
use modules with charge 1, or abandon this module decomposition.

## 3. Exact one-root ordinal response

For a tournament M and a new root x, let O be the vertices beaten by x.
For S subset V(M), let

- h(S)=H(M[S]), with h(empty)=1;
- e(S)=number of Hamilton paths of M[S] ending outside O;
- s(S)=number of Hamilton paths of M[S] starting inside O.

Set e(empty)=s(empty)=0, and define

    h=h(V), E=e(V), S=s(V),
    D=sum_(empty<L<V) e(L)s(V\L),
    A=sum_(empty<L<=V) e(L)h(V\L),
    B=sum_(empty<R<=V) s(R)h(V\R),
    Q=sum_(L subset V) h(L)h(V\L).

For the ordinal core M_1=>...=>M_k with this root, write
`P(J)=product_(j in J) h_j`. Empty products are 1.

**Theorem 2 (one-root boundary formula).** The exact Hamilton-path count is

    S_1 P(2..k) + E_k P(1..k-1)
    + sum_i D_i P([k]\{i})
    + sum_(i<k) E_i S_(i+1) P([k]\{i,i+1})
    + sum_(r<l) B_r A_l P([k]\[r,l]) product_(r<j<l) Q_j.  (3)

**Proof.** Delete x from a path; the remaining prefix and suffix each move
forward through the ordered blocks. If either is empty, the first two
terms apply. Otherwise let l be the last prefix block and r the first
suffix block. If r=l, both pieces occupy one block and give its D term.
If r>l, necessarily r=l+1, giving the adjacent E*S term. If r<l, the suffix
starts in block r and the prefix ends in block l. Their boundary partitions
contribute B_r and A_l; each intermediate block may be divided arbitrarily
between the two ordered paths, giving Q_j. Earlier/later blocks lie wholly
on one side and give h_j. This is a disjoint exhaustive classification. QED.

For a binary choice of two modules with the same h and Q, all higher mixed
finite differences of (3) vanish beyond order two. Its pair coefficient,
for r<l, is

    (Delta B_r)(Delta A_l) P([k]\[r,l]) product_(r<j<l)Q_j,

plus, when l=r+1, `(Delta E_r)(Delta S_l) P([k]\{r,l})`.
This specifies exactly which pointed responses an invisible switch must
change. It is stronger than saying that a root might make it observable.

## 4. An affine physical cube of every dimension

In every A/B module let x beat exactly vertices {0,2}, mask 5. The marked
states are

    A: (h,Q,E,S,D,A,B)=(15,216,9,2,30,119,38),
    B: (h,Q,E,S,D,A,B)=(15,216,8,2,33,119,41).

Thus Delta S=Delta A=0, while Delta E=-1, Delta D=3, Delta B=3.
Every pair coefficient in (3) vanishes. Consequently

    H(T_b)=C_k + sum_(i=1)^k g_(k,i)b_i, b_i in {0,1},     (4)

where `C_k` is (3) evaluated on all A modules. Explicitly,

    g_(k,k)=2*15^(k-1),
    g_(k,i)=43*15^(k-2)
             +357 sum_(d=1)^(k-i) 15^(k-d-1)216^(d-1), i<k.  (5)

For k=1 the second line is unused. Every coefficient is positive and even.
All binary vertices are literal tournaments on 6k+1 vertices. Adjacent bits
reverse the same actual designated internal arc. The tournaments are strong:
each block is strong and x has both in- and out-neighbors in every block,
so every block communicates with x in both directions.

This is an H-affine marked realization of Session170's formerly constant
cube, proved from its boundary profile. No profile equality after marking is
claimed. Distinct bit words need not be distinguished after forgetting all
block and root labels; the theorem is on physical labelled tournaments.

Examples:

    k=1: H=41,43.
    k=2: C_2=5605, g=(400,30),
         H in {5605,5635,6005,6035}.

The large gaps are real. Theorem 1 implies a factor 3^(k-2) for k>=3, so
this particular family cannot be the all-target answer. Positive affine
selection and prime-target freedom are distinct requirements. The useful
new search prescription is to seek equal-profile flip pairs with charge 1
and marked identities Delta A=Delta S=0, then control their actual g_i and
C_k. Those conditions are explicit, finite response identities; neither
gcd one nor their existence alone asserts completeness.

## 5. Connection ledger

- Source: Session170 constant-profile flip cubes. Map: append one marked
  root. Preserved: literal designated arc flips and full physical vertices.
  Lost: unmarked module invariance. Restored sidecar: E,S,D,A,B and Q.
- Source: THM-3121 factorial path-cover walk kernel. Map: cut paths into
  module runs while retaining exceptional-vertex incidences. Preserved:
  exact multiplicity. New consequence: bounded-interface valuation survival.
- Source: old coefficient-tiling/fixed-divisor reflections. Map: slot count
  F(s) to its Newton coefficients c!pc(c). Preserved: fixed divisor, not
  monomial content. New consequence: a primitive slot polynomial can carry
  an unavoidable connector-budget charge.
- Source: Session246 affine exchange selection. Map: replace disjoint cut
  swaps by invisible internal module flips with an asymmetric observer.
  Preserved: additive H-cube. Lost: small weights and baseline control.
  New positive realization (4); exact arithmetic obstruction (1).

Sources:
`05-knowledge/results/session170_holonomy_homotopy_20261103/metagraph_topology/NOTE.md`;
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`;
`07-reflections/the-metagraph-needs-a-stalk-canonical-addresses-tiling-fibres-and-the-4d-lrc-object-codex-S4.md`;
`07-reflections/coefficient-tiling-and-prime-irreducible-addresses.md`;
`07-reflections/irreducibility-prime-values-and-tournament-recombination.md`;
`01-canon/theorems/THM-5122-shared-core-exchange-intervals-and-negative-inertia-capacity.md`.
