# Fixed nontransitive prefixes and the odd part of the slope

2026-11-03. **PROVED reduction; FINITE-EXACT complete small core inventory.**
This rules out a specified prefix escape for the uniform 44 target. It does
not rule out consuming prefix vertices into a changing suffix or a
nonuniform attachment.

## Inheritance and the connection

The attachment identity and ordered-sum law are inherited from
`session182_physical_returns_20261104/induction/PROOF.md`. The recovered older
mechanism is [THM-4744, interval chains and bounded-count inventories](../../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md).
Its classical order bound is checked against Busch's primary paper,
[Theorem 1 and Corollary 1, pages 3–4](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/).
The bound gives minimum 45 at strong order eight and at least 125 at strong
orders ten and above. These are **CITED** inputs, not new results here.

The map considered here replaces a variable transitive prefix by a fixed
tournament Q in an ordered sum with a shorter variable transitive prefix.
It preserves the exact attachment count formula, but carries H(Q) and the
ordered two-path-cover count P_Q(2) separately. Equal H(Q) does not permit
discarding the odd part of P_Q(2). Names and native guards are additional
obligations; the arithmetic exclusion below does not need to forget them.

## The necessary coefficient equations

Consider a fixed nonempty Q, a fixed suffix B with actual port set W meeting
its terminal SCC, and a fixed integer k. The output prefix is

    A_t = TT_(t-k) -> Q,

at every sufficiently large integer t. Placing Q before that transitive
block gives the same two scalar equations, although not the same graph.
Write a=H(Q), p=P_Q(2), f=F_W(B), and c=H(B)-J_W(B). The output count is

    a*c + 2^(t-k)*p*f.

The desired second-level families have counts

    63*2^t - 35       (word 24),
    63*2^t - 29       (word 44).

Equality at two distinct legal clocks forces respectively

    a*c = -35 or -29,       p*f = 63*2^k.                 (P1)

In particular a divides 35 or 29, and the odd part of p divides 63. Since
p and f are positive integers, negative k cannot solve the latter equation.
These are necessary conditions; they do not produce a physical rewrite.

## A complete prefix-core inventory

Hamilton counts multiply across the ordered strongly connected components
of Q, and P_Q(2) also multiplies. A singleton contributes (H,P2)=(1,2).
For H(Q)=29, exactly one component has H=29; every other component is a
singleton. For H(Q)=35, the same holds with 35, because the other possible
nontrivial factorization would require H=7. The all-order exclusion of 7
is inherited from
[THM-1370, h-spectrum-omits-7-21-all-n](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
The case H(Q)=5 likewise has one component of count 5.

Busch's minimum at order eight is 45, so every such nontrivial component has
at most seven vertices. An exact isomorphism-class generation through seven
and independent permutation counts on all selected components give:

| H of the nontrivial strong component | Possible odd parts of P2 |
|---:|---|
| 5 | 1 |
| 29 | 5, 145 |
| 35 | 169 |

Arbitrarily many singleton components multiply P2 by powers of two, so the
table controls Q of **every order** with these Hamilton counts. The table is
a complete finite-exact result with an explicit all-order reduction, not an
extrapolation from small examples.

It follows from (P1) that the **44 family admits no nontransitive fixed Q**
in this prefix model: a=29 is its only nontransitive option, and neither
5 nor 145 divides 63. For 24, a=35 is excluded by 169, and a=7 is already
impossible. The only nontransitive possibility left by this test is a=5.

Such Q consists of one strong four-vertex component and singleton components.
Consequently p=32*2^z for some z>=0, and (P1) requires

    c=-7,       f=63*2^(k-5-z),       k>=5+z.            (P2)

Thus this prefix route to 24 must consume at least five transitive clock
positions. Equation (P2) is a residual target, not a constructed repair.
If a=1, Q is transitive and merely changes the transitive clock; that is the
separate clock-consumption lane.

## Reproduction and exact scope

`probe.cpp` generates every tournament class inductively by all one-vertex
birth cuts and exact degree-respecting canonicalization. The counts of classes
at orders one through seven are 1,1,2,4,12,56,456. This is a class-coverage
control, not the only completeness argument: deleting the newly born vertex
shows inductively that every class is visited. Canonicalization explores every
permutation inside each degree class, which is an exact isomorphism key.

The C++ counter uses forward subset paths. `check.py` independently enumerates
all path permutations on every subset of each of the nine selected strong
representatives, then reconstructs P2 from ordered complementary subsets.
The receipt retains every actual adjacency matrix and its two counts.

Run from the repository root:

    python 05-knowledge/results/c4_port_returns_20261103/prefix/check.py
    python -O 05-knowledge/results/c4_port_returns_20261103/prefix/check.py

No assertion about arbitrary prefix families follows. In particular Q may
depend on t, prefix vertices may be promoted into the suffix, or incidences
between the prefix and suffix may cease to be uniform. Those operations have
different coefficient equations and remain available.
