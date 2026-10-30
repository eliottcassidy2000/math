# Ramified symmetry meets spectral growth and native complements

2026-10-30. **PROVED scoped deductions / CITED inputs / FINITE-EXACT
controls / OPEN full carrier.** Same-agent proof review; no independent
referee or claim of literature priority.

The ramified prime-tournament norm theorem now connects to actual
vertex attachments and Hamilton-path complements. We obtain rigid prime
children at the same spectral scale, an all-depth affine cycle inventory,
and a quantitative obstruction to short metagraph routes.

## Rigidity can be selected within a uniform spectral budget

The inherited simply ramified family supplies prime regular tournaments
T on N=p^(2k) vertices, k>=2, with

    sqrt(pN) <= ||K_T|| <= 200 sqrt(pN).

The source is [Session110's logarithm-free theorem](../session110_overnight_20261030/LOGFREE_PRIME_TOURNAMENTS.md),
with its balance constraints and discrepancy input retained. We prove
that EVERY one-vertex attachment obeys

    ||K_T||^2 <= ||K_U||^2 <= ||K_T||^2+N.

Combining the regular-parent root theorem, the older prime-extension
classification and the classical Gluck theorem gives at least two
nonisomorphic **rigid prime** children at each of these sizes N+1, with

    sqrt(pN) <= ||K_U|| <= sqrt((40000p+1)N).

Thus the arithmetic lower scale persists in a graph with no nontrivial
automorphisms. It survives through the induced source matrix; it is not
evidence that a symmetry must survive in the final graph. The cut can be
chosen for rigidity because the spectral bound holds for every cut.
The existence statement does not impose compatibility of cuts across
successive ring reductions.

There is also an orientation result: **every prime one-vertex child of
a regular tournament is nonself-converse**. Its regular-deletion root
is unique, while a converse symmetry cannot fix a vertex at even order.
The root and parity recover the orientation distinction. Clone children
are the explicit failure boundary when uniqueness is removed.

The finite certificate contains an actual625-vertex prime source and
an explicit626-vertex rigid prime child. All624 root-pair module closures
certify the parent; a unique regular card identifies the new root; two
refinement rounds distinguish all626 vertices. The source norm is
25 sqrt(5+2sqrt5), and the child norm is at most25(1+sqrt5).

## A complete cycle recurrence for the named arithmetic symmetry

Let A_k be the actual subgroup x -> g^t x+b on R_k, of order p^(3k).
For every m=1,...,k its order-p^m elements have exactly two cycle types:

| type | multiplicity | total cycles |
|---|---:|---:|
|has fixed points|(p-1)p^(3m-2)|2N/p^m-N/p^(2m)|
|fixed-point-free|(p-1)(p^2+1)p^(3m-3)|N/p^m|

[S2 in PROOFS.md](PROOFS.md) records every cycle length, not only its
total. If Phi_k(Q) sums Q^(cycle count) over A_k, then

    Phi_0(Q)=Q,
    Phi_k(Q)=Phi_(k-1)(Q^(p^2))
      +(p-1)p^(3k-2)Q^(2p^k-1)
      +(p-1)(p^2+1)p^(3k-3)Q^(p^k).

Divide by p^(3k) for complete Q-color orbit counts under A_k. Taking
Q=2^r retains r joint named binary ports. These are specified-action
counts, not full graph-isomorphism counts unless that full action has
separately been certified.

For uniform Q colors the exact mean residual subgroup order is

    1+sum_m [L_m Q^(-N(1-p^-m)^2)+U_m Q^(-N(1-p^-m))].

Its leading correction is p(p-1)Q^(-N(1-1/p)^2), with an explicit
exponentially smaller remainder. Even one binary port removes A_k with
probability tending to1. Every fixed positive residual-order moment
tends to1, uniformly for measured IID palettes with maximum color
probability bounded below1. This is a different mechanism from the
previous wreath towers, where local symmetries remain in small blocks.

The full-group guard is essential. At the25-vertex base, the named
arithmetic group has order125, but the full tournament group has
order15625. Session110's higher-depth prime scalar hostile remains
in force. The probability statement concerns A_k; the separate Gluck
argument and explicit626 certificate concern full graph rigidity.

## An operation law for the actual metagraph

For a legal Hamilton path P, retain its ordered vertices and use the
actual relation T^P=T^op xor E(P). Its skew matrix satisfies

    K_(T^P)=-K_T+2J_P,
    | ||K_(T^P)||-||K_T|| | <=4 cos(pi/(N+1))<4.

Re-marking by vertex permutation costs zero. More generally, a word
with r insertions and c legal complements, starting at size N and
norm M, obeys

    max(0,M-4c) <= M_final
      <= sqrt(M^2+rN+binom(r,2))+4c.

This is a genuine composition bound for those operations. It neither
reconstructs their joint-port responses nor preserves primeness.
At fixed size it gives a linear asymptotic lower bound on complement
distance from the ramified quasirandom sources to transitive tournaments.
For the certified625-vertex source, **at least81 actual complements**
are required along any such route. Existence of a route is not assumed.

## Evidence and continuation

The exact checker covers22,860 affine permutations in four complete
finite universes,65 larger-depth probes, joint palette recurrences,
1,098 insertion cases and7,886 legal path-complement cases. It verifies
the native626 example, three further literal626-vertex complements,
and replays the81/27 hostile with equal67 labelled cyclic triples but
1620/1728 regular binary masks. Normal and optimized runs are compared;
see [AUDIT.md](AUDIT.md) for counts, exactness and source pins.

From the repository root:

    python -X utf8 -B 05-knowledge/results/ramified_symmetry_spectral_20261030/check.py
    python -X utf8 -O -B 05-knowledge/results/ramified_symmetry_spectral_20261030/check.py --output /path/to/separate-receipt.json
    python -X utf8 -B 05-knowledge/results/ramified_symmetry_spectral_20261030/seal.py

The next targets are compatible distinguishing cuts under the real ring
reductions; full graph groups beyond the named affine subgroup; and
joint response carriers that retain the selected paths while respecting
the spectral budget. General complement/re-marking closure and new full
A000568 terms remain OPEN.
