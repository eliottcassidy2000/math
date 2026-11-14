# Exact modular tests on a sparse contact grid

The one-apex gcd in [WINDING_LOCK](WINDING_LOCK.md) extends to
several independent apices. This offers a rigorous way to test
for an arithmetic obstruction without mistaking it for an image
or surjectivity theorem.

Fix an n-vertex tournament B, add k distinguished vertices
with prescribed arcs among them, and assign Boolean contacts
u_(a,i) for apex a against base vertex i. Then the exact
Hamilton path count is an integer multilinear polynomial
in these kn variables, of degree at most two in each apex's
row. Consequently:

* To recover the entire count function, it is sufficient to
  test every mask with at most two outcontacts **per apex**.
  The number of masks needed is
  ((1+n+\binom{n}{2})^k), rather than (2^{kn}).
* The criterion is valid modulo **every** modulus, with no
  assumption of prime modulus, and can simultaneously test
  exact count equality of two different fixed bases.
* In particular, a congruence that holds for all those sparse
  masks holds for every mask in the full chart. Its converse
  is obvious. This detects congruence locks, not gaps inside
  the allowed classes.

Proof. The odd-cycle packing formula for H sums over
vertex-disjoint directed odd cycles with weight 2 per cycle.
Write each required variable-arc orientation as u or 1−u
against its bit. In any pack, each apex is on at most one
cycle, whose two incident arcs use at most two bits from
that apex's row. Multiplying and expanding proves the degree
claim. If a possible monomial is indexed by S, its coefficient
is the inclusion-exclusion sum of evaluations at masks
R⊆S, namely ∑_(R⊆S)(−1)^(|S|−|R|) H(1_R).
All these evaluations are in the asserted sparse grid; no
division is needed. The same formula works for differences
of functions and over Z/mZ.

More generally, allow arbitrary tournament arcs to have
their own bits. The support of every nonzero monomial is a
subgraph of a vertex-disjoint union of odd cycles: it can
have maximum degree at most two and no even cycle. Indeed
every monomial comes from the expansion of one of the
literal packs and is a subset of its arcs. Thus even for
arbitrary arc mutations one can recover the function from
masks supported on disjoint paths and odd cycles. This uses
at most (sum_{j\le n}\binom{\binom{n}{2}}{j}=\exp(O(n\log n)))
masks on an n-vertex tournament, compared to the whole
(2^{\binom{n}{2}})-cube. There is no claim that every
allowed support has a nonzero coefficient: cancellation
may create further constraints.

This tells us exactly which multi-contact interactions can
carry a hidden obstruction. A single new vertex cannot
support a genuinely cubic interaction in its own contact
row. A broader chart can acquire genuine interactions by
distributing the endpoints over several vertices, but the
sharp congruence scan still says nothing about realizing
every allowed absolute count.

Run `python replay_sparse_contacts.py` for a direct
two-apex finite difference replay.
