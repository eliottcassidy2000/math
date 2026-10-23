# Even every transitive motif and every outneighborhood transitive profile can miss H

23 October 2026. **FINITE-EXACT:** the explicit order-fifteen pair and all
counts below. **PROVED:** preservation of the complete transitive polynomial
and the multiset of rooted outneighborhood polynomials under source padding
and strong cyclic lifting. The inherited exact-E/H propagation then gives
an every-order obstruction. **NOT CLAIMED:** agreement of all induced motifs,
a counterexample within the doubly regular family, or failure of complete
authenticated port reconstruction.

## 0. Incoming priority and the precise new boundary

Incoming commit
[`f2cac6db1bb8ef780b9211d1cf3f7adc109c5fca`](https://github.com/eliottcassidy2000/math/commit/f2cac6db1bb8ef780b9211d1cf3f7adc109c5fca)
contains the [minimal-signature packet](https://github.com/eliottcassidy2000/math/blob/f2cac6db1bb8ef780b9211d1cf3f7adc109c5fca/05-knowledge/results/metagraph_minimal_signatures_20261023/README.md)
and its independent
[`neighborhood_defect.py`](https://github.com/eliottcassidy2000/math/blob/f2cac6db1bb8ef780b9211d1cf3f7adc109c5fca/05-knowledge/results/metagraph_minimal_signatures_20261023/neighborhood_defect.py).
That work already proves and tests that a transitive-five motif repairs its
specific 31-vertex doubly regular hostile: the two total counts are 10,695
and 10,395. It also identifies the rooted quantity

    B_v = number of transitive four-subsets in N^+(v)
        = number of transitive five-subsets with source v.

Those identities and that repair have incoming priority. The packet
explicitly says the repair is scoped to the named pair, not universal.
Nothing here contradicts it.

Our seed is the order-fifteen exact-E collision already recorded in
`INCIDENCE_COLLISION_PROPAGATION.md`, not a newly discovered pair. The new
calculation shows that adding **all** transitive motif counts, and even
**all** outneighborhood transitive polynomials, still does not repair that
pair. The propagation below carries this strengthened collision to larger
orders. These regular circulants are not doubly regular, so this does not
settle completeness inside the doubly regular subclass.

## 1. Definitions and the complete finite witness

Let t_k(T) count induced transitive k-subsets, with t_0=1, and define

    F_T(z) = sum_(k=0)^n t_k(T) z^k,
    R(T) = multiset_(v in T) F_(T[N^+(v)])(z).

Each nonempty transitive set has one source, so the all-size identity is

    t_(k+1)(T) = sum_v t_k(T[N^+(v)]),
    F_T(z) = 1 + z sum_v F_(T[N^+(v)])(z).              (1)

Thus R determines F, but the full multiset R usually contains more
information than their sum. We retain both in the propagation formulas
because the distinction is useful.

On Z/15 define T_i by x -> x+d for d in the following sets:

    S_1 = {14,13,12,11,5,6,8},
    S_2 = {14,13,12,4,10,9,8}.

The translations act transitively on vertices. Both are regular of score
seven and strong, since step 14 gives a directed Hamiltonian cycle. The
exact edge-incidence signature

    E(T) = multiset_v multiset_(v->w) (c_w,c_vw)

is identical: every c_v is 28, and each outgoing row has edge incidences
3,3,3,4,4,5,6. Here c_v counts cyclic triangles through v and c_vw counts
those through the unordered pair {v,w}.

Independent enumeration of **every subset** gives, for both seeds,

    F_T = 1 + 15z + 105z^2 + 315z^3 + 375z^4 + 150z^5 + 15z^6.

All coefficients of degrees seven through fifteen are zero. Independently
counting every subset of each of the thirty outneighborhoods gives

    F_(N^+(v)) = 1 + 7z + 21z^2 + 25z^3 + 10z^4 + z^5

for every v in either seed. In particular B_v=10 everywhere, and the total
transitive-five count is 150 in both. Equation (1) supplies a separate
cross-check of all these coefficients.

Nevertheless their Hamiltonian counts are different:

| Quantity | T_1 | T_2 |
|---|---:|---:|
| H(T) | 196264905 | 196277235 |
| H(T[N^+(v)]), every v | 101 | 89 |
| Ordered spanning two-path covers | 4857720210 | 4857619470 |

These are exact subset-DP counts, not estimates. The H difference proves
that the seeds are nonisomorphic. Both therefore give the same complete
signature (E,F,R), but not the same H.

The common-outneighbor histogram inside each seven-vertex outneighborhood
is also identical:

    {0:5, 1:8, 2:7, 3:1}.

Consequently sum_pairs binom(common_out,2)=10, independently agreeing with
the literal transitive-four count. All four induced four-vertex motif counts
inside these outneighborhoods agree as well. We do not apply the incoming
DRT-specific variance formula to these non-DRT neighborhoods.

## 2. Source padding preserves the stronger signature

Let S(T) be T with a new dominating source. A subset of S(T) is transitive
exactly when its old vertices induce a transitive subset of T, whether or
not the new source is included. Therefore

    F_(S(T)) = (1+z) F_T.

Old vertices have unchanged outneighborhoods, while the new source has
outneighborhood T. Hence

    R(S(T)) = R(T) multiset-union {F_T}.                (2)

Simultaneous equality of F and R is preserved. The inherited incidence
propagation proves that exact E is preserved too, and H is unchanged.
Iterating (2) therefore gives different-H collisions of (E,F,R) at **every
order n>=15**. This padding need not preserve strong connectivity.

A sink variant also holds: adding one sink multiplies every old rooted
outneighborhood polynomial by 1+z and adds the polynomial 1 for the new
sink. Its global F again multiplies by 1+z.

## 3. Strong cyclic lifting preserves every rooted profile

For r>=1 form

    L_r(T) = C3[T,{a},P_r],    T -> a -> P_r -> T,

where P_r is a transitive r-vertex tournament. This is strong for every
nonempty T. Write P(z)=(1+z)^r.

A transitive set cannot meet all three cyclic blocks. Conversely, a subset
meeting at most two blocks is transitive exactly when each occupied block
is transitive. Inclusion-exclusion over the three allowed block-pairs gives

    F_(L_r(T)) = F_T F_a + F_a P + P F_T - F_T - F_a - P + 1,
    F_a = 1+z.

Equivalently,

    F_(L_r(T)) = (P+z) F_T + z(P-1).                  (3)

This is exact at every order, including the empty-subset coefficient.
It does not count merely maximal transitive sets.

The rooted outneighborhood polynomials transform just as explicitly:

- For an old vertex v, its outneighborhood is N_T^+(v) followed by the sink
  a. Its polynomial is (1+z)F_(N_T^+(v)).
- The outneighborhood of a is P_r, giving P(z).
- A vertex of P_r with j later vertices, j=0,...,r-1, sees the ordinal sum
  P_j -> T. Its polynomial is (1+z)^j F_T.

Thus

    R(L_r(T))
      = multiset_v {(1+z)F_(N_T^+(v))}
        union {P}
        union {(1+z)^j F_T : 0<=j<r}.                 (4)

Equations (3)-(4) prove that simultaneous equality of F and R survives
strong cyclic lifting. This preserves the entire unbounded transitive
family, not just a fixed maximum motif size.

The existing exact-E propagation proves E(L_r(T_1))=E(L_r(T_2)). Its
Hamiltonian calculation, independently rechecked here at the seed via
ordered two-path covers, is

    H(L_r(T)) = (2^r+1)H(T) + (2^r-1)F2(T),

where F2 is the ordered spanning two-path-cover count. In the T_2-T_1
order, the seed differences are 12330 and -100740. Consequently

    H(L_r(T_2)) - H(L_r(T_1)) = -88410*2^r + 113070 < 0

for every r>=1. There is no exceptional positive r. We obtain **strong**
(E,F,R) collisions with different H at order fifteen and every order
n>=17. The strong order-sixteen case remains unspecified by this argument.

## 4. A small nontransitive motif does repair this particular pair

The negative result is about transitive motifs, not all motifs. In each
seven-vertex outneighborhood, count five-subsets whose five internal
outdegrees are all two: the regular five-tournament R_5. The exact counts
are

    T_1 outneighborhood: 0 copies,
    T_2 outneighborhood: 1 copy.

The checker counts this directly from degrees without relying on an
isomorphism classifier. Thus a regular-five local motif distinguishes the
pair even though its entire local transitive polynomial does not.

Adjoining the unique source gives the six-vertex motif K_1 -> R_5.
A direct enumeration of all six-subsets confirms global counts **0 versus
15**. These are also the sum of the rooted counts above because that
six-tournament has a unique source. Several other nontransitive five-motif
counts, both locally and globally, also differ; full lists are retained in
the receipt. In particular we make no claim that every five-motif agrees.

Since all induced motif counts through size four agree inside the two
seven-vertex outneighborhoods, five is the smallest size of an unrooted
induced motif **inside those neighborhoods** that separates them. This
restricted local minimality is not minimality among global motifs, rooted
records, or arbitrary invariants.

The research lesson is precise: repeatedly adding larger transitive motifs
can miss a stable ambiguity. A different motif family can repair it. A
successful repair on one hostile does not certify universal realization,
and none of these compressed records replaces the exact signed port test.

## 5. Independent finite checks

`check_transitive_profile.py` reconstructs both seeds from the difference
sets, imports no research producer, and verifies:

- every subset of each seed and each root outneighborhood;
- the unique-source coefficient identity (1);
- exact E, full H, neighborhood H and ordered two-path-cover counts;
- the common-out histograms and all local four-/five-motif classes;
- the regular-five repair by direct degree tests and all global six-subsets;
- formulas (2)-(4) on every labelled tournament of orders one through four,
  for r=1,2,3;
- literal order-sixteen source and order-seventeen cyclic extensions of
  both fifteen-vertex seeds, including every rooted transitive profile.

Run

    python check_transitive_profile.py
    python -O check_transitive_profile.py --out transitive_profile_receipt.optimized.json

Only the Python standard library is required. Explicit checks survive
optimized mode. The finite computations test the identities on substantial
controls; the proofs, not their finite range, supply the all-size claims.
