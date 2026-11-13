# Projective observations become exact after two affine response refinements

Status: **PROVED IN NOTES, G1–G2, after independent audit.** The finite
controls below are **FINITE-EXACT**. This note generalizes the separately
frozen [64-state norm example](NORM_CLOSURE.md); it does not change its source,
proof, receipts, or hashes. The formulation was proposed by the root agent
and independently checked here, including the small-field boundaries.

## G1. Projective-to-affine observer theorem

Let q be a prime power, V a d-dimensional vector space over F_q with d≥2,
and v a nonzero vector. Let B be an invertible F_q-linear map such that
Bv is not in the line L=F_q v. Put A(x)=x+v. Let P0 be the partition

    {0}, and F_q* x for each projective point [x].

Thus P0 has 1+(q^d−1)/(q−1) cells. Given a partition P, define its simultaneous
response refinement by the label

    R(P)(x) = (P(x), P(Bx), P(Ax)),

identifying equal labels. Set P_(k+1)=R(P_k). The exact strict-refinement
sequences are:

| Field size | Cell counts | Number of strict rounds |
|---|---|---:|
| q=2 |2^d, already discrete|0|
| q=3 |1+(3^d−1)/2 →3^d|1|
| q≥4 |1+(q^d−1)/(q−1) →q^d−q+3 →q^d|2|

In particular, for q≥4 the second round is both sufficient and necessary.
For q≥3, the first refined partition has every point outside L as a singleton,
together with {0}, {−v}, and L\{0,−v}. At q=3 the last set is a singleton;
at q=2 it is empty and is **not** counted as a cell. The middle formula
q^d−q+3 therefore applies only for q≥3.

### Proof candidate

Linearity and invertibility imply that B sends each nonzero scalar orbit to
a nonzero scalar orbit. Thus P0(Bx) is a function of P0(x), and B contributes
no extra distinction to the first refinement.

Fix an initial cell F_q* x outside L. The vectors x,v are linearly independent.
If cx+v and ex+v have the same P0 observation for c,e in F_q*, their nonzero
values differ by some scalar lambda in F_q*. Hence

    cx+v = lambda(ex+v).

Independence gives lambda=1 and c=e. Thus A splits every such cell into
singletons. There are q^d−q vectors outside L.

On L, write a state as av. Its initial observation distinguishes a=0 from
all other a, while the observation after A distinguishes a=−1 from all
other a. Consequently {0}, {−v}, and the remaining q−2 points are exactly
the first-round cells on L, omitting the last cell when empty. This proves
the stated P1 description and counts. For q=2, P0 was already discrete;
for q=3, P1 is discrete.

Suppose q≥4. The only nonsingleton cell after one round is
L\{0,−v}, of size q−2≥2. Under B its elements av go to aBv, which are
distinct points outside L. Every one is already a singleton of P1. The B
response therefore separates the surviving cell in the next round. All
other points stay separate, so P2 is discrete and stable. Both rounds are
strict: the first splits a nonzero scalar orbit outside L, and the second
splits the surviving q−2 points. This proves the exact depth. QED candidate.

### Minimal observation meaning

For any deterministic operations on a finite set, equality in P_k means
equality of the P0 observation after **every word of length at most k** in
those operations. This follows by induction: equality at k+1 retains equality
at k and adds equality at k after each allowed first operation. Therefore the
stable limit is the coarsest partition refining P0 through which the declared
operations descend. G1 proves that this particular limit is discrete, with
the stated sharp depth. It does not impose a two-round bound on arbitrary
operator refinement, sheafification, or a tournament metagraph.

The source is an affine vector space with scalar-orbit observations. The map
is the observation quotient followed by operation-response refinement. P0
retains zero and projective direction, losing the nonzero scalar coordinate.
Translation fixes a distinguished affine unit v; the line-moving linear map
exposes the remaining scalar ambiguity. The theorem concerns this exact
interface. The earlier norm example is its q=4,d=3 instance: F64 is a
three-dimensional F4-space, the norm fibres are its nonzero F4-scalar orbits,
v=1, and multiplication by b∉F4 moves the line F4·1.

## G2. Why the line-moving hypothesis matters

This hypothesis is sufficient, not claimed necessary. If B preserves L,
the displayed proof no longer exposes the residual cell, and either failure
of complete recovery or a longer repair depth can occur.

- Over F4² with B=identity and v=(1,0), counts are **6 →15**, stable. The two
  points av with a in F4\F2 stay together forever: adding1 exchanges them.
- Over F5² with the same choices, counts are **7 →23 →24 →25**. The final
  observer is exact, but needs **three** strict rounds.
- Over F7² the analogous sequence is **9 →45 →46 →47 →48 →49**, requiring
  **five** strict rounds.

These claims have a short direct explanation. On L, a word in A and the
identity simply adds jv for some integer j. The initial observation is zero
exactly when a+j=0 in the field; otherwise it is the common nonzero-line
label. Let p be the characteristic. Through depth k, the partition on L
consists of the individual cells

    −jv, for 0≤j≤k,

with repeats modulo p, and one residual cell containing all other states
(omitted if empty, and already a singleton if it has size one). Outside
L every point was already singled after one round. If q=p≥3, the final
remaining cell has one element at k=p−2, so precisely p−2 strict rounds
suffice. These examples give arbitrarily large repair depth even for affine
permutations and scalar-orbit starting observations, once the line-moving
hypothesis is removed.

If q>p, after k=p−1 the p points F_p v are singled but the q−p points of
L\F_p v remain one nonsingleton cell forever. The stable total is
q^d−q+p+1. For q=4 this is15 in dimension2. This boundary also distinguishes
full recovery from a merely stable observer. The q=2 initial-discrete case
is separate; the prime formula p−2 then correctly gives zero strict rounds.

## Finite controls and reproducibility

The independent standard-library checker [check_general_closure.py](check_general_closure.py)
uses a direct Cartesian vector-space representation. F4 is
F2[X]/(X²+X+1); it imports no code from the old sextic checker.

It enumerates **every** invertible2×2 matrix and nonzero v over q=2,3,4,5,
checks that B alone descends to P0, and checks the exact theorem sequence
whenever Bv leaves F_q v. Matrix invertibility is tested by the complete
image set, independently of a determinant formula. The counts of matrices
are checked against (q²−1)(q²−q). It also identifies the unique first-round
residual cell, exercises the three line-preserving hostiles above and the
q=2/3 boundaries, and compares refinement to direct bounded-word observation.
A separate coordinate-cycle matrix over F4³ gives22 →63 →64 without using
multiplication in the older F64 encoding.

From the repository root:

```sh
python3 05-knowledge/results/integer32_twoedge_20261113/census/check_general_closure.py
python3 -O 05-knowledge/results/integer32_twoedge_20261113/census/check_general_closure.py
```

Normal and optimized outputs are stored separately. Their receipt records
its source hash, complete finite universe, cell-count sequences, and active
check count. Finite verification supports G1 in those universes; the
all-q,d claim rests on the displayed proof candidate and its independent
audit. No external literature-priority claim is made.
