# Session138: actual-kernel structure for deleted-triangle zero rigidity

**Status:** PROVED identities and strict comparison; FINITE-EXACT implementation
controls; a stronger FORMAL-NETWORK firewall. **The all-order zero-rigidity
conjecture is still OPEN. No actual tournament counterexample is claimed.**

Baseline: `/tmp/session137-source`, commit
`252dd4a839b945f0b2be77b4f181dc6241aa6868`. Inputs are the Session137
`coverage/{README,SYMMETRY_EDIT,TRIANGLE_CARD_RIGIDITY}.md`,
`zero_closure_controls.py`, and the independent endpoint implementation in
`explore.py` (now frozen locally as `reference_paths.py` for portable audits). No Git mutation, fetch, push, or expanded quotient census was done.

## Inheritance / board

* Anchor: scalar unique-zero reconstruction for `Q[C3^k]`, Q strong. Closest
  proved mechanisms are marked-boundary divisibility and unique homogeneous
  pair/canonical family completion. Lost coordinate: actual path-kernel
  realizability, not positivity or symmetry.
* Niche: exact homogeneous-pair lift of an ear polynomial. Cheap hostile:
  Session137's formal positive-kernel perturbation.
* Wildcard: quotient run-word edge derivatives. Sidecar: ordered boundary
  block identities, not just local path-cover scalars.
* Hostile retained: nonstrong `C3 => C3` has nonoriginal zero completions;
  generic positive formal kernels can have them even when the quotient is
  strong. No finite check is extrapolated to all orders.

## 1. Homogeneous-pair lift and a strict cyclic comparison

Let D be **any** finite tournament, possibly empty. For a set R of its
vertices, adjoin a homogeneous pair u->w, each beating exactly R in D; call
this parent B. For an arbitrary S subset of V(D), compare its two completions

    C(R,S): x->u->w->x,  x beats exactly S in D,
    M(R,S): u->x->w,     x beats exactly S in D.

Write the actual ear polynomial of D as

    F_D(z) = H(D) + sum_(i->j in D) c^D_ij z_i(1-z_j),

where `c^D_ij=c^D_ji>0`, with `A^D_i=sum_(i->j)c^D_ij`.
The gradient means the gradient of this specified multilinear quadratic,
not a finite-difference convention at a Boolean point.

### Theorem 1 (all orders)

With r=1_R and s=1_S,

    H(C(R,S))-H(M(R,S))
      = 2F_D(r)+(s-r)^T grad F_D(r)
      = 2H(D) + sum_(i->j) c^D_ij
                      [r_i(1-s_j)+s_i(1-r_j)]
      >= 2H(D) > 0.                                      (1)

Equality in the displayed lower bound is exact precisely when every bracket
is zero. This includes D empty and singleton, when the edge sum is empty.
For strong D of order at least three, equality occurs only at
`R=S=empty` or `R=S=V(D)`; this last classification is not needed below.

### Elementary proof (no new literature input)

Swap u,w in the cyclic completion. It now agrees with the middle completion
except that u->w has been reversed. A Hamiltonian word of the middle
completion that uses adjacent u,w cannot have x immediately before u or
after w. Reverse the adjacent pair: all remaining boundary arcs are unchanged
because u,w have identical orientations to D. This injects those words into
words of the cyclic completion after the swap. Words not using the changed
edge agree. The unmatched words are exactly those containing one of the two
consecutive triples

    x,w,u       or       w,u,x.

They cannot both occur in one word. Contract either triple. If `sigma_j` and
`epsilon_i` count D-paths starting at j and ending at i, and `Q_ij` counts
D-words with the single exposed adjacent gap i,j, their total is

    sum_j (r_j+s_j) sigma_j
    + sum_i (2-r_i-s_i) epsilon_i
    + sum_(i,j) [(1-s_i)r_j+(1-r_i)s_j] Q_ij.             (2)

For nonempty D, the standard insertion formula is

    F_D(z) = sum_j z_j sigma_j + sum_i(1-z_i)epsilon_i
             + sum_(i,j)(1-z_i)z_j Q_ij.

Polarizing this quadratic gives exactly (2) as
`2F_D(r)+(s-r)^T grad F_D(r)`. For empty D there are two unmatched words,
so the same identity holds with F_D=1. Expansion of the orientation-supported
quadratic gives the nonnegative expression in (1).

For the optional strong-D equality classification, put
A=R intersect S, B0=R minus S, C0=S minus R, E=outside(R union S).
Zero brackets imply that A receives every edge from its complement,
`|B0|,|C0|<=1`, and E beats both B0,C0. Strongness excludes all nonconstant
possibilities. This argument is unnecessary for strict positivity itself.

### Corollary: exact capacity recursion

For the actual capacity kernel of B,

    c^B_uw = 2F_D(r),
    c^B_ui-c^B_wi = sum_j c^D_ij r_j-A^D_i
                  = -partial_i F_D(r).                   (3)

Indeed, directly expanding the B-ear polynomial, cyclic minus middle is

    c^B_uw + sum_i(r_i-s_i)(c^B_ui-c^B_wi).

Compare its coefficients in s with (1); put s=r to get the pair entry.
This is a realizability identity connecting **two actual induced kernels**,
not a free constraint on an arbitrary positive matrix.

## 2. Internal homogeneous triangle identity

If a,b,c form a homogeneous directed triangle in any tournament B, then

    c^B_ab = c^B_bc = c^B_ca = 2H(B-a).                  (4)

To prove it, contract the marked free adjacency a,b in both orders in the
exposed-word interpretation of c_ab. With the third vertex c in a separate
block run, both orders of a,b are allowed, giving twice the analogous
count in B-a. If c is consecutive to the contracted pair, one order of a,b
allows both boundary orientations and the other allows neither; again the
sum is twice the single permitted orientation in B-a. External arcs agree
by homogeneity. Triangle rotation gives the other equalities. Equivalently,
the local ordered run weights of the exposed pair are twice those of TT2.

## 3. Actual quotient-kernel compiler

For a quotient Q, let a block run be a maximal consecutive segment of a
Hamiltonian word lying in one block. A run-word has no consecutive equal
block labels, and every transition follows the quotient. If block i occurs
r_i times, the local number of ordered internal path covers is

    C3:   f(1),f(2),f(3) = 3,6,6;
    TT2:  f(1),f(2)      = 1,2;
    K1:   f(1)           = 1.

These follow by enumerating run-length compositions, or cutting the directed
triangle into ordered path components. In particular C3's weights are the
falling factorials `(3)_r`; their successive ratios are `(3,2,1)`.

Give each ordered quotient edge i,j an independent variable a_ij. Let Z_Q(a)
be the sum over all allowed multiplicity ranges of all adjacent-distinct
block words, weighting transitions by a_ij and block i by f_i(r_i).
At the 0/1 adjacency matrix of Q,

    Z_Q = H(Q[blocks]).                                  (5)

This is an exact run contraction, not a substitution-product assumption.
For all-C3 blocks, equivalently

    H(Q[C3^k]) = 3^k sum_r N_Q(r) product_i (1,2,2)_(r_i).

For a deleted-triangle parent B, the root block is TT2 and every other
block is C3. Let `G_ij = partial Z_Q / partial a_ij` at Q. Then:

* For distinct intact blocks A,B and any a in A,b in B,

      c^B_ab = (G_AB+G_BA)/9.                            (6)

* For the broken root block r={u,w} and a vertex a in intact block A,

      c^B_ua+c^B_wa = (G_rA+G_Ar)/3.                     (7)

* Their difference is (3), computed from the smaller all-C3 quotient Q-r.
  Equations (3), (4), (6), (7) recover the entire actual parent kernel.

Proof of (6),(7): differentiating marks one free transition between the two
blocks. All other transitions remain valid; this also correctly counts a
single forbidden transition when the differentiated variable is zero.
Independent C3 rotations divide the marked endpoint choices equally, by9
for two intact blocks and by3 for one. The TT2 endpoints are summed, not
incorrectly declared symmetric. Internal entries use (3),(4).

`quotient_kernel.py` implements (5)--(7) by a forward/backward weighted-walk
DP. The state is (run counts,last block), with at most 4^k*k entries, and all
first derivatives follow by one forward/backward transition summation. It
uses Python arbitrary-precision integers. It does not enumerate 2^(3k)
vertex subsets. This is a new implementation of an exact reduction, not
an expanded zero-rigidity census.

## 4. What the new identities exclude—and what they do not

The old Session137 formal perturbation violates (1). For outside cut224,
its cyclic cut225 has formal value198963, while its middle cut226 has
formal value199053. The gap is -90, whereas the required lower bound is54.
The actual gap for these same cuts is1404. The minimum actual twin gap over
all outside cuts is432.

But imposing the new identities still does **not** prove zero closure.
`formal_twin_firewall.py` starts with the same actual order-eleven parent
and instead changes just two intact-block cross-capacity orbits:

    block1 versus block3: 13770 -> 14268  (all nine entries),
    block2 versus block3: 10002 ->  9504  (all nine entries).

It preserves positivity, independent C3 symmetries, every actual pair-deletion
lower bound, exact ear gcd6, total capacity582822, every cardinality-layer
mean, parent H, original child H, **all twin capacity differences and the twin
pair capacity**, **every internal homogeneous-C3 identity (4)**, and all the
strict cyclic-versus-middle inequalities (1). The minimum twin gap remains432.
Nevertheless it has formal target cuts

    480,736,1248,1793

at H=241137, where only1793 is the actual canonical zero. The actual value
at480 is239643; its block profile is (0,3,1), so the changed cross orbits
have coefficients3 and0 there, producing the required increase1494. At the
canonical profile both orbits have coefficient9, so their changes cancel.
Their equal orbit sizes preserve total capacity and hence every layer mean.

This is still **not a tournament counterexample**. It violates the exact
quotient derivative values in (6). It locates the remaining missing sidecar
more sharply: the jointly realizable **cross-block** derivative tensor,
not merely the homogeneous-pair gradient and internal block entries.

## 5. Controls and replay

`python3 -B twin_audit.py` checks (1) literally by full vertex-subset path DP
for every labelled D through order4, every homogeneous-pair attachment R,
and every outside attachment S: 16,933 cases. It also prints the old formal
countercontrol's explicit negative twin gap.

`python3 -B quotient_kernel.py` compares (5) and every entry reconstructed
by (3),(4),(6),(7) with independent full-vertex endpoint DP, for every labelled
quotient through order4: 75 quotients and285 deletion rows, including strong
and nonstrong controls. This is a coefficient audit, not a new zero census.

`python3 -B formal_twin_firewall.py` verifies every listed invariant of the
new stronger formal countercontrol and prints the actual/formal zero cuts.

Matching `.out` files retain semantic SHA256 hashes. The independent endpoint
DP is imported only by the optional audit, not by the quotient compiler's
functions. Remaining target: use the exact cross-block derivative identities
to prove scalar zero uniqueness or construct an actual strong-quotient
counterexample. This packet does neither.

## 6. Later proof extensions and exact boundaries

* **`ORDINAL_ZERO_CLASSIFICATION.md`** proves every H-zero completion of
  `TT_k[C3]` for all k: canonical, or a middle-pair flip plus one wrong
  incidence in an adjacent triangle. Counts are1 at k=1,4 at end roots,
  and7 at interior roots. The bounded audit through k=8 is of this closed
  nonstrong theorem, not an expansion of the strong census.
* **`PROFILE_DICHOTOMY.md`** strengthens it: every star completion of this
  ordinal module has its whole ordered path-cover profile entrywise below
  or above the original, never equal unless canonical. An all-order bounded
  core reduction uses15 explicitly listed cores; a separate C++ literal
  verifier checked all3,633,840 leaf permutations. Consequently
  `C3[TT_k[C3],A,B]` is zero-rigid under star edits confined to the ordinal
  module, for arbitrary nonempty A,B. Outside attachment changes remain open.
* **Actual hostile to extrapolation:** `strong_profile_hostile.py/.out`
  gives a strong Q4 card with cut996 whose H falls by3726 while all higher
  nonconstant cover coordinates rise. Thus the profile dichotomy really
  needs its ordinal hypothesis. Cuts867 and995 have the opposite eventual
  profile sign, supplying further genuine mixed-profile controls.

Exploratory parameter scans and scratch precursors are excluded from the curated
public manifest. The actual strong-Q4 mixed-profile hostile is retained because
it is the explicit boundary of the profile theorem. `public_manifest.json` and
`SHA256SUMS` identify the frozen, self-contained public packet. The all-order
ordinal classification and profile/restricted-strong-rigidity proofs were
accepted by the independent Session138 reconciliation audit; arbitrary strong
quotient zero rigidity remains OPEN.
