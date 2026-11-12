# One physical scalar as a lossless substitution-response code

**Status: PROVED in this note; finite controls are separate. Chan–Pak remains OPEN.**

## Inheritance, novelty, and question

The full ordered spanning path-cover profile
`F_T(r)=r! pc_T(r)` is already a substitution congruence, with the exact
quotient-walk contraction in [THM-3121](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
Its ordinal product/Newton coordinates are [THM-3166](../../../../01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md).
Actual physical tomography and contextual minimality are already Session176,
`session176_tournament_information_20261104/information/PROOF.md`, and the
specific symmetric triangle cap below is already Session215,
`session215_exotic_composition_20261107/fixed-caps/PROOF.md` §3. Session268,
`session268_composable_tournaments_20261111/modular/PROOF.md`, carries response
buckets with actual class multiplicities. None of these is rediscovered here.

The new strengthening is a **single-observation, constructive digit decoder**:
at known order n, one sufficiently long *fixed* positive tournament context
separates every possible path-cover profile simultaneously. No adaptive probe,
signed formal tournament, candidate catalogue, or Pisot seed restriction is
needed. This changes the many-probe tomography into a superincreasing physical
code. Its cost is a larger context and integer; no computational speedup is
claimed.

Live concepts: module substitution / positive path-cover coefficients /
superincreasing coding / exact signed correction / nonuniform attachment.
The connection is literal: cap coefficients are digit place values, and the
factorial path-cover bound supplies a universal digit range. Unlike mere
base-expansion analogy, every encoded integer is H of an actual strong tournament.

## 1. The inherited positive cap weights

Let TT_m be the transitive m-vertex tournament, and define

    C_m(T) = C3[T, TT_m, TT_m].

All blocks are nonempty. Write `g_s(m)=s! S(m,s)`, the number of surjections
from an m-element set onto an s-element set; put g_0=0. The inherited C3 kernel
says

    H(C_m(T)) = sum_(r=1)^n K_r(m) F_T(r),                  (1)
    K_r = 3g_r² + 2g_r g_(r+1) + g_(r+1)²
                 + g_(r-1)² + 2g_(r-1)g_r.              (2)

This is Session215's formula, now used in a regime with m much larger than n,
not m=r-1 for triangular tomography. Every coefficient is nonnegative.
`C_m(T)` is strong: between any two modules the directed triangle supplies a
route, and within a module a detour through the other two supplies a route.

Every ordered r-path cover concatenates to a vertex permutation with r-1
cuts; hence

    0 <= F_T(r) <= n! binom(n-1,r-1) <= B_n,
    B_n = n! 2^(n-1).                                    (3)

The concatenation is bijective with permutations and cuts whose uncut edges
are directed forward, so the bound does not count covers with a missing
factorial.

## 2. Universal one-count theorem

For n>=1 define the integer L_n by

    L_n = ceil(log2(72 (n+1) B_n)),   m_n = n L_n.          (4)

The logarithm need not be computed numerically: L is the least nonnegative
integer such that `2^L >= 72(n+1)B_n`.

**Theorem.** For every n>=1, every m>=m_n, and every 1<=r<=n,

    K_r(m) > B_n sum_(s<r) K_s(m).                       (5)

Consequently the single integer `Z_n(T)=H(C_(m_n)(T))`, together with n,
recovers the complete ordered path-cover profile of every n-vertex tournament
by ordinary greedy division. In particular it separates precisely the full
substitution-response classes, not tournament isomorphism classes.

**Proof of the quantitative inequality.** For s>=1,

    g_s <= s^m,
    g_s >= s^m - s(s-1)^m.                               (6)

The first inequality counts all functions; the second is the union bound
for a missing image. For 1<=r<=n, Bernoulli's inequality gives
`((r+1)/r)^n >= 1+n/r >= 2`. Thus for m>=nL,

    (r/(r+1))^m <= 2^(-L),
    (r+1) 2^(-L) <= 1/(72B_n) < 1/2.

It follows that `g_(r+1) > (r+1)^m/2` and therefore
`K_r > (r+1)^(2m)/4`. For s<r, the sum of the five displayed numerical
coefficients in (2) is 9, so `K_s <= 9(s+1)^(2m) <= 9r^(2m)`. There are
at most n such indices. Finally

    ((r+1)/r)^(2m) >= 4^L > 36 n B_n.

Combining these inequalities gives
`K_r > (r+1)^(2m)/4 > 9n B_n r^(2m) >= B_n sum_(s<r)K_s`.
For r=1 the lower sum is zero and positivity suffices. QED.

**Decoder.** Start R=Z. For r=n,n-1,...,1 set

    F(r)=floor(R/K_r),   R <- R-F(r)K_r.                 (7)

Inductively the remainder is the sum of lower places. It is nonnegative and,
by (3),(5), strictly smaller than K_r. Thus (7) recovers exactly the next
digit, without floating arithmetic or a subtraction of nearly equal approximations.
The final remainder is zero. The same argument works for ANY integer digit
vector in `[0,B_n]^n`; tournament realizability is not assumed by the decoder.
It must NOT be inferred from successful decoding of an arbitrary integer.

**Order separation.** If profiles F and G differ, let r be their largest
differing coordinate. The absolute contribution from all smaller differences
is at most `B_n sum_(s<r)K_s < K_r`. Hence

    sign(Z(F)-Z(G)) = sign(F(r)-G(r)).                   (8)

So reverse-lexicographic profile order is literally realized as numerical
Hamiltonian-count order in one common positive context.

**Size/bit cost.** Formula (4) has `L_n=O(n log(n+1))`, so the output graph has
`n+2m_n=O(n² log(n+1))` vertices. Computing g_1,...,g_(n+1) by the Stirling
recurrence costs O(m_n n) integer operations; there is no graph enumeration
in this encoder once F is supplied. Equation (2) yields
`Z <= 9n B_n (n+1)^(2m_n)`, so its bit length is
O(n² log²(n+1)). The underlying F remains expensive on arbitrary inputs.
The bounds are deliberately conservative, not minimal-cap or optimal-bit claims.

### Certified shorter caps, without changing the proof obligation

The crude common bound B_n can be replaced computationally by the individual
bounds `b_s=n! binom(n-1,s-1)`. Search m from1 through m_n and accept the first
one satisfying all exact integer inequalities

    K_r(m) > sum_(s<r) b_s K_s(m),  r=1,...,n.            (9)

Termination follows from (5), since b_s<=B_n. Every accepted cap has the same
greedy decoder. This is an exact finite compiler for a sufficient digit-box
condition, not a conjectured asymptotic approximation. First sufficient caps
at n=1,2,3,4,5,6,8,10,16,20 are respectively

    1,3,6,11,18,27,54,91,276,465.

For example the universal formula gives m_5=100, but the certified digit-box
compiler gives m=18, so41 vertices suffice for a single observation separating
all five-vertex profiles. These are NOT asserted minimal physical separators:
realizable profiles obey additional dependencies, while (9) separates every
vector in the much larger rectangular digit box.

## 3. An exact coded grammar, and its unavoidable sidecars

For any known orders, replace a full profile by the pair `(n,Z_n)`. Decode,
apply the inherited quotient-walk profile contraction, then encode at the new
order. This is an exact change of coordinates on the substitution-response
algebra, valid for every finite quotient and arbitrary nesting. It is NOT a
new family-generation theorem: prime quotient supply remains whatever the
underlying grammar actually provides.

For A000568, distinct ordinary classes sharing a code must retain their actual
class multiplicity and the action of quotient automorphisms, as in Session268.
Using one code as one isomorphism class undercounts. For metagraph growth,
nonuniform attachment may split a code fibre. The next concrete witness shows
why a complete growth grammar needs more than even this lossless substitution
code.

## 4. Sharp physical loss: equal complete profile, different next-count menu

Take the two five-vertex tournaments with outgoing bit rows

    A = (14,28,24,16,1),    B = (10,28,25,16,1).

They are inherited Session224 H9 first-appearance witnesses, not new graphs.
The new control explicitly compares their substitution and growth responses.
Both have

    F = (9,66,186,240,120).

They are nonisomorphic: their unique outdegree-two vertices have outgoing
neighbor degree multisets `{1,1}` and `{1,3}`. They therefore receive the SAME
one-count code at every m, not only at the safe bound. All uniform substitution
contexts preserve this equality.

Nevertheless the complete one-new-vertex H menus are

    A: {9,15,17,19,23,25,27,29,33,37,41},
    B: {9,15,17,19,23,29,31,33,37,43}.

For example, a new vertex beating only vertex 2 of B gives H=31; no attachment
to A gives31. The same bit cut applied to A gives23. This is an entire menu
comparison over all32 cuts, not only different values for one labelling.
Thus substitution congruence is **strictly weaker than the equivalence needed
for unrestricted vertex growth**, even for equal-order strong first creators.
The distinction is operational, not a failure of the new decoder.

## 5. Chan–Pak boundary and useful next operation

The code provides exact recursively composable scalar state labels, ordered
by an explicit dominant-coordinate law. It may be used to identify, cache,
or compare entire families without storing a list of physical probe outcomes.
It does not provide a tournament for an arbitrary odd code, or an arbitrary
odd original H. A universal one-count observer is a forward map and inverse
*on its image*, not a surjectivity theorem. The H9 witness further shows why
optimizing arbitrary attachment using only substitution codes can lose real
construction options. A next target compiler should retain code plus its
rooted attachment response, or an exact set of those responses.

No completeness result for Chan–Pak, no finite prime catalogue, no exact large
A000568 value, and no asymptotic enumeration speedup is claimed by this lane.

## Reproduction

Run `python3 verify.py` and `python3 -O verify.py` in this directory. The verifier
has explicit exceptions instead of assertions. It checks every labelled
tournament through order5, independent permutation-cut profiles on the small
orders, literal small-cap counts, the H9 menus, and safe-bound coding/decoding.
These finite controls support the displayed all-order argument; they are not
its proof. See the machine-readable receipt for exact test cardinalities.
