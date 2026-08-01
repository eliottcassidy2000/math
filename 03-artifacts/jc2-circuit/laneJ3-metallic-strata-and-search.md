---
lane: J3
title: "Metallic and maximal-alternation strata for JC(2), and an actual search: the metallic stratum is NOT distinguished, and the two readings of 'metallic' answer the reciprocal-closure question oppositely"
status: mixed -- see the per-claim status table (headline result is a clean NEGATIVE)
script: 04-computation/jc2_metallic_search_laneJ3.py
output: 05-knowledge/results/jc2_metallic_search_laneJ3.out
script_sha256: a57e574df5d4177efa8e49b354a575ce318720186a951d0bf74b772cad872ea9
output_sha256: ab5ad93cf5a353096a6bbe125c58f7929fbf69a84b2278783ca7e2c7a515c407
hash_basis: LF-normalized bytes
checks: 26/26 pass
depends_on:
  - THM-3000 (circuit convention, sec 1)
  - THM-3003 (antipodal rigidity: R_k = R_(d-k) <=> reciprocal-closed roots)
  - THM-3004 (sign-change law 2K-3 for K distinct roots; doubled-geometric attainment 3b(ii))
  - THM-3010 (metallic recurrences attain maximal circuit alternation, via Simson)
  - 03-artifacts/jc2-circuit/laneJ2-leading-form-circuit-gate.md (J2-1..J2-9)
  - L0/L1/L2 (upstream this session: leading form P_n = cH^a, Q_m = c'H^b; Jung--van der Kulk)
scope_guard: >
  MISTAKE-237.  NOTHING below is a reduction of JC(2) to a circuit problem, an
  equivalence between programs, or a bridge from NC2/GMC(2)/THM-3010.  Every
  statement is either (i) an exact identity about binary forms, (ii) a
  finite-exact search result with its window written out, or (iii) a
  STRATIFICATION of the hypothetical counterexample locus with the direction of
  implication stated.  In particular the appearance of the golden ratio in the
  exponent pair (a,b) (Lame / continued fractions) and its appearance in a root
  multiset (metallic orbits) are TWO DIFFERENT OBJECTS; no map between them is
  claimed or found.
---

# Lane J3 -- metallic strata and a real search

Setting is lane J2's.  `(P,Q)` is a Jacobian pair in `C[x,y]`, `n=deg P`,
`m=deg Q`, `g=gcd(n,m)`, `a=n/g`, `b=m/g`, `gcd(a,b)=1`, and by **L0**

    P_n = c H^a,   Q_m = c' H^b,   deg H = g   (the leading form).

By **L2** (Jung--van der Kulk) an automorphism has `a=1` or `b=1`, so a
counterexample needs `a,b >= 2` coprime.  `H` gets the repo circuit
(THM-3000 sec 1) through `A_k`, `h_k = A_k/C(g,k)`, `R_k = h_k^2/(h_(k-1)h_(k+1))`,
`c_k = log(R_k/R_(k-1))`.

---

## 0. Claim table

| # | claim | status | scope |
|---|---|---|---|
| J3-0 | **The leading-form circuit is EMPTY unless `g >= 3` and can change sign only from `g >= 4`.**  `c_k` runs over `k=2..g-1`, so it has `g-2` entries.  Combined with `a,b >= 2` coprime this means the circuit can say *nothing at all* about a counterexample until `(n,m) = (8,12)` | **PROVED** (index count) | all Jacobian pairs |
| J3-1 | the metallic quadratic `H_q = x^2 - qxy - y^2` has roots `{lam_q, -1/lam_q}`, product `-1`, and IS reciprocal-closed with `mu = -1` | **VERIFIED-EXACT** (`q=1..5`, exact in `Z[lam]`) | `g=2` |
| J3-2 | but at `g=2` the circuit is **empty**, `R_k=R_(g-k)` is **vacuous**, and `H_q` is `PGL_2(C)`-equivalent to `xy`: the metallic quadratic carries **no moduli whatsoever** | **PROVED** + **VERIFIED-EXACT** | `g=2` |
| J3-3 | root-metallic geometric orbits `{1,lam,...,lam^(g-1)}` are reciprocal-closed and have **exactly one** circuit sign change -- the **MINIMUM**, not the maximum | **VERIFIED-EXACT** (`q=1..4`, `g=4..9`) | positive geometric root sets |
| J3-4 | that sign word is **independent of the ratio**: `T=2,3,10,1000,101/100,1/3,7/5,10^6` all give the same word at each `g`.  Reason (exact): `e_k = T^(k(k-1)/2) gauss(g,k)_T` and `k(k-1)/2` is quadratic in `k`, hence killed by the third difference | **VERIFIED-EXACT** + structural proof of the quadratic-kill step | geometric root sets |
| J3-5 | THM-3010's maximal-alternation object is the **coefficient**-metallic form `h_k = a_(k+1)`; it attains the maximum `g-3` sign changes for every `q=1..5`, `g=5..10`, has all `A_k>0` and `K=g` distinct roots -- and is **NEVER palindromic**, hence lies strictly **OUTSIDE** the reciprocal-closed stratum | **VERIFIED-EXACT** | `q=1..5`, `g=5..10` |
| J3-6 | **the two readings of "metallic" answer lane J2's question oppositely.**  Root-metallic is inside the reciprocal-closed stratum and minimally alternating; coefficient-metallic is maximally alternating and outside it.  THM-3010's maximal alternation therefore does **not** transfer to the reciprocal-closed stratum through the word "metallic" | **PROVED** (from J3-3 and J3-5) | -- |
| J3-7 | the family that *is* both reciprocal-closed and maximally alternating is THM-3004 3b(ii)'s **doubled geometric** family, and maximality switches on at explicit algebraic thresholds in the ratio: `T* = 2.4931270138811760006...` (`K=3`, minimal polynomial of degree 10) and `T* = 3.6952124348324745171...` (`K=4`, degree 16).  **Neither is a metallic number**, and metallic ratios fall on both sides | **VERIFIED-EXACT** (exact real roots of the exact numerators) | doubled geometric orbits, `K=3,4,5` |
| J3-8 | **scope finding for THM-3004.**  THM-3004's `2K-3` cap is a POSITIVE-coefficient statement.  For `K=2` with roots of **opposite sign** -- which is exactly `P_n = cH_q^a` for the metallic quadratic -- the observed sign-change count reaches `9` against the cap `1`, and repeatedly attains the absolute maximum `d-3`.  Same-sign `K=2` obeys the cap exactly (`1`) in all `28` instances | **VERIFIED-EXACT** (`a=2..8`; metallic `q=1..6` and rational `r=2,3,5,3/2,7/3,10`) | binary forms over `R`, some `A_k` of mixed sign |
| J3-9 | **Frobenius kernel lemma (char `p`).**  `Jac(P,R)=0` for every `R` in `F_p[x^p,y^p]`, and `ker(L_P) = {Q : Jac(P,Q)=0}` is a **ring** containing `F_p[x^p,y^p]` and `P`.  Hence over `F_p` the set of mate degrees is unbounded, and "`(n,m)` admits a Jacobian pair" is a *different question* in char `p` | **PROVED** + **VERIFIED-EXACT** | any `p` |
| J3-7b | within the doubled-geometric family the **golden** ratio is the metallic value that fails maximality hardest (sub-maximal at `K=3,4,5`), silver next, and only `q >= 5` is maximal at `K=5`.  THM-3010 sec 3's "`phi` is the extreme case" is a statement about `h`-sequences and **inverts** when read on roots | **VERIFIED-EXACT** (`q=1..6,8,10`; `K=3,4,5`) | doubled geometric orbits |
| J3-10 | the honest char-`p` invariants are therefore `d_min(P)` / `d_max(P)` = minimal / maximal degree of a Jacobian mate, and the *primitive* pair is `(n,d_min)`.  Search results are reported in these terms | **method** | -- |
| J3-11 | **the realisable-degree law.**  Over `F_p` the set of mate degrees of `P` is `{d_min(P)} u (<n,p> n (d_min,inf))`, `<n,p>` the numerical semigroup on `n=deg P` and `p`.  The inclusion `>=` is PROVED (`ker(L_P)` is a ring containing `F_p[x^p,y^p][P]`); equality is observed on every row of the search.  Over `C` this specialises to `{d_min} u nZ`, i.e. **Jung--van der Kulk's divisibility recovered as an output**, and it makes the target explicit: a counterexample is a `P` with `d_min(P) = m` and `n` not dividing `m` | **PROVED (one inclusion)** + **VERIFIED-EXACT** (16 search rows, 0 discrepancies) | see sec 4.4 |
| J3-11a | no `P` in any char-0 scan has `d_min > 1`: every char-0 hit is a **coordinate** with a linear partner.  `(2,3)`, `(2,5)`, `(3,4)`, `(3,5)`, `(4,6)` admit **no** Jacobian pair in the scanned boxes; the `(4,6)` char-0 scan is re-verified exactly over `Q` on all 54 hits | **VERIFIED-EXACT within the boxes** | boxes in sec 4.5 |
| J3-11b | a **primitive characteristic-3** Jacobian pair exists at `(n,m) = (4,5)`, `(a,b) = (4,5)` -- both `>= 2` and coprime: `P = y+y^4+xy^2`, `Q = 2y^2+y^5+2x+2xy^3+x^2y`, `Jac = 1` in `F_3[x,y]`, `d_min(P) = 5`, not an automorphism.  It is a char-`p` phenomenon and says **nothing** about JC(2) | **VERIFIED-EXACT** (and by hand) | `char 3` only |
| J3-11c | no Jacobian pair with `K >= 2` was produced by any search in this lane, in any characteristic, at any degree.  Lane J2's closing caveat is unchanged | **VERIFIED-EXACT within the windows** | see sec 4.7 |
| J3-12 | **(c) verdict: the metallic idea gives NO exclusion and NO prioritisation.**  What controls circuit alternation is (i) the sign pattern of the roots and (ii) the size of the root ratio against explicit **non-metallic** thresholds.  Searching by degree remains strictly better | **NEGATIVE, established** | -- |

---

## 1. (a) There are two different "metallic strata", and they answer oppositely

The task's premise -- *"metallic pairs are reciprocal-closed, hence they sit
inside lane J2's stratum; confirm and push"* -- is **half right, and the half
that is right is the half that carries no alternation.**

### 1.1 Root-metallic: inside the reciprocal-closed stratum, minimally alternating

`lam_q = (q+sqrt(q^2+4))/2` satisfies `lam^2 = q lam + 1`.  Its Galois conjugate
is `-1/lam`, **not** `1/lam`, so the `Q`-rational metallic pair is

    {lam_q, -1/lam_q},   sum = q,  product = -1,   H_q = x^2 - q x y - y^2.

Exactly verified for `q=1..5` in `Z[lam]`.  It is reciprocal-closed in
THM-3003's sense with `mu = -1` (`{r_i} = {mu/r_i}`), so **J3-1 confirms the
premise**.

But three facts empty it out:

1. **`g=2` has no circuit.**  `c_k` is indexed by `k=2..g-1`; at `g=2` that range
   is empty and at `g=3` it has one entry (no sign change possible).  So the
   circuit is silent about the metallic quadratic.
2. **`R_k = R_(g-k)` is vacuous at `g=2`** -- there is a single ratio `R_1`.
   "Reciprocal-closed" is not a restriction here; it is automatic.
3. **No moduli.**  `H_q` has distinct roots (discriminant `q^2+4>0`), and any two
   distinct points of `P^1` are `PGL_2(C)`-equivalent to `{0,infinity}`.  Since a
   linear change of `(x,y)` carries Jacobian pairs to Jacobian pairs (Jac scales
   by `det`), the metallic-quadratic stratum **is** the whole `g=2`, `K=2`
   stratum, i.e. `H = xy` after a coordinate change.  The parameter `q` is a
   coordinate artefact.

Raising the degree does not help.  For a positive geometric metallic orbit
`{1, lam, ..., lam^(g-1)}` the circuit is palindromic (as THM-3003 requires) and
has **exactly one sign change** for every `q=1..4` and `g=4..9` -- the word is
`+-`, `++--`, `+++---`, ... , a single hump.  That is the *minimum* a nonzero
circuit can have, and

> the same word is produced by **every** ratio: `T = 2, 3, 10, 1000, 101/100,
> 1/3, 7/5, 10^6` all give the identical word at each `g`.

The reason is exact and structural.  For roots `{T^j}`,

    e_k(1,T,...,T^(g-1)) = T^(k(k-1)/2) * gaussian_binom(g,k)_T     (VERIFIED-EXACT)

and `k(k-1)/2` is **quadratic in `k`**, so it is annihilated by the third
difference `-Delta^3(log h)` that defines the circuit.  The circuit of a
geometric orbit sees only `gauss(g,k)_T / C(g,k)`, and that turns out to be
unimodal for every `T`.  **Metallicity of the ratio is invisible to the circuit.**

### 1.2 Two structural degeneracies of the rational metallic root strata

* `H_q^a` -- which is exactly the leading form `P_n = c H^a` of a hypothetical
  Jacobian pair -- has roots of **opposite sign** (`lam > 0 > -1/lam`), so its
  coefficient sequence is not sign-definite and `h_k` can **vanish** (e.g.
  `q=1, a=3`, where `h_2 = h_4 = 0`).  The repo circuit needs `A_k != 0`
  (THM-3000 sec 1, J2-2), so the circuit is simply **undefined** there.
* The smallest `Q`-rational reciprocal-closed metallic root set is the Galois
  closure `{lam, 1/lam, -lam, -1/lam}`, i.e.
  `H = x^4 - (q^2+2) x^2 y^2 + y^4`.  It has `A_1 = A_3 = 0`: the circuit is
  **undefined** for every `q`.

So the "metallic quartic" -- the first place where the root-metallic idea could
have had circuit content -- is outside the circuit's domain of definition.

### 1.3 Coefficient-metallic: OUTSIDE the reciprocal-closed stratum, maximally alternating

THM-3010 sec 3's object is **not** a root condition.  It is the `h`-sequence
condition `h_k = a_k` with `a_k = q a_(k-1) + a_(k-2)`, whose Simson identity
`a_(k-1)a_(k+1) - a_k^2 = (-1)^k` forces `R_k = a_k^2/(a_k^2+(-1)^k)` to
alternate about `1` at every index.  Realised as a binary form
`H = sum_k C(g,k) a_(k+1) x^(g-k) y^k` (shifted so `h_0 != 0`):

| | |
|---|---|
| Simson prediction `R_k = a_(k+1)^2/(a_(k+1)^2+(-1)^(k+1))` | exact, `q=1..5`, `g=5..10` |
| sign word | `+-+-+-...`, **`g-3` changes = the maximum** |
| all `A_k > 0` | yes -- circuit well defined |
| squarefree, so `K = g` distinct roots | yes |
| real roots | `0` (even `g`) or `1` (odd `g`) -- **not real-rooted** |
| `R_k = R_(g-k)` palindromic | **NO**, for every `q` and every `g` |

The last row is decisive.  `R_k = a_(k+1)^2/(a_(k+1)^2+(-1)^(k+1))` tends to `1`
as `k` grows, so it cannot equal `R_(g-k)`; by THM-3003 sec 1 (extended over `C`
by J2-4) the root multiset is **not** reciprocal-closed.  Hence

> **J3-6.  The metallic object that attains maximal circuit alternation lies
> strictly OUTSIDE the reciprocal-closed stratum, while the metallic object that
> lies inside it has minimal alternation.  "Metallic" does not transport
> THM-3010's maximality into lane J2's stratum.**

This is the correct answer to (a), and it is the opposite of the expected one for
the reading that matters.  (A second useful by-product: the circuit is well
defined for these forms **without real-rootedness** -- positivity of the `A_k`
is all THM-3000's convention needs.)

---

## 2. Is "metallic" distinguished at all?  Exact thresholds say no

The family that is simultaneously reciprocal-closed *and* maximally alternating
is THM-3004 sec 3b(ii)'s **doubled geometric** family `prod_j (x - T^j y)^2`
(equal multiplicities on a geometric progression are log-symmetric, hence
palindromic).  THM-3004 proves attainment "for `T` large".  Exactly how large?

Computing `R_k(T)` as exact rational functions and taking exact real roots of the
numerators of `R_k - R_(k-1)`:

| `K` | `d=2K` | max alternation `d-3` | threshold `T*` | `deg` of `T*` over `Q` |
|---|---|---|---|---|
| 3 | 6 | 3 | `2.4931270138811760006...` | 10 |
| 4 | 8 | 5 | `3.6952124348324745171...` | 16 |
| 5 | 10 | 7 | `2.8975831738658444...` and `4.2852336357059...` | -- |

with (`K=3`) minimal polynomial

    8z^10 + 17z^9 + 10z^8 - 94z^7 - 226z^6 - 330z^5 - 226z^4 - 94z^3 + 10z^2 + 17z + 8

(a reciprocal polynomial, as the `T -> 1/T` symmetry demands).  Below the
threshold the word is the single hump `++--` (1 change); above it, `+-+-`
(3 changes = maximal).  For `K=4`: `+++---` below, `+-+-+-` above.

The metallic ratios are `lam_1 = 1.6180`, `lam_2 = 2.4142`, `lam_3 = 3.3028`,
`lam_4 = 4.2361`, `lam_5 = 5.1926`.  Therefore

* **no threshold is a metallic number** (checked exactly);
* at `K=3` the **silver** doubled orbit (`2.4142 < 2.4931`) is sub-maximal while
  the **bronze** one (`3.3028`) is maximal;
* at `K=4` **bronze** (`3.3028 < 3.6952`) is sub-maximal while `lam_4 = 4.2361`
  is maximal.

Metallic ratios straddle the thresholds with no pattern.  **The operative
quantity is the SIZE of the ratio (separation), exactly as THM-3004 sec 3a says;
metallicity is an irrelevant arithmetic decoration on that ray.**

### 2.1 The golden ratio is the WORST metallic ratio here, not the best

Exact alternation of the doubled metallic orbit `prod_(j<K) (x - lam_q^j y)^2`:

| `q` | `lam_q` | `K=3` (`max 3`) | `K=4` (`max 5`) | `K=5` (`max 7`) |
|---|---|---|---|---|
| 1 (golden) | 1.61803 | 1 | 1 | 1 |
| 2 (silver) | 2.41421 | 1 | 1 | 1 |
| 3 (bronze) | 3.30278 | **3 = max** | 1 | 3 |
| 4 | 4.23607 | **3 = max** | **5 = max** | 3 |
| 5 | 5.19258 | **3** | **5** | **7 = max** |
| 6, 8, 10 | 6.16, 8.12, 10.10 | **3** | **5** | **7** |

Because `lam_q` is increasing in `q` and maximality is governed by exceeding a
separation threshold, the **golden** ratio -- the smallest metallic number -- is
the one that *fails* hardest, at every `K`.  THM-3010 sec 3's remark that
"`phi` is the smallest metallic parameter, hence the extreme case of maximal
Newton-circuit alternation" is a statement about the **`h`-sequence** `a_k`,
where the norm `-1` does the work and the size of `lam` is irrelevant.  Read on
**roots**, the same ordering runs the other way.  The two faces of THM-3010
sec 3 -- "interleaved bands whose dominance alternates" and "a quadratic order
of norm `-1`" -- come apart here: only the second one survives when the object
is a binary form's root multiset.

---

## 3. A sharp scope finding for THM-3004, produced by the metallic form

The one place the metallic quadratic *does* produce something new is as a
counterexample to an over-read of THM-3004.

`P_n = c H_q^a` has `K=2` distinct roots `{lam_q, -1/lam_q}`, each of
multiplicity `a`, and degree `d=2a`.  THM-3004 sec 3b(i)'s cap for `K=2` is
`2K-3 = 1`.  Observed exact sign-change counts (`a = 2..8`):

| roots (each with multiplicity `a`) | `a=2` | 3 | 4 | 5 | 6 | 7 | 8 | cap `2K-3` | max `d-3` at `a=8` |
|---|---|---|---|---|---|---|---|---|---|
| `{lam_1,-1/lam_1}` (golden) | 1 | undef | 3 | 3 | **9** | 7 | undef | 1 | 13 |
| `{lam_2,-1/lam_2}` (silver) | 1 | 3 | 3 | **7** | 5 | 5 | 7 | 1 | 13 |
| `{lam_3,-1/lam_3}` (bronze) | 1 | 3 | **5** | 3 | 3 | 7 | 5 | 1 | 13 |
| `{lam_4,-1/lam_4}` | 1 | 3 | 3 | **7** | 5 | 3 | 7 | 1 | 13 |
| `{2,-1/2}` (not metallic) | 1 | 1 | **5** | 5 | 5 | 7 | **9** | 1 | 13 |
| `{3/2,-2/3}` (not metallic) | 1 | 3 | 3 | **7** | 5 | 7 | **9** | 1 | 13 |
| `{7/3,-3/7}` (not metallic) | 1 | 3 | 3 | **7** | 7 | 5 | 7 | 1 | 13 |
| `{10,-1/10}` (not metallic) | 1 | 3 | 3 | 3 | 3 | 3 | 3 | 1 | 13 |
| `{r,s}` both **positive**, any `r,s,a` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 13 |

(`undef` = some `h_k = 0`; boldface marks a count equal to the absolute maximum
`d-3` for that `a`.)  Two readings:

1. **THM-3004's `2K-3` law is a positivity statement and fails badly without it.**
   The last row shows the cap is exact and attained for same-sign roots (`1`, in
   every one of the `28` instances tested); the rows above show `K=2` forms
   reaching `9` sign changes against a cap of `1`, and repeatedly hitting `d-3`,
   the absolute maximum a length-`(d-2)` circuit admits -- as soon as the two
   roots have opposite signs.  This is consistent with (and sharpens) lane J2's
   scope flag 1 in sec 1.3, and is worth recording as a boundary note on
   THM-3004 sec 3b(i): the bound needs "all `A_k` of one sign", not merely
   "`K` distinct roots".
2. **Metallic is not distinguished among these either.**  The metallic counts
   (`1,3,5,7,9`) are interspersed with the counts of arbitrary rational
   reciprocal pairs `{r,-1/r}`; `{2,-1/2}` reaches `9` at `a=8` while the golden
   pair reaches `9` at `a=6`.  There is no metallic signature.

All rows are palindromic (reciprocal-closed), as THM-3003 sec 1 / J2-4 require.

---

## 4. (b) The search

### 4.1 Method

For a fixed degree pair `(n,m)` and field `k`, the map

    L_P : V_m -> V_(n+m-2),   Q |-> Jac(P,Q)

is **linear in `Q`**.  So the search enumerates `P` only, and answers by linear
algebra.  Define the **deficiency**

    delta(P) = min { total degree of a nonzero element of image(L_P) }.

`delta(P) = 0` iff a nonzero constant is in the image, i.e. iff `P` has a
Jacobian mate at all.  When `delta(P)=0` the mates form the coset
`Q_0 + ker(L_P)`, and we record

    d_min(P) = minimal degree of a mate,   d_max(P) = maximal degree of a mate.

`(n,m)` **admits a Jacobian pair** iff some `P` of degree exactly `n` has
`d_max(P) >= m` with a mate of degree exactly `m`; the pair is a *primitive*
`(n,d_min)` inflated inside `ker(L_P)` otherwise.  Every found pair is verified
by recomputing `Jac(P,Q)` exactly, and classified by an **automorphism test**
that is exact in both directions: `F=(P,Q)` is an automorphism of the affine
plane iff `k[P,Q] = k[x,y]`, i.e. iff `x` and `y` lie in the span of the
monomials `P^i Q^j`.  In dimension `2`, Jung--van der Kulk gives
`deg F^(-1) = deg F = D := max(n,m)`, so `x = G(P,Q)` with `deg G <= D` and every
monomial `P^i Q^j` that can occur has degree at most `D^2` in `(x,y)`.  Taking
the span over all `P^i Q^j` of degree `<= D^2` therefore neither misses an
automorphism nor admits a non-automorphism (`x, y in k[P,Q]` forces
`k[P,Q] = k[x,y]`).  An earlier draft used the bound `D+2`, which reports every
degree-`(n,m)` automorphism with `n,m >= 2` as a non-automorphism; the numbers in
sec 4.3 are from the corrected bound.

Reductions used, all rigorous because they preserve `Jac` up to a nonzero
constant: `P -> P + const`, `P -> alpha P`, `(x,y) -> A(x,y)` with `A` in `GL_2`
(applied to the leading form, giving `GL_2`-class representatives), and
`(x,y) -> (x+u,y+v)` (applied by keeping only lexicographically minimal
representatives of each translation orbit).

### 4.2 The characteristic-`p` obstruction, stated before the numbers

**J3-9 (PROVED).**  In characteristic `p`, every `R` in `F_p[x^p,y^p]` has
`R_x = R_y = 0`, so `Jac(P,R) = 0` for **every** `P`.  Moreover
`Jac(P,RS) = R Jac(P,S) + S Jac(P,R)`, so `ker(L_P)` is a **ring**; it contains
`F_p[x^p,y^p]` and `P`, hence `F_p[x^p,y^p][P]`.

Consequences, and they are severe:

* mate degrees over `F_p` are **unbounded**: if `(n,m_0)` admits a pair then so
  does `(n, m_0 + 5)` over `F_5` (add `x^5`), `(n, m_0+7)` over `F_7`, and so on;
* the explicit witness is `P = x^2+y`, `Q = x^5 - x` over `F_5`:
  `Jac(P,Q) = 1`, `deg Q = 5`, and `(P,Q)` is **not** an automorphism
  (Artin--Schreier), while `d_min(P) = 1` with the automorphism mate `Q = -x`.
  `(n,m) = (2,5)` has `(a,b) = (2,5)`, both `>= 2` and coprime -- **precisely the
  corner a char-0 counterexample would have to live in**;
* therefore **a char-`p` hit in the `a,b >= 2` corner is not evidence about
  JC(2)**.  It is the Jacobian conjecture being false in char `p`, which is
  classical.  This is exactly the confound the lane brief warned about, and the
  right response is to report `d_min`, not "(n,m) admits a pair".

### 4.3 Results over `F_p`

Exhaustive over all `GL_2(F_p)`-classes of leading forms and all
translation-canonical lower parts, with `P` of degree exactly `n` and the
constant term of `P` set to `0`.  `hist` = deficiency histogram; `REALISE m` =
number of `P` carrying a mate of degree exactly `m`; `auto/non-auto` is the
classification of those mates (sampled, capped at 40 per row).

| `(n,m)` | `p` | `(a,b)` | `g` | `P` scanned | deficiency hist | mates | `d_min` | `d_max` | REALISE `m` | auto / non-auto |
|---|---|---|---|---|---|---|---|---|---|---|
| (2,3) | 5 | (2,3) | 1 | 7 | `{0:4, 1:3}` | 4 | `[1]` | `[2]` | **0** | -- |
| (2,5) | 5 | (2,5) | 1 | 7 | `{0:4, 1:3}` | 4 | `[1]` | `[5]` | 4 | 0 / 4 |
| (2,7) | 5 | (2,7) | 1 | 7 | `{0:4, 1:3}` | 4 | `[1]` | `[7]` | 4 | 0 / 4 |
| (2,2) | 5 | (1,1) | 2 | 7 | `{0:4, 1:3}` | 4 | `[1]` | `[2]` | 4 | 4 / 0 |
| (2,4) | 5 | (1,2) | 2 | 7 | `{0:4, 1:3}` | 4 | `[1]` | `[4]` | 4 | 4 / 0 |
| (3,4) | 5 | (3,4) | 1 | 645 | `{0:20, 1:124, 2:501}` | 20 | `[1]` | `[3]` | **0** | -- |
| (3,5) | 5 | (3,5) | 1 | 645 | `{0:20, 1:124, 2:501}` | 20 | `[1]` | `[5]` | 20 | 0 / 20 |
| (3,3) | 5 | (1,1) | 3 | 645 | `{0:20, 1:120, 2:505}` | 20 | `[1]` | `[3]` | 20 | 20 / 0 |
| (3,6) | 5 | (1,2) | 3 | 645 | `{0:20, 1:124, 2:501}` | 20 | `[1]` | `[6]` | 20 | 8 / 12 |
| (2,3) | 7 | (2,3) | 1 | 9 | `{0:6, 1:3}` | 6 | `[1]` | `[2]` | **0** | -- |
| (2,5) | 7 | (2,5) | 1 | 9 | `{0:6, 1:3}` | 6 | `[1]` | `[4]` | **0** | -- |
| (2,7) | 7 | (2,7) | 1 | 9 | `{0:6, 1:3}` | 6 | `[1]` | `[7]` | 6 | 0 / 6 |
| (3,4) | 7 | (3,4) | 1 | 1757 | `{0:42, 1:336, 2:1379}` | 42 | `[1]` | `[3]` | **0** | -- |
| (3,5) | 7 | (3,5) | 1 | 1757 | `{0:42, 1:336, 2:1379}` | 42 | `[1]` | `[3]` | **0** | -- |
| **(4,6)** | **3** | **(2,3)** | **2** | 30726 | `{0:156, 1:1410, 2:7632, 3:21528}` | 156 | `[1,2,5]` | `[6]` | 156 | 0 / 40 |
| **(4,6)** | **5** | **(2,3)** | **2** | 234475 | `{0:116, 1:1320, 2:83004, 3:150035}` | 116 | `[1,2]` | `[5]` | **0** | -- |

The `(4,6)` row at `p=5` restricts leading forms to `P_4 = H^2`
(3 `GL_2(F_5)`-classes: `y^4`, `x^2y^2`, `x^4+4x^2y^2+4y^4`), which is L0 used as
an ansatz over `F_5`; the `p=3` row is **unrestricted** over all 12
`GL_2(F_3)`-classes of quartic leading forms.

### 4.4 One law explains every row

**Observed law (VERIFIED-EXACT over the whole window above).**  The realisable
mate degrees of a given `P` are

    { d_min(P) }  union  ( <n, p>  intersect  (d_min, infinity) ),

where `<n,p>` is the **numerical semigroup generated by `n = deg P` and `p`**.

Half of this is PROVED: `ker(L_P)` is a ring (from
`Jac(P,RS) = R Jac(P,S) + S Jac(P,R)`) containing `P` (degree `n`) and
`F_p[x^p,y^p]` (degrees `p Z_(>=0)`), so `F_p[x^p,y^p][P]` is inside it and
contributes exactly the degrees `<n,p>`.  The law asserts that in this window
nothing else is there.  It reproduces every single `REALISE` entry:

| row | `<n,p>` | is `m` in it? | REALISE |
|---|---|---|---|
| (2,3) `p=5` | `<2,5> = {0,2,4,5,6,7,...}` | `3` no | 0 |
| (2,5) `p=5` | same | `5` yes | 4 |
| (2,7) `p=5` | same | `7` yes | 4 |
| (3,4) `p=5` | `<3,5> = {0,3,5,6,8,...}` | `4` no | 0 |
| (3,5) `p=5` | same | `5` yes | 20 |
| (2,5) `p=7` | `<2,7> = {0,2,4,6,7,...}` | `5` no | 0 |
| (2,7) `p=7` | same | `7` yes | 6 |
| (3,4),(3,5) `p=7` | `<3,7> = {0,3,6,7,9,...}` | `4`,`5` no | 0, 0 |
| **(4,6) `p=5`** | `<4,5> = {0,4,5,8,9,10,...}` | **`6` no** | **0** |
| (4,6) `p=3` | `<4,3> = {0,3,4,6,7,...}` | `6` yes | 156 |

**Char-0 specialisation.**  Over `C` there is no Frobenius, so the same ring is
just `C[P]` and the law reads

    realisable mate degrees  =  { d_min(P) }  union  { n, 2n, 3n, ... }.

That is **Jung--van der Kulk's divisibility, recovered as an output of the linear
algebra**, and it makes the shape of a counterexample completely explicit:

> a char-0 counterexample at `(n,m)` with `a,b >= 2` is exactly a `P` of degree
> `n` with `d_min(P) = m` and `n` not dividing `m`.

Every `P` in the entire char-0 search (E1, E2, and both large-prime proxies) has
`d_min(P) = 1`: it is a **coordinate**, one half of an automorphism whose other
half is linear.  No `P` with `d_min > 1` and `a,b >= 2` was found in char 0.

### 4.5 Characteristic 0 (PART E)

Exact rational arithmetic, leading form `P_n = x^n` (WLOG over `Q` when `g=1`,
since a linear form has a rational root), lower coefficients in `{-1,0,1}`:

| `(n,m)` | `P` scanned | deficiency hist | mates | `d_min` | `d_max` | mate of degree `m`? |
|---|---|---|---|---|---|---|
| (2,3) | 9 | `{0:6, 1:3}` | 6 | `[1]` | `[2]` | **no** |
| (2,5) | 9 | `{0:6, 1:3}` | 6 | `[1]` | `[4]` | **no** |
| (3,4) | 243 | `{0:18, 1:216, 2:9}` | 18 | `[1]` | `[3]` | **no** |
| (3,5) | 243 | `{0:18, 1:216, 2:9}` | 18 | `[1]` | `[3]` | **no** |

For `(4,6)` (the `g=2`, `(a,b)=(2,3)` corner) the box `{-1,0,1}^9` on the lower
part with `P_4 in {x^4, x^2y^2}` (the two `PGL_2(C)`-classes of `H^2`) is scanned
at two large primes, `p = 1000003` and `p = 999983`, where Frobenius is absent
because `p >> deg`:

    39366 P scanned at each prime;  identical histograms {0:54, 1:672, 2:21114, 3:17526};
    54 hits at each; ALL 54 re-verified EXACTLY over Q (0 mod-p artefacts);
    the only (d_min,d_max) that occurs is (1,4)  ->  no mate of degree 6.

Consistent with the law: `6` is not in `<4> = {0,4,8,...}` and `d_min = 1 != 6`.

### 4.6 The characteristic-`p` hit, stated exactly

The `(4,6)` scan at `p=3` produces `d_min = 5`, i.e. a **primitive** char-3
Jacobian pair at `(n,m) = (4,5)`, `(a,b) = (4,5)` -- **both `>= 2` and coprime,
exactly the corner a char-0 counterexample would need**.  Verified exactly (and
by hand):

    P = y + y^4 + x y^2                                   deg 4
    Q = 2y^2 + y^5 + 2x + 2x y^3 + x^2 y                  deg 5
    Jac(P,Q) = 1   in  F_3[x,y]
    deficiency(P) at m = 1,2,3,4 is 2  ->  no mate of degree < 5, so d_min = 5
    (P,Q) is NOT an automorphism (x,y are not in F_3[P,Q])

This is **not** evidence about JC(2).  The Jacobian conjecture is false in
characteristic `p` -- `x - x^p` already breaks it -- and the object above is a
char-3 phenomenon.  It is recorded because it is the concrete demonstration that
**a positive search result in the `a,b >= 2` corner over `F_p` proves nothing
about `C`**, which is precisely the confound the lane brief flagged.

### 4.7 Which `(n,m)` admit a Jacobian pair -- summary answer

Over `C`, within the scanned boxes: **none** of `(2,3)`, `(2,5)`, `(3,4)`,
`(3,5)`, `(4,6)`.  Over `F_p`: `(n,m)` is admitted iff `m` lies in
`{d_min} u <n,p>`, and every admitted pair in the window is either an
automorphism or a Frobenius inflation of a `d_min = 1` coordinate -- with the
single exception of the primitive char-3 pair of sec 4.6.

**Observed `H`, `K` and circuit data for the `P` that do have mates.**  Every
`P` in the char-0 part of the search has `d_min = 1`, so its partner is linear,
`m = 1`, `g = gcd(n,1) = 1`, `H` is a linear form and `K = 1` -- lane J2's
automorphism stratum, with `R(H) == 1` and circuit `c(H) == 0` identically
(J2-2).  **No Jacobian pair with `K >= 2` was produced by any search in this
lane, in any characteristic, at any degree.**  There is therefore still no
object on which the leading-form circuit can be evaluated non-trivially; J2's
closing caveat ("no Jacobian pair with `K >= 2` exists to test the coupled law
on") is unchanged by this lane.

---

## 5. (c) Honest verdict: a clean negative

**The metallic / maximal-alternation idea produces no exclusion and no
prioritisation for JC(2), and searching by degree remains strictly better.**
The reasons, in decreasing order of decisiveness:

1. **Index arithmetic kills it at the bottom (J3-0).**  The leading-form circuit
   has `g-2` entries.  It is empty for `g <= 2`, has one entry at `g=3`, and can
   change sign only from `g >= 4`.  With `a,b >= 2` coprime, the first degree
   pair at which the circuit can carry *any* alternation information is
   `(n,m) = (8,12)`.  Everything the repo currently calls a "corner" with `g=2`
   -- including the golden corner `(178,288)` singled out in J2 sec 4.3, which
   has `g = 2` -- is **circuit-invisible**.
2. **At `g=2` the metallic form has no moduli (J3-2).**  `x^2-qxy-y^2` is
   `PGL_2(C)`-equivalent to `xy`.  Any statement about "the metallic stratum" at
   `g=2` is a statement about the whole `K=2` stratum, which lane J2 already
   showed cannot be excluded.
3. **The two metallic readings disagree (J3-6),** so the phrase does not name a
   single stratum.  Root-metallic is reciprocal-closed and *minimally*
   alternating; coefficient-metallic is maximally alternating and *not*
   reciprocal-closed.
4. **Where alternation really lives is separation, not arithmetic (J3-7).**  The
   thresholds are explicit algebraic numbers of degree 10 and 16 that are not
   metallic, and metallic ratios sit on both sides of them.
5. **The rational metallic root strata are outside the circuit's domain (sec 1.2).**

What the lane *did* produce that is worth keeping:

* **J3-0**, the index bound `g >= 4` for any circuit alternation content, which
  is a sharp and cheap scoping statement for lane J2's gate;
* **J3-8**, the scope correction that THM-3004's `2K-3` cap needs positivity, with
  explicit `K=2` witnesses reaching `9` and the absolute maximum `d-3`;
* **J3-9**, the Frobenius kernel lemma, which explains why char-`p` searches
  cannot address the `a,b >= 2` corner and gives the right invariants
  (`d_min`, `d_max`) for reporting them;
* **J3-11**, the realisable-degree law `{d_min} u <n,p>` (char `p`) /
  `{d_min} u nZ` (char `0`).  The char-0 form *is* Jung--van der Kulk's
  divisibility, produced as an output of pure linear algebra rather than assumed;
  it turns "does `(n,m)` admit a pair?" into the single scalar question
  "is there a `P` of degree `n` with `d_min(P) = m`?", which is what the search
  engine computes;
* **J3-7b**, the observation that the golden ratio is the *least* alternating
  metallic root ratio, which is the exact inversion of the intuition this lane
  was asked to push;
* a reusable exact search engine for Jacobian mates at fixed `(n,m)`, in any
  characteristic, reporting `(deficiency, d_min, d_max)` per `P`.

### 5.1 One prioritisation that survives, stated as a criterion (not a reduction)

In the minimal counterexample corner `(a,b) = (2,3)` we have `P_n = c H^2`: the
leading form of `P` is a **perfect square**, so every direction at infinity has
**even multiplicity**.  That is exactly the `m=2` row of THM-3004 sec 3b(ii)'s
attainment family.  Hence:

> **If** a counterexample has `(a,b) = (2,3)` and its `g` directions at infinity
> are real, of one sign, and separated above the threshold of sec 2, **then**
> `P_n` has the maximal circuit alternation `n-3`.

The direction of implication is one-way and the hypotheses are strong (they are
not known to hold, and by J3-0 they are empty unless `g >= 4`, i.e.
`(n,m) >= (8,12)`).  It is recorded as a computable criterion on a hypothetical
object, **not** as a reduction of JC(2), and **not** as an equivalence with the
Lame/Zaremba continued-fraction program on the exponent pair `(a,b)`.

### 5.2 Explicit non-result: the two golden ratios do not meet

The repo's PROBLEM-LEDGER section C names a "golden/worst-approximable degree
corner (Lame-for-polygons)".  That golden ratio is the limit of the **exponent**
ratio `a/b` along Fibonacci pairs -- a statement about the continued fraction of
`a/b`, i.e. about the depth of the L1 subtractive-Euclid tower.  The golden ratio
of this lane is a **root** `lam_1` of a leading form.  These live in different
coordinates (`(a,b)` versus the direction divisor `D` in `P^1`), and this lane
found **no map between them**.  Worse, J2 sec 4.3's lattice observation pushes
them apart: the deepest continued fractions in range come with the smallest `g`
(the golden corner `(178,288)` has `g=2`), and by J3-0 small `g` is exactly where
the circuit is empty.  Anyone tempted to fuse the two golden ratios should read
MISTAKE-237 first: this is the shape of error it records.

---

## 6. Boundaries -- what is NOT claimed

* No arrow from THM-3010, THM-3004 or any circuit statement **to** JC(2).
* J3-3/J3-4 are verified over the scanned `(q,g)` and ratio windows, plus the
  proved structural step (the `T^(k(k-1)/2)` factor is circuit-invisible); the
  claim "one sign change for **all** ratios and all `g`" is not proved.
* J3-5 is verified for `q=1..5`, `g=5..10`; the general-`g` statement is
  THM-3010's, which is about `h`-sequences, and is cited, not re-proved.
* J3-7's thresholds are exact for `K=3,4,5`; nothing is claimed for larger `K`.
* J3-8 is a **counterexample to an over-read** of THM-3004, not a refutation of
  THM-3004, whose statement is made under positivity.  No THM-3004 claim as
  written is contradicted.
* The searches are exhaustive **only** over the windows written in sec 4.3.  In
  particular the `p=5` scan at `(4,6)` restricts the leading form to
  `P_4 = H^2`, which is L0 -- a **char-0** theorem, used over `F_5` as an
  *ansatz*, flagged as such.  The unrestricted `(4,6)` scan is done at `p=3`.
  (L0's binary-form step `Jac(F,G)=0 => F=cH^a` divides by `n` and `m`, and its
  char-`p` analogue genuinely fails -- `Jac(x^p, G) = 0` for every `G` -- which
  is why the ansatz must be flagged rather than assumed.)
* J3-11's equality (as opposed to the proved inclusion `ker(L_P) >= F_p[x^p,y^p][P]`)
  is an **observation on 16 rows**, not a theorem.  A `P` with a larger kernel
  would break it, and such a `P` is exactly what a counterexample would be.
* The char-0 boxes are small (`{-1,0,1}` on the lower coefficients, one or two
  leading-form classes).  "No pair at `(4,6)` over `C`" is a statement about that
  box, not a proof; the corresponding theorem for `gcd(n,m) = 2` is classical and
  is cited, not re-proved here.
* J3-11b is char 3 only and is reported precisely so that it is **not** mistaken
  for progress on JC(2).
* Nothing here bears on the repo's degree-26 frontier, the Faber valuation
  obstruction (THM-2823), or any of the closed strata.

---

## 7. Reproduction

    python3 04-computation/jc2_metallic_search_laneJ3.py

Deterministic: the only randomness is a fixed-seed `RandomState(20260731)` used
for the D0 lemma's random instances.  Two independent full runs were diffed and
agree on **every** number; the only differing bytes are the wall-clock timings
printed at the end of each search row, so the `output_sha256` above certifies one
specific run rather than a timing-independent normal form.  Runtime `~6.5` min,
dominated by the `234475`-point `(4,6)` scan at `p=5` (`~4.5` min).  Output in
`05-knowledge/results/jc2_metallic_search_laneJ3.out`; `26/26` checks report `OK`.

Total exact work: `A`/`B`/`C` in `Z[lam]` and `Q` (no floating point in any
asserted quantity; `float` appears only in printed labels), `D`/`E` over `F_p`
with `numpy int64` and over `Q` with `Fraction`.
