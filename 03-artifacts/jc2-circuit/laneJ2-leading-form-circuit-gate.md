---
lane: J2
title: "The circuit of the leading form as a JC(2) gate: K=1 is forced for automorphisms, the P<->Q swap IS the repo reversal, and the antipalindromic locus is a PGL_2 involution"
status: mixed -- see the per-claim status table
script: 04-computation/jc2_leading_form_circuit_laneJ2.py
output: 05-knowledge/results/jc2_leading_form_circuit_laneJ2.out
script_sha256: 879aefbd15572f2cc36e07aa98514644667d55473180d3a98f3175a864c146b5
output_sha256: e7aabb3ac4e205f32727266b3c5942f22219fa8f2d822937f0617404207142d1
hash_basis: LF-normalized bytes
depends_on:
  - THM-3000 (circuit convention, sec 1)
  - THM-3003 (antipodal circuit rigidity: R_k = R_(d-k) <=> reciprocal-closed roots)
  - THM-3004 (sign-change bound 2K-3 for K distinct roots; two-cluster scope)
  - THM-1300 (the dim-3 JC counterexample)
  - 04-computation/jc2_degree_lattice_hyp9070.py (E1/E2/E3 = L0/L2/Euclid depth)
  - 04-computation/jc2_resonance_dictionary_leadingform_deathstar_S107.py
scope_guard: >
  MISTAKE-237.  NOTHING below is a reduction of JC(2) to a circuit problem, an
  equivalence between programs, or a bridge from NC2/GMC(2).  Every statement is
  either (i) a theorem about automorphisms, (ii) an exact identity about binary
  forms, or (iii) a STRATIFICATION of the hypothetical counterexample locus by a
  computable invariant, with the direction of implication stated explicitly.
---

# Lane J2 -- the circuit of the leading form as a JC(2) gate

Setting throughout.  `(P,Q)` is a **Jacobian pair** in `C[x,y]`: `Jac(P,Q)` is a
nonzero constant.  `n = deg P`, `m = deg Q`, `n+m >= 3`.  By **L0** (upstream
this session; independently the repo's `jc2_degree_lattice_hyp9070.py` E1) the
degree `n+m-2` part of the Jacobian gives `Jac(P_n,Q_m)=0`, hence

    P_n = c H^a,   Q_m = c' H^b,   g = deg H = gcd(n,m),  a = n/g,  b = m/g,
    gcd(a,b) = 1.

`H` is the **leading form**; its zeros in `P^1` are the **directions at
infinity** shared by the generic fibres of `P` and of `Q`.  Write

    D = sum_i m_i [p_i]   (the direction divisor),   K = #distinct p_i,
    sum_i m_i = g,        h(t) = H(1,t),

and give `H` the repo circuit (THM-3000 sec 1) via its coefficient sequence
`H = sum_k A_k x^(g-k) y^k`, `h_k = A_k/binom(g,k)`, `R_k = h_k^2/(h_(k-1)h_(k+1))`,
`c_k = log(R_k/R_(k-1))`.

---

## 0. Claim table

| # | claim | status | scope |
|---|---|---|---|
| J2-1 | every polynomial automorphism of `C^2` of degree `> 1` has `K = 1`, i.e. `H = c L^g` for a single linear form `L`; multiplicity vector `(g)` | **PROVED** (induction, using classical Jung--van der Kulk) + **VERIFIED-EXACT** (119 random automorphisms, 0 violations) | automorphisms only |
| J2-2 | for a binary form with all `A_k != 0`: `R_k = 1` for every `k` **iff** `K = 1`.  Hence an automorphism has `R(H) == 1` and leading-form circuit `c(H) == 0` | **PROVED** + **VERIFIED-EXACT** (383 forms; 87 automorphism leading forms) | any field, `A_k != 0` |
| J2-3 | contrapositive gate: `c(H) != 0` (equivalently `K >= 2`) **implies** `(P,Q)` is not an automorphism.  The converse is **NOT** proved here | **PROVED** (one direction) / **OPEN** (converse) | see sec 2.3 -- this is the honest form of the "gate" |
| J2-4 | THM-3003 sec 1 rigidity holds verbatim over `C` (no positivity, no real-rootedness), by replacing the log-affine step with "geometric sequence in `C*`" | **PROVED** + **VERIFIED-EXACT** (complex `mu`, complex roots; exhaustive rational converse at `g=4,5`) | all `e_k != 0` |
| J2-5 | the coordinate swap `sigma(x,y)=(y,x)`, i.e. the Jacobian-pair involution `tau(P,Q) = (Q o sigma, P o sigma)` which swaps `(n,m)` and `(a,b)`, acts on the leading-form circuit **exactly as THM-3001's reversal**: `R_k(H o sigma) = R_(g-k)(H)` | **PROVED** + **VERIFIED-EXACT** | all `A_k != 0` |
| J2-6 | "some coordinate frame makes `c(H)` antipalindromic" **iff** "`D` admits a nontrivial involutive `PGL_2(C)` symmetry".  Always true for `K <= 4` with reduced `D`; **never** true when `K >= 3` and the multiplicities `m_i` are pairwise distinct | **PROVED** + **VERIFIED-EXACT** (exact Mobius search; 110 SL_2(Z) controls; `K=3..6` census) | `D` = direction divisor of any binary form |
| J2-7 | (J2-LAW) at each direction: `e_i - f_i = m_i(a-b) + eps_i - delta_i`, where `e_i = ord_(L_i)P_(n-1)`, `f_i = ord_(L_i)Q_(m-1)`, and `delta_i = 0` unless `e_i g = m_i(n-1)`, `eps_i = 0` unless `f_i g = m_i(m-1)` | **PROVED** from L1 + an exact order lemma; **VERIFIED-EXACT** (81/81 pair instances; 1004/1004 lemma instances at `K = 2..6`) | `n,m >= 2`, both Jacobians in L1 nonzero |
| J2-8 | corollary of J2-7 with `n > m`: `e_i >= 1` for **every** `i`, i.e. `rad(H)` divides `P_(n-1)` -- every direction at infinity reappears in the sub-leading form of the higher-degree component | **PROVED** (same scope) | `n > m >= 2` |
| J2-9 | the only known JC counterexample (THM-1300, dim 3) has all three leading forms equal to `C[y]`-multiples of the single form `x^3 z`; its direction divisor is `x^3 z`, so `K = 2` with multiplicity vector `(3,1)` -- exactly one step above the automorphism stratum `K = 1` | **VERIFIED-EXACT** | dim 3, one map; NOT a statement about dim 2 |

---

## 1. (a) What is forced about `K` and the multiplicity vector

### 1.1 Theorem J2-1: automorphisms have `K = 1`.  PROVED.

**Statement.**  Let `F=(P,Q)` be a polynomial automorphism of `C^2` with
`deg F = max(n,m) >= 2`.  Then there is a linear form `L` and constants `c_1,c_2`
with `P_n = c_1 L^n` and `Q_m = c_2 L^m`.  Consequently `H = c L^g`, the
direction divisor is `g[p]` for a single point `p in P^1`, `K = 1`, and the
multiplicity vector is the one-part partition `(g)`.

**Proof.**  Induction on `N = n+m`.

*Base, and the case `min(n,m)=1`.*  If `n=1` then `H` has degree `g=gcd(1,m)=1`,
so `H` is itself linear; L0 gives `P_1 = c_1 H` and `Q_m = c_2 H^m`, and `L=H`
works.  Symmetrically for `m=1`.  (This covers `N=3`.)

*Step.*  Let `n,m >= 2`, wlog `n <= m`.  By Jung--van der Kulk (L2) `n | m`; put
`k = m/n >= 1`.  Then `g = n`, `a = 1`, `b = k`, so `P_n = c_1 H` and
`Q_m = c_2 H^k`, hence the leading form of `Q - lam P^k` cancels for
`lam = c_2/c_1^k`.  Set `Q' = Q - lam P^k`.  Then

* `(P,Q')` is again an automorphism -- it is `F` post-composed with the
  triangular map `(u,v) -> (u, v - lam u^k)`;
* `deg Q' =: m' < m`, and `m' >= 1` (if `Q'` were constant then
  `Jac(P,Q) = Jac(P, lam P^k) = 0`);
* `n + m' < N` and `max(n,m') >= n >= 2`.

The induction hypothesis applied to `(P,Q')` yields `L` with `P_n = c_1 L^n`;
and then `Q_m = lam (P_n)^k = lam c_1^k L^m`.  QED.

The two ingredients used are: **L0** (proved upstream from the top graded piece
of the Jacobian) and **Jung--van der Kulk's divisibility** (classical, cited).
The reduction step is checked mechanically in PART A of the script (degree
strictly drops, `Jac` preserved, every trial).

**Census (VERIFIED-EXACT).**  119 random automorphisms built as alternating
compositions of `SL_2(Z)`-affine maps and de Jonquieres maps `(x, y+p(x))` /
`(x+q(y), y)` with `deg p, deg q in {2,3}`, 22 distinct degree pairs up to
`(27,27)`:

    Jac(P_n,Q_m)=0 (L0)          119/119
    n|m or m|n (L2)              119/119
    P_n=cH^a, Q_m=c'H^b          119/119
    L1 identity (n+m>=4)         104/104
    K = 1                        119/119
    R(H) == 1 in general position 87/87
    observed multiplicity vectors: (1),(2),(3),(4),(6),(8),(9),(12),(18),(27)
    observed (a,b): (1,1),(1,2),(1,3),(2,1),(3,1)   -- always a=1 or b=1

**Caveat found and fixed during the run.**  L1 is the degree `n+m-3` graded
piece, so it is the statement "`= 0`" only when `n+m >= 4`.  For `(n,m)=(1,2)`
the piece `n+m-3=0` is the **constant Jacobian itself** and equals `1`, not `0`.
Fifteen census rows initially "failed L1" for exactly this reason.  Any use of
L1 must carry the hypothesis `n+m >= 4`.

### 1.2 Theorem J2-2: the circuit detects `K=1` exactly.  PROVED.

For a binary form `H` of degree `g >= 2` with all `A_k != 0`:

    R_k = 1 for every k = 1..g-1
      <=> h_k h_(k-1)^(-1) is constant in k, i.e. (h_k) is geometric
      <=> A_k = binom(g,k) r^k A_0
      <=> H = A_0 (x + r y)^g
      <=> K = 1.

(The argument is THM-3003 sec 1's second-difference step run **multiplicatively**
in `C*`, so no positivity and no logarithms are needed.)  Verified on 383 random
forms, `K` ranging over `1..7`, zero violations.

Combining with J2-1:

> **An automorphism of `C^2` of degree `> 1` has trivial leading-form Newton
> ratio sequence, `R(H) == 1`, hence identically vanishing circuit `c(H) == 0`.**

### 1.3 The gate, stated honestly

The correct logical form is a **one-way criterion**, not an equivalence:

    c(H) != 0   (equivalently K >= 2)   ==>   (P,Q) is NOT an automorphism.

The converse -- "`K = 1` forces an automorphism" -- is **OPEN** and is *not*
supplied by anything in this lane.  `K=1` says `H = cL^g`, i.e. after a linear
change of coordinates `P_n = c x^n` and `Q_m = c' x^m`: the projective closures
of the generic fibres of `P` and `Q` meet the line at infinity in the **single**
point `[0:1:0]`.  That is "one *point* at infinity", which is strictly weaker
than "one *place* at infinity" (a point may carry several branches).  The
classical Abhyankar--Moh--Suzuki circle of results lives here, and THM-1345 sec 5
already names the propagation of leading-form dependence as "the AMS-hard open
content".  So:

* **what the lane adds:** `K` and the multiplicity vector `(m_1,...,m_K)` are
  *computable, coordinate-free* invariants of a hypothetical counterexample, and
  the entire automorphism locus is confined to the single stratum `K=1`,
  `lambda = (g)`;
* **what the lane does not claim:** that `K >= 2` is *necessary* for a
  counterexample, nor any reduction of JC(2) to a circuit statement.

**Consequence, correctly scoped.**  *If* a counterexample has `K >= 2`, then by
THM-3004 sec 3b(i) its leading-form circuit has at most `2K-3` sign changes, and
by THM-3004 sec 3a the reversal sites localise at the partial sums of the
multiplicities `m_i` sorted by root modulus.  Two scope flags on that use:

1. THM-3004's sign-change count is a statement about **real** `c_k`, hence needs
   `h` to have real coefficients of one sign (e.g. all directions real and of one
   sign after a real coordinate choice).  For general complex `H` the substitute
   is the winding of `R_k/R_(k-1)` in `C*`; the `2K-3` bound is **not** asserted
   over `C` here.
2. `K >= 2` bounds `2K-3 <= 2g-3`, and `g = gcd(n,m) >= 2` for a counterexample
   (`a,b >= 2` coprime), so the bound is informative only when `K` is small
   relative to `g` -- i.e. exactly in the **high-multiplicity** regime.

### 1.4 What L1 forces on the multiplicity vector: the order law

See sec 3 (J2-7, J2-8).  The headline: the multiplicities `m_i` enter an exact
integer law that couples them to `a-b`, and for `n > m` every direction is forced
to divide `P_(n-1)`.

---

## 2. (b) THM-3003 applied: the reciprocal-closed stratum

### 2.1 The swap is the reversal.  PROVED + VERIFIED-EXACT.

Let `sigma(x,y) = (y,x)`.  Then `tau(P,Q) := (Q o sigma, P o sigma)` is again a
Jacobian pair (`det D sigma = -1` twice, so `Jac` is preserved), with

    (n,m) -> (m,n),   (a,b) -> (b,a),   H -> H o sigma.

Writing `H = prod_i (x + r_i y)`, `H o sigma = (prod_i r_i) prod_i (x + r_i^(-1) y)`:
the direction parameter `t = y/x` is inverted, the coefficient sequence `A_k` is
**reversed**, `h_k -> h_(g-k)`, and therefore

    R_k(H o sigma) = R_(g-k)(H)        (VERIFIED-EXACT, 60 forms, g=3..7)

which is precisely THM-3001's reversal involution `N -> N*`.  So the repo's
antipodal map on log-root space (THM-3003 sec 2) **is** the `P<->Q` swap of the
Jacobian problem, read at infinity.  This is a dictionary entry, not a reduction.

### 2.2 THM-3003 rigidity over `C`.  PROVED + VERIFIED-EXACT.

THM-3003 sec 1 is stated for positive coefficients.  It extends verbatim to
arbitrary complex roots with all `e_k != 0`, because the only analytic step --
"two sequences with equal second differences differ by an affine function" --
is the multiplicative statement "a sequence in `C*` with trivial second
multiplicative difference is geometric".  Hence

    R_k = R_(g-k) for all k   <=>   {r_i} = {mu/r_i} as multisets,  mu^g = e_g^2
                              <=>   D is invariant under t -> mu/t on P^1.

Verified: (i) 40 reciprocal-closed multisets with **complex** `mu` (including
`mu = (1+i)^2/2`) and complex roots, degrees 2--7, all palindromic; (ii)
exhaustive converse over a 9-element rational pool at `g=4,5`: 103 palindromic
multisets out of 1698, **every one** reciprocal-closed.

### 2.3 The geometric meaning, and whether the stratum can be excluded

`t -> mu/t` is an involution of `P^1` with fixed points `+-sqrt(mu)`, swapping
`0` and `infinity`.  **Every** non-identity involution of `PGL_2(C)` is conjugate
to such a map, and a linear change of `(x,y)` -- which preserves the Jacobian-pair
property up to a constant -- acts on directions by the **full** `PGL_2(C)`.
Therefore:

> **Theorem J2-6.**  There exists a linear coordinate frame in which the
> leading-form circuit of `(P,Q)` is antipalindromic **iff** the direction
> divisor `D` admits a nontrivial involutive `PGL_2(C)` symmetry (an involution
> of `P^1` preserving `D` **with multiplicities**).

The stratum is therefore *not* coordinate-invariant, but the criterion is.  Both
directions of the answer to "can it be excluded?" are now computable, and the
answer is a clean dichotomy:

| regime | involutive symmetry | conclusion |
|---|---|---|
| `K = 1` | yes (trivially; `c(H)=0` is both palindromic and antipalindromic) | automorphism stratum is **inside** the antipalindromic locus |
| `K = 2`, any `(m_1,m_2)` | **always** (the unique involution fixing both points) | cannot be excluded |
| `K = 3`, some two `m_i` equal | yes (the involution swapping the equal pair) | cannot be excluded |
| `K = 3`, `m_1,m_2,m_3` pairwise distinct | **never** | antipalindromic in **no** coordinate frame |
| `K = 4`, reduced (`m_i = 1`) | **always** (`t -> mu/t` pairing `0<->inf`, `1<->lambda`) | cannot be excluded |
| `K >= 5`, reduced | generically **no** (codimension `>= 1`) | a genuine closed condition |

Proof of the "never" row: an involution preserving `D` permutes the `p_i`
preserving multiplicities; with all `m_i` distinct that permutation is the
identity; a Mobius map fixing `>= 3` points is the identity, which is not a
nontrivial involution.  Exact control in PART D: `H` with `K=3`,
`lambda = (3,2,1)`, `g=6`; 103 random `SL_2(Z)` coordinate changes produced
**zero** palindromic circuits.  Mirror control: `roots {2, 1/2, 1,1,1}`
(`K=3`, `lambda=(3,1,1)`) gives the exactly palindromic
`R = (121/115, 23/22, 23/22, 121/115)`.

Census over reduced configurations (60 random rational configurations each):

    K=3:  60/60 symmetric      K=4:  60/60 symmetric
    K=5:   2/60 symmetric      K=6:   5/60 symmetric

matching the dimension count `dim(symmetric locus) ~ floor(K/2)-1` against
`dim(moduli) = K-3`.

**Verdict for (b).**  The reciprocal-closed stratum **cannot** be excluded --
it *contains* the whole automorphism locus and, more sharply, it is
**unavoidable** whenever `K <= 4` with reduced `D` or whenever `K = 2` (any
multiplicities).  Since `K = 2` is the first stratum above the automorphisms,
and since (J2-9) the only known JC counterexample in any dimension sits at
`K = 2`, "the counterexample lives in the reciprocal-closed stratum" is the
**default expectation**, not an exotic possibility.  The stratum becomes a real
constraint only for direction divisors with `K >= 3` and unequal multiplicities,
or `K >= 5` reduced -- and there it is a genuine, checkable exclusion.

This is the honest inversion of the naive hope: the antipalindromic locus is not
where a counterexample is trapped; it is where a counterexample is **cheap**,
and the informative statement is its complement.

---

## 3. The order law extracted from L1

**Order lemma (exact, any binary forms).**  Let `H = L^m U`, `F = L^e V` with
`L` linear and `L` coprime to `U,V`, `g = deg H`, `d = deg F`.  Then

    Jac(F,H) = L^(e+m-1) [ L Jac(V,U) + m U Jac(V,L) + e V Jac(L,U) ],

and restricting the bracket to the line `L = 0` (where forms restrict to
monomials in the line parameter) gives bracket `= u v t^(...) (e(g-m) - m(d-e))`.
Hence

    ord_L Jac(F,H) = e + m - 1 + extra,   extra >= 1  <=>  e * g = m * d.   (*)

Verified exactly: 1004/1004 instances on forms with `K = 2..6` distinct roots
and mixed multiplicities, plus 162/162 on automorphism leading forms.

**J2-LAW.**  Apply `(*)` to both terms of L1
(`Jac(P_n,Q_(m-1)) + Jac(P_(n-1),Q_m) = 0`, valid for `n+m >= 4`) after
substituting L0.  With `e_i = ord_(L_i) P_(n-1)`, `f_i = ord_(L_i) Q_(m-1)`:

    m_i(a-1) + [f_i + m_i - 1 + eps_i]  =  m_i(b-1) + [e_i + m_i - 1 + delta_i]

    ==>    e_i - f_i  =  m_i (a - b)  +  eps_i - delta_i        (J2-LAW)

with `delta_i = 0` unless `e_i g = m_i (n-1)` and `eps_i = 0` unless
`f_i g = m_i (m-1)`.  Scope: `n, m >= 2` and both Jacobians in L1 nonzero (the
degenerate branch `Jac(P_(n-1),H) = 0` means `P_(n-1)` is a power-form of a
common root form of `H`, and must be treated separately).

**J2-8 (corollary).**  Suppose `n > m >= 2`, so `a > b`.  If some `e_i = 0` then
`e_i g = 0 != m_i(n-1)`, so `delta_i = 0`, and J2-LAW gives
`-f_i = m_i(a-b) + eps_i >= m_i >= 1`, impossible since `f_i >= 0`.  Hence

> **every direction at infinity divides `P_(n-1)`: `rad(H) | P_(n-1)`.**

Quantitatively, in the non-resonant case `delta_i = eps_i = 0` the law reads
`e_i = m_i(a-b) + f_i`, so `L_i^(m_i(a-b))` divides `P_(n-1)` for every `i`, i.e.
`H^(a-b) | P_(n-1)` -- a **`(n-m)`-fold** vanishing of the sub-leading form along
the directions at infinity.

Verified exactly on automorphism data (`K=1`, so a single `i`): 81/81 instances
satisfy J2-LAW with the measured `delta,eps`; 81/81 non-resonant instances
satisfy `e - f = g(a-b)` on the nose.  Sample rows:

    (n,m)=(4,2) a,b=(2,1) g=2  e=2 f=0  delta=eps=0   e-f = 2 = m_1(a-b)
    (n,m)=(2,4) a,b=(1,2) g=2  e=0 f=2  delta=eps=0   e-f = -2 = m_1(a-b)
    (n,m)=(6,6) a,b=(1,1) g=6  e=4 f=4  delta=eps=0   e-f = 0

Note the law is exactly the *first* step of the L1 subtractive-Euclid tower of
the upstream L1 note: it converts the exponent-pair difference `a-b` into a
**vanishing order at infinity**, which is why the continued fraction of `a/b`
controls a depth.  That correspondence is an observation about the same
bookkeeping, **not** an equivalence of programs (MISTAKE-237).

---

## 4. (c) Observed `(K, alternation)` on the repo's actual objects

Greps: `04-computation/jc2_*` and `keller_*` were searched for explicit leading
forms.  The repo's degree-18/22/26 strata work (THM-2411/2778/2784/2796/2823)
runs in a *response/Faber* chart whose symbol `H` is a **different object**
(`H = V z^2 + B z + C_0`, a source-data polynomial) and carries no leading form
at infinity; those files supply no `(K, alternation)` data.  The two files that
do are below, together with THM-1300.

### 4.1 THM-1300, the only known JC counterexample (dim 3).  VERIFIED-EXACT.

With `u = 1+xy`,

    F1 = u^3 z + y^2 u (4+3xy),   deg 7,   top form   x^3 y^3 z
    F2 = y + 3x u^2 z + 3x y^2(4+3xy),  deg 6,  top form  3 x^3 y^2 z
    F3 = 2x - 3x^2 y - x^3 z,     deg 4,   top form   -x^3 z

    det JF = -2 (re-verified);  det of the top-form Jacobian = 0  (L0 analogue)
    gcd(top forms) = x^3 z,   cofactors  y^3,  3y^2,  -1

So all three leading forms are `C[y]`-multiples of the **single** form `x^3 z`:
the three-variable shadow of `P_n = cH^a`, `Q_m = c'H^b`.  The invariant being
read is the **gcd of the leading forms** -- in dim 2 that gcd is exactly
`H^min(a,b)`, so it carries the same distinct-component count `K` and the same
multiplicity ratios as `H` itself.  Here it is `x^3 z`, the divisor
`3[{x=0}] + 1[{z=0}]` in `P^2`:

> **`K = 2` distinct components, multiplicity vector `(3,1)`.**

Restricting to a generic line (the two-cluster profile `(t-p)^3(t-q)`), the
circuit sign word over a ratio scan is `--` or `++` for real positive ratios
(0 sign changes) and `-+` for `q/p = -2` (1 sign change); the observed maximum
is `1`, exactly THM-3004's cap `2K-3 = 1` and consistent with THM-3004 sec 4's
exhaustive two-cluster scope.

**Pattern.**  The single known counterexample in any dimension sits at
`K = 2` -- the *first* stratum above the automorphism stratum `K = 1` -- with an
**unequal** multiplicity vector `(3,1)`, and with maximal-for-its-`K` alternation
`1`.  Scope: this is one map, in dimension three; it is not a Jacobian pair in
two variables and it does not constrain JC(2).  It is recorded because it is the
only real data the invariant has ever been evaluated on.

### 4.2 The repo's homogeneous-Keller family.  VERIFIED-EXACT.

`jc2_resonance_dictionary_leadingform_deathstar_S107.py` proves: for `F = I + H`
with `H` homogeneous of degree `d >= 2`, Keller `<=>` `H = c (bx-ay)^d (a,b)`,
the shear.  Then `P_d = ca(bx-ay)^d`, `Q_d = cb(bx-ay)^d`, so
`g = d`, `(a,b)=(1,1)`, `H = (bx-ay)^d`:

    (a,b,c,d)=(1,2,1,3): K=1, mults [3]
    (a,b,c,d)=(2,-1,3,4): K=1, mults [4]
    (a,b,c,d)=(1,1,1,5): K=1, mults [5]

`K = 1` throughout -- as J2-1 requires, since shears are automorphisms.  The
homogeneous stratum therefore gives **no** `K >= 2` data.

### 4.3 The degree lattice

`jc2_degree_lattice_hyp9070.py` (E1--E3) plus J2-1 give the bracketing
`1 <= K <= g = gcd(n,m)`, with `K = 1` the automorphism stratum.  On the repo's
named degree corners:

    (n,m)=(178,288) = 2*(F_11,F_12): g=2,  (a,b)=(89,144),  cf = [1,1,...,1,2],  K in {1,2}
    (n,m)=(16,24):   g=8,  (a,b)=(2,3),  K in {1..8}
    (n,m)=(26,39):   g=13, (a,b)=(2,3),  K in {1..13}
    (n,m)=(64,96):   g=32, (a,b)=(2,3),  K in {1..32}
    (n,m)=(55,89):   g=1  -- excluded by gcd=1 (Magnus), and K is forced to 1

**Observation worth flagging (CONJECTURAL, no proof offered).**  The golden
corner `(178,288)` that the repo's `jc2_golden_corner_census...` singles out has
`g = 2`, hence `K <= 2` and `deg H = 2`: the deepest continued fraction in range
comes with the **shallowest** direction divisor the counterexample condition
allows.  The two "complexity" axes -- Euclid depth of `a/b` and the multiplicity
structure of `H` -- pull in opposite directions along `n+m = const`, since
`a*b = nm/g^2` grows as `g` shrinks.  Whether this trade-off is an artefact of
the census window or a real constraint is **open**; it is stated as an
observation about the lattice, not as a result about JC(2).

---

## 5. Boundaries -- what is NOT claimed

* No arrow from NC2, GMC(2) or any circuit theorem **to** JC(2), in either
  direction.  J2-3 is a one-way criterion on automorphisms; MISTAKE-237's
  quarantine is respected throughout.
* `K = 1 => automorphism` is **not** proved; it is the AMS-hard content named in
  THM-1345 sec 5, and this lane does not touch it.
* THM-3004's `2K-3` sign-change bound is a **real-sign** statement; it is applied
  here only in the scoped form of sec 1.3, and is **not** asserted for complex
  direction divisors.
* J2-7/J2-8 assume `n, m >= 2`, `n+m >= 4`, and that neither Jacobian in L1
  vanishes; the degenerate branch (`P_(n-1)` a power-form of `H`'s root form) is
  untreated.
* The J2-LAW census is on automorphism data only, hence `K = 1`; the order lemma
  behind it is verified independently at `K = 2..6`, but no Jacobian pair with
  `K >= 2` exists to test the *coupled* law on -- that is the whole point of the
  conjecture.
* sec 4.1 is a dim-3 datum.  It says nothing about dim 2.

---

## 6. Reproduction

    python3 04-computation/jc2_leading_form_circuit_laneJ2.py

The run is **deterministic** (a private `random.Random(20260731)` is used, because
sympy perturbs the global random stream -- two consecutive runs are byte-identical,
which is what the output hash above certifies).

Seven parts, all reporting `True`:
A automorphism census `K=1` (+ mechanical induction-step check),
B `R == 1 <=> K == 1`,
C THM-3003 over `C` and swap `=` reversal,
D involutive-symmetry criterion with the `SL_2(Z)` control and the `K=3..6`
census, E the order lemma and J2-LAW, F the THM-1300 / shear / degree-lattice
anatomy.
