# Posets, DAGs and the zeta(5) shape: the base-3/2 tree behind the Syracuse map, the branching of the residual DAG, and the triple Diophantine criticality of a Collatz orbit (Apéry exponent 1, Ridout exponent 2, Borel–Dwork product 1)

**Session:** opus, `collatz-posets-zeta5-20260927`, 2026-09-27 (S15 of the
Collatz thread; a second opus session, `collatz-poset-dag-20260927`, worked
the same directive in parallel and is cited, not repeated).
**Owner's directive:** "consider deeply posets and their directed acyclic
graphs as you push these final steps towards a creative Collatz proof
finishing move. consider recent proofs made public about zeta(5) being
irrational."
**Inherits (read, cited, not re-derived):**
[THM-4503](../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md)
(no-descent words of a cell = linear extensions of a width-2 poset),
[THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
(coefficient-descent cylinders, thresholds `N(w)`, exact modulus `2^(A+1)`),
[THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md)
(divergent orbits: `sum 1/m_l < infinity`, `C_inf` finite, the Bernstein
series at two places),
[THM-4495](../../01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md)
(`|W_k|` of order `2^(h* k) k^(-3/2)`),
the parallel session's note
[`collatz_posets_dags_20260927_spine_descent_tree.md`](collatz_posets_dags_20260927_spine_descent_tree.md)
and its provisional [THM-4514](../../01-canon/theorems/THM-4514-collatz-value-time-poset-spine-blocks-and-two-place-descent-tree.md)
(value-time poset, excursion forest, spine = orbit ∩ `E_inf`, tight spine
blocks, the descent tree with both places, and the typing of every DAG
finishing move as the transversality statement),
the S13 note [`collatz_directions_20260926.md`](collatz_directions_20260926.md)
(`E_inf`; Proposition 4: no prefix rank),
the S14 note [`collatz_shadow_error_flp_20260927.md`](collatz_shadow_error_flp_20260927.md)
(shadow error `eta_l`, the `3x-1` copy, `N = (m+1)/2` coordinates),
the HARD-class note [`collatz_procgen_20260922_hard_class.md`](collatz_procgen_20260922_hard_class.md)
(Theorems S, D, Y: 2-adic Liouville and Tschakaloff–Padé arguments that
exclude words as parity vectors of rationals),
the repo's audit of the p-adic zeta draft
[`.scratch/padic_repo_audit_20260825/REPORT.md`](../../.scratch/padic_repo_audit_20260825/REPORT.md),
and the Mahler core-paper sheet
[`CORE-PAPERS-MAHLER-THREE-HALVES.md`](../reference/CORE-PAPERS-MAHLER-THREE-HALVES.md)
(Akiyama–Frougny–Sakarovitch rational-base numeration).

**Status: PROVED elementary propositions (Propositions 1–10; the genuinely
new exact objects are the base-3/2 tree dictionary of Proposition 1, the
Apéry-form recurrence of Proposition 3, and the triple criticality of
Propositions 5, 6, 8 with the Borel–Dwork dichotomy of Proposition 10) +
FINITE-EXACT (identity checks on 2·10^6 odd integers, the AFS expansion to
10^5, the residual-tree branching to length 60, the exact criticality
identities on three orbits, the recurrence on 4·10^6 steps) + CITED (the
2025–2026 p-adic zeta(5) record, verified in the browser today) + DIRECTION
(section 4, marked) + AUDIT pending (section 6). Collatz OPEN. No finishing
move: the honest content is that each of the three Diophantine shapes behind
the zeta(5)-type proofs is exactly critical on a Collatz orbit, and the note
says by which factor.** Scripts
`04-computation/experiments/collatz_posets_dags_zeta5_20260927.py` and
`collatz_posets_dags_zeta5_20260927_criticality.py`, outputs `.out` beside
them.

Notation (Syracuse map): `U(m) = (3m+1)/2^v` on odd `m`; along the orbit
`m_0 = n, m_1, ...` the valuations are `v_l`, `d_l = v_1 + ... + v_l`,
`S_(L-1) = sum_(k<L) 3^(L-1-k) 2^(d_k)`, `C_L = prod_(j<L) (1 + 1/(3 m_j))`,
and the basic identity is `2^(d_L) m_L = 3^L n + S_(L-1) = 3^L n C_L`. The
valuation word of a 2-adic integer `x` is its sequence `(v_l(x))`; `E_inf`
is the set of `x in Z_2` with `3^l < 2^(d_l)` for no `l` (no coefficient
descent ever; S13). The shadow error of S14 is `eta_l = sum_(k>=0)
2^(d_(l+k) - d_l)/3^(k+1)` (finite on divergent orbits, THM-4476), with
`xi = n C_inf = n + eta_0` and `xi 3^l/2^(d_l) = m_l + eta_l`.

## 0. The answer in one paragraph

Posets first. The Syracuse step in the coordinate `N = (m+1)/2` is "append
the digit 0 to the base-3/2 expansion of `N`" when `N` is even (that is the
`v = 1` step, the left-child edge of the Akiyama–Frougny–Sakarovitch tree,
and also Mahler's even branch), and "apply the `3x-1` Syracuse map to `N`
and shift back" when `N` is odd (Proposition 1): the `3x+1` map on odd
integers is the `3x-1` map on the shifted integers, interleaved with the
digit-append edges of the base-3/2 tree, and the three negative cycles sit at
`N`-values `{0}`, `{-2,-3}`, `{-8,-12,-18,-27,-20,-30,-45}`. The residual
DAG (no-descent prefixes, THM-4503's width-2 posets) branches with the
exact next-letter law `P(odd) = |W_j|/|W_(j+1)| in [1/2, 1]`, equal to
`1/2` exactly at the steps where the coefficient threshold does not move
(Proposition 2): the DAG is asymptotically balanced, with a forced-odd
boundary of vanishing weight. Every rank-type finishing move on these DAGs
is the one the parallel session typed (its section 6), and the
label-precision state graph of the reset lane's construction is not acyclic
(section 2.3), so no rank factors through it. Then the zeta(5) shape. The
public record today is p-adic: `zeta_2(5)` is irrational (Calegari–Dimitrov–
Tang by arithmetic holonomy bounds, and Lai–Sprang–Zudilin by explicit
2-adic Apéry-type approximants with `0 < |zeta_2(5) - p_n/q_n|_2 <
max(|p_n|,|q_n|)^(-1-delta)` and `mu(zeta_2(5)) <= 20.342`, IMRN 2026);
real `zeta(5)` has no confirmed proof that today's search found. The three
mechanisms behind such proofs are (A) approximants with exponent above 1,
(B) a nonzero integer bounded below (Baker), (C) integrality plus
convergence radii at several places (Borel–Dwork, Pólya–Bertrandias,
holonomy bounds). On a Collatz orbit each is exactly critical: (A) the
orbit's own approximants `-S_(L-1)/3^L -> n` have 2-adic exponent `nu_L =
d_L/(L log_2 3)` with `nu_L - 1 = log_2(n C_L/m_L)/(L log_2 3)` exactly
(Proposition 4), so "exponent above 1" *is* coefficient descent and the
target is the rational `n`; (B) is the cycle half and is already the
classical Steiner–Simons–de Weger–Hercher route; (C) the series `F_n(z) =
sum 2^(d_k) z^k in Z[[z]]` has `R_oo R_2 = 2^(liminf d_k/k - limsup d_k/k)
<= 1`, the exact boundary of the Borel–Dwork criterion (Proposition 8), so
the periodicity conjecture for `n` is equivalent to a meromorphic
continuation of `F_n` past the circle `|z| = 2^(-liminf d_k/k)` (Proposition
10), and the real and 2-adic errors of the orbit's approximant multiply to
exactly `eta_L/3^L` (Proposition 5), the three-place product to `eta_L (n
C_L)^2 H_L^(-2)`, Roth–Ridout's exponent 2 with the shadow error as the
constant (Proposition 6). The factor by which each criterion misses is the
quantity the conjecture itself controls: the descent depth, the shadow
error, the oscillation of the mean valuation. That is the sharp form of
"the word does not know the integer", and it is why the finishing move is
not here.

## 1. What is new against the parallel note

The parallel session built the value-time poset, the excursion forest with
its spine, the tight spine blocks and the two-place descent tree, and typed
the finishing-move shapes (rank, well-quasi-order, balanced pairs, König,
self-consistency, `D`-coarsened chains). This note does not touch those
objects. It adds (i) a third tree, the base-3/2 numeration tree, in which
the `v = 1` steps are edges (section 2.1); (ii) the exact branching law of
the residual DAG (2.2); (iii) the observation that the reset lane's
label-precision state graph has cycles (2.3); and (iv) the whole of section
3, the zeta(5) shape, which the parallel note does not treat.

## 2. Three DAGs (PROVED, elementary)

### 2.1 The base-3/2 tree behind the Syracuse map

Write `N = (m+1)/2` for odd `m`, so `m = 2N - 1` and `3m + 1 = 2(3N - 1)`.

**Proposition 1 (N-coordinates and the AFS tree).**
(a) If `N` is even, then `v = 1` and the next `N` is `3N/2`. If `N` is odd,
then `v >= 2` and the next `N` is `(U^-(N) + 1)/2`, where `U^-(x) =
(3x-1)/2^(v_2(3x-1))` is the `3x-1` Syracuse map. In particular `v = 1` iff
`N` is even.
(b) Every positive integer has a unique finite least-significant-first
base-3/2 expansion `N = sum_i (a_i/2)(3/2)^i` with digits `a_i in {0,1,2}`,
obtained by `a = 2N mod 3`, `N -> (2N - a)/3 = floor(2N/3)`
(Akiyama–Frougny–Sakarovitch). The tree with parent `floor(2N/3)` has
children `3N/2, 3N/2 + 1` (N even) and `(3N+1)/2` (N odd): appending the
digit 0, 2, or 1.
(c) Hence a Syracuse step with `v = 1` appends the digit 0 to the base-3/2
expansion of `N` (the left-child edge of an even node), and Mahler's map
`g -> ceil(3g/2)` (the recursion of the integer parts `floor(xi (3/2)^n)` of
a Z-number) is the walk that appends the digit `g mod 2` at every step. The
`v >= 2` steps are not tree edges: for `v = 2` the step is `N -> (3N+1)/4`,
the digit-1 child halved once; for `v >= 3` it is `N -> (3N - 1 +
2^(v-1))/2^v`, which is not a child of `N` halved.
(d) The plus-sheet negative cycles `{-1}`, `{-5,-7}`, `{-17,-25,-37,-55,
-41,-61,-91}` sit at `N = {0}`, `{-2,-3}`, `{-8,-12,-18,-27,-20,-30,-45}`.

*Proof.* (a) `3m+1 = 2(3N-1)`. For `N` even `3N-1` is odd, so `v = 1`,
`m' = 3N-1`, `N' = 3N/2`. For `N` odd `3N-1` is even; with `2^(v') || 3N-1`
we get `v = v'+1 >= 2`, `m' = (3N-1)/2^(v') = U^-(N)` and `N' = (m'+1)/2`.
(b) is the cited theorem; the digit-append reading of the children is the
identity `sum_i (a_i/2)(3/2)^(i+1) + a/2 = (3N + a)/2` for a new lowest
digit `a`. (c) For a Z-number `xi (3/2)^n = g_n + f_n` with `0 <= f_n < 1/2`;
`g_(n+1) + f_(n+1) = 3g_n/2 + 3f_n/2`, and `f_(n+1) < 1/2` forces `g_(n+1) =
3g_n/2` (then `f_n < 1/3`) when `g_n` is even and `g_(n+1) = (3g_n+1)/2`
(then `f_n >= 1/3`) when `g_n` is odd: `g_(n+1) = ceil(3g_n/2)`, which is the
digit `g_n mod 2` appended. The `v = 2` formula: `N' = ((3N-1)/2 + 1)/2 =
(3N+1)/4`. (d) Direct. Checked for all odd `m < 2·10^6` (part 1 of the
script) and for `N <= 10^5` (part 2). ∎

What this buys: the `v = 1` runs of an orbit, which are the `N = 0 mod 2^k`
classes (`m = -1 mod 2^(k+1)`), are runs of appended zeros in base 3/2; a
Z-number's whole orbit is such a digit walk with digits in `{0,1}`; a
Collatz orbit is the same walk broken by halvings at the odd nodes. This is
the exact tree behind S14's "Mahler's `3N/2` branch", and it places the
Mahler frontier (THM-3848/4072/4074/4077/4082) and the Collatz thread in one
tree. It is a dictionary, not a mechanism: the halvings are exactly what the
tree does not see.

### 2.2 The branching law of the residual DAG

Let `W_j` be the no-descent T-words of length `j` (every prefix with `i`
odd and `j'` total letters satisfies `3^i > 2^(j')`; THM-4503, THM-4495).
The prefix DAG (a tree) of these words is the residual DAG of the
certification problem: its infinite paths are `E_inf`, its leaves are the
certified cylinders of THM-4512.

**Proposition 2.** Every `w in W_j` extends by the letter 1; it extends by
the letter 0 iff `3^(o(w)) > 2^(j+1)`. Hence, with `o_min(j) = min{o :
3^o > 2^j}`,

    |W_(j+1)| = |W_j| + #{w in W_j : 3^(o(w)) > 2^(j+1)},
    P_j := (fraction of W_(j+1) ending in 1) = |W_j|/|W_(j+1)| in [1/2, 1],

with `P_j = 1/2` iff `o_min(j+1) = o_min(j)` (the threshold does not move),
and `P_j > 1/2` iff it moves, in which case the words of `W_j` at the
minimal number of ones are forced odd.

*Proof.* Adding a 1 keeps every prefix inequality and gives `3^(o+1) >
3 · 2^j > 2^(j+1)`. Adding a 0 keeps the old prefixes and requires the new
one. The rest is counting. ∎

FINITE-EXACT (script part 3): `P_j` for `j = 1, 2, 3, 5, 8, 13, 21, 34, 55,
60` is `1, 1, 0.5, 0.75, 0.684, 0.616, 0.586, 0.544, 0.5, 0.5`; the range
over `j <= 60` is `[0.5, 1]`. So the residual DAG is asymptotically
balanced, its only bias being the forced-odd boundary whose weight
`#{o = o_min}/|W_j|` tends to zero (THM-4495's `k^(-3/2)` against the
boundary count). This is a branching statement about the DAG of words. It is
not the 1/3–2/3 statement (Linial's theorem for the width-2 cell posets),
which concerns a pair of positions, and the parallel note's section 6.3
already records why balanced pairs do not transfer to one integer.

### 2.3 Ranks, and the label-precision state graph

The parallel note's section 6.1 is the rank statement: a well-foundedness
proof is a function on spine points decreasing along the spine, S13's
Proposition 4 forbids it to depend on the word prefix, and the only known
such function is the stopping time itself. One addition. The reset lane's
construction (S13's directive: "regenerates arbitrarily deep precision
around the same negative cycle") is a walk on the state graph whose nodes
are pairs (shadowed negative cycle `w`, precision `K`) and whose edges are
completed excursions. In that graph the precision can increase along an
edge (an orbit re-enters a `2^(K')`-shadow of the same `x_w` with `K' >
K`), so the graph has directed cycles, and a rank cannot factor through it.
Read as a DAG statement this is the whole obstruction: the natural quotient
of the orbit by its shadow labels is not acyclic, and the acyclic object
(the descent tree `m -> D(m)`, parallel note section 4) is acyclic only
because it carries the integer.

### 2.4 The reachability order

The Syracuse functional graph, oriented `m -> U(m)`, is acyclic away from
the loop at 1 iff there is no other cycle, and well-founded iff there is no
divergent orbit; the Collatz conjecture says it is a tree with root 1. Its
finite truncations are exact (script part 5): all odd `m <= X` reach 1 for
`X = 10^3, ..., 10^6`, against the proved lower bound `X^0.84` for the
tree of 1 (Krasikov–Lagarias 2003). Nothing new here; the point of recording
it is that every linear extension of this order of type omega (an
enumeration in which `U(m)` precedes `m`) is a rank, and conversely, so
"find a topological sort" is the conjecture restated.

## 3. The zeta(5) shape (CITED record; PROVED criticalities)

### 3.1 The public record (browser, 2026-09-27)

* **Lai, Sprang, Zudilin, "A note on the irrationality of zeta_2(5)"**,
  arXiv:2505.05005 (v2, 26 May 2026), Int. Math. Res. Notices 2026, no. 16,
  rnag180. Abstract: in the spirit of Apéry, a sequence `p_n/q_n` of
  rational approximations to the 2-adic zeta value with `0 < |zeta_2(5) -
  p_n/q_n|_2 < max{|p_n|,|q_n|}^(-1-delta)` for an explicit `delta > 0`,
  giving a new proof of the irrationality of `zeta_2(5)`, "the result
  established recently by Calegari, Dimitrov and Tang using a different
  method", and the measure `mu(zeta_2(5)) <= 16 log 2/(8 log 2 - 5) =
  20.342...`.
* **Calegari–Dimitrov–Tang**: arithmetic holonomy bounds (the linear
  independence of `1, zeta(2), L(2, chi_-3)`, arXiv:2408.15403; "Arithmetic
  holonomy bounds and effective Diophantine approximation",
  arXiv:2510.04156), with the `zeta_2(5)` irrationality as reported in the
  abstract above. The repo's own p-adic audit cites both.
* **Lai–Sprang, "Many p-adic odd zeta values are irrational"** (arXiv,
  2024–2025) and **Zudilin 2018** ("one of the odd zeta values from
  `zeta(5)` to `zeta(25)` is irrational, by elementary means"): the
  "one-of / many-of" linear-form statements.
* **An unrefereed draft** (C. Long, GitHub `octonion/p-adic-zeta-irrationality`,
  August 2026) claims 22 Kubota–Leopoldt irrationalities including
  `zeta_2(5)`, `zeta_3(5)`, `zeta_5(5)` by a "hybrid arithmetic holonomy"
  method; the repo's audit of 2026-08-25 classifies it AUTHOR-CLAIMED /
  UNREFEREED and this note keeps that classification.
* **Real `zeta(5)`**: today's arXiv search found no confirmed proof of its
  irrationality; a 2024 claim (Suman) is marked "found incorrect" on
  arXiv. So "recent proofs made public about zeta(5) being irrational" means,
  as of today, the 2-adic (and other p-adic) `zeta(5)`, refereed for `p = 2`.

Three mechanisms carry these proofs. **(A) Approximants with exponent above
1.** If `alpha in Z_p` is rational, `|q alpha - p|_p >= c/max(|p|,|q|)` for
every nonzero linear form with integer `p, q`; so infinitely many nonzero
forms with `|q alpha - p|_p <= H^(-1-delta)` force irrationality. This is
Apéry's shape, and LSZ's theorem is exactly its 2-adic instance. **(B) A
nonzero integer bounded below.** The hardest step of an Apéry proof is
nonvanishing; Baker's theory supplies lower bounds `|a log 2 - b log 3| >
b^(-C)`, hence `|2^a - 3^b| > 3^b b^(-C)`. **(C) Integrality plus radii.**
Borel–Dwork: a series in `Z[[z]]` meromorphic in a complex disc of radius
`R_oo` and in a `p`-adic disc of radius `R_p` with `R_oo R_p > 1` is
rational; Pólya–Bertrandias replaces discs by capacity; the
Calegari–Dimitrov–Tang holonomy bounds are dimension bounds for spaces of
*holonomic* functions under such arithmetic constraints, and the p-adic
zeta proofs run by assuming rationality, building a holonomic function with
too-good arithmetic, and contradicting the bound.

The repo already uses (A) in the correct direction: Theorems S, D, Y of the
HARD-class note are 2-adic Liouville and Tschakaloff–Padé arguments (the
latter a 2-adic transcription of a Zudilin construction) proving that the
2-adic point `Phi(w)` of a given word `w` is *irrational*, hence that `w`
is the parity vector of no rational. What follows is what (A), (B), (C) look
like when the target is one integer `n` rather than one word.

### 3.2 The Apéry form of an orbit and the 2-adic exponent

**Proposition 3 (the Apéry-form recurrence).** Put `A_L = 3^L n + S_(L-1)`.
Then `A_0 = n`, `A_(L+1) = 3 A_L + 2^(v_2(A_L))`, `v_2(A_L) = d_L` and the
odd part of `A_L` is `m_L`; equivalently `A_L = 2^(d_L) m_L = 3^L n C_L`.

*Proof.* `A_(L+1) = 3^(L+1) n + 3 S_(L-1) + 2^(d_L) = 3 A_L + 2^(d_L)`, and
`3 · 2^(d_L) m_L + 2^(d_L) = 2^(d_L)(3 m_L + 1) = 2^(d_L + v_(L+1)) m_(L+1)`.
Induction. Checked on every orbit of odd `n < 2·10^5` (4,041,077 steps) and
against `3^L n C_L` exactly for `n < 2000` (criticality script, part 7). ∎

So the whole Syracuse dynamics is one integer recurrence, "multiply by 3 and
add the 2-part", and the linear forms `A_L = 3^L · n + S_(L-1)` in the
integer `n`, with integer coefficients `(3^L, S_(L-1))`, are the orbit's
Apéry forms: their 2-adic size is `|A_L|_2 = 2^(-d_L)`. In LSZ the forms
`q_n zeta_2(5) - p_n` come from a linear recurrence with polynomial
coefficients whose 2-adic valuations are computed by explicit
hypergeometric formulas; here the coefficient recurrence `S_L = 3 S_(L-1) +
2^(v_2(3^L n + S_(L-1)))` feeds the valuation of the form back into the next
coefficient. That feedback is the difference between the two problems.

**Proposition 4 (the 2-adic exponent is the descent ratio).** The
approximants `-S_(L-1)/3^L` converge 2-adically to `n` with `|n +
S_(L-1)/3^L|_2 = 2^(-d_L)`. With respect to the denominator `3^L` define
`nu_L = d_L/(L log_2 3)`. Then

    nu_L - 1 = log_2(n C_L/m_L)/(L log_2 3)   exactly,

so `nu_L > 1` iff `m_L < n C_L` iff the class of the first `L` valuations
has coefficient descent `3^L < 2^(d_L)` (THM-4512's notion); `E_inf` is the
set of 2-adic points whose own approximants never reach exponent 1; and
the divergence half of the conjecture says every positive integer's
approximants reach exponent above 1 at some `L`.

*Proof.* `d_L = log_2 A_L - log_2 m_L = L log_2 3 + log_2(n C_L) - log_2
m_L`. ∎

Two consequences fix the shape. First, the target `n` is rational, so
exponent above 1 is possible only inside the window `nu_L <= 1 + log_2(n
C_L)/(L log_2 3)`, of width `log_2(n C_L)` bits in `d_L`; on divergent
orbits `C_L` is bounded (THM-4476) and the window is `log_2 n + O(1)` bits
wide for ever. An irrationality proof needs exponent above 1 infinitely
often for an irrational target; the Collatz statement needs exponent above
1 *once*, for a rational target, and the approximants are not chosen but
dictated. The direction is inverted, and the quantity that must be
controlled (the valuations `v_2(A_L)` of a feedback recurrence) is exactly
the unknown. Numerically (script part 4 and the criticality script): along
the orbit of 27, `nu_L = 0.757, 0.883, 0.883, 0.883, 1.041` at `L = 5, 10,
20, 30, 40`; for 703, `nu_61 = 1.076`; for 6171, `nu_40 = 1.025`. The
exponent crosses 1 exactly where the orbit falls below `n C_L`.

### 3.3 The real side, the product identity and Ridout's exponent

The orbit's approximant is one rational number seen at two places:
`r_L := 2^(d_L) m_L/3^L = n + S_(L-1)/3^L` is 2-adically small (`|r_L|_2 =
2^(-d_L)`) and, on a divergent orbit, real-close to `xi = n C_inf`.

**Proposition 5 (product identity).** On a divergent orbit, for every `L`,

    |xi - r_L|_oo = eta_L 2^(d_L)/3^L,     |r_L|_2 = 2^(-d_L),
    |xi - r_L|_oo · |r_L|_2 = eta_L / 3^L,

and with `mu_L = -log_3 |xi - r_L|/L` (the real exponent against the
denominator `3^L`): `mu_L + nu_L = 1 - log_3(eta_L)/L <= 1`, with equality
iff `eta_L = 1` (never on the plus sheet, S14).

*Proof.* `xi 3^L/2^(d_L) = m_L + eta_L` (S14), so `xi - r_L = eta_L
2^(d_L)/3^L`; multiply by `|r_L|_2`. ∎

**Proposition 6 (Ridout criticality).** Let `H_L = max(2^(d_L) m_L, 3^L) =
2^(d_L) m_L = 3^L n C_L` be the height of `r_L`. Then `3` does not divide
`m_L` for `L >= 1`, `|1/r_L|_3 = 3^(-L)`, and

    |xi - r_L|_oo · |r_L|_2 · |1/r_L|_3 = eta_L/9^L = eta_L (n C_L)^2 · H_L^(-2).

So the three-place product of the local distances of `r_L` to the point
`(xi, 0, oo) in R × Q_2 × Q_3` is `H_L^(-2)` up to the factor `eta_L (n
C_L)^2 >= 1`: the exponent is exactly 2, the Roth–Ridout exponent. Ridout's
theorem (the S-adic Roth theorem) gives finitely many solutions of `product
< H^(-2-eps)` when the targets are algebraic; the orbit supplies exponent 2
with a constant at least 1, for a `xi` not known to be algebraic. The orbit
is Ridout-critical, not Ridout-violating, and the slack is the shadow error
times the square of the drift.

*Proof.* `m_L = (3 m_(L-1) + 1)/2^v = 2^(-v) mod 3`; the product is
Proposition 5 times `3^(-L)`; `H_L^2 = 9^L (n C_L)^2`. ∎

FINITE-EXACT (criticality script, part 6): all identities of Propositions
5–6 hold exactly (rational arithmetic) for `L = 1, ..., K-1` on the
truncated words of 27 (`K = 41`), 703 (`K = 62`) and 6171 (`K = 96`), where
`xi = n + eta_0` is the truncated shadow; the table prints `nu_L`, `mu_L`,
`eta_L`, `n C_L` and `eta_L (n C_L)^2`. For 27: `(L, m_L, nu_L, mu_L, eta_L)
= (10, 103, 0.883, -0.121, 13.66)`, `(30, 1367, 0.883, -0.035, 148.3)`,
`(40, 5, 1.041, -0.016, 0.333)`. The real exponents are negative: the
orbit's approximants do not even approach `xi` at the Dirichlet rate, since
`|xi - r_L| = xi eta_L/m_L` is the shadow error divided by the current
value.

### 3.4 The generating function, Borel–Dwork, and the periodicity conjecture

For `x in Z_2` with valuation word `(v_k)` put `F_x(z) = sum_(k>=0) 2^(d_k)
z^k in Z[[z]]`.

**Proposition 7 (PC in generating-function form).** The valuation word of
`x` is eventually periodic iff `F_x` is a rational function. Hence
Lagarias's periodicity conjecture for positive integers (no divergent
orbit) says: `F_n` is rational for every positive integer `n`.

*Proof.* If `v_(k+p) = v_k` for `k >= k_0` then `2^(d_(k+p)) = 2^A 2^(d_k)`
for `k >= k_0` with `A` the period sum, so `(1 - 2^A z^p) F_x` is a
polynomial. Conversely, if `F_x` is rational, its coefficients are all
powers of 2, and Pólya's theorem (1921: a rational series whose
coefficients have their prime factors in a finite set is, on each residue
class of some modulus `N` and from some index on, a geometric progression)
gives `d_(k+N) - d_k` constant for large `k`, i.e. an eventually periodic
word. ∎

**Proposition 8 (Borel–Dwork criticality).** `F_x` has complex radius `R_oo
= 2^(-limsup d_k/k)` and 2-adic radius `R_2 = 2^(liminf d_k/k)`, so

    R_oo · R_2 = 2^(liminf d_k/k - limsup d_k/k) <= 1,

with equality iff `d_k/k` converges to a finite limit (the mean valuation). The Borel–Dwork
criterion (`R_oo R_2 > 1` forces rationality) therefore never applies to a
Collatz series, and its boundary `R_oo R_2 = 1` contains irrational Collatz
series: the Sturmian valuation word `d_k = ceil(k log_2 3)` (letters 1, 2;
mean `log_2 3`; in `E_inf`) has `R_oo = 1/3`, `R_2 = 3`, and its `F` is not
rational (Proposition 7), while Theorem S says its 2-adic point is
irrational as well.

*Proof.* `|2^(d_k)|_oo^(1/k) = 2^(d_k/k)` and `|2^(d_k)|_2^(1/k) =
2^(-d_k/k)`; Cauchy–Hadamard at both places. ∎

**Proposition 9 (two values of one series).** `F_x(1/3)` converges in `Q_2`
for every `x` (`|1/3|_2 = 1 < R_2`) and equals `-3x`. For a positive
integer `n` with divergent orbit `F_n(1/3)` also converges in `R` (THM-4476)
and equals `3 n (C_inf - 1) = 3 eta_0 > 0`; for an orbit reaching 1 the real
series diverges (`d_k ~ 2k`).

*Proof.* `n = -S_(L-1)/3^L mod 2^(d_L)` (Proposition 3) and `S_(L-1)/3^L =
sum_(k<L) 2^(d_k)/3^(k+1) = F_n(1/3)/3` truncated; in `R`, `eta_0 = sum
2^(d_k)/3^(k+1)`. ∎

**Proposition 10 (the Borel–Dwork dichotomy for one integer).** Let `n` be a
positive integer with divergent orbit. Then `F_n` is holomorphic in `|z| <
2^(-limsup d_k/k)`, a disc of radius at least `1/3`, and is *not*
meromorphic in any disc `|z| < rho` with `rho > 2^(-liminf d_k/k)`. In
particular, if the mean valuation `v̄ = lim d_k/k` exists, the circle of
convergence `|z| = 2^(-v̄)` admits no meromorphic continuation of `F_n` to
any larger disc. Equivalently: for any positive integer `n`, a meromorphic
continuation of `F_n` to a disc of radius above `2^(-liminf d_k/k)` implies
that the orbit of `n` is eventually periodic.

*Proof.* `d_k <= k log_2 3 + log_2(n C_k)` with `C_k` bounded on a divergent
orbit (THM-4476) gives `limsup d_k/k <= log_2 3` and `R_oo >= 1/3`. If `F_n`
were meromorphic in `|z| < rho` with `rho · R_2 > 1`, Borel–Dwork would make
it rational and Proposition 7 the orbit eventually periodic, which a
divergent orbit is not. ∎

This is the shape (C) statement for Collatz, and it is exact: the
periodicity conjecture for `n` is a continuation statement about one
integer power series whose coefficients are the 2-parts of the orbit's
Apéry forms. It is also where the shape stops. The holonomy-bound method
needs a *holonomic* carrier; `F_n` is holonomic when it is rational and
(expectedly, section 4) only then, so the carrier is holonomic exactly when
there is nothing to prove. And the boundary `R_oo R_2 = 1` is populated by
irrational series (Proposition 8), so no capacity refinement that ignores
the integer can separate `n` from the Sturmian point. The one input that
distinguishes the integer, the finiteness of its binary expansion (S13
section 2), enters `F_n` nowhere: `F_n` is a function of the word alone.

### 3.5 The cycle half is already shape (B)

For a cycle of `L` odd steps with `d = d_L`, `n(2^d - 3^L) = S_(L-1) > 0`,
so `2^d > 3^L` and the nonzero integer `2^d - 3^L` must divide `S_(L-1)`;
Baker-type lower bounds `|2^d - 3^L| > 3^L L^(-C)` (equivalently the
irrationality measure of `log_2 3`) bound the number of "circuits" a cycle
can have. This is Steiner's theorem (no 1-cycles, 1977), Simons–de Weger
(no `m`-cycles for `m <= 68`, 2005) and Hercher ("There are no Collatz
`m`-cycles with `m <= 91`", J. Integer Seq. 26 (2023), Art. 23.3.5,
arXiv:2201.00406; a Lean formalization project's README, seen in today's
search, reports a machine-verified `m <= 82` and an error it found in the
paper's Theorem 21, which this note records as reported, not verified). So
the "nonzero integer forced small" half of the zeta(5) shape has been in
Collatz for fifty years, on the cycle half only; the divergence half has no
nonzero integer to bound, because the real quantity that would play its
role is the shadow error, which is not small (HYP-9164 says it is
unbounded).

### 3.6 Verdict

The three shapes are each exactly critical on a Collatz orbit:

| shape | the Collatz object | critical value | the slack |
|---|---|---|---|
| (A) 2-adic Apéry exponent | `nu_L = d_L/(L log_2 3)` of `-S_(L-1)/3^L -> n` | 1 | `log_2(n C_L/m_L)`: the descent depth |
| (C') Roth–Ridout, three places | `r_L = 2^(d_L) m_L/3^L` against `(xi, 0, oo)` | exponent 2 | `eta_L (n C_L)^2`: the shadow error |
| (C) Borel–Dwork radii | `F_n(z) = sum 2^(d_k) z^k` | `R_oo R_2 = 1` | `2^(liminf - limsup) d_k/k`: the oscillation of the mean valuation |

and the sum rule `mu_L + nu_L = 1 - log_3(eta_L)/L` ties (A) to the real
side. Every slack is a quantity the conjecture controls and nothing else
does. An Apéry-shaped Collatz proof would have to bound from below the
cumulative valuation of the feedback recurrence `A_(L+1) = 3 A_L +
2^(v_2(A_L))` for a single positive start (its negative starts `-1, -5,
-17` have `nu_L < 1` for ever, so positivity must enter); a holonomy-shaped
one would have to find a holonomic carrier other than `F_n`. Neither is a
finishing move, and this note claims none.

## 4. Directions (DIRECTION; none proved)

* **D1. A holonomic carrier from the functional equation.** The only
  Collatz object with a functional equation over `Q(z)` is Berg–Meinardus's
  generating function of the tree of 1 (synthesis row 22, typed only). Its
  equation is of Mahler type (`z -> z^3` with cube roots of unity), and
  Mahler functions have their own arithmetic rigidity theory
  (Mahler's method; Adamczewski–Bell: an algebraic Mahler function is
  rational). Whether the tree-of-1 series is a genuine Mahler function, and
  whether the Pólya–Carlson dichotomy (integer coefficients, radius 1:
  rational or natural boundary) can be combined with the functional
  equation, is the one place where shape (C) has a carrier with an
  equation. Cheapest test: write the Berg–Meinardus equation for the
  Syracuse tree and check whether it is linear over `Q(z)` in `h(z), h(z^3),
  h(z^9)`.
* **D2. `F_n` holonomic implies rational.** Expected: a P-recursive sequence
  whose terms are powers of 2 with `d_k = O(k)` has asymptotics `rho^k
  k^alpha (log k)^beta` times periodic factors (Birkhoff–Trjitzinsky), and
  the ratio `2^(v_(k+1))` lies in a discrete set, so the valuation word is
  eventually periodic. A proof would make precise that the holonomy-bound
  method has no carrier in `F_n`. Not attempted.
* **D3. Digit-append runs and Z-numbers.** Proposition 1 puts Z-numbers
  (digit walks in `{0,1}`) and Collatz orbits (digit walks broken by
  halvings) in one tree. The Mahler frontier's drift control `5/2` and
  Dubickas–Mossinghoff's bounds concern the unbroken walk; the question
  whether a broken walk can stay in `E_inf` is the conjecture; the
  sub-question "can the unbroken walk from an integer `N` stay in `{0,1}`
  digits for ever" is Mahler's. The dictionary suggests transporting the
  repo's Z-number lower bounds to the "`v = 1` runs" of THM-4512's
  classes; whether that says anything not already in the class thresholds
  is untested.
* **D4. The branching law as an entropy statement.** Proposition 2 gives
  `|W_(j+1)|/|W_j| = 1/P_j`, so `h* = lim (1/j) sum_i log_2(1/P_i)`; the
  boundary weights `1 - P_i` at the threshold steps (density `1 - log_3 2 =
  0.369` of steps, the parallel note's Beatty lengths) are the whole
  deficit `1 - h* = 0.05`. A closed form for the boundary weight would
  give `h*` exactly.

## 5. Reproduction

    cd 04-computation/experiments
    python3 collatz_posets_dags_zeta5_20260927.py > collatz_posets_dags_zeta5_20260927.out
    python3 collatz_posets_dags_zeta5_20260927_criticality.py > collatz_posets_dags_zeta5_20260927_criticality.out

Both run in under two minutes with the standard library only; the
identities are checked in exact rational arithmetic.

## 6. Independent audit

Pending; to be appended.
