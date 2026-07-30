# Lane B: C = 1 is impossible — no fair extractor meets the deadline T(n) = n + D, nor even T(n) = n + o(n)

Lane B deliverable for the AMM 12592 frontier (HYP-9061 context; canon THM-2160, THM-2225).
Referee: `laneB_referee.py` (sha256 `3000852d…bfd08`), output `laneB_referee.out`
(final line `LANE-B REFEREE: ALL CHECKS PASS`), both in this directory.

---

## 0. Results

**Setting.** Bits `X_1, X_2, … ∈ {0,1}` are i.i.d. with `P(X_i = 0) = p`,
`P(X_i = 1) = q = 1 − p`, `p ∈ (0,1)` unknown. For a nonconstant stream the
*critical value* is `n = min{k ≥ 1 : X_{k+1} ≠ X_1}` (length of the initial
constant run; matches THM-2160). A *(deterministic) extractor* is a pair
`(Σ, o)` where `Σ ⊆ {0,1}^*` is a stop set — the rule stops at the first
prefix lying in `Σ` — and `o : Σ → {H, T}`. Write `τ` for the stop time.
The extractor is *fair* if

> **(Fair)**  for every `p ∈ (0,1)`:  `P_p(τ < ∞ and output = H) = 1/2`.

(Nothing about `P_p(T)` is assumed; only (Fair) is used.)
Given `T : ℤ_{≥1} → ℤ_{≥1}`, the extractor *meets the deadline `T`* if for
every nonconstant stream `x` (pathwise, not just a.s.):

> **(Deadline)**  `τ(x) ≤ T(n(x))`.

**Theorem B1 (lane target, PROVED).**
*For every integer `D ≥ 0` there is no fair extractor meeting the deadline
`T(n) = n + D` (for all `n ≥ 1` simultaneously).*

**Theorem B2 (strengthening, PROVED).**
*Let `T : ℤ_{≥1} → ℤ_{≥1}` satisfy `(T(n) − n)/n → 0` as `n → ∞` (i.e.
`T(n) = n + o(n)`). Then there is no fair extractor meeting the deadline `T`.*

**Corollary B3.** *Every fair extractor which meets some pathwise deadline
function `T` (finite for each `n`) satisfies `limsup_{n→∞} (T(n) − n)/n > 0`.*
(If the limsup were `0`, then since `T(n) − n ≥ 0` may be assumed — see
Lemma 1 / Remark 6.4 — the limit would be `0`, contradicting B2.)

**Proof status, brutally honest.**
The proofs below are complete and self-contained **except** for four classical
theorems, which are quoted with precise statements and standard references in
§3–§4 and are *not* reproved here: the Pólya–Carlson theorem (Carlson 1921),
Szegő's theorem (1922, used only as an alternative route for B1), Fatou's
lemma on rational power series (1906), and Kronecker's theorem (1857). Every
finite computation used anywhere in the argument is machine-verified in exact
arithmetic by `laneB_referee.py` (checks C1–C9, listed in §7). The spine
reformulation supplied to this lane was **verified**, not assumed: Lemmas 1–3
below prove exactly the parts of it that the impossibility needs (the
realizability converse of the "decided-tree" fact is *not* needed and is only
referee-corroborated). A positive control (check C6) validates the telescoped
identity on the classical von Neumann rule.

**What is *not* proved.** No quantitative lower bound `limsup (T(n) − n)/n ≥ c`
with a universal `c > 0`; hence the frontier constant `C* = inf{C : ∃D,
T(n) ≤ Cn + D achievable}` is *not* shown to exceed 1 (the infimum could
a priori still be 1 via schemes with `T(n) ≤ (1+ε)n + D_ε` for every ε).
B1/B2 exactly kill the `C = 1` lane: no bounded — indeed no sublinear —
additive excess over `n` is possible. See §6.

---

## 1. The spine reformulation, verified

Throughout, fix a fair extractor `(Σ, o)`. For a word `x ∈ {0,1}^*` not
having a prefix in `Σ` ("undecided"), and `p ∈ (0,1)`, define

```
G_x(p) := P_p( output = H | X_1…X_|x| = x )  ∈ [0,1].
```

Since the flips are i.i.d., the conditional law of the future given the prefix
`x` is again i.i.d.(p), so for undecided `x`:

```
G_x(p) = p·G_{x0}(p) + q·G_{x1}(p),                                (R)
```

where for a decided child `y` (i.e. `y` has a — necessarily unique, first-hit —
prefix in `Σ`) we set `G_y := 1` if its stop prefix has output H, else `0`.
Fairness says `G_∅(p) = 1/2` for all `p ∈ (0,1)`.

**Lemma 1 (no stopping on constant prefixes).** *A fair extractor's stop set
`Σ` contains no word of the form `0^k` or `1^k`, `k ≥ 0` (including the empty
word).*

*Proof.* Suppose `0^k ∈ Σ` with `k` minimal (so the rule genuinely stops
there). If `o(0^k) = H` then `P_p(H) ≥ p^k` for every `p`; choosing `p` with
`p^k > 1/2` contradicts (Fair). If `o(0^k) = T`, then `P_p(H) ≤ 1 − p^k < 1/2`
for such `p`, again contradicting (Fair). Same for `1^k` with `q^k`. For
`k = 0` the output is deterministic and `P_p(H) ∈ {0,1}`. ∎

Consequently every `0^m` and `1^m` (`m ≥ 0`) is undecided, and the nodes
`0^m 1` and `1^m 0` (`m ≥ 1`) are reached with the rule still running. Define

```
W_m(p) := G_{0^m 1}(p),      V_m(p) := G_{1^m 0}(p)      (m ≥ 1).
```

**Lemma 2 (spine telescoping; identity (S)).** *For every fair extractor and
every `p ∈ (0,1)`:*

```
Σ_{m≥1} p^m q · W_m(p)  +  Σ_{m≥1} q^m p · V_m(p)  =  1/2.        (S)
```

*Proof.* By (R) at `0^m` (undecided by Lemma 1): `G_{0^m} = p G_{0^{m+1}} +
q G_{0^m 1}`. Iterating from `m = 1` to `M`:

```
p·G_0(p) = Σ_{m=1}^{M} p^m q W_m(p) + p^{M+1} G_{0^{M+1}}(p),
```

and `0 ≤ p^{M+1} G_{0^{M+1}} ≤ p^{M+1} → 0` as `M → ∞` (the constant rays
carry vanishing mass). Hence `p G_0 = Σ_{m≥1} p^m q W_m`, and symmetrically
`q G_1 = Σ_{m≥1} q^m p V_m`. Add and use `1/2 = G_∅ = p G_0 + q G_1`. ∎

**Lemma 3 (decided-tree polynomial form under the deadline).** *Suppose the
extractor meets the deadline `T` and set `d_m := T(m) − m − 1`. Then `d_m ≥ 0`
for all `m`, and for each `m ≥ 1`:*

```
W_m(p) = Σ_{k=0}^{d_m} w_{m,k} · p^{d_m−k} q^k,   w_{m,k} ∈ ℤ,  0 ≤ w_{m,k} ≤ C(d_m, k),
```

*and likewise `V_m(p) = Σ_{k=0}^{d_m} v_{m,k} p^k q^{d_m−k}` with integers
`0 ≤ v_{m,k} ≤ C(d_m, k)`.*

*Proof.* Every infinite stream extending `0^m 1` has critical value exactly
`m` (its initial run is the `m` zeros), so by (Deadline) it stops by absolute
time `m + d_m + 1`; the stopping prefix strictly extends `0^m` (Lemma 1), so
the whole subtree below the node `r = 0^m 1` (depth `m+1`) is decided within
relative depth `≤ d_m`. If `T(m) ≤ m` for some `m` the stream `0^m 1 …` would
have to stop at a constant prefix, contradicting Lemma 1 — so `d_m ≥ 0`.

Now condition on reaching `r` and read the next `d_m` flips `y ∈ {0,1}^{d_m}`.
Each `y` has a unique initial segment at which the rule stops (at relative
depth ≤ `d_m`), with output `out(y)`. Partitioning by `y`:

```
W_m(p) = Σ_{y ∈ {0,1}^{d_m}, out(y)=H} p^{#0(y)} q^{#1(y)} = Σ_k w_{m,k} p^{d_m−k} q^k,
```

with `w_{m,k} := #{y : y has k ones, out(y) = H} ∈ {0, 1, …, C(d_m,k)}`. The
`V_m` case is the mirror image (`#0 ↔ #1`). ∎

*Remarks.* (i) The converse — every such integer vector is realized by some
decided tree — is true (choose any `w_k` of the `C(d,k)` leaves) but is **not
needed** for impossibility; referee check C1 verifies both directions
exhaustively for `d ≤ 3` (e.g. 1446 rules at `d = 3` produce exactly the 64
target polynomials). (ii) A depth-`d_m` form can be refitted to any depth
`d ≥ d_m` by multiplying by `(p+q)^{d−d_m} = 1`; the new weights again satisfy
`0 ≤ w'_j ≤ C(d,j)` by Vandermonde (check C2). This normal form is convenient
but not load-bearing below. (iii) Early stopping *on* a spine, were it
possible, would merely truncate the sums in (S) at that depth with a constant
final term; the analytic argument below would be unchanged. Lemma 1 rules it
out anyway.

---

## 2. Two integer power series

Fix a fair extractor meeting the deadline `T`, and let `d_m = T(m) − m − 1 ≥ 0`.
Write everything in the single variable. Let

```
W̄_m(p) := Σ_k w_{m,k} p^{d_m−k} (1−p)^k        (polynomial, degree ≤ d_m),
P_m(p)  := (1−p)·W̄_m(p) = Σ_{j=0}^{d_m+1} a_m(j) p^j,   a_m(j) ∈ ℤ,
```

and symmetrically, in the variable `u` (which will play the role of `q`),

```
Ṽ_m(u) := Σ_k v_{m,k} (1−u)^k u^{d_m−k},   Q_m(u) := (1−u)·Ṽ_m(u) = Σ_j b_m(j) u^j.
```

**Lemma 4 (integer window coefficients).**
*(a) `|a_m(j)| ≤ 2·3^{d_m}` and `|b_m(j)| ≤ 2·3^{d_m}` for all `j`.*
*(b) The series*

```
F(p) := Σ_{m≥1} p^m P_m(p) = Σ_{n≥0} f_n p^n,     G(u) := Σ_{m≥1} u^m Q_m(u) = Σ_{n≥0} g_n u^n
```

*have integer coefficients `f_n = Σ_{m=1}^{n} a_m(n−m)`, `g_n = Σ_m b_m(n−m)`
(finitely many nonzero terms per `n`), satisfy `|f_n|, |g_n| ≤
2(n+1)·3^{max_{m≤n} d_m}`, converge absolutely on `|p| < 1`, `|u| < 1`, and*

```
F(p) + G(1−p) = 1/2      for all p ∈ (0,1).                        (S')
```

*(c) If `d_m ≤ d` for all m (Theorem B1's regime), then `|f_n|, |g_n| ≤
2(d+2)·3^d`: the coefficients take finitely many values. If `d_m = o(m)`
(Theorem B2's regime), then `limsup |f_n|^{1/n} ≤ 1` and likewise for `g`, so
both series have radius of convergence ≥ 1.*

*Proof.* (a) The coefficient of `p^i` in `W̄_m` is
`Σ_k w_{m,k} (−1)^{i−d_m+k} C(k, i−d_m+k)`, bounded in absolute value by
`Σ_k C(d_m,k) 2^k = 3^{d_m}`; multiplying by `(1−p)` at most doubles the
bound. (Referee check C3 confirms the bound at all `2^{d+1}` vertices of the
weight box for `d ≤ 6` — the coefficients are linear in `w`, so the box
maximum is attained at a vertex; empirically the max ratio to `3^d` is 1.)

(b) The double series `Σ_m Σ_j a_m(j) p^{m+j}` converges absolutely for
`|p| < 1`: with `d_m ≤ εm + C_ε` for every `ε > 0` (which holds in both
regimes), `Σ_m (d_m + 2)·2·3^{d_m} |p|^m < ∞`. Fubini gives the stated `f_n`
(the contributing `m` satisfy `m ≤ n ≤ m + d_m + 1`, at most `n+1` of them,
each `|a_m| ≤ 2·3^{max_{m≤n} d_m}`). By Lemma 2 and Lemma 3,
`F(p) = p·G_0(p)` and `G(u)|_{u=1−p} = (1−p)·G_1(p)` for `p ∈ (0,1)`, so (S)
becomes (S').

(c) Bounded case: only `m ∈ [n−d−1, n]` contribute, i.e. ≤ `d+2` terms of
size ≤ `2·3^d`. Sublinear case: `max_{m≤n} d_m ≤ εn + C_ε` for every ε, so
`|f_n|^{1/n} ≤ (2(n+1))^{1/n} 3^{(εn+C_ε)/n} → 3^ε` for every `ε > 0`. ∎

---

## 3. Analytic continuation and rationality

**Lemma 5 (continuation across p = 1).** *`F` extends to a function analytic
on `D(0,1) ∪ D(1,1)` (open unit disks about 0 and 1). In particular `F` is
analytic at the boundary point `p = 1` of its disk of convergence, and at the
whole open arc `{e^{iθ} : |θ| < π/3}`. Symmetrically for `G`.*

*Proof.* `G` has radius ≥ 1 (Lemma 4), so `h(p) := 1/2 − G(1−p)` is analytic
on `D(1,1)`. `F` is analytic on `D(0,1)`. They agree on `(0,1)` by (S'),
hence on the connected (convex) open lens `D(0,1) ∩ D(1,1)` by the identity
theorem; gluing gives the extension. A point `e^{iθ}` lies in `D(1,1)` iff
`|1 − e^{iθ}| = 2|sin(θ/2)| < 1` iff `|θ| < π/3`. ∎

We now quote the classical inputs. All are stated for a power series
`f(z) = Σ_{n≥0} a_n z^n`.

> **Theorem PC (Pólya–Carlson; F. Carlson, "Über Potenzreihen mit
> ganzzahligen Koeffizienten", Math. Z. 9 (1921), 1–13).** If `a_n ∈ ℤ` for
> all `n` and the radius of convergence equals 1, then either `f` is a
> rational function, or the unit circle is a natural boundary of `f` (every
> boundary point is singular).

> **Theorem Sz (Szegő, "Über Potenzreihen mit endlich vielen verschiedenen
> Koeffizienten", Sitzungsber. Preuß. Akad. Wiss. (1922), 88–91; see also
> Pólya–Szegő, *Problems and Theorems in Analysis* II, Part 8, and Dienes,
> *The Taylor Series*).** If the coefficients `a_n` take only finitely many
> distinct values, then either `f(z) = P(z)/(1 − z^k)` for some polynomial
> `P` and some `k ≥ 1` — equivalently `(a_n)` is eventually periodic — or the
> unit circle is a natural boundary of `f`.

> **Lemma Ft (Fatou's lemma; Fatou, "Séries trigonométriques et séries de
> Taylor", Acta Math. 30 (1906); see Salem, *Algebraic Numbers and Fourier
> Analysis*, or Stanley, *Enumerative Combinatorics* I, §4).** If `a_n ∈ ℤ`
> and `f` is a rational function, then `f = P/Q` with `P, Q ∈ ℤ[z]`,
> `gcd(P,Q) = 1` and `Q(0) = 1`.

> **Theorem K (Kronecker, J. reine angew. Math. 53 (1857)).** A monic
> polynomial in `ℤ[z]` with all roots on the unit circle has all roots equal
> to roots of unity.

**Lemma 6 (pigeonhole periodicity).** *If `f = Σ a_n z^n` is rational and
`(a_n)` takes finitely many values, then `(a_n)` is eventually periodic.*

*Proof.* Write `f = P/Q`, `Q(0) ≠ 0`, `Q(z) = q_0 + … + q_r z^r`. From
`Qf = P`: `q_0 a_n = −Σ_{j=1}^{r} q_j a_{n−j}` for all `n > deg P`, so the
window `(a_{n−r+1},…,a_n)` determines `a_{n+1}`. Windows live in a finite
set; two coincide at indices `N_1 < N_2` beyond `deg P`; forward determinism
then forces period `N_2 − N_1` from `N_1` on. ∎

**Proposition 7 (rationality with root-of-unity poles).** *In the setting of
Lemma 4:*

*(a) (B1 regime, `d_m ≤ d`.) The sequences `(f_n)`, `(g_n)` are eventually
periodic; `F` and `G` are rational with all poles simple and contained in the
roots of unity.*

*(b) (B2 regime, `d_m = o(m)`.) `F` and `G` are rational, and all their poles
(of any multiplicity) are roots of unity.*

*Proof.* If `(f_n)` is eventually zero, `F` is a polynomial with integer
coefficients and (a),(b) hold trivially; likewise if the radius of
convergence exceeds 1 (then `f_n → 0` by Cauchy–Hadamard and `f_n ∈ ℤ`
forces eventual vanishing). Otherwise the radius is exactly 1: `f_n ∈ ℤ \ {0}`
infinitely often gives `limsup |f_n|^{1/n} ≥ 1`, and Lemma 4(c) gives `≤ 1`.
By Lemma 5, `F` is analytic at `p = 1`, so the circle is not a natural
boundary. Theorem PC then makes `F` rational. (For (a) one can quote Theorem
Sz instead and skip PC.)

(a) By Lemma 6 (or Sz directly), `(f_n)` is eventually periodic, say
`f_{n+k} = f_n` for `n ≥ N`. Then

```
F(p) = Σ_{n<N} f_n p^n + p^N C(p)/(1 − p^k),   C ∈ ℤ[p], deg C < k,
```

whose poles are among the simple roots of `1 − p^k`.

(b) By Lemma Ft, `F = P/Q`, `P, Q ∈ ℤ[z]` coprime, `Q(0) = 1`, `deg Q = r`,
leading coefficient `q_r ≠ 0`. The power series converges on `|z| < 1`, so
every pole `ρ` (root of `Q`; genuine pole by coprimality) has `|ρ| ≥ 1`. The
reversed polynomial `Q*(z) := z^r Q(1/z) = z^r + q_1 z^{r−1} + … + q_r` is
monic in `ℤ[z]` with roots `1/ρ_i`, all of modulus ≤ 1, and
`Π |1/ρ_i| = |q_r| ≥ 1` (a nonzero integer). Hence every `|1/ρ_i| = 1`, and
Theorem K makes every `1/ρ_i` — hence every pole `ρ_i` — a root of unity.
Same for `G`. ∎

---

## 4. Pole localization at the primitive sixth roots

**Lemma 8 (the ζ₆ lemma).** *`{z ∈ ℂ : |z| = 1 and |1 − z| = 1} =
{ζ₆, ζ̄₆}` where `ζ₆ = e^{iπ/3} = (1 + i√3)/2`. These are the primitive 6th
roots of unity, the roots of `Φ₆(p) = p² − p + 1`, and `1 − ζ₆ = ζ̄₆ = ζ₆⁵`.
Consequently `{z : z and 1 − z are both roots of unity} = {ζ₆, ζ̄₆}`.*

*Proof.* `|1 − e^{iθ}|² = 2 − 2cos θ = 1` iff `cos θ = 1/2` iff
`θ = ±π/3`. Both are roots of unity, so the second set equals the first.
The remaining identities are direct; referee check C4 verifies all of them
exactly (sympy, exact arithmetic). ∎

**Proposition 9.** *In either regime, all poles of `F` lie in `{ζ₆, ζ̄₆}`,
and `F` is analytic at `p = 0` and `p = 1`; likewise for `G`.*

*Proof.* (S') holds on `(0,1)`; both sides are rational (Prop. 7), so
`F(p) = 1/2 − G(1−p)` holds as an identity of rational functions. Therefore
the pole set of `F` equals `1 − (pole set of G)` (with multiplicities).
Poles of `F` are roots of unity; poles of `G` are roots of unity; so every
pole `z` of `F` has `z` and `1 − z` both roots of unity, and Lemma 8 applies.
`p = 0` is regular (power series); `p = 1` is regular because `G` is regular
at `u = 0`. ∎

---

## 5. The integrality endgame

**Lemma 10 (integer polynomial part).** *Suppose a power series
`F = Σ f_n p^n` with `f_n ∈ ℤ` is rational with all poles in `{ζ₆, ζ̄₆}`
(any multiplicities). Then*

```
F(p) = A(p) + h(p),   A ∈ ℤ[p],   h(p) = N(p)/Φ₆(p)^e,  N ∈ ℝ[p],  deg N < 2e
```

*(`h = 0` allowed). The same holds for `G` in the variable `u`.*

*Proof.* As a rational function with poles only at the two conjugate points
`ζ₆, ζ̄₆`, `F = A + h` with `A ∈ ℝ[p]` a polynomial and `h` a proper rational
function (`deg N < 2e`) with denominator `Φ₆^e`, since
`(1 − p/ζ₆)(1 − p/ζ̄₆) = Φ₆(p)` (Euclidean division of the Fatou fraction
`P` by `Q = Φ₆^e`; that `Q = Φ₆^e` for real `F` with these poles is forced:
`Q(0) = 1` and conjugate-symmetric multiplicities).

Coefficients of `h`: partial fractions over ℂ give
`h = Σ_{j=1}^{e} [α_j/(1 − p/ζ₆)^j + ᾱ_j/(1 − p/ζ̄₆)^j]`, so the `n`-th
coefficient is

```
c_n = R(n)·ζ₆^{−n} + conj(R(n))·ζ₆^{n},   R ∈ ℂ[n], deg R ≤ e − 1
```

(using `[p^n] (1 − p/ζ)^{−j} = C(n+j−1, j−1) ζ^{−n}`). Restricted to a fixed
residue class `n ≡ ρ (mod 6)`, `ζ₆^{−n} = ζ₆^{−ρ}` is constant, so
`n ↦ c_n` is a real polynomial of degree ≤ `e − 1` on each class (referee
check C5 verifies this classwise quasi-polynomial structure exactly for
`e = 1, 2, 3` to order 72; for `e = 1` the coefficient pattern of `1/Φ₆` is
the period-6 word `1,1,0,−1,−1,0`).

Since `A = F − h` is a polynomial, `c_n = f_n ∈ ℤ` for all `n > deg A`. A
real polynomial of degree ≤ `e−1` taking integer values at all sufficiently
large integers of an arithmetic progression takes integer values on the whole
progression: writing `g(t)` for its restriction and `T` large,
`g(t) = Σ_{i<e} Δ^i g(T)·C(t−T, i)` (Newton; referee check C9 verifies the
identity symbolically for degrees ≤ 5), where `Δ^i g(T) ∈ ℤ` and `C(x,i) ∈ ℤ`
for **all** integers `x`, including negative. Hence `c_n ∈ ℤ` for all
`n ≥ 0`, and `A`'s coefficients `f_n − c_n` are integers: `A ∈ ℤ[p]`. ∎

**Theorem B2 (hence B1).** *No fair extractor meets a deadline
`T(n) = n + o(n)`.*

*Proof.* Suppose one exists. Build `F, G` (Lemma 4; `d_m = T(m) − m − 1 ≥ 0`,
`d_m = o(m)`). By Prop. 7(b), Prop. 9 and Lemma 10:

```
F(p) = A(p) + N(p)/Φ₆(p)^e,      G(u) = B(u) + Ñ(u)/Φ₆(u)^{ẽ},
A ∈ ℤ[p],  B ∈ ℤ[u],  deg N < 2e,  deg Ñ < 2ẽ.
```

Substitute into (S') and use the key involution identity

```
Φ₆(1 − p) = (1−p)² − (1−p) + 1 = p² − p + 1 = Φ₆(p)
```

(referee check C5). Then

```
1/2 − A(p) − B(1−p)  =  N(p)/Φ₆(p)^e + Ñ(1−p)/Φ₆(p)^{ẽ}.          (†)
```

The right side is a *proper* rational function (numerator degree <
denominator degree after clearing to the common denominator `Φ₆^{max(e,ẽ)}`),
so it tends to 0 as `p → ∞`; the left side is a polynomial. A polynomial
with limit 0 at infinity is identically 0, so both sides of (†) vanish:

```
A(p) + B(1 − p) = 1/2       identically.
```

Evaluate at `p = 0`: `A(0) + B(1) = 1/2`. But `A(0) ∈ ℤ` and `B(1) ∈ ℤ`
(integer polynomial at an integer). An integer equals 1/2 — contradiction. ∎

*Proof of B1 as a stand-alone (simpler chain, `d_m ≤ D − 1`).* Coefficients
of `F, G` take finitely many values (Lemma 4(c)); continuation across `p = 1`
(Lemma 5) plus Theorem Sz (or PC + Lemma 6) makes `(f_n), (g_n)` eventually
periodic; poles are then simple and lie among roots of `1 − p^k`
(Prop. 7(a)); pole matching localizes them to `{ζ₆, ζ̄₆}` (Prop. 9); with
simple poles, `e ≤ 1` in Lemma 10 and the tail of `(f_n)` is the exact
period-6 pattern `c_n = b₀α_n + b₁α_{n−1}` with `α = (1,1,0,−1,−1,0)^∞`
(check C5), whose integrality at two large indices already gives
`b₀, b₁ ∈ ℤ`; the endgame (†) is identical. ∎

*Step-1 corollary in the assigned formulation (depth 0, `T(n) = n + 1`).*
There are no subsets `S, T ⊆ {1, 2, …}` with `q·Σ_{m∈S} p^m + p·Σ_{m∈T} q^m
= 1/2` on `(0,1)` — this is B1 with `D = 1`, where `W_m, V_m ∈ {0,1}`.
Referee check C7 independently corroborates this exhaustively for all
eventually periodic `S` with pre-period ≤ 2 and period ≤ 8 (3570 configs:
the forced companion `T`-series never has `{0,1}` coefficients through order
60), and check C8 does the analogous corroboration at depth `d = 1` (period
≤ 5 for `(W_m)`, arbitrary `(V_m)` via an exact sliding-window DP; all 1364
configs die before depth 60).

---

## 6. Remarks, scope, and frontier placement

**6.1 Why the argument stops at `o(n)`.** The only place the deadline enters
is Lemma 4: `|f_n| ≲ n·3^{max_{m≤n} d_m}` must be subexponential to give
radius 1 for the Pólya–Carlson dichotomy. At any linear depth `d_m ~ γm`
(γ > 0) the natural radius drops below 1 and the entire rigidity mechanism
(integer series on the *unit* circle) disappears. This is consistent with the
existence of the canon scheme at `T(n) ≤ max(2, 2n−1)` (`d_m ~ m`): check C6
verifies the telescoped identity (S) exactly on the classical von Neumann
rule (`W_m = 1` for `m` odd, `q/2` for `m` even; `V_m = 0` odd, `q + p/2`
even; restart value `pq/(1 − p² − q²) = 1/2`), a positive control showing the
reformulation's bookkeeping is correct where solutions do exist.

**6.2 What B1/B2 mean for HYP-9061.** The frontier question is the minimal
`C` with `T(n) ≤ Cn + D` achievable. Canon gives `C = 2` (`T ≤ max(2,2n−1)`,
THM-2225) and `T(n) ≤ 2n − 2` for `n ≥ 3` (THM-2160 §5). Lane B settles the
bottom of the scale: **the pair (C, D) = (1, D) is infeasible for every D**,
and more strongly no `n + o(n)` envelope is feasible. Lane B does **not**
decide whether `C* > 1` as an infimum (nor anything about `C* = 9049/6592`
or the decoded certificate (27), which belongs to the upper-construction /
lower-bound-rate discussion at linear rates; nothing here contradicts it).
The natural next target is a quantitative version: a lower bound
`limsup (T(n) − n)/n ≥ c₀ > 0` uniform over schemes, which requires a
genuinely different (quantitative, e.g. deficit-flow/carry-transport) method:
the present proof is an all-or-nothing analytic rigidity argument.

**6.3 Robustness.** The proof uses of the weights only that they are
*integers* with `|w_{m,k}| ≤ C(d_m,k)`-type size (via Lemma 4(a)); their
nonnegativity and the exact binomial cap are never used beyond the size
bound. Fairness is used only as `P_p(H) = 1/2` for all `p` (Lemma 1 needs it
for `p` near 1; Lemma 2 pointwise). The n = 2 causality obstruction of
THM-2160 §4 (envelope `max(2, 2n−2)`) is independent of and not implied by
Lane B: that envelope is not of the form `n + o(n)`.

**6.4 WLOG notes.** If `T(m) ≤ m` for some `m`, Lemma 3 already yields a
contradiction (a constant prefix would have to be terminal), so `T(n) ≥ n+1`
i.e. `d_m ≥ 0` is automatic, and Corollary B3's `T(n) − n ≥ 0` normalization
is harmless. The extractor need not be assumed to stop a.s.: only (Fair) for
H is used.

**6.5 Classical inputs — exact dependence list.**
- Theorem PC (Carlson 1921): used in Prop. 7 (both regimes; avoidable in the
  B1 regime via Theorem Sz).
- Theorem Sz (Szegő 1922): optional alternative for B1.
- Lemma Ft (Fatou 1906): used only in the B2 regime (Prop. 7(b)).
- Theorem K (Kronecker 1857): used only in the B2 regime (Prop. 7(b)).
Everything else is proved above; every finite computation appears in
`laneB_referee.py` (§7).

---

## 7. Referee coverage map

| Check | Verifies | Used in |
|---|---|---|
| C1 | decided-tree fact, both directions, `d ≤ 3` (exhaustive over all 2, 6, 38, 1446 rules) | Lemma 3 (corroboration; forward direction proved) |
| C2 | Vandermonde refinement bound, `d ≤ 8` | Lemma 3 remark (ii) |
| C3 | `|coeff((1−p)W̄)| ≤ 2·3^d`, `d ≤ 6`, all weight-box vertices | Lemma 4(a) |
| C4 | `{|z|=1} ∩ {|1−z|=1} = {ζ₆, ζ̄₆}`, primitivity, `1−ζ₆ = ζ̄₆ = ζ₆⁵`, `Φ₆ = p²−p+1` | Lemma 8 |
| C5 | `Φ₆(1−p) = Φ₆(p)`; `1/Φ₆` period-6 pattern; linear/Φ₆ coefficient formula; classwise quasi-polynomials for `1/Φ₆^j`, `j ≤ 3` | Lemma 10, Theorem proof (†) |
| C6 | von Neumann positive control: exact identity (S) | Lemma 2 sanity (§6.1) |
| C7 | depth 0: 3570 eventually periodic `S`-configs, forced companion fails by order 60 | Step-1 corollary (corroboration) |
| C8 | depth 1: 1364 periodic `(W_m)`-configs, exact DP over all companions dies by depth 60 | B1 at `d=1` (corroboration) |
| C9 | Newton backward-propagation identity (deg ≤ 5), `C(x,i) ∈ ℤ` for all `x ∈ ℤ` | Lemma 10 |

Output: `laneB_referee.out`, final line `LANE-B REFEREE: ALL CHECKS PASS`
(exit 0). Script sha256 `3000852d5403792bcfe77ece03b99de6c9519393cae7e419b018eb5aa90bfd08`,
output sha256 `61e6bb3eb7d150539fe041691d9a974a941e7962f0ce5877672c06dd56eb34e7`.
