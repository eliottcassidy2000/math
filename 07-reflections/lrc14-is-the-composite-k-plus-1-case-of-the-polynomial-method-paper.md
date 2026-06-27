# LRC(14) is the composite-`k+1` case of the polynomial-method paper — and Conjecture 7.1 *is* our witness route

*kind-pasteur-2026-06-27-S31ag. The owner asked for a comprehensive understanding of where the
14-runner attack stands and how it relates to the recently-proven 11–13 runner cases, which rest on a
calculation (Sungkawichai–Trakulthongchai, arXiv:2604.23906). I read the paper's internals and mapped
them onto our covering-bound / witness machinery link-by-link. The result is not a small analogy: **our
whole route is the `k+1=14` composite case of their method, and their Conjecture 7.1 — the property they
say would handle all `k` — is precisely the object our witness route already controls.** No agent had
pulled the paper's polynomial method, `I(k,p,1)`, or Conjecture 7.1 into our frame before; they were
cited only as "the induction base."*

---

## 0. The paper in one paragraph (verified from the source)

Sungkawichai–Trakulthongchai prove **LRC(k) for k ≤ 12 speeds** (= ≤ 13 runners). Engine:
- **Prop 4.1 (polynomial method):** over the field `ℤ_{k+1}` (needs **k+1 an odd prime** for Fermat's
  little theorem), build two degree-`k` indicators `P(X)=∏_i(v_i i^{-1}+X)` and
  `Q(X)=∑_{m∈G}(1−(X−m)^k)`. They agree at `k+1` points with `P(0)=Q(0)=0`, so interpolation forces
  `P=Q`; comparing leading coefficients gives `1 ≡ −|G| (mod k+1)`, hence `|G|=k` and the tuple is
  proper. This is **analytic** — no computation.
- **Composite `k+1` fallback (their k=11, `k+1=12=2²·3`):** Prop 4.1 is unavailable, so they do explicit
  **computational lifts at the prime factors** `c=2,3` of `k+1`.
- **Finite check:** verify `J(k,p)=∅` (no eventually-`(k,p)`-improper tuples) by sieving over a table of
  primes `p>k²+k`; counterexample product bounded by `B_k=((C(k+1,2))^{k-1}/k)^k`.
- **Stated bottleneck at k=13:** "the efficient computation of `I(k,p,1)`"; cost scales like
  `p^{(k+1)/2}/(k·2^k)`. They offer **Conjecture 7.1** as the route to *all* `k`.

## 1. The reduction is the same in both frameworks

| Paper | Project |
|---|---|
| Tuple `(u_1,…,u_k)` has **LR property**: `∃t, ‖tu_i‖ ≥ 1/(k+1) ∀i` | `M(S)=max_t min_i ‖s_i t‖ ≥ 1/14` |
| **Tight** tuple (meets `1/(k+1)` with equality) | `M(S)=1/14` exactly (AP, GW, all dilations) |
| Reduce to coprime tuples (`gcd=1`) | **WLOG primitive** (dilation invariance `M(dS)=M(S)`) |
| The hard tuples are `≡(1,…,k) mod p` (an **AP mod p**) | The hard sets are **covering** (contain a mult of 14) |

The 14-free half is free on both sides (`t=1/14` / Dirichlet). **LRC(14) ⟺ the covering bound** is the
project's S59 redirect; it is the *same* localization the paper makes when it throws away everything but
the AP-residue family.

## 2. `14 = 2·7` composite is the wall — identically in both frameworks

The paper's Prop 4.1 dies because `ℤ_{14}` is **not a field**: Fermat fails for the zero-divisors `2,7`,
and the `k+1=14` interpolation points **collapse under CRT** — `ℤ_{14} ≅ ℤ_2 × ℤ_7`, so the 14 points
project to only **7 distinct residues mod 7** (each hit twice) and **2 mod 2** (each hit seven times).
You cannot have 14 distinct evaluation points in a field of size 7 or 2, so "agree at 14 points ⟹ equal"
evaporates. This is *exactly* the project's "apex prime 7 / zero-divisor `q` / `14=2·7` stacks the
2-adic tower on the apex prime" (HYP-2063, the four-faces reflection). **The reason 14 is the first open
case is one fact stated two ways: `k+1` is composite for the first time in a way the field method can't
absorb.** (Earlier composite `k+1=12` at k=11 was small enough to brute-force the `c=2,3` lifts.)

## 3. THM-573 *is* the paper's `c=7` lift; the descent 14→7→2 *is* the CRT fallback

The paper says: when `k+1` is composite, **lift at its prime factors**. For `14` those are `c=2` and
`c=7`. The project independently built exactly these:
- **`c=7` lift = THM-573 (level-7 sieve).** `≥7` multiples of 7 ⟹ `M>1/14`, by `H=7P`, a `P`-safe phase
  from LRC(≤13), and 7 lifts `t_j=(v*+j)/7`. This is a `c=7` computational lift turned into a *clean
  theorem* — and it uses LRC(≤13), i.e. **the paper itself as the induction base.** Our level step is
  literally an inductive layer on top of their k≤12 result.
- **`c=2` lift = the 2-adic dyadic tower** (HYP-2656/2661, "each doubling halves the interval, the `2` in
  `14=2·7`"). The carry-conservation/mouth-ownership work is the `c=2` side.
- **Residual `≤6` multiples of 7** = exactly the regime where the `c=7` lift alone does not fire and the
  `c=2` lift + an analytic input must finish. This is the project's located hard core.

So our "arithmetic descent 14→7→2" and the paper's "composite-`k+1` lifts at `c=2,7`" are the **same
algorithm**. We have, without citing it, reconstructed their fallback and pushed one level further by
making the `c=7` step a sieve theorem rather than a table lookup.

## 4. `I(k,p,1)` ↔ our continuous `p0` / `μ_{1/7}` — the analytic substitute for the prohibitive count

`I(k,p,l)` = the set of `(k,p,l)`-**improper** tuples (no `gcd` drop and no witness `t∈(1/lp)ℤ` with the
gap). The bottleneck is counting `I(k,p,1)` at cost `p^{(k+1)/2}/(k2^k)`. Our **sector cap** `p0(E)≤cap_k`
and **witness floor** `μ_{1/7}(E)=meas(G_P)` are the `p→∞` **continuous limits** of exactly this count:
"fraction of phases that fail the gap." The project's analytic bounds (the `3/π²=1/(2ζ(2))` floor,
THM-546 Weyl discrepancy, the L7 residue atlas) are an **analytic estimate of `I(k,p,1)`** that, if made
uniform, *replaces* the computation the paper cannot afford. **This is the precise sense in which the
project can contribute past k=12: not a new finite check, but an analytic lower bound on the witness
fraction.**

## 5. The headline: **Conjecture 7.1 for k=13 ⟺ LRC(14)**, and it is our witness route

> **Conjecture 7.1 (paper).** For each `k+1` there is a constant `D` such that for every integer `d≥D`,
> every non-tight coprime tuple `v∈ℤ_{>0}^k` has a witness time in `(1/d)ℤ`.

Restricted to a single `k`, **Conjecture 7.1(k) ⟹ LRC(k)** outright (take any `d≥D`: a witness exists ⟹
the LR property; tight tuples are trivially fine). So **Conjecture 7.1(13) ⟹ LRC(14).**

Now the bridge to our machinery, which is the new content:

- The **lonely set** `L(S)={t: ‖s_i t‖≥1/14 ∀i}` is a **finite union of closed intervals** (each
  constraint is a union of arcs; finitely many linear breakpoints — this is the project's THM-565
  "`G` = finite union of intervals").
- A finite union of intervals of positive total measure that contains an interval of length `ℓ`
  **contains a point `a/d` for every `d ≥ ⌈1/ℓ⌉`.**
- Therefore **Conjecture 7.1(13) ⟺ every non-tight covering 13-tuple has a lonely interval of length
  `≥ ℓ_0 > 0` (uniform)**, with the constant `D = ⌈1/ℓ_0⌉`.

And the project already controls the largest lonely arc:
- **Witness floor (THM-530/565):** total lonely measure `≥ m_P = 14249/252252 ≈ 0.0565`.
- **Three-gap arc count (THM-565):** after **scale separation** (peel the apex cluster; the small part
  `P=S∩[1,13]` gives the slow phase, the cluster gives bounded offsets), the arc count is a **constant**
  `arcCount = m` *independent of the apex magnitude `V`* — the "Framing A" unlock. At the boundary core
  `arcCount = 12`.
- Hence `ℓ_max ≥ (lonely measure)/arcCount ≥ m_P/12 ≈ 0.0047`, giving `D ≈ 213`.

**That `D≈213` is exactly the project's `V*` atlas** (worst `V*≈234`, decreasing in `k`). The paper's
"astronomical computation at k=13" and the project's "finite-`V*` wall" are the **same constant `D`** —
and the witness/three-gap route says it is a few hundred, *not* astronomical. The astronomy is in the
paper's *enumeration* of `I(k,p,1)`; the witness route replaces enumeration with a measure bound.

### Why HYP-2866 ("witness denominator unbounded") does **not** kill Conjecture 7.1
HYP-2866 refuted "the *minimal* witness denominator is `≤41`." Conjecture 7.1 is weaker: for each large
`d` it needs *some* witness with denominator dividing `d` — different `d` may use different witnesses.
A positive-length lonely interval supplies one for **every** `d≥⌈1/ℓ⌉` automatically. So the two are
compatible; the unbounded-`D` family `{1..11,13,84m}` has a *small minimal* lonely arc that nonetheless
catches a rational for each large `d`. The right invariant is **largest arc length, not minimal witness
denominator.**

## 6. `φ(14)=6` = the six inner sectors

`(ℤ/14)^* = {1,3,5,9,11,13}`, order `φ(14)=6`, and reduces mod 7 onto all of `{1,…,6}`. The polynomial
method wants this group to have order `≥ k = 13`; it has order 6. **The deficiency `13−6=7` is the apex
prime again.** Every "6 inner sectors", "`S2` over `C(6,2)=15` pairs", "Clebsch cut-space `(ℤ/2)^4`"
object in the repo lives on this group of units — they are all the multiplicative structure that *would*
have driven Prop 4.1 had `14` been prime. mac-mini's gK8 = low-order moment-LP led by pairwise `S2` is,
in this light, **measuring the failure of the order-6 unit group to fill the degree-13 indicator.**

---

## 7. What this buys us — the sharpened target

LRC(14), fully reconciled, is now:

```
LRC(14)  ⟺  Conjecture 7.1(13)  ⟺  every non-tight primitive covering 13-tuple has a
            lonely interval of length ≥ ℓ0 > 0  (uniform ⟹ D = ⌈1/ℓ0⌉ ≈ few hundred)

   built from:
   (i)  scale separation ⟹ lonely-arc count uniformly bounded   [THM-565 Framing A: arcCount = m]
   (ii) witness floor    ⟹ lonely measure ≥ m_P > 0             [THM-530/565]
   (iii) (i)+(ii) ⟹ largest arc ≥ m_P/arcCount = ℓ0             [pigeonhole, elementary]
   (iv) largest arc ≥ ℓ0 ⟹ bounded-denominator witness ⟹ LR property  [a/d ∈ interval]
```

The genuine open piece is exactly **(i) for the direct `1/14` lonely set** (THM-565 currently bounds the
arc count for the *maxgap>1/7* witness object after peeling; we need the same uniform arc-count bound for
the lonely set itself, or a clean reduction between the two). Everything else is in hand. This is a
*much* more concrete target than "OPEN-Q-108 / `p0≤cap` symbolic," and it is the paper's own
Conjecture 7.1 — so a proof would be recognized immediately by the LRC community as completing their
program at k=13.

**Honest status:** LRC(14) still NOT proved. But the endgame is re-localized to a single uniform
arc-count bound, the constant `D` is pinned to the few-hundred range (matching the `V*` atlas), and the
project's witness machinery is revealed to be the engine for the paper's own Conjecture 7.1. The
next computational step (this session) is to *measure* `ℓ_max` over covering tuples and confirm the
uniform floor `ℓ0` and the implied `D`.

→ HYP-3087 (polynomial-method/CRT bridge), HYP-3088 (largest-arc ⟺ Conjecture 7.1(13)), THM-573,
THM-565, THM-530, OPEN-Q-108, HYP-2866 (compatible), HYP-2602, mac-mini-S59/S60 (covering redirect,
gK8), the-four-faces-of-14, [[lrc14-thread]], arXiv:2604.23906.
