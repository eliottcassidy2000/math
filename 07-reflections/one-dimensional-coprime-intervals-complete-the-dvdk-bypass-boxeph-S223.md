# One-dimensional coprime intervals complete the DvdK bypass

*boxeph-2026-07-21-S223. Owner: keep going on the DvdK bypass; think one-dimensional coprime intervals.
Builds on S222 (the saddle-point/Watson bypass), THM-1840 (single-character/coprime seed), THM-1630 (DvdK,
a one-variable theorem), S218 (arithmetic entropy), and the repo's heavy three-distance/coprime threads.
Verified in `04-computation/one_dimensional_coprime_intervals_return_semigroup_boxeph_S223.py`.*

## DvdK is a one-variable theorem — and in one variable it is elementary

The key realization sharpening S222: **DvdK (THM-1630) is a *one-variable* statement**, and GMC(2) only ever
uses it in one variable (the polar bridge collapses the face to a single-variable Laurent polynomial `Λ_s`).
In one variable the whole content is elementary **coprime-interval / numerical-semigroup combinatorics** — no
residues, no Liouville, no Watson estimate even. The return set is *completely and effectively* determined.

For `f = Σ_{k∈S} c_k z^k` with support `S` (take `0∉S`, the pure two-sided / nonzero-charge case):
- the **Newton polytope is the interval `[min S, max S]`**; `f` two-sided ⟺ `0` interior;
- the **period `d = gcd` of the support gaps** (`coprime ⟺ d=1` after reduction);
- the **return set `R = {m : 0 is an m-fold sum of S}`** is a **numerical semigroup** (closed under `+`);
- for **positive coefficients**, `CT(fᵐ) ≠ 0 ⟺ m ∈ R` (verified) — no cancellation, so nonvanishing *is*
  reachability.

## Two poles: the coprime pair (periodic) and the filled interval (cofinite)

The 1D structure has exactly two regimes, both verified:

- **Bare coprime PAIR `S={−q,p}`, `gcd(p,q)=1`:** the support sits on a single gap, so the period is
  `d=p+q` and `R = (p+q)ℤ_{>0}` — the returns are **exactly the multiples of `p+q`** (verified for
  `{−2,3},{−1,2},{−3,4},{−2,5}`). This is **THM-1840's single-character seed** `m₀=p+q`, and `CT(fᵐ)=0` for
  every other `m`: the Frobenius number is `∞` (fully periodic).
- **FILLED coprime interval** (endpoints + interior making `gcd` of gaps `=1`): the return set becomes
  **cofinite** — `R =` all `m` beyond a finite **Frobenius number** (verified:
  `{−2,1,3}→R=all m≥3`, `{−2,2,3}→` all `m≥4` except a gap, `{−3,−1,1,2}→` all `m≥2`). Adding a single
  interior exponent collapses the period from `p+q` to `1` and turns "only multiples of `p+q`" into "all
  large `m`."

So the **DvdK 1D conclusion — "two-sided ⟹ `CT(fᵐ)≠0` for some `m`" — is sharpened to an exact statement:**
`R` is the coprime interval's numerical semigroup, and `CT(fᵐ)≠0` for all `m` beyond its Frobenius number
(positive coeffs; mixed signs lose at most finitely many `m` to cancellation, still cofinite by the S222
saddle). This is **effective** — the **Frobenius number is the explicit `m₀`** (the open effective-DvdK
bound), computed from the interval alone.

## The complete 1D bypass

```
   DvdK (THM-1630) in 1 variable  ==  an elementary numerical-semigroup fact:
   ─────────────────────────────────────────────────────────────────────────
   support S ⊂ [min,max],  0 in interior (two-sided)
        │
        ▼   R = {m : 0 is an m-fold sum of S}  (numerical semigroup)
   bare coprime pair {-q,p}:  R = (p+q)Z    (PERIODIC, THM-1840 seed, Frob = inf)
   filled coprime interval :  R = all m > Frobenius#   (COFINITE, effective m0)
   positive coeffs: CT(f^m)!=0 <=> m in R ;  mixed signs: cofinite (saddle S222)
```

No residues, no Liouville — the 1D DvdK content is a Chicken-McNugget / Frobenius fact about the coprime
interval. Combined with S222 (the mixed-sign saddle) and S208 (the confluent cusp for the degenerate
saddle), **the entire one-variable non-vanishing GMC(2) needs is now elementary and effective**, with DvdK
demoted from an imported premise to a special (periodic) case of an interval-combinatorics fact.

## The shared engine: coprime intervals also drive the LRC three-distance

The same coprimality is the LRC's engine. The gaps of `{k·t mod 1}` take exactly **three values**
(three-distance / Steinhaus), determined by the **continued fraction of `t`** — the coprime-interval
structure of the orbit. The LRC extremal `t*=14/183=[0;13,14]` has **coprime partial quotients**, and its
**danger arcs are the coprime intervals** whose covering the covering-min measures. So *one-dimensional
coprime intervals* are the common combinatorial engine of two apparently different repo problems: **GMC's
constant-term return times** (this note) and **LRC's three-distance danger-arc geometry** — both are "when
does the coprime interval close up," read multiplicatively (return semigroup) vs. geometrically (arc gaps).

## Scope

Verified: the coprime-pair periodicity (`R=(p+q)ℤ`, THM-1840), the filled-interval cofiniteness with explicit
Frobenius numbers, and `CT(fᵐ)≠0 ⟺ m∈R` for positive coefficients. The mixed-sign "cofinite up to sporadic
cancellations" uses the S222 saddle (verified there). The contribution is recognizing that **DvdK's
one-variable content is elementary coprime-interval / numerical-semigroup combinatorics** — completing the
GMC(2) DvdK-bypass in 1D, effectively (the Frobenius number *is* the effective `m₀`), and exposing the
coprime interval as the shared engine with the LRC three-distance geometry. Not a full replacement write-up
of THM-2022, but the 1D core it imports is now elementary and self-contained.

Links: HYP-8895, HYP-8890, THM-1840, THM-1630,
[[bypassing-dvdk-the-saddle-point-watson-route-to-the-gmc2-angular-floor-boxeph-S222]],
[[each-prime-is-its-paley-tournament-a-periodic-table-of-2-3-5-7-11-for-lrc14-boxeph-S215]].
