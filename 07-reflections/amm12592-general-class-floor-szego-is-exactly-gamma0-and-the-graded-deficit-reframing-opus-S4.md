---
source: opus-2026-07-31-S4 (the AMM 12592 GENERAL-class lower bound, beyond balanced-block)
status: >
  SUPERSEDED IN ITS INFIMUM CLAIM by THM-3342 and MISTAKE-368; the analytic
  delimitation and graded-deficit reframing survive.  The general result is
  nonattainment of every n+o(n) deadline, not C*>1.  THM-3009's golden floor
  remains balanced-block only.  The analytic/Szego rigidity method is confined to
  EXACTLY gamma=0 -- the two-circle point e^{+-i pi/3} that drives the gamma=0 proof lies ON the boundary of
  the spine's analyticity domain at gamma=0 and EXITS it for every gamma>0 (domain value 2^gamma>1), and the
  radius of convergence drops below 1; so no quantitative general floor follows from this method. (3) Records
  the resummation trap (von Neumann): term-by-term/truncation/p->0 arguments are INVALID because the identity
  is an infinite resummation. (4) Reframes the general floor via klein's LRC graded-atom insight: the general
  (non-block) obstruction is a GRADED family of forced deficits over all critical values m, of which the block
  case (single middle word) is the collapsed special case. Open: the general floor value; two routes named.
tags: [amm12592, minimal-C, lower-bound, general-class, szego, analyticity-domain, resummation, graded-deficit, golden, open]
related: [THM-3009, HYP-9061, amm12592-C-eq-1-szego-rigidity, amm12592-two-ray-threshold-is-golden]
---

# AMM 12592, the general-class floor: Szego rigidity is exactly gamma=0, and the graded-deficit reframing

## 1. The problem, and what is actually proved for the GENERAL class

`C* = inf{C : some D and an exactly-fair extractor with pathwise deadline T <= Cn+D on critical value n}`,
`C* = 1+gamma*`. The general characterization (HYP-9061, no block assumption): pick a prefix-free decided
tree, label `a_{h,t}` of the `N_{h,t}` leaves at composition `(h,t)` with output `1`, such that

```
   sum_{h,t} a_{h,t} p^h (1-p)^t = 1/2   for all p in (0,1),   0 <= a_{h,t} <= N_{h,t} in Z,   (F)
```

with the deadline bounding leaf depths. Equivalently the spine identity
`sum_m p^m q W_m(p) + sum_m q^m p V_m(p) = 1/2`, `deg W_m,V_m <= gamma m + D`, integer Bernstein coefficients.

**Slope one is unattainable in the general class.** THM-3342 works from (F)/(S) directly and proves the
stronger statement that no fixed extractor has `T(n)=n+o(n)`.  In the bounded case,
`gamma=0`, `d_m=D-1` is constant, so the spine coefficients take finitely many values, Szego forces the
generating function rational, the two circles `|p|=1`, `|p-1|=1` meet only at `e^{+-i pi/3}`, and integrality
kills every period. No block-balancedness is used.  Endpoint nonattainment does not separate an infimum, so
the honest general state is

```
   1 <= C*_general <= C*_block,       slope 1 unattained,
   C*_block >= log_5(5 phi^2) = 1.59799                    (THM-3009, within its transfer scope).
```

Whether `C*_general = golden` or is strictly smaller is OPEN: the general class allows UNBALANCED splittings
(`a_{h,t} != N_{h,t}/2`), strictly more freedom than the balanced-block schemes THM-3009 bounds, so a priori
the general floor could be lower.

## 2. NEW (rigorous): the analytic/Szego method is confined to EXACTLY gamma=0

Why doesn't the `gamma=0` proof extend to a quantitative general floor? Precisely because its engine leaves
the stage. Bounding `|W_m(p)| <= (|p|+|1-p|)^{d_m}` gives `g(p)=sum p^m W_m(p)` analytic on
`D_gamma = {|p|(|p|+|1-p|)^gamma < 1}`, and by the `p<->1-p` mirror `g` continues to

```
   D_gamma U D'_gamma = { min(|p|,|1-p|) (|p|+|1-p|)^gamma < 1 }.
```

The whole `gamma=0` argument runs at the two-circle point `w=e^{i pi/3}` (`|w|=|1-w|=1`). Its domain value is

```
   min(|w|,|1-w|) (|w|+|1-w|)^gamma = 1 * 2^gamma = 2^gamma:
      gamma = 0 : value 1  -> w is ON the boundary (rigidity applies, two circles are tangent to the domain),
      gamma > 0 : value 2^gamma > 1 -> w is OUTSIDE the domain (no analytic control there at all).
```

Simultaneously the radius of convergence `R(gamma) = min|p|` on the singular set, attained on the negative
axis at `r(1+2r)^gamma=1`, is `R(0)=1` but `R(gamma)<1` for every `gamma>0` (`R(0.1)=0.902`, `R(1/2)=0.657`,
`R(1)=1/2`), so the spine coefficients are UNbounded and Szego (finite value set) cannot even be invoked.
**Conclusion: the Carlson-Szego-two-circle route proves exactly bounded `gamma=0` impossible and, provably,
nothing quantitative beyond the endpoint.** THM-3342 reaches all sublinear excess by the stronger
Polya--Carlson/Fatou/Kronecker route, but likewise supplies no uniform gap.  A quantitative general floor
must use a new capacity mechanism. (Verified:
`04-computation/amm12592_general_floor_szego_delimitation_opus_S4.py`.)

## 3. The resummation trap (a MISTAKES-style guard)

Naive lower-bound attempts on (F) are INVALID and must be flagged, because (F) holds only after an infinite
resummation. Von Neumann's extractor: label-`1` leaves are `HT`, `(HH|TT)HT`, `(HH|TT)^2 HT, ...`, and

```
   sum_{label 1} p^h q^t = pq + (p^2+q^2)pq + ... = pq/(1-(p^2+q^2)) = pq/(2pq) = 1/2.
```

Each individual term `-> 0` as `p->0`, yet the sum is `1/2` for all `p`. Therefore:
(i) truncating the `m`-sum in (S) DESTROYS the low-degree coefficients (they arise from the tail resummation,
   e.g. `[p^0]=1/2` is not carried by any finite set of spine cells);
(ii) `p->0` / `p->1` limits term-by-term are meaningless (`sum_t a_{0,t}` would be a divergent-or-integer,
   never `1/2`);
(iii) LP/ILP feasibility on a degree-`N` truncation is NOT a valid necessary condition.
Any general-floor argument must respect the resummation -- work with the closed rational form of an
eventually-periodic tree, or with a genuinely convergent functional. (This is why THM-3009 works inside a
finite dyadic block, where the sum is finite and the archimedean expansion at `u=-1` is legitimate.)

## 4. The reframing that survives: a GRADED deficit family (from klein's LRC atom)

klein (LRC covering, S428) found the LRC lower bound's obstruction is a central atom shared by all speeds at
`t=0`, and sharpened it: the atom is **graded by speed** -- arc length `2h/v` -- so the correct peel is
WEIGHTED, not a set subtraction. That is exactly the shape of the general AMM floor. In the balanced-block
case the entire forced deficit collapses to a SINGLE object (the middle word `z=1^m`, THM-3009 sec 1); in the
GENERAL case the forced parity deficits (`binom(d_m,k)` odd, Lucas) sit at EVERY critical value `m`, a graded
family with per-cell budget `binom(d_m,k)` growing like `2^{gamma m}`. So:

> **General floor = the archimedean capacity of a GRADED forced-deficit family across all `m`; the block floor
> is the collapse of that family onto a single middle word.** The two problems (AMM-general, LRC-covering)
> share this exact shape: peel the graded central deficit, bound the remainder by capacity. In both, the
> grading is the source of difficulty and the reason the single-object (block / removable-atom) treatment does
> not transfer.

**The exact general-vs-block distinction (clean).** Every fair extractor equals the balanced baseline plus a
DEFICIT FIELD `b_{h,t} = a_{h,t} - N_{h,t}/2`: since `sum N_{h,t} p^h q^t = sum_leaves p^h q^t = 1`,

```
   sum_{h,t} b_{h,t} p^h (1-p)^t = 0,   b_{h,t} in (1/2)Z,  b_{h,t} half-integer exactly when N_{h,t} is odd,
   |b_{h,t}| <= N_{h,t}/2.
```

The balanced extractor is `b == 0` (needs every `N_{h,t}` even -> dyadic blocks, THM-3007). A general extractor
carries a nonzero `b` -- the forced half-integer deficits (odd `N_{h,t}`, Lucas) -- which must cancel under
`sum b p^h q^t = 0`. **THM-3009's `(ARCH)` is exactly the capacity of this cancellation WITHIN one dyadic shell;
the general class differs only by allowing the deficit to route ACROSS shells (anywhere in the cone
`t <= gamma h + O(D)`).** So the entire "general vs block" gap is one question:

> **Does cross-shell deficit routing have more capacity than within-shell, asymptotically?** If no, then
> `C*_general = C*_block = log_5(5 phi^2)`. If yes, the general floor is strictly smaller and cross-shell
> routing is the mechanism.

This says the missing lemma is one of two:
(a) **Reduction / rearrangement:** every exactly-fair extractor can be symmetrized to a balanced-block one
   without increasing `C` (then `C*_general = golden`). Plausible -- a balanced tree is fair for free
   (`a_{h,t}=N_{h,t}/2` gives `sum = (1/2) sum_leaves p^h q^t = 1/2`) -- but symmetrization must preserve the
   deadline, and the graded deficit is exactly what obstructs a naive averaging.
(b) **Direct graded capacity:** a resummation-respecting archimedean bound on the graded deficit family,
   generalizing THM-3009's `(ARCH)` off the single block. The negative real axis is the right place (it sets
   `R(gamma)` and is where THM-3009 expands), but on the infinite spine it becomes an alternating-sum /
   Abel-summation statement, not a point evaluation.

## 5. Honest status and conjecture

Rigorous: `1 <= C*_general <= C*_block`, slope one and every `n+o(n)` envelope are unattainable; the bounded
Szego method is exactly `gamma=0` (sec 2). Open: the value of
`C*_general`. Conjecture (weakly supported): `C*_general = log_5(5 phi^2)` via route (a) -- because all three
converging lines (my two-ray threshold, THM-3009's block floor, death-star's certified `rho(m)`) are the
golden constant, and no mechanism is known by which unbalanced freedom would beat a balanced-block scheme on
the deadline. Falsifier: an explicit unbalanced extractor with `C < golden`. The concrete next step I can
attempt is route (a)'s rearrangement lemma; route (b) needs the resummed graded-capacity machinery that
death-star / the concurrent opus writeup are building block-locally.
