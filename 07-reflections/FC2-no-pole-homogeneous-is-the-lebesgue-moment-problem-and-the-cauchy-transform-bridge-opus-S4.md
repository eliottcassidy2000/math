---
source: opus-2026-07-31-S4 (FC(2) via the no-pole / exponential-period approach; tournament analogy)
status: >
  RESULT + reduction + verified-through-degree-4 + methodological guard. FC(2) is the HOLOMORPHIC moment
  problem for the complex random variable Z=f(X,Y), X,Y~Exp(1); the no-pole reformulation is exactly the
  Cauchy-transform condition C_mu(w)=1/w for mu=f_*(Exp^2). For HOMOGENEOUS f the radial substitution
  collapses it EXACTLY: L(f^m)=(dm+1)! int_0^1 phi^m da with phi(a)=f(a,1-a) ranging over ALL polynomials
  (Bernstein basis). So FC(2)-homogeneous <=> the LEBESGUE polynomial moment problem
  (int_0^1 phi^m=0 forall m => phi=0), i.e. Pakovich's PMP with weight Q'=1, Q=a. Simple-arc images force
  mu=delta_0 (Cauchy transform); a counterexample would need a self-intersecting/area-enclosing image, which
  Pakovich's composition condition rules out for an injective Q (no W with W(0)=W(1)). EXACT confirmation:
  no non-constant solution for degree <= 4. GUARD: double-precision int phi^m via coefficient sums suffers
  catastrophic cancellation and manufactured false "solutions"; high-precision/exact is mandatory. BRIDGE: the
  Cauchy/R-transform here is the SAME free-probability object behind the AMM golden floor and the Paley
  tournament cluster -- FC no-pole, AMM edge, Paley R->e are three faces of the Cauchy transform.
tags: [factorial-conjecture, FC2, no-pole, cauchy-transform, moment-problem, pakovich, polynomial-moment-problem, exponential-periods, free-probability, paley, tournament, homogeneous]
related: [factorial-conjecture-exponential-periods-and-the-repo-state, cross-shell-floor-is-golden-and-the-amm-paley-semicircle-bridge, THM-2022, THM-1510]
---

# FC(2) via no-pole: the homogeneous case IS the Lebesgue moment problem (Cauchy transform = 1/w)

## 1. The no-pole reformulation, concretely

`FC(2)`: for `f in C[x,y]`, `L(f^m)=0 forall m>=1 => f=0`, `L(g)=int_{[0,inf)^2} g e^{-x-y}dxdy = E[g(X,Y)]`,
`X,Y ~ Exp(1)`. So `L(f^m)=E[Z^m]`, `Z=f(X,Y)`: FC(2) is the **holomorphic moment problem** for the complex
random variable `Z` (only the moments `int z^m dmu`, never `int zbar^n`, are constrained -- which is exactly
why real moment-determinacy fails and cancellation is possible). The exponential generating function is the
two-sided transform `E[e^{tZ}] = int e^{tz} dmu(z)`, `mu = f_*(e^{-x-y}dxdy)`; the repo's **no-pole**
criterion (deathstar-S61e) -- FC holds iff `E[Q e^{tf}]` is entire -- is the statement that the **Cauchy
transform** `C_mu(w) = int dmu(z)/(w-z) = 1/w` (all moments `>=1` vanish `<=>` `C_mu(w)=1/w` for large `w`).
The failure mode is a pole: Long's/repo's counterexample has `E[Q e^{tP}] = 1/(1-t)`.

## 2. The homogeneous case collapses EXACTLY to the Lebesgue polynomial moment problem

For `f` homogeneous of degree `d`, the radial substitution `x=ra, y=r(1-a)` (`a in [0,1]`, `dxdy = r dr da`,
`f = r^d phi(a)`, `phi(a):=f(a,1-a)`) gives, using `int_0^inf r^{dm+1}e^{-r}dr = (dm+1)!`,

```
   L(f^m) = (dm+1)! int_0^1 phi(a)^m da          (verified exactly, m=1,2,3).
```

As `f` runs over homogeneous polynomials of degree `d`, `phi = f(a,1-a) = sum_j a_j a^j(1-a)^{d-j}` runs over
ALL polynomials of degree `<= d` (Bernstein basis), and `phi == 0 <=> f == 0` (since
`f(x,y)=(x+y)^d phi(x/(x+y))`). Hence

> **FC(2)-homogeneous `<=>` the LEBESGUE polynomial moment problem:**
> `int_0^1 phi(a)^m da = 0 for all m >= 1  =>  phi == 0`.

This is Pakovich's **Polynomial Moment Problem** with weight `Q'(a)=1`, `Q(a)=a` -- the PMP that arose in the
center problem for Abel differential equations. The no-pole/Cauchy criterion is transparent here: with
`mu = phi_*(Leb_{[0,1]})` supported on the arc `Gamma = phi([0,1])`, `C_mu(w) = 1/w` on the unbounded
component of `C\Gamma`. If `Gamma` is a SIMPLE arc, `C\Gamma` is connected, `C_mu = 1/w` everywhere, so
`mu = (1/pi) d_wbar (1/w) = delta_0`, forcing `phi == 0`. A counterexample therefore REQUIRES `phi([0,1])` to
self-intersect / enclose area -- and Pakovich's composition condition rules this out for the Lebesgue weight:
a solution needs `phi = P o W`, `Q = R o W` with `W(0)=W(1)`, but `Q(a)=a` is injective, forcing `W` linear,
so `W(0) != W(1)`. **No composition solution exists**, hence FC(2)-homogeneous holds.

**Exact confirmation:** direct symbolic solution of `int_0^1 phi^m = 0` (`m=1..d+2`) for `deg phi = d` returns
NO non-constant solution for **`d = 2, 3, 4`** (extending the known quadratic case).

## 3. Methodological guard (a real trap on this problem)

A NUMERICAL search (double precision) for PMP solutions reported "solutions" at `d = 3,4,5` with residual
`~1e-13`, and even the UNFITTED higher moments looked `~1e-15` -- convincingly fake. High-precision recompute
exposed it: the same `phi` has `|int phi^1| = 5.5e-6` (nonzero), not `1e-12`. The cause is **catastrophic
cancellation**: `int_0^1 phi^m = sum_j [a^j](phi^m)/(j+1)` sums large near-cancelling coefficients, so double
precision manufactures spurious near-zeros that grow with degree (`d=2` stayed at `2.8e-4`, no cancellation;
higher `d` "converged" as cancellation worsened). **Moment computations on this problem must be exact or
high-precision** -- a MISTAKE-class warning for the FC/GMC cluster, since the whole conjecture is a statement
about moments being exactly zero.

## 4. The free-probability bridge (creative tournament analogy)

The Cauchy transform is the same object in three of the repo's frontiers, which lets ideas cross:

- **FC no-pole:** `C_mu(w) = 1/w` (all holomorphic moments vanish) `<=>` `f=0`.
- **Paley tournament (THM-438):** the signed cluster GF `F(x)=sum(-1)^k C_k x^k` solves the loop equation
  `xF^2+F-1=0` -- the **R-transform** of the semicircle, i.e. the free-cumulant inverse of the Cauchy
  transform; `A_L = 0` for odd `L` (from `chi(-1)=-1`) is a SYMMETRY-forced moment vanishing, the tournament
  analog of a PMP symmetry solution.
- **AMM golden floor:** `F(1) = C(-1) = 1/phi` is the edge value of that same Cauchy/R-transform.

So the tournament analogy for FC is precise: **`L(f^m)` is a moment sequence, and FC asks whether its Cauchy
transform can be `1/w` non-trivially -- the exact question the Paley cluster answers with "only the semicircle
survives" (`A_odd=0`, cherry-only).** A creative next probe: build a "tournament" whose arc-character is the
pushforward `phi` of a candidate FC counterexample and read its cluster expansion -- a nonzero surviving
cluster would be a nonzero moment (`f != 0` detector), turning FC(2) into a tournament cluster-survival
statement. The symmetry that kills Paley's odd clusters (`chi(-1)=-1`) is exactly the kind of
measure-preserving involution `sigma` with `f o sigma = -f` that would give FC vanishing for odd `m` -- and
the Paley result that this is NOT enough to kill ALL clusters is the tournament shadow of "no single
root-of-unity symmetry solves the PMP."

## 5. Outlook

Homogeneous FC(2) is the graded core; this reduces it to the (essentially solved) Lebesgue PMP with the
Cauchy-transform / no-pole criterion, exact through degree 4. General FC(2) needs the non-homogeneous mixing
(leading-form + descent, as in the repo's JC(2) rigidity cage), and the Kontsevich-Zagier exponential-period
control of when `int e^{tf}e^{-x-y}` acquires a pole is the tool for the inhomogeneous poles. The clean
takeaway for the friend's approach: **for homogeneous data the K-Z / no-pole route is already a solved moment
problem (Pakovich); the open difficulty is entirely in the inhomogeneous coupling.**
