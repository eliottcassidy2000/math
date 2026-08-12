---
source: kind-pasteur-2026-07-24-S163 (Opus 4.8)
status: RESULT (honest frontier) + SYNTHESIS + a confabulation flag. (A) FC(3) exact D=3: characterized the
  wall precisely -- the Nullstellensatz certificate degree is high (>12 even for the PROVED-empty D=2), so the
  Macaulay-certificate route is infeasible and the exact test genuinely needs an F4/F5 engine (none installed);
  the strong evidence (finite-field 0 over 5 good primes + D=2 exact + transversality D<=5) is the achievable
  state. (B) AMM 12592 in tandem: FLAGS that `C_arch = log_5(5 phi^2)` is a CONFABULATION (not in the repo; real
  constant is 9049/6592); shows my artanh decode (kps-S129/S130) IS the AMM 12592 rate certificate; and identifies
  the AMM zeta_6 lemma as the SAME Conway-Jones roots-of-unity governor as FC(3)/series/LRC.
tags: [factorial-conjecture, amm12592, roots-of-unity, conway-jones, confabulation, honest-frontier, synthesis]
related: [kps-S129, kps-S159, kps-S162, THM-2966, THM-3342, HYP-9061, HYP-9023]
---

# (A) FC(3) exact D=3: the wall is tooling, and it is precise

Tried the Nullstellensatz certificate `1 = sum h_j leak_j` via Macaulay linear algebra (flint `nmod_mat`, fast
mod-`p` rank; "is the constant in the row space of `{m*leak_j : deg m <= d}`?"). **Validation:** trivial empty
ideal `(x,x-1)` -> certificate at `d=0` (search correct). **But** the PROVED-empty D=2 ideal (`Groebner=(1)`,
re-confirmed; `{L1,L2}=0` solutions have `|L3|~2.3e5 != 0`) has **certificate degree > 12** -- the ideal has high
regularity. So the certificate matrix explodes before it succeeds, and for D=3 (5 vars) it is hopeless.
> **Honest frontier: exact D=3 requires a real F4/F5 Groebner engine (msolve / Magma / Singular `std`), none
> installed.** sympy Buchberger does not finish (deg-18 in 5 vars); flint exposes fast `nmod_mat` but no Groebner;
> the low-degree Macaulay certificate cannot reach the high regularity. The achievable rigorous state is
> **D=2 exact-empty (kps-S161) + D=3 no-counterexample over 5 good primes (kps-S162) + transversality D<=5
> (kps-S159)** -- strong, converging, not a proof. This is a computational wall, not a mathematical gap.

# (B) AMM 12592 in tandem

## B1. Confabulation flag (honest)
The premise `opus AMM-12592: C_arch = log_5(5 phi^2)` is **not in the repository** (whole-repo grep empty) and is
numerically wrong for this lane: `log_5(5 phi^2) = 1.5980`, whereas the lane's actual constant is
**`C* candidate = 9049/6592 = 1 + 2457/6592 = 1.37272`** (HYP-9061). There is no golden ratio and no log base 5 in
the biased-coin lane. Treat the `log_5(5 phi^2)` claim as confabulated; do not build on it.

## B2. My artanh decode IS the AMM 12592 certificate
The AMM 12592 rate certificate (HYP-9061/HYP-9023) is exactly
`(2457/6592) log(8847357/2974400) - log(1285/896) = 0.045725 > 1/25` (margin 0.00572) -- **verbatim my kps-S129/
S130 "equation (27)"** artanh decode. So the snippet I decoded at the start of this thread (then filed as
"LRC-wider-gap OR irrationality-measure") is the **biased-coin minimal-rate gate**, and my decode work feeds this
lane directly. The "13" link is real: `2457 = 3*sum(k^2, 1..13)` (`= 3*S_2(AP{1..13})`), the same `13` as LRC(14)'s
13 speeds -- the arithmetic is shared, the problem is AMM 12592.

## B3. The zeta_6 lemma = the Conway-Jones governor (shared with FC(3)/series/LRC)
THM-3342 (sublinear excess impossible) turns on the lemma `{z : z and 1-z are both roots of unity} = {zeta_6, zeta_6-bar}`
(roots of `Phi_6(p)=p^2-p+1`, `1-zeta_6 = conj(zeta_6)`; verified). The `p <-> q=1-p` (biased-coin) symmetry forces
the pole set of the window generating function `F(p)=sum_m p^m q W_m(p)` to satisfy `z, 1-z` both roots of unity,
pinning all poles to `zeta_6`, and integrality then contradicts `= 1/2`.
> This is **the same vanishing-sums-of-roots-of-unity (Conway-Jones / Lam-Leung) governor** that runs the FC(3)
> seed (`L((sum zeta^m X_m)^k)=k![m|k]`, kps-S156), the series closed-form locus (kps-S149), and the LRC dichotomy
> (THM-415). Different symmetry, same engine: `p<->1-p` forces `zeta_6` for AMM; the `C_3`-cyclic symmetry forces
> `zeta_3` for FC(3); the modulus forces the Lam-Leung vanishing for LRC.

## B4. The frontier through my lens: rigidity-until-a-threshold
THM-2966: `C* = 1 + gamma*`, `gamma*` = minimal degree-growth of the box polynomials `W_m` (degree `d_m=(C-1)m+D-1`)
solving the spine identity `sum_m p^m q W_m + sum_m q^m p V_m = 1/2`. The `C=1` lane (`d_m` **bounded**) is rigid:
`F(p)` is radius-1 with bounded-degree numerators -> Polya-Carlson forces rational -> poles at roots of unity ->
`zeta_6` -> contradiction (THM-3342). Linear excess (`d_m` **growing**) unlocks flexibility (no forced unit-circle rationality).
> **This is the exact "capacity holds until a threshold" shape** seen across the cluster (FC(2)-true / FC(3)-open;
> the series elementary-until-`k=P(n)`; the `n=2 | n>=3` wall). AMM's threshold is at the growth rate `gamma=0`
> (bounded degree = rigid); the open question is the *minimal* `gamma*>0` -- the same "how much capacity does the
> constraint force" question as the FC(3) transversality capacity.

## B5. Open frontier (from the lane) + a concrete handle
Most concrete open sub-problem: **is `T(4)=5` attainable?** (`T(4) in {5,6}`, the first open envelope value) --
a finite rational-LP/Farkas feasibility on the truncated spine identity (THM-2966 box constraints, length `<=24`).
The certificate rate `9049/6592` is claimed but its *status* (lower bound vs construction) is **open/unverified**.
My roots-of-unity lens suggests the C=1->C>1 rigidity transition (B4) is where a lower bound on `gamma*` should
come from: an Erdos-Turan / cyclotomic-nonvanishing bound on how fast `W_m` degree must grow to break the
`zeta_6` pole-rigidity -- the analogue of the FC(3) transversality "capacity=P" forcing.

Files: `/tmp/{debug,calib2}.py` (FC3 wall), AMM verify inline. Builds on kps-S129/S159/S161/S162; engages
THM-2966/2967, HYP-9061/9023.
