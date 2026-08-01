---
source: opus-2026-08-01-S4 (general reformulation of the exp-integral transcendence conjecture, all degrees)
status: >
  STRUCTURAL REFORMULATION + a SELF-CORRECTION (degree<=2 closed independently; the general "easy step" was
  NOT easy). The fixed integral I = int_0^1 e^{Q(t)} dt is the value at z=1 of Phi_Q(z) := int_0^1 e^{z Q(t)} dt,
  and Phi_Q is an E-FUNCTION in the parameter z for EVERY Q in Qbar[t] (Taylor coefficients mu_m = int_0^1 Q^m dt:
  algebraic, geometric growth |mu_m|<=max|Q|^m, geometric denominators) -- this part is correct and unconditional.
  CORRECTION (thanks to death-star THM-3031, 2026-08-01): my first draft claimed "Phi_Q not in Qbar(z) for
  nonconstant Q" as a free lemma via an indicator/Laplace argument. THAT IS WRONG. Phi_Q == 1 IDENTICALLY exactly
  when all moments int_0^1 Q^m=0 (m>=1), i.e. exactly when Q is an FC(2) counterexample -- so "Phi_Q non-constant"
  is EQUIVALENT to (homogeneous) FC(2), not a free lemma. The indicator only bounds |Phi_Q| ABOVE; g(t)=e^{2 pi i t}
  is nonconstant with Phi_g==1 by total cancellation (verified). What remains valid and useful: (i) Phi_Q is an
  E-function; (ii) I=Phi_Q(1); (iii) IF Phi_Q is non-rational then (Beukers 2006, refined LINEAR Siegel-Shidlovskii,
  no exceptional point) I is transcendental as soon as {1,Phi_Q,...,Phi_Q^{(r-1)}} are Qbar(z)-linearly independent
  = the horizontal-endomorphism RIGIDITY (1 not a horizontal section of the differential module). So the E-function
  device converts BOTH difficulties (non-rationality AND rigidity) into the differential module of Phi_Q -- it does
  NOT make either free. Degree<=2 is closed by the INDEPENDENT explicit 1F1 proof (companion file), which does not
  use the flawed indicator step and is unaffected. death-star's operative point stands: the FC(2) bridge needs
  "value != 1" (a counterexample pins the period to 1), and transcendence delivers that.
tags: [factorial-conjecture, FC2, exponential-integral, transcendence, e-functions, beukers, siegel-shidlovskii, rigidity, critical-values, saddle, parameter-e-function, all-degrees]
related: [degree-2-exp-integral-is-transcendental-rigorous-via-linear-beukers-shidlovskii, HYP-9078, exp-integral-radial-bridge-gives-fc-homog-and-the-weighted-generalization-for-full-fc2, fc3-decisive-max-ridge-is-the-fold-equals-kps-rank-drop-locus]
external: ["Beukers, Ann. of Math. 163 (2006) 369-379 (refined linear Siegel-Shidlovskii)", "Shidlovskii 1989", "Andre, 'Series Gevrey de type arithmetique' (E-functions)"]
---

# The exp-integral is an E-function in the parameter: a general reformulation, and rigidity = critical-value resonance

The degree-2 theorem (companion file) proves `int_0^1 e^{alpha t^2+beta t+gamma} dt` transcendental via the
confluent E-function `1F1(1/2;3/2)`. That route is degree-specific. Here is the **degree-free** reformulation it
generalizes to, with two of its three ingredients proved *unconditionally for all degrees*, isolating the whole
remaining difficulty as one clean rigidity statement.

## 1. The parameter device: `I = Phi_Q(1)` with `Phi_Q` an E-function

For `Q in Qbar[t]` set
```
        Phi_Q(z) := int_0^1 e^{z Q(t)} dt = sum_{m>=0} mu_m z^m/m!,     mu_m := int_0^1 Q(t)^m dt.
```
Then `I = int_0^1 e^{Q} dt = Phi_Q(1)`. The point is that **`Phi_Q` is an E-function in `z`**:

- **Algebraic coefficients.** `mu_m = int_0^1 Q^m dt` is a `Qbar`-linear combination of `int_0^1 t^k dt = 1/(k+1)`,
  hence lies in the fixed number field `K = Q(coeffs of Q)`. (Verified: for `Q=(7/10)t^2+(9/10)t-1/5`,
  `mu_0..4 = 1, 29/60, 1349/3000, 9123/20000, 452467/900000` -- all rational.)
- **Geometric growth.** `|mu_m| <= int_0^1 |Q|^m <= (max_{[0,1]}|Q|)^m`; every Galois conjugate obeys the same
  bound with `Q` replaced by `Q^sigma`, so the house is `O(C^m)`.
- **Geometric denominators.** `den(mu_0,...,mu_m) | lcm(1,...,dm+1) * den(Q)^m = e^{O(m)}` (prime number theorem
  for the `1/(k+1)`; `den(Q)^m` for the coefficients of `Q^m`).

So `Phi_Q` satisfies Siegel's E-function conditions for every `Q`. The fixed integral `I` is thereby an
**E-function value at the algebraic point `z=1`** -- the setting where Siegel-Shidlovskii theory speaks.

## 2. `Phi_Q` is not rational (nonconstant `Q`) -- this is NOT free; it IS the FC(2) content [CORRECTED]

**My first draft was wrong here, and the error is instructive.** I argued by Laplace that
`log|Phi_Q(R e^{i theta})| = R h(theta) + o(R)` with `h(theta)=max_t Re(e^{i theta}Q(t))`, concluding `Phi_Q`
grows (non-polynomial) whenever `Q` is nonconstant. **The Laplace equality is false in general**: for a complex
integrand the oscillation can cause total cancellation, and the indicator only gives the UPPER bound
`|Phi_Q(z)| <= int_0^1 e^{Re(z Q)} dt = e^{R h(theta)+o(R)}`, never a lower bound.

**Decisive counterexample (death-star THM-3031 phenomenon).** `Phi_Q` is CONSTANT `== 1` precisely when
`mu_m = int_0^1 Q^m dt = 0` for all `m >= 1` -- i.e. precisely when `Q` is an FC(2) moment-problem counterexample.
Take `g(t) = e^{2 pi i t}` (nonconstant): `int_0^1 g^m dt = 0` for every `m != 0`, so `Phi_g(z) = int_0^1 e^{z g} dt
= sum_m z^m/m! * 0 = 1` IDENTICALLY -- verified numerically to 24 digits at several `z`, even while
`int_0^1 |e^{z g}| dt` grows like `e^{|z|}`. (`g` is not a polynomial, but it exhibits exactly the cancellation my
indicator argument ignored; whether a POLYNOMIAL `g` can do the same is FC(2) itself.)

**Consequence.** `Phi_Q notin Qbar(z)` (equivalently `Phi_Q` non-constant, equivalently `mu_m != 0` for some `m>=1`)
is **EQUIVALENT to the homogeneous FC(2) moment problem**, not a free lemma. The E-function device (S1) is real,
but it does **not** dispose of non-rationality for free: both non-rationality AND the rigidity of S3 are packaged
into the differential module of `Phi_Q`. For REAL `Q` non-constancy is still trivial (`mu_{2m}=int Q^{2m}>0`), and
for `deg Q <= 2` it holds because a degree-2 polynomial maps `[0,1]` to a non-closed parabolic ARC (its pushforward
cannot be an all-moments-zero measure) -- but for general complex `Q` it is exactly the open problem. `//`

## 3. Transcendence <= rigidity (Beukers, no exceptional point)

Let `r` be the order of the minimal linear ODE of `Phi_Q` over `Qbar(z)` and consider the E-functions
`1, Phi_Q, Phi_Q', ..., Phi_Q^{(r-1)}` (all E-functions: derivatives and constants of E-functions are
E-functions).

**Refined linear Siegel-Shidlovskii (Beukers 2006).** If E-functions `f_0,...,f_n` are `Qbar(z)`-linearly
independent, then `f_0(xi),...,f_n(xi)` are `Qbar`-linearly independent for *every* nonzero algebraic `xi` that
is not a singularity of the system -- with **no exceptional set**.

Apply it at `xi = 1` (the only finite singularity of the natural system is where `Q'` or the ODE degenerates;
`z=1` is regular). If `{1, Phi_Q, ..., Phi_Q^{(r-1)}}` are `Qbar(z)`-linearly independent, then
`{1, Phi_Q(1), ...}` are `Qbar`-linearly independent; in particular `Phi_Q(1) = I notin Qbar`, i.e. **`I` is
transcendental**. The hypothesis

```
   (RIGIDITY)   1, Phi_Q, Phi_Q', ..., Phi_Q^{(r-1)}  are linearly independent over Qbar(z)
              <=>  Phi_Q satisfies NO inhomogeneous ODE  L Phi_Q = c_0  (c_0 in Qbar(z), constant term present)
              <=>  1 is not a horizontal section of the differential module M(Phi_Q).
```

is the **horizontal-endomorphism rigidity** -- but note (per S2's correction) that `{1, Phi_Q, ...}` linear
independence subsumes BOTH `Phi_Q notin Qbar(z)` (the FC(2)-equivalent non-rationality) AND the inhomogeneous
refinement that `1` is not dragged in. Neither half is free in general; the device localizes both into the one
differential module of `Phi_Q`.

## 4. What rigidity means geometrically: critical-value non-resonance

Integrating `d/dt e^{zQ} = z Q'(t) e^{zQ}` by parts repeatedly expresses the module `M(Phi_Q)` through the
**boundary and saddle exponentials**
```
        e^{z Q(0)},   e^{z Q(1)},   and   e^{z Q(c)}  for interior critical points  Q'(c)=0.
```
(For degree 2 the relation is exactly `e^{zQ(1)} - e^{zQ(0)} = z(2 alpha Phi_1 + beta Phi_0)`, coupling `Phi`
to the two endpoint exponentials with rates `Q(0)=gamma`, `Q(1)=alpha+beta+gamma`, and the single saddle
`Q(-beta/2alpha)=delta`.) The constant function `1 = e^{z * 0}` is the exponential of **rate 0**. It becomes a
horizontal section -- RIGIDITY fails -- precisely when `0` is *resonant* with the critical-value set: when a
`Qbar(z)`-relation among the `e^{z Q(c)}` produces the rate-0 line. Generically the critical values
`{Q(0), Q(1)} cup {Q(c): Q'(c)=0}` avoid `0` and are `Q`-linearly unrelated, and RIGIDITY holds.

So the transcendence conjecture becomes, cleanly (modulo the S2 correction -- this describes the RIGIDITY half;
non-rationality must also hold, and by S2 the two are entangled in the same module):

> **`int_0^1 e^{Q} dt` is transcendental once `{1, Phi_Q, ..., Phi_Q^{(r-1)}}` are `Qbar(z)`-independent**, and
> the obstruction is a resonance of rate `0` with the critical values `{Q(c)}` of `Q`.

This is the precise home of the friend's "decisive horizontal-endomorphism rigidity theorem": that the
endpoint/saddle exponential connection has no horizontal `1`, i.e. minimal endomorphism algebra. For `deg Q <= 2`
we *prove* it (companion file: the two critical-value exponentials `e^{z eta_i}` have distinct rates, which
simultaneously gives non-rationality AND kills any relation yielding rate 0 -- so BOTH halves close at once). For
`deg Q >= 3` there are up to `d-1` interior saddles and the resonance lattice among `d+1` critical values is the
object to control -- the genuine remaining work, and (per S2) it is not lighter than FC(2) itself.

## 5. Status ledger (honest, after the S2 correction)

- **Proved, all degrees:** `Phi_Q(z)=int_0^1 e^{zQ}dt` is an E-function (S1). This is the one genuinely free,
  degree-free ingredient -- the useful device.
- **NOT free (retracted):** `Phi_Q notin Qbar(z)` for nonconstant `Q`. My S2 indicator "proof" was wrong (upper
  bound only; `g=e^{2 pi i t}` gives `Phi_g==1`). This statement is EQUIVALENT to homogeneous FC(2). It IS free for
  real `Q` (`mu_{2m}>0`) and for `deg<=2` (parabolic arc), but not in general.
- **Proved, `deg Q <= 2`:** `I` transcendental (companion file; `1F1(1/2;3/2)` distinct-rate argument -- does
  NOT use S2's flawed step, so unaffected). Degree 1 = Lindemann-Weierstrass (HYP-9078); degree 2 closed.
- **Reduced, all degrees:** the full conjecture `<=` "`{1, Phi_Q, ..., Phi_Q^{(r-1)}}` linearly independent over
  `Qbar(z)`" (S3), which packages BOTH non-rationality AND rigidity. Not a shortcut past FC(2) -- a relocation of
  it into one differential module, which is still the honest content.
- **death-star THM-3031 (adopted):** the FC(2) bridge needs "value `!= 1`", NOT "`!= 0`" -- a counterexample pins
  the period to `1`. Transcendence (deg<=2, proved) gives `!= 1` for free since `1` is rational; that is why the
  transcendence half is the operative one.
- **Not claimed:** RIGIDITY / non-rationality for `deg Q >= 3`; and the inhomogeneous `poly+rational` exponent
  that FULL FC(2) needs (saddle weight `exp(phi_{D-1}/D phi_D)` -- a rational term in the exponent is outside the
  pure E-function frame of S1).
- **Cross-checks:** `I=Phi_Q(1)` matches `sum mu_m/m!` to 15 digits; the golden constant
  `C=1+2 log(phi)/log5=1.59798...` (klein THM-3027 = opus two-ray) is the AMM neighbor, not used here.

## 6. Why this is the right frame

The device "`fixed integral = E-function value in a scaling parameter`" converts a hard transcendence question
into the Siegel-Shidlovskii world *without* needing a hypergeometric closed form. It works for `int_0^1 e^{zQ}`,
and the same move (`Phi(z) = int_gamma e^{z Q}` over any fixed algebraic cycle `gamma`) turns a whole class of
exponential periods into E-function values. What the device does NOT do -- and my first draft wrongly said it did
-- is make non-rationality free: by S2 that half is FC(2) itself. So the honest summary is that the reformulation
**localizes the entire difficulty (non-rationality + rigidity) into the single differential module of `Phi_Q`**,
and the degree-2 proof shows how, for a genuine sub-family, distinct critical-value rates close both halves at
once. The lesson (again): a complex oscillatory integral can cancel to a constant while its absolute value grows
-- `|int| != int|.|` -- so "it obviously grows, hence isn't rational" is never a proof. death-star THM-3031 caught
exactly this, and it is the same trap the moment-problem catastrophic-cancellation guard flags on the FC side.
