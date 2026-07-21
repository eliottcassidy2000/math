---
id: THM-2017
title: "GMC(2) THREE-WEIGHT DEGREE-GAP + BOUNDARY THEOREM: for P=Z^p a(s)+b(s)+Zbar^q c(s), the exact primitive-return channels have affine factorial degree D(k)=dm+(e-rd)k. If |e-rd|>=r+1, one endpoint dominates uniformly over all channels. At the sharp boundary |e-rd|=r, the complete boundary layer instead converges to an explicit generalized hyper-Bessel function; NC2 holds whenever that value is nonzero, hence generically on the boundary. The proof controls channels proportional to m and permits first-return cancellation."
status: >
  PROVED. The charge-balance identity is exact. The asymptotic proof uses a uniform
  mixed-factorial lemma plus a small-channel/linear-channel split; it does not infer
  scalar vanishing term by term. Dominated convergence at |e-rd|=r gives the stated
  hyper-Bessel limits. Direct Wick expansion and the channel formula agree for m=1..8
  in six controls; strict and sharp boundary limits are sampled. This proves a genuine
  NC2 stratum, not full NC2. The inner resonance band and the discrete zero loci of the
  two boundary functions remain open except for subfamilies such as THM-2014.
source: codex-2026-07-21-gmc2-degree-gap
depends_on:
  - THM-1510  # one-variable factorial asymptotic / two-weight NC2
  - THM-1540  # NC2 implies GMC(2)
  - THM-2014  # constant-endpoint p=q=1 model, including its resonant d=0,1 closure
related:
  - THM-1775  # weighted toral/radial channel formulation
  - THM-1765  # two-straddle tower
  - HYP-8765  # radial-channel tower from THM-2014
  - HYP-8766  # close the finite degree-resonance band
script: 04-computation/gmc2_three_weight_degree_gap_codex_20260721.py (+ .out)
---

# THM-2017 — the three-weight degree-gap theorem

Let `Z` be a circular complex Gaussian, write `Zbar` for its conjugate and
`s=Z Zbar`. Thus

```text
E[Z^A Zbar^B] = A! if A=B, and 0 otherwise,
L(f) := E[f(s)],                 L(s^N)=N!.
```

Fix positive integers `p,q` and nonzero polynomials `a,b,c in C[s]`. Set

```text
g  = gcd(p,q),       p0=p/g, q0=q/g,       r=p0+q0=(p+q)/g,
h  = s^(pq/g) a^q0 c^p0,
d  = deg b,          e=deg h.
```

The polynomial under study is

```text
P = Z^p a(s) + b(s) + Zbar^q c(s).
```

## 1. Exact primitive-return decomposition

In a term of `P^m`, let `i,j,l` count the charge `+p`, charge `-q`, and
charge-zero factors. Wick survival requires `pi=qj`. Since `(p0,q0)=1`,

```text
i=q0 k,        j=p0 k,        l=m-rk
```

for a unique `0<=k<=floor(m/r)`. The product of the charged factors is

```text
(Z^p a)^(q0 k) (Zbar^q c)^(p0 k)
  = [s^(pq/g) a^q0 c^p0]^k = h^k.
```

Therefore the moment has the exact all-channel expansion

```text
E[P^m] = sum_{0<=k<=m/r}
          m! / ((q0 k)! (p0 k)! (m-rk)!) * L(h^k b^(m-rk)).       (1)
```

This is a sum of scalar complex numbers. Different `k`, and even different
primitive-return atoms inside a fixed moment, can cancel. Nothing below treats
the channels as algebraically independent.

## 2. The uniform mixed-factorial lemma

We use two standard consequences of the factorial functional.

**One-polynomial asymptotic.** If
`f(s)=alpha s^D+alpha_1 s^(D-1)+...` with `D>=1`, then

```text
L(f^n) = alpha^n (Dn)! exp(alpha_1/(D alpha)) (1+O(1/n)).          (2)
```

In particular this is nonzero for all sufficiently large `n`. This is the EMP
asymptotic of THM-1510 (and the analytic input used uniformly in THM-2014).

**Mixed upper bound.** Fix finitely many nonzero polynomials `f_i`, of degrees
`D_i` and leading coefficients `alpha_i`. For nonnegative exponent vectors
`n_i`, put `N=sum D_i n_i`. If `sum n_i=O(M)` and `N>=eta M` for a fixed
`eta>0`, then

```text
|L(prod_i f_i^n_i)| <= C (prod_i |alpha_i|^n_i) N!,              (3)
```

where `C` is independent of the exponent vector.

Here is a proof, included because uniformity in channels is load-bearing. Write

```text
f_i(s)=alpha_i s^D_i F_i(1/s),       F_i(0)=1,
prod F_i(t)^n_i = sum_j A_j t^j.
```

After dividing by `(prod alpha_i^n_i)N!`, the absolute factorial sum is at
most `sum_j |A_j|/(N)_j`. For `j<=N/2`, replace `(N)_j` by `(N/2)^j`; the
sum is bounded by the coefficient-majorant product
`prod |F_i|_maj(2/N)^n_i=exp(O(M/N))=O(1)`. For `j>N/2`, the total coefficient
mass is at most `C_0^M`, while `(N-j)!/N!` is at most
`(ceil(N/2))!/N!=exp(-Omega(N log N))`. This proves (3). If `N` is not
linear in `M`, the cruder bound `|L(prod f_i^n_i)|<=C_0^M N!` will suffice.

## 3. The `b`-dominant half

Assume

```text
Delta := rd-e >= r+1.                                                   (4)
```

The degree in channel `k` is

```text
N_k = e k + d(m-rk) = dm-Delta k.                                      (5)
```

The `k=0` term of (1) is `L(b^m)`, which by (2) is a nonzero constant times
`beta^m(dm)!`, where `beta` is the leading coefficient of `b`.

Condition (4) forces `d>=1`. Fix
`0<epsilon<d/(2 Delta)`, so `N_k>=dm/2` throughout the small-channel
range. For `1<=k<=epsilon m`, (3), (2),
`m!/(m-rk)!<=m^(rk)`, and
`(dm-Delta k)!/(dm)!=O((c m)^(-Delta k))` give

```text
|T_k/L(b^m)|
  <= C [A m^(r-Delta)]^k / ((q0 k)!(p0 k)!)
  <= C [A/m]^k / ((q0 k)!(p0 k)!).                                    (6)
```

The sum of (6) over all small positive channels is `O(1/m)`. For
`k>epsilon m`, use the crude coefficient bound following (3). The trinomial
coefficient is at most `3^m`; the coefficient `l1` norm of
`h^k b^(m-rk)` is at most `C_0^m`; and normalizing by the fixed nonzero
leading coefficients adds only another `C_1^m`. Thus every non-factorial
quantity costs `exp(O(m))`. Equation (5) loses at least
`Delta epsilon m` units of factorial degree, giving

```text
sum_{k>epsilon m} |T_k/L(b^m)|
  <= exp(-c m log m).                                                    (7)
```

Equations (6)-(7) prove the uniform, all-channel limit

```text
E[P^m] / L(b^m) -> 1.                                                    (8)
```

The gap `r+1` is the exact threshold for a ratio-to-one argument: the first
charged channel has multinomial gain `m^r` and factorial loss `m^Delta`.

### The sharp boundary `Delta=r`

The same proof gives a nontrivial limit at equality. Let `alpha,beta` be the
leading coefficients of `h,b` and define

```text
Phi_(p0,q0)(x) = sum_{k>=0} x^k/((q0 k)!(p0 k)!),
xi = alpha/(beta^r d^r).
```

For each fixed `k`, the `k`-channel divided by `L(b^m)` tends to
`xi^k/((q0k)!(p0k)!)`: the multinomial gain `m^(rk)` exactly cancels the
factorial loss `(dm)^(-rk)`. Estimate (6), now without its factor `m^-k`, is
a summable majorant independent of `m`; (7) still removes linear channels.
Dominated convergence therefore proves

```text
E[P^m]/L(b^m) -> Phi_(p0,q0)(xi)             when rd-e=r.              (9)
```

The entire function `Phi` has value one at zero, so it is not identically
zero and its zeros are discrete. Thus the boundary is NC2-clear whenever
`Phi_(p0,q0)(xi)!=0`, in particular for generic leading coefficients.

## 4. The all-return-dominant half

Assume instead

```text
Gamma := e-rd >= r+1.                                                   (10)
```

Condition (10) forces `e>=1`. Restrict to `m=rn` and put `j=n-k`, the
number of primitive returns removed from the endpoint channel. Its degree is

```text
N_j = en-Gamma j.
```

The `j=0` term is

```text
A_n = (rn)! / ((q0 n)!(p0 n)!) * L(h^n),                               (11)
```

and is nonzero for large `n` by (2). Fix
`0<epsilon<e/(2 Gamma)`, so `N_j>=en/2` in the small-channel range. For
`j<=epsilon n`, the multinomial
ratio to (11) satisfies

```text
(q0 n)_(q0 j) (p0 n)_(p0 j) / (rj)! <= A^j n^(rj)/(rj)!.
```

Combining this with the factorial-degree loss and (3) gives

```text
|T_(n-j)/A_n| <= C [A n^(r-Gamma)]^j/(rj)!
                <= C [A/n]^j/(rj)!.                                   (12)
```

The positive-`j` sum is `O(1/n)`. For `j>epsilon n`, all falling-factorial
multinomial ratios and fixed polynomial coefficient norms are `exp(O(n))`,
while the loss of at least `Gamma epsilon n` factorial degrees is
`exp(-Omega(n log n))`. Thus the tail contributes `exp(-c n log n)`. Hence

```text
E[P^(rn)] / A_n -> 1.                                                   (13)
```

### The sharp boundary `Gamma=r`

At equality, define

```text
Psi_r(y) = sum_{j>=0} y^j/(rj)!,
eta = (q0^q0 p0^p0/e^r) * (beta^r/alpha).
```

For fixed `j`, the falling-factorial multinomial ratio contributes
`q0^(q0j)p0^(p0j)n^(rj)/(rj)!`, the factorial degree contributes
`(en)^(-rj)`, and the leading coefficients contribute
`beta^(rj)/alpha^j`. The bound (12) at equality is again summable in `j`,
and linear `j` remain superfactorially small. Hence

```text
E[P^(rn)]/A_n -> Psi_r(eta)                  when e-rd=r.              (14)
```

Again `Psi_r(0)=1`, so its zero set is discrete and the boundary is
NC2-clear whenever `Psi_r(eta)!=0`. Equivalently,
`Psi_r(y)=r^-1 sum_(omega^r=1) exp(omega y^(1/r))`; the displayed power
series makes this identity branch-independent.

### The exceptional boundary zeros do not survive in the symmetric monomial model

There is one useful case where the discrete exception can be removed entirely.
Assume `p0=q0=1` (so `r=2`) and both `b=beta s^d` and `h=alpha s^e`
are monomials.

On the neutral boundary `e=2d-2`, the exact normalized `k`-term is

```text
xi^k/(k!)^2 * [1-(1-1/d)(2k)(2k-1)/(2m)+O_k(m^-2)],
xi=alpha/(beta^2 d^2).
```

The factorial majorant above permits termwise summation of this expansion.
With `theta=xi d/dxi` and `Phi(xi)=sum xi^k/(k!)^2`,

```text
E[P^m]/L(b^m)
 = Phi(xi) -(1-1/d)(4 theta^2-2 theta)Phi(xi)/(2m)+O(m^-2).           (15)
```

The Bessel equation is `theta^2 Phi=xi Phi`. If `Phi(xi)=0`, then
`xi!=0` and (15) reduces to

```text
E[P^m]/L(b^m) = ((d-1)/d) theta Phi(xi)/m+O(m^-2).
```

Here `d>=2` (because `h` is a nonzero polynomial on this boundary), and
`theta Phi(xi)!=0`: otherwise the second-order ODE with zero value and
derivative at the ordinary point `xi` would force `Phi` to be identically
zero. Thus the supposedly exceptional leading zero is detected one order
later.

On the all-return boundary `e=2d+2`, write
`eta=beta^2/(alpha e^2)` and `Psi(eta)=sum eta^j/(2j)!=cosh(sqrt(eta))`.
The exact normalized removed-return term gives

```text
E[P^(2n)]/A_n
 = Psi(eta) + [(-1+2/e)theta^2+(1-1/e)theta]Psi(eta)/n+O(n^-2).       (16)
```

Now `theta^2 Psi=(theta Psi)/2+(eta Psi)/4`. At a zero of `Psi`, the
`1/n` coefficient in (16) is exactly `(theta Psi)/2`, which is nonzero by
the same ODE-uniqueness argument. Consequently **both sharp boundaries are
fully NC2-clear for symmetric primitive charges with monomial `b,h`**, even
when the leading hyper-Bessel limit vanishes.

## 5. NC2 on the degree-gap stratum

The degree gate is asserted only when `a,b,c` are all nonzero, so `d,e` are
defined. In that case, if `|e-rd|>=r+1`, (8) or (13) supplies infinitely
many nonzero moments. The degenerate cases are already exact:

- if `a=0` or `c=0` while `b!=0`, charge balance leaves `L(b^m)`. This is
  eventually nonzero by EMP when `deg b>=1`, and is exactly `b^m!=0` when
  `b` is a nonzero constant;
- if `b=0` and `a,c` are both nonzero, the two-weight theorem (THM-1510)
  supplies a nonzero moment;
- if `b=0` and exactly one of `a,c` is nonzero, every moment vanishes by
  strict charge; these, together with `P=0`, are precisely the one-sided NC2
  members.

These `a=0`, `b=0`, or `c=0` cases are boundary cases of the slice, not
instances of a degree inequality with an artificially assigned degree to the
zero polynomial.

Thus NC2, and therefore GMC(2), holds throughout the strict degree-gap slice.
Equations (9) and (14) also prove it generically on both sharp boundaries;
(15)-(16) close both boundaries completely in the symmetric monomial model.

The unresolved degrees form the finite band

```text
-r <= e-rd <= r.                                                       (17)
```

The endpoints in (17) are now reduced further to the discrete exceptional
sets `Phi(xi)=0` and `Psi(eta)=0`; only the inner `2r-1` offsets require a
new saddle regime for every leading coefficient.

THM-2014 is the important constant-endpoint model inside and across this
picture: for `p=q=1` and constant `a,c`, it proves the full slice for arbitrary
`b`, treating `deg b>=2` by a uniform channel estimate and the resonant degrees
zero and one by exact generating functions. The present theorem extends the
uniform mechanism to arbitrary charges and arbitrary radial endpoints, while
making the residual resonance band explicit.

## Verification and challenged assumptions

`04-computation/gmc2_three_weight_degree_gap_codex_20260721.py` verifies (1)
against direct Wick expansion for `m=1..8`, samples both boundary gaps, and
records the required Tournament Analysis. Candidate vertex sets were
monomials, charges, primitive atoms, radial-deficit channels, and proof
obligations. The chosen vertices are channels `k`; their pairwise observable is
`D(k)-D(l)` for `D(k)=dm+(e-rd)k`, gauged by its sign, with increasing `k` as
the tie Hamiltonian path. These tournaments are transitive (zero directed
3-cycles, singleton SCCs, one Hamiltonian path). The quotient preserves the
factorial slope used in the proof but destroys coefficient phases and
within-channel cancellations. That destruction is why it is a diagnostic,
not a proof by itself.

The challenged assumption is explicit: a scalar moment does **not** remember
which first-return atom produced it. Cancellation between different atoms is
possible, so the proof compares the complete channel sum to one nonzero
endpoint rather than declaring its summands separately zero.
