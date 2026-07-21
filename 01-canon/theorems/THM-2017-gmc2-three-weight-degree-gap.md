---
id: THM-2017
title: "GMC(2) THREE-WEIGHT DEGREE-GAP THEOREM: for P=Z^p a(s)+b(s)+Zbar^q c(s), s=Z Zbar, put g=gcd(p,q), r=(p+q)/g and h=s^(pq/g)a^(q/g)c^(p/g). The exact radial-channel formula is E[P^m]=sum_k m!/((q/g k)!(p/g k)!(m-rk)!) L(h^k b^(m-rk)). If d=deg b and e=deg h satisfy |e-rd|>=r+1, one endpoint channel dominates UNIFORMLY over every k: E[P^m]/L(b^m)->1 when rd-e>=r+1, while on m=rn the all-return channel gives ratio ->1 when e-rd>=r+1. Hence NC2 (and GMC(2)) holds on this infinite-dimensional slice outside the finite resonance band -r<=e-rd<=r. The proof controls k proportional to m and does not separate first-return atoms, which can cancel."
status: >
  PROVED. The charge-balance identity is exact. The asymptotic proof uses a uniform
  mixed-factorial lemma plus a small-channel/linear-channel split; it does not infer
  scalar vanishing term by term. Direct Wick expansion and the channel formula agree
  for m=1..8 in four controls. Boundary cases |e-rd|=r+1 show the asserted ratios
  numerically. This proves a genuine NC2 stratum, not full NC2; the resonance band
  remains open except for previously closed subfamilies such as THM-2014.
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

Choose a small fixed `epsilon>0`. For `1<=k<=epsilon m`, (3), (2),
`m!/(m-rk)!<=m^(rk)`, and
`(dm-Delta k)!/(dm)!=O((c m)^(-Delta k))` give

```text
|T_k/L(b^m)|
  <= C [A m^(r-Delta)]^k / ((q0 k)!(p0 k)!)
  <= C [A/m]^k / ((q0 k)!(p0 k)!).                                    (6)
```

The sum of (6) over all small positive channels is `O(1/m)`. For
`k>epsilon m`, use the crude coefficient bound following (3). Equation (5)
loses at least `Delta epsilon m` units of factorial degree, whereas the
multinomial and all fixed coefficient norms cost only `C_1^m`. Consequently

```text
sum_{k>epsilon m} |T_k/L(b^m)|
  <= exp(-c m log m).                                                    (7)
```

Equations (6)-(7) prove the uniform, all-channel limit

```text
E[P^m] / L(b^m) -> 1.                                                    (8)
```

The gap `r+1` is the exact threshold for this ratio argument: the first
charged channel has multinomial gain `m^r` and factorial loss `m^Delta`.

## 4. The all-return-dominant half

Assume instead

```text
Gamma := e-rd >= r+1.                                                    (9)
```

Restrict to `m=rn` and put `j=n-k`, the number of primitive returns removed
from the endpoint channel. Its degree is

```text
N_j = en-Gamma j.
```

The `j=0` term is

```text
A_n = (rn)! / ((q0 n)!(p0 n)!) * L(h^n),                               (10)
```

and is nonzero for large `n` by (2). For `j<=epsilon n`, the multinomial
ratio to (10) satisfies

```text
(q0 n)_(q0 j) (p0 n)_(p0 j) / (rj)! <= A^j n^(rj)/(rj)!.
```

Combining this with the factorial-degree loss and (3) gives

```text
|T_(n-j)/A_n| <= C [A n^(r-Gamma)]^j/(rj)!
                <= C [A/n]^j/(rj)!.                                    (11)
```

The positive-`j` sum is `O(1/n)`. As before, `j>epsilon n` loses a linear
fraction of `en` and contributes `exp(-c n log n)`. Hence

```text
E[P^(rn)] / A_n -> 1.                                                    (12)
```

## 5. NC2 on the degree-gap stratum

If `a,b,c` are all nonzero and `|e-rd|>=r+1`, (8) or (12) supplies infinitely
many nonzero moments. The degenerate cases are already exact:

- if `a=0` or `c=0` while `b!=0`, charge balance leaves `L(b^m)`, eventually
  nonzero by EMP;
- if `b=0` and `a,c` are both nonzero, the two-weight theorem (THM-1510)
  supplies a nonzero moment;
- if `b=0` and exactly one of `a,c` is nonzero, every moment vanishes by
  strict charge; these, together with `P=0`, are precisely the one-sided NC2
  members.

Thus NC2, and therefore GMC(2), holds throughout this degree-gap slice.

The unresolved degrees form the finite band

```text
-r <= e-rd <= r.                                                        (13)
```

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
