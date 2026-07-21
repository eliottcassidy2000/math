---
id: THM-2018
title: "GMC(2) EXACT RESONANCE CLOSURES: the affine/quadratic p=q=1 band has a Catalan-algebraic EGF whose nullity would make the exponential of a nonconstant algebraic germ algebraic; for arbitrary charges, the central hypersurface h=kappa*b^r factors into one radial power times a hyper-Bessel coefficient sequence that is not eventually zero"
status: >
  PROVED. Part A is an exact coefficient identity, valid for arbitrary complex
  coefficients: if p=q=1, deg b<=1, and deg(a*c)<=1, then all Gaussian moments
  vanish exactly on the one-sided locus b=0 and a*c=0. Part B is an
  arbitrary-charge, arbitrary-radial-degree theorem: on h=kappa*b^r, all moments
  vanish exactly when b=h=0. The only imported asymptotic fact is the one-variable
  EMP theorem that L(b^m) is eventually nonzero. Exact direct-Wick, channel,
  Catalan-coefficient, closed-EGF, and proportional-factor checks pass. These are
  genuine resonance-band strata. Full NC2 is now supplied separately by
  THM-2022's Frobenius lowest-face theorem.
source: codex-2026-07-21-NC2-followup
depends_on:
  - THM-1510  # one-variable EMP / eventual nonvanishing of L(b^m)
  - THM-1540  # NC2 implies GMC(2)
  - THM-2017  # three-weight primitive-return notation and channel identity
related:
  - THM-2014
  - THM-2022
  - HYP-8766
  - HYP-8769
script: 04-computation/gmc2_resonance_closures_codex_20260721.py (+ stored .out)
---

# THM-2018 -- two exact closures inside the three-weight resonance band

Let `Z` be a circular complex Gaussian, put `W=Zbar` and `s=ZW`, and write

```text
E[Z^A W^B] = A! if A=B, and 0 otherwise,
L(f) := E[f(s)],                         L(s^N)=N!.
```

This theorem closes two loci on which the degree tournament of THM-2017 is tied.
The mechanisms are exact resummation and symmetry, not endpoint domination.

## 1. A small lemma: an algebraic germ cannot have algebraic exponential

We will use the following standard fact in a form that includes all hypotheses.

**Algebraic-exponential lemma.** Let `A(t)` be a complex algebraic function germ,
holomorphic at `t=0`. If `exp(A(t))` is also algebraic over `C(t)`, then `A` is
constant.

**Proof.** Take a connected compact Riemann surface `X` whose function field
contains both algebraic germs, and continue the germ identity to `X` away from
the finitely many poles and branch points. The pullback of `A` is meromorphic on
`X`. If it had a pole, then in a local parameter `u` it would have a nonzero
principal part `c u^(-d)+...`; its exponential would have an essential
singularity there. An algebraic function on `X` is meromorphic and cannot have
an essential singularity. Hence `A` has no poles. A holomorphic function on the
compact connected surface `X` is constant. QED.

The same proof applies if `exp(A)` is merely a nonzero rational function. This
is the exact no-common-zero replacement used in both parts below.

## 2. Part A: the full affine/quadratic `p=q=1` resonance slice

Consider

```text
P = Z a(s) + (b0+b1 s) + W c(s),
a(s)c(s) = delta + alpha s,             b0,b1,delta,alpha in C.       (1)
```

Thus `b=b0+b1s` and the primitive-return polynomial is

```text
h=s a(s)c(s)=delta s+alpha s^2.
```

The conclusion is

```text
E[P^m]=0 for every m>=1
       iff b0=b1=delta=alpha=0
       iff b=0 and a*c=0.                                      (2)
```

Because `C[s]` is a domain, `a*c=0` means `a=0` or `c=0`; the right side of
(2) is exactly the one-sided locus in this slice.

### 2.1 Exact Wick and EGF identities

Charge balance gives, without separating scalar summands,

```text
M_m:=E[P^m]
   = sum_(0<=k<=m/2) m!/(k!^2 (m-2k)!) L(h^k b^(m-2k)).          (3)
```

For sufficiently small `|t|`, the exponential generating function converges
absolutely: `b(s)=O(s)` and `sqrt(h(s))=O(s)`, so the defining radial integrand
is bounded by `exp(-(1-C|t|)s)`. Termwise integration in the Bessel series gives

```text
F(t):=sum_(m>=0) M_m t^m/m!
     =L(exp(t b) I_0(2t sqrt(h))).                               (4)
```

We now derive (4) coefficient by coefficient, so no contour or branch identity
is being assumed. Put

```text
p=1-b1 t.
```

Expanding the `k`-th Bessel term and choosing `j` of its `k` factors from
`alpha s^2` gives

```text
L(exp(b1 t s) h^k)
 =sum_(j=0)^k binom(k,j) alpha^j delta^(k-j)
                    (k+j)!/p^(k+j+1).                            (5)
```

Indeed the chosen monomial is `s^(k+j)`, and
`integral_0^infinity exp(-p s)s^(k+j)ds=(k+j)!/p^(k+j+1)` near
`t=0`. Set `n=k-j` and

```text
x=alpha t^2/p^2,                 y=delta t^2/p.
```

Substitution of (5) into (4) gives every factorial explicitly:

```text
exp(-b0 t)F(t)
 =1/p sum_(n,j>=0) (n+2j)!/((n+j)!j!n!) x^j y^n.                 (6)
```

Let

```text
C(x)=(1-sqrt(1-4x))/(2x)=1+x+2x^2+... ,                         (7)
```

with `C(0)=1`. The elementary Catalan coefficient identity

```text
[x^j] C(x)^n/sqrt(1-4x)=binom(n+2j,j)                            (8)
```

follows directly from Lagrange--Buermann applied to
`C=1+xC^2` (and includes `n=0`, the central-binomial series). Since
`(n+2j)!/((n+j)!j!n!)=binom(n+2j,j)/n!`, summing first in `j` and then in `n`
turns (6) into

```text
F(t)= exp(b0 t+y C(x))/(p sqrt(1-4x)).                            (9)
```

Define the algebraic germ

```text
q=p sqrt(1-4x)=sqrt((1-b1 t)^2-4 alpha t^2),      q(0)=1,         (10)
A(t)=b0 t+y C(x).                                                 (11)
```

Then the exact closed form is simply

```text
F(t)=exp(A(t))/q(t).                                              (12)
```

When `alpha!=0`, (11) may equivalently be written

```text
A(t)=b0 t + delta(p-q)/(2 alpha).                                (13)
```

When `alpha=0`, no division by `alpha` is made: `x=0`, `C(x)=1`,
`q=p`, and

```text
A(t)=b0 t+delta t^2/p,
F(t)=exp(b0 t+delta t^2/(1-b1t))/(1-b1t).                        (14)
```

Thus (9)--(12), including the zero-polynomial and degree-drop cases, are one
uniform identity; (13) is only shorthand off `alpha=0`.

### 2.2 Nullity is impossible off the one-sided locus

If every positive moment vanishes, then `F(t)=1` near zero. Equation (12)
therefore says

```text
exp(A(t))=q(t).                                                   (15)
```

Both `A` and `q` are algebraic germs holomorphic at zero. The lemma forces
`A` to be constant, and `A(0)=0` forces `A=0`. Hence `q=1`. Squaring (10) and
comparing the coefficients of `t` and `t^2` gives

```text
b1=0,                         alpha=0.                            (16)
```

Now the division-free formula (11) has `p=1`, `x=0`, and hence

```text
A(t)=b0t+delta t^2=0,
```

so `b0=delta=0`. This proves the forward implication in (2). Conversely,
`b=0` and `a*c=0` leave only one nonzero charge, so every positive moment is
zero by charge. Part A is proved.

This closes arbitrary complex coefficients at the central degree pair
`deg b=1, deg h=2`, its degree-drop inner offsets, and the low-degree sharp
boundary `deg b=0, deg h=2`. In particular no exceptional hyper-Bessel leading
zero survives on that boundary.

## 3. Part B: arbitrary charges and arbitrary degree on `h=kappa b^r`

Fix positive charges `p,q` and nonzero radial polynomials when degrees are
mentioned. Put

```text
g=gcd(p,q),       p0=p/g,       q0=q/g,       r=p0+q0,
h=s^(pq/g) a(s)^q0 c(s)^p0.                                      (17)
```

Assume, for some `kappa in C`, the exact polynomial identity

```text
h=kappa b^r.                                                       (18)
```

There is no bound on the degrees or coefficients of `a,b,c`. We prove

```text
E[(Z^p a+b+W^q c)^m]=0 for every m>=1
       iff b=h=0.                                                  (19)
```

Under (18), the degree pair is central, `deg h=r deg b`, whenever both sides
are nonzero.

### 3.1 Exact radial-shadow factorization

The primitive-return formula of THM-2017 is

```text
M_m=sum_(0<=k<=m/r)
       m!/((q0k)!(p0k)!(m-rk)!) L(h^k b^(m-rk)).                  (20)
```

Using (18), every channel has the same complete radial shadow:

```text
h^k b^(m-rk)=kappa^k b^m.
```

Consequently

```text
M_m=H_m(kappa)L(b^m),                                             (21)
H_m(kappa)=sum_(0<=k<=m/r)
 m! kappa^k/((q0k)!(p0k)!(m-rk)!).                               (22)
```

No termwise vanishing has been asserted: all channel cancellation is retained
in the single scalar `H_m`.

### 3.2 The hyper-Bessel coefficient sequence is not eventually zero

Let

```text
Phi_(p0,q0)(z)=sum_(k>=0) z^k/((q0k)!(p0k)!).                     (23)
```

Reindexing `m=rk+l` in (22) yields the exact entire EGF

```text
G(t):=sum_(m>=0) H_m(kappa)t^m/m!
     =exp(t) Phi_(p0,q0)(kappa t^r).                              (24)
```

Suppose `H_m(kappa)=0` for all sufficiently large `m`. Then `G` is a
polynomial `R(t)` with `R(0)=1`. Since `r>=2`, choose a nontrivial `r`-th root
of unity `omega`. The hyper-Bessel factor in (24) is invariant under
`t -> omega t`, so

```text
R(omega t)=exp((omega-1)t)R(t).                                  (25)
```

The left side divided by `R(t)` is rational, while the right side is the
exponential of the nonconstant algebraic germ `(omega-1)t`. This contradicts
the algebraic-exponential lemma. Therefore

```text
H_m(kappa)!=0 for arbitrarily large m.                            (26)
```

This also covers `kappa=0`, when in fact `H_m(0)=1`.

### 3.3 Finish of Part B

If `b` is a nonzero constant, `L(b^m)=b^m` never vanishes. If `deg b>=1`,
the one-variable EMP theorem (THM-1510) gives an `m0` such that

```text
L(b^m)!=0 for every m>=m0.                                       (27)
```

Choose in (26) an `m>=m0`. Equations (21), (26), and (27) give `M_m!=0`.
Thus all-moment nullity forces `b=0`, and (18) then forces `h=0`.

Conversely, if `b=h=0`, the domain property in (17) forces `a=0` or `c=0`.
The remaining polynomial has only one strict charge, so all its positive
moments vanish. This proves (19) and Part B.

## 4. Method boundary and assumption challenge

Part A is the maximal ordinary-EGF range suggested by degree: `b=O(s)` and
`sqrt(h)=O(s)` keep the radial integral convergent near `t=0`. At higher
radial degree the moment EGF is generally Gevrey-divergent, so (12) cannot be
claimed without a Borel--Leroy replacement. Part B escapes that wall because
the exact relation (18) collapses every channel to `b^m` before any summation.

For the required Tournament Analysis, candidate vertices were monomials,
charges, Catalan lattice points `(n,j)`, primitive-return channels `k`, and
proof obligations. The useful vertices in Part B are the channels `k`. Their
pairwise observable is the radial exponent left after (18):

```text
(rk+m-rk)-(rl+m-rl)=0.
```

Thus every edge is a tie, gauged by increasing `k`; the resulting tournament
is transitive, with score histogram `0,1,...,floor(m/r)`, zero directed cycles,
singleton SCCs, and the increasing-`k` Hamiltonian path. This quotient preserves
the *entire* common radial factor `b^m`, but destroys the scalar phase cancellation
inside `H_m`. Equation (24) and root-of-unity observability restore exactly that
destroyed coordinate.

The challenged assumption is the same one exposed by MISTAKE-211: return
channels are not independent equations. In Part A they are resummed into one
algebraic-exponential germ; in Part B they are resummed into `H_m`. Neither proof
declares an individual channel zero. At the time of this local theorem the
outside loci were open; THM-2022 now closes them globally by preserving a
complete lowest-face sum modulo a good prime.
