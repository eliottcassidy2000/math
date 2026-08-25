---
id: THM-4071
title: "Stern prime-power stationary phase and all-odd apex balance"
status: >
  PROVED + CITED PRIME WEIL INPUT + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED. A self-contained p-adic stationary formula gives
  |K_(p^a)(h,k)| <= 2 sqrt(p^a gcd(h,k,p^a)) for every odd prime power;
  no prime-power exponential-sum estimate is imported. Fourier completion
  yields S(p^a),B(p^a)=O(sqrt(p^a) log^2(p^a)) uniformly. Combining the
  local bound with THM-4068's CRT kernel bookkeeping proves
  |S(q)|,|B(q)| <= q^(1/2+o(1)) uniformly over all odd q, hence normalized
  balance of every odd tournament apex family.
source: codex-frontier-synthesis-creative-20260825c / Stern prime-power lane
audit: >
  PASS. The primary companion checks 93 prime-power packets (68
  nonsquarefree) through q=390625, 34 complete lift layers and 47,872
  fibers, all 40,672 frequencies on eight prime-power universes, 90 exact
  Gauss norms, 56 hostile parity plaquettes, and all 62,244 frequencies in
  four mixed prime-power CRT decompositions. The independent companion
  reconstructs signs by Euclidean continued-fraction depth, computes lower
  stars directly, uses extended Euclid, generic polynomial division by
  Phi_(p^a), the z-z^(-1) p-adic Morse coordinate, and tuplewise CRT. It
  independently checks 49 packets, 20 lift layers, the same 40,672 spectral
  frequencies and 62,244 mixed frequencies. Both normal/optimized pairs
  byte-match; neither script uses Python assertions or floating-point
  literals.
depends_on:
  - THM-4059-stern-brocot-depth-packet-character-and-divisor-star-convolution
  - THM-4061-stern-depth-modular-hyperbola-and-prime-apex-balance
  - THM-4068-squarefree-stern-packet-and-tournament-apex-balance
related:
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
script: 04-computation/stern_prime_power_stationary_all_odd_thm4071.py
output: 05-knowledge/results/stern_prime_power_stationary_all_odd_thm4071.out
independent_audit_script: 04-computation/stern_prime_power_stationary_all_odd_thm4071_independent_audit.py
independent_audit_output: 05-knowledge/results/stern_prime_power_stationary_all_odd_thm4071_independent_audit.out
script_sha256: 2ec61a857c6139ad3beed446a71167aa12a3475d52fd714751d9f0d0777217b9
output_sha256: db160ed6d0644070f6ae06ef21891abbc703ffd7568df24a262bfce6246c2df3
independent_audit_script_sha256: cdbdd0e00c989a2cdca1a709ec2971c400a118f6e941024b970cbbe71cfe0eb9
independent_audit_output_sha256: 7f28a5566477b3484db2f0e874aa1618c4705a41c573ebf5dfbf2dfda2079ae6
hash_basis: raw LF bytes
---

# THM-4071 -- Stern prime-power stationary phase and all-odd apex balance

**PROVED in the stated odd-modulus scopes.** THM-4068 stopped at squarefree
denominators because it imported no prime-power Kloosterman estimate. The
missing local estimate is proved here from valuation reduction, exact
derivative cancellation, and a two-critical-class p-adic Morse coordinate.
Only the residual exponent-one case uses the prime Weil estimate already
pinned by THM-4061. No prime-power estimate is cited or imported.

## 1. Statement and notation

For every odd `q>1`, retain THM-4059's canonical-representative packet and
signed lower-star imbalance

```text
eta_q(x)=(-1)^[x]_q,
S(q)=sum_(x in U_q)eta_q(x)eta_q(x^(-1)),
B(q)=sum_(d|q,d>1)S(d).                                      (1)
```

The maximal vertex `q` of the initial depth tournament has

```text
indeg(q) =(q-1+B(q))/2,
outdeg(q)=(q-1-B(q))/2.                                      (2)
```

Write `e_q(z)=exp(2*pi*i*z/q)` and define the complete kernel

```text
K_q(h,k)=sum_(x in U_q)e_q(hx+kx^(-1)).                       (3)
```

The local theorem is the following uniform prime-power estimate.

```text
boxed: for q=p^a, p odd prime,
       |K_q(h,k)| <= 2 sqrt(q gcd(h,k,q)).                     (4)
```

With `gcd(0,0,q)=q`, equation `(4)` includes the origin and every degenerate
axis. Fourier completion then gives, uniformly over odd prime powers,

```text
boxed: |S(p^a)|=O(sqrt(p^a)log^2(p^a)),
       |B(p^a)|=O(sqrt(p^a)log^2(p^a)).                        (5)
```

More strongly, combining `(4)` over the prime-power CRT blocks proves,
uniformly as odd `q` tends to infinity,

```text
boxed: |S(q)|<=q^(1/2+o(1)),
       |B(q)|<=q^(1/2+o(1)).                                  (6)
```

Consequently every odd apex family is normalized-balanced:

```text
|indeg(q)-outdeg(q)|/(q-1)=|B(q)|/(q-1)->0.                   (7)
```

## 2. Exact lift recursion and its noncontractive boundary

The THM-4056 lift clock has an exact packet-level description. Let `q=p^n`,
`Q=pq`, and choose `x in U_q` with canonical inverse `u`. Put

```text
kappa=(xu-1)/q.                                                (8)
```

Every lift is `X=x+tq`, `0<=t<p`. Write its canonical inverse modulo `Q` as
`U=u+sq`, `0<=s<p`. Expanding `XU=1 mod Q` and using `xu=1+kappa q` gives

```text
kappa+x s+u t=0 mod p,
s=-u kappa-u^2 t mod p.                                       (9)
```

Because `q` is odd,

```text
(-1)^(X+U)=(-1)^(x+u)(-1)^(t+s).                              (10)
```

Define the affine parity correlation

```text
C_p(alpha,beta)
 =sum_(t=0)^(p-1)(-1)^(t+[alpha+beta t]_p).                    (11)
```

Equations `(9)` and `(10)` prove the exact recursion

```text
S(p^(n+1))
 =sum_(x in U_(p^n))(-1)^(x+x^(-1))
    C_p(-u kappa,-u^2).                                        (12)
```

The carry `kappa mod p` is the sidecar lost by merely reducing a lifted unit
and inverse. Equation `(12)` does not itself contract. For `x=u=1`, one has

```text
C_p(0,-1)=2-p,                                                 (13)
```

which is asymptotically full fiber size. Full correlation occurs in the
actual tower: at the first checked witness `p=3`, `3^3 -> 3^4`, the lower
unit data `(x,u,kappa)=(4,7,1)` give `C_3(2,2)=3`. Thus a fiberwise triangle
inequality in `(12)` is a stopping obstruction, not a proof of `(5)`.

At the mandatory first nonsquarefree hostile,

```text
S(3)=2,       S(9)=-2,       B(9)=0.                           (14)
```

The lifts of `1 mod 3` have signs `+--`, while the lifts of `2 mod 3` have
signs `--+`; both fiber sums are `-1`.

## 3. Prime-power Kloosterman stationary phase

Let `q=p^a`, and use the convention `v_p(0)=a`. Put

```text
v=min(v_p(h),v_p(k),a).                                       (15)
```

If `v<a`, write `h=p^v H`, `k=p^v K`, and `m=a-v`. Reduction of
`U_(p^a)` onto `U_(p^m)` has fibers of size `p^v`, commutes with inversion,
and reduces the additive character. Hence

```text
boxed: K_(p^a)(h,k)=p^v K_(p^m)(H,K).                         (16)
```

We now evaluate every primitive residual pair for `m>=2`. Split the unit sum
by residue modulo `p`. Inside one residue class write

```text
x=y+t p^(m-1),       0<=t<p.
```

Then

```text
x^(-1)=y^(-1)-t p^(m-1)y^(-2) mod p^m,                        (17)
```

so the top-digit sum is

```text
sum_(t mod p)e_p(t(H-Ky^(-2))).                                (18)
```

It vanishes unless

```text
H y^2=K mod p.                                                 (19)
```

Therefore, when exactly one of `H,K` is a unit, the residual kernel is zero.
When both are units but `K/H` is a quadratic nonresidue modulo `p`, it is
also zero.

Suppose both are units and `(19)` is soluble. Hensel lifting gives exactly
two roots `+s,-s mod p^m` with

```text
s^2=K/H mod p^m.                                               (20)
```

On the `+s` critical class every unit is uniquely `x=s z^2` with
`z in 1+p Z/(p^m)`. Squaring is an automorphism of the principal units since
`p` is odd. Moreover

```text
z |-> y=z-z^(-1)                                               (21)
```

is a bijection from the principal units to `p Z/(p^m)`: for each such `y`,
the equation `z^2-yz-1=0` has exactly one root `z=1 mod p`, because its
derivative is `2 mod p`. The identity

```text
H(x+s^2x^(-1))
 =Hs(z^2+z^(-2))
 =2Hs+Hs(z-z^(-1))^2                                          (22)
```

is exact. The `-s` class gives both signs reversed.

For a unit `c`, define the quadratic Gauss sum

```text
G_j(c)=sum_(z mod p^j)e_(p^j)(c z^2),       G_0(c)=1.          (23)
```

Writing the coordinate in `(21)` as `y=pz` and counting its `p` lifts gives
the exact stationary formula

```text
boxed:
K_(p^m)(H,K)
 =p e_(p^m)(2Hs)G_(m-2)(Hs)
  +p e_(p^m)(-2Hs)G_(m-2)(-Hs).                               (24)
```

The Gauss magnitude is elementary and imported from nowhere. Splitting off
the last base-`p` digit gives

```text
G_j(c)=p G_(j-2)(c),          j>=2.                            (25)
```

For `j=1`, the change of variables `(r,s)->(r-s,r+s)` is bijective over
`F_p` and proves `|G_1(c)|^2=p`. Thus

```text
|G_j(c)|=p^(j/2),
|K_(p^m)(H,K)|<=2p^(m/2).                                     (26)
```

If `m=1`, the elementary prime axes are

```text
K_p(0,0)=p-1,
K_p(H,0)=K_p(0,H)=-1,       H!=0,                              (27)
```

and the two-unit case uses exactly THM-4061's pinned prime Weil bound
`|K_p(H,K)|<=2sqrt(p)`. Combining `(16)`, `(26)`, and `(27)` proves `(4)`.
At the origin, `phi(q)<=q<=2sqrt(q*q)`.

The degenerate `q=9` controls display every boundary:

```text
K_9(0,0)=6,       K_9(1,0)=0,
K_9(3,0)=3K_3(1,0)=-3,       K_9(1,2)=0,
K_9(1,1)=3(e_9(2)+e_9(-2)).                                  (28)
```

These are respectively the origin, primitive axis, residual-prime axis,
nonsquare pair, and two-critical-class pair.

## 4. Fourier completion at a prime power

For every odd `q`, not only primes, the canonical parity Fourier coefficient
is

```text
c_q(h)=2/[q(1+e_q(-h))],                                      (29)
```

and THM-4068's exact completion is

```text
S(q)=sum_(h,k mod q)c_q(h)c_q(k)K_q(h,k).                      (30)
```

Put

```text
Hstar(n)=2H_(n-1)-H_((n-1)/2),
L_q=sum_h|c_q(h)|,
A_q=sum_h|c_q(h)|sqrt(gcd(h,q)).                               (31)
```

For `q=p^a`, equations `(4)` and `(30)`, together with
`sqrt(gcd(h,k,q))<=sqrt(gcd(h,q))`, give

```text
|S(q)|<=2sqrt(q)L_q A_q.                                      (32)
```

The primality-free harmonic calculation of THM-4068 gives

```text
L_q<=1/q+Hstar(q).                                             (33)
```

Partitioning nonzero frequencies by `d=gcd(h,q)` and using
`c_q(dj)=c_(q/d)(j)/d` gives

```text
A_q<=1/sqrt(q)
     +sum_(r=0)^(a-1)Hstar(p^(a-r))/p^(r/2).                   (34)
```

Since `Hstar(n)<=2(1+log n)` and

```text
sum_(r>=0)p^(-r/2)<=1/(1-3^(-1/2)),                            (35)
```

equations `(32)`--`(35)` prove the packet estimate in `(5)` with an absolute
constant. THM-4059 gives

```text
B(p^a)=sum_(j=1)^a S(p^j).                                    (36)
```

The geometric bound

```text
sum_(j=1)^a sqrt(p^j)
 <=sqrt(p^a)/(1-p^(-1/2))
 <=sqrt(p^a)/(1-3^(-1/2))                                    (37)
```

proves the star estimate in `(5)` and therefore prime-power apex balance.

## 5. CRT closure for every odd denominator

Let

```text
q=prod_(i=1)^r q_i,       q_i=p_i^(a_i),       r=omega(q),
Q_i=q/q_i,                t_i=Q_i^(-1) mod q_i.                (38)
```

CRT for a unit and its inverse gives the exact kernel factorization

```text
boxed: K_q(h,k)=prod_(i=1)^r K_(q_i)(t_i h,t_i k).             (39)
```

If `g=gcd(h,k,q)`, the local bound `(4)` therefore yields

```text
|K_q(h,k)|<=2^r sqrt(qg).                                      (40)
```

This factors the complete kernel, not canonical representative parity.
Indeed, for every nontrivial coprime odd factorization `q=uv`, the canonical
CRT idempotents satisfy `e_u+e_v=q+1`; hence the representative-parity flux
on the `{0,1}x{0,1}` plaquette is `-1`, exactly as in THM-4068.

Applying the Fourier separation from `(32)` now gives

```text
|S(q)|<=2^r sqrt(q)L_q A_q,
A_q<=1/sqrt(q)
     +sum_(d|q,d<q)Hstar(q/d)/sqrt(d).                          (41)
```

Put

```text
C=(1-3^(-1/2))^(-1).
```

Then

```text
sum_(d|q)d^(-1/2)
 =prod_(p^a||q)(1+p^(-1/2)+...+p^(-a/2))<=C^r.                 (42)
```

Thus `(41)` gives

```text
|S(q)|<=sqrt(q)(2C)^r O(log^2(2q)).                            (43)
```

The product of the `r` distinct odd primes dividing `q` is at least
`(r+2)!/2`; hence `r=o(log q)` uniformly. Every fixed-base exponential in
`r`, as well as every fixed power of `log q`, is `q^o(1)`. This proves the
packet half of `(6)`.

The divisor star requires retaining the local exponential factor in `(43)`;
an unweighted estimate for `sum sqrt(d)` alone is insufficient. Choose a
fixed `A_0>2C` large enough that, for every odd `d>1`,

```text
|S(d)|<=A_0^omega(d) sqrt(d) O(log^2(2d)).                      (44)
```

For `q=prod p^a`, the weighted divisor sum factors exactly:

```text
sum_(d|q) A_0^omega(d)sqrt(d)
 =prod_(p^a||q)(1+A_0 sum_(b=1)^a p^(b/2)).                     (45)
```

After factoring out `sqrt(q)`, each local block in `(45)` is

```text
p^(-a/2)+A_0 sum_(b=1)^a p^((b-a)/2)
 <=1+A_0/(1-3^(-1/2))=:D.                                     (46)
```

Therefore THM-4059's divisor convolution gives

```text
|B(q)|
 <=O(log^2(2q))sum_(d|q)A_0^omega(d)sqrt(d)
 <=sqrt(q)D^r O(log^2(2q))
 =q^(1/2+o(1)).                                                 (47)
```

This proves the star half of `(6)` without dropping the
`A_0^omega(d)` factor and completes `(7)`.

## 6. Exact audits and boundary

The primary companion checks:

- 93 prime-power packets, including 68 nonsquarefree powers, through
  `q=390625`, together with every half-box identity;
- 34 complete affine lift layers and 47,872 lower fibers, including `(14)`
  and the first full-correlation witness;
- all 40,672 frequency pairs for
  `(p,a)=(3,1..4),(5,2..3),(7,2),(11,2)` by exact equality in
  `Z[zeta_(p^a)]`, with every degenerate category separated;
- 90 exact quadratic Gauss norm identities;
- 56 non-tensorial parity plaquettes, `S(3)S(5)=0` versus `S(15)=8`, and all
  62,244 frequency pairs in `45=9*5`, `63=9*7`, `75=3*25`, and `225=9*25`.

The independent companion shares none of the primary packet, inverse, star,
normal-form, or stationary-histogram routines. It reconstructs packet signs
from Euclidean continued-fraction digit sums, computes every lower star
directly, uses extended Euclid, divides exponent polynomials generically by

```text
Phi_(p^a)(X)=sum_(t=0)^(p-1)X^(t p^(a-1)),                     (48)
```

audits the `z-z^(-1)` Morse bijection, and enumerates local CRT tuples rather
than convolving the primary histograms. It checks 49 packets, 20 lift layers,
the same 40,672 spectral frequencies, and the same 62,244 mixed frequencies.
Both companions reproduce byte-identically in normal and optimized modes.

Run from the repository root:

```text
python3 -B 04-computation/stern_prime_power_stationary_all_odd_thm4071.py
python3 -B -O 04-computation/stern_prime_power_stationary_all_odd_thm4071.py
python3 -B 04-computation/stern_prime_power_stationary_all_odd_thm4071_independent_audit.py
python3 -B -O 04-computation/stern_prime_power_stationary_all_odd_thm4071_independent_audit.py
```

### Scope ledger

- **THM-4059:** supplies the exact Stern packet and divisor-star identity.
- **THM-4061:** supplies only the pinned prime Weil input at residual exponent
  one, plus the odd-modulus Fourier coefficient already reused by THM-4068.
- **THM-4068:** supplies the global Fourier/CRT architecture, harmonic
  separation, and canonical parity plaquette. Its squarefree boundary is
  crossed only after the self-contained prime-power lemma `(4)`.
- **THM-4056:** is a related denominator-depth/lift-clock sidecar; it is not an
  analytic proof dependency.

The theorem proves no Euler product or multiplicativity for `S`, no
zero-packet classification, no improvement beyond the two logarithms from
pointwise completion, no even-denominator asymptotic, and no LRC consequence.
In particular, **LRC(14) remains OPEN. QED.**
