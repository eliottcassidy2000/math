---
id: THM-4068
title: "Squarefree Stern-packet and tournament-apex balance"
status: >
  PROVED + CITED PRIME WEIL INPUT + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED. For every odd squarefree denominator q, the THM-4059
  Stern packet and the full lower star of the apex q have absolute imbalance
  at most q^(1/2+o(1)), uniformly over the squarefree family. The complete
  Kloosterman kernel factors over the prime CRT coordinates, but canonical
  representative parity has cross-ratio -1 on the canonical `{0,1}x{0,1}`
  plaquette for every nontrivial CRT factorization and therefore is not a
  rank-one local tensor. Prime powers and
  odd nonsquarefree denominators, and hence the general unrestricted odd-
  composite family, remain OPEN.
source: codex-frontier-synthesis-creative-20260825b / squarefree Stern lane
audit: >
  PASS. The primary exact companion checks all 2,026 odd squarefree q through
  5,000, 388 exact cyclotomic-exponent CRT factorizations, every one of the
  1,358 canonical composite parity plaquettes, the half-box identity, divisor stars,
  and rational-square corollaries of the analytic bounds. The independent
  companion reconstructs continued-fraction depth by Euclid, computes lower
  tournament stars directly, uses an extended-Euclid half-box counter, and
  checks all 28,705 Fourier pairs on nine hostile composite moduli by iterated
  local-histogram convolution. Both normal/optimized pairs byte-match; neither
  script uses Python assertions or floating-point constants.
depends_on:
  - THM-4059-stern-brocot-depth-packet-character-and-divisor-star-convolution
  - THM-4061-stern-depth-modular-hyperbola-and-prime-apex-balance
related:
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
script: 04-computation/stern_squarefree_packet_apex_balance_thm4068.py
output: 05-knowledge/results/stern_squarefree_packet_apex_balance_thm4068.out
independent_audit_script: 04-computation/stern_squarefree_packet_apex_balance_thm4068_independent_audit.py
independent_audit_output: 05-knowledge/results/stern_squarefree_packet_apex_balance_thm4068_independent_audit.out
script_sha256: 99d055bd73ca8f47663829e954e33e1751c18280d1fef0e2e424674679df0584
output_sha256: d2b6b42769bc5e07bfaeef75771b1b6d8744fa22c7dcb7ae77dc7a381ae6bf45
independent_audit_script_sha256: 97339258924ffd6c9bd7c14f585d8f0a66e187e4348cbf7c2df79f7620b681ce
independent_audit_output_sha256: 22b5fe52045b91d4483bbae6e7b4a3cd978a3e1730548ace8203e74497f7858d
hash_basis: raw LF bytes
---

# THM-4068 -- squarefree CRT cancellation without parity factorization

**PROVED + CITED PRIME WEIL INPUT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.** The result is uniform over odd squarefree denominators.
Its only non-elementary input is the prime Kloosterman bound already pinned
and consumed by THM-4061; no prime-power estimate is imported here.

## 1. Exact completion and squarefree CRT factorization

Let `q>1` be odd and squarefree, put `r=omega(q)`, and use canonical
representatives `0<=[x]_q<q`. Retain THM-4059's packet

```text
eta_q(x)=(-1)^[x]_q,
S(q)=sum_(a in U_q) eta_q(a)eta_q(a^(-1)).             (1)
```

For `h,k mod q`, write

```text
c_q(h)=2/[q(1+exp(-2*pi*i*h/q))],
K_q(h,k)=sum_(x in U_q) exp(2*pi*i(hx+kx^(-1))/q).     (2)
```

The finite geometric series giving `c_q` uses only that `q` is odd. Thus
Fourier inversion in the two parity factors gives the exact identity

```text
S(q)=sum_(h,k mod q)c_q(h)c_q(k)K_q(h,k).              (3)
```

For every prime `p|q`, put

```text
q_p=q/p,                    t_p=q_p^(-1) mod p.        (4)
```

CRT writes each unit and its inverse as

```text
x=sum_(p|q)q_p t_p x_p mod q,
x^(-1)=sum_(p|q)q_p t_p x_p^(-1) mod q.               (5)
```

Consequently the complete kernel factors exactly:

```text
boxed: K_q(h,k)=prod_(p|q)K_p(t_p h,t_p k).           (6)
```

This is an identity of cyclotomic sums, not merely an absolute-value bound.
At a prime, the local cases are

```text
h=k=0 mod p:       K_p(h,k)=p-1,
exactly one zero:  K_p(h,k)=-1,
hk!=0 mod p:       |K_p(h,k)|<=2sqrt(p).               (7)
```

The final line is exactly the classical prime Weil estimate pinned in
[the THM-4061 reference sidecar](../../05-knowledge/reference/stern-depth-kloosterman-weil-pin.md).
The first two lines are elementary complete character sums.

If `g=gcd(h,k,q)`, multiplying the local estimates in `(7)` gives

```text
|K_q(h,k)| <= 2^omega(q/g)sqrt(qg)
             <= 2^r sqrt(qg).                         (8)
```

Primes dividing `g` contribute at most `p`; a prime not dividing `g`
contributes either `1` or at most `2sqrt(p)`. This proves `(8)` with every
degenerate local case included.

## 2. An explicit packet bound

For odd `n>1`, define

```text
Hstar(n)=2H_(n-1)-H_((n-1)/2),
L_q=sum_(h mod q)|c_q(h)|,
A_q=sum_(h mod q)|c_q(h)|sqrt(gcd(h,q)).               (9)
```

Since `sqrt(gcd(h,k,q))<=sqrt(gcd(h,q))`, equations `(3)` and `(8)` separate
the two Fourier variables:

```text
|S(q)| <= 2^r sqrt(q)L_q A_q.                         (10)
```

The THM-4061 harmonic estimate is primality-free, so

```text
L_q <= 1/q+Hstar(q).                                  (11)
```

There is also a divisor-sensitive estimate for `A_q`. The `h=0` term is
`q^(-1/2)`. For `h!=0`, set `d=gcd(h,q)`, write `h=dj` and `n=q/d`. Directly
from `(2)`,

```text
|c_q(dj)|=|c_n(j)|/d.                                 (12)
```

Therefore the exact-gcd `d` block contributes at most

```text
d^(-1/2)sum_(1<=j<n)|c_n(j)|
 <= d^(-1/2)Hstar(n).                                 (13)
```

Summing the blocks proves

```text
A_q <= 1/sqrt(q)
       +sum_(d|q,d<q) Hstar(q/d)/sqrt(d).             (14)
```

Combining `(10)`, `(11)`, and `(14)` yields the explicit uniform bound

```text
boxed:
|S(q)| <= 2^r sqrt(q)(1/q+Hstar(q))
          * (1/sqrt(q)
             +sum_(d|q,d<q)Hstar(q/d)/sqrt(d)).       (15)
```

No cancellation between different complete sums is assumed in `(15)`.

## 3. Uniform `q^(1/2+o(1))` and the tournament apex

The elementary estimate

```text
Hstar(n)<=2H_(n-1)<=2(1+log n)                        (16)
```

and `tau(q)=2^r` give, with an absolute implied constant,

```text
L_q=O(log q),              A_q<=2^r O(log q),
|S(q)|<=sqrt(q)4^r O(log^2 q).                        (17)
```

This is uniform in `q`. To absorb the factor depending on `r`, list the
distinct odd prime factors increasingly. The `j`-th is at least `j+2`, hence

```text
q>=prod_(j=1)^r(j+2)=(r+2)!/2.                        (18)
```

For `r>=4`, the last half of this product alone gives
`log q >= (r/2)log(r/2)`; for bounded `r`, the conclusion is immediate.
Thus, uniformly over odd squarefree `q -> infinity`, `r=o(log q)`. Every
fixed-base exponential in `r`, as well as `log^2 q`, is consequently
`q^o(1)`. Equation `(17)` proves

```text
boxed: |S(q)|<=q^(1/2+o(1)).                          (19)
```

Moreover,

```text
phi(q)=q prod_(p|q)(1-1/p)>=q(2/3)^r=q^(1-o(1)),      (20)
```

so `|S(q)|/phi(q)->0` uniformly on this family.

Now let `B(q)` be THM-4059's signed lower-star imbalance of the maximal
vertex `q` in the intrinsic depth tournament on `{1,...,q}`. That theorem
proves the exact divisor convolution

```text
B(q)=sum_(d|q,d>1)S(d).                               (21)
```

Every `d|q` is again odd and squarefree. Bounding each term in `(21)` by the
corresponding `q`-quantities from `(17)` and using at most `2^r` divisors gives

```text
|B(q)|<=sqrt(q)8^r O(log^2 q)=q^(1/2+o(1)).           (22)
```

Therefore the complete star of the apex satisfies

```text
indeg(q) =(q-1+B(q))/2,
outdeg(q)=(q-1-B(q))/2,
indeg(q),outdeg(q)=(q-1)/2+O(q^(1/2+o(1))),           (23)
```

and `|indeg(q)-outdeg(q)|/(q-1)->0`. This is the full star because `q` is
the maximal vertex of the initial tournament; it says nothing about the
degree of the same vertex inside a larger tournament.

## 4. The kernel factors, but canonical parity does not

The CRT estimate above is not an Euler product for `S`. The obstruction is
uniform and already visible on a two-by-two plaquette.

Let `q=uv` with coprime odd `u,v>1`, and define the canonical CRT idempotents

```text
e_u=v(v^(-1) mod u),          e_v=u(u^(-1) mod v).     (24)
```

Both lie strictly between `0` and `q`. Since their sum is congruent to `1`
modulo `q` and lies strictly between `1` and `2q`,

```text
e_u+e_v=q+1.                                           (25)
```

On the CRT plaquette `{0,1} x {0,1}`, the four canonical lifts are
`0,e_u,e_v,1`. Hence the cross-ratio of the global representative-parity
function is

```text
[eta_q(0)eta_q(1)]/[eta_q(e_u)eta_q(e_v)]
 =(-1)^(e_u+e_v+1)=-1.                                (26)
```

Every nonzero rank-one tensor `f_u(x)f_v(y)` has cross-ratio `+1`.
Therefore `eta_q` is not a tensor of local parity functions, or even a
rank-one tensor of arbitrary nonzero local functions. The separate local
Fourier transforms are invertible row and column operations, so the Fourier
coefficient array also cannot have rank one in CRT coordinates.

The exact connection ledger is

| source | target | preserved | destroyed / needed sidecar |
|---|---|---|---|
| unit hyperbola modulo `q` | prime CRT factors | additive character and complete Kloosterman kernel `(6)` | canonical order and representative parity; retain the dyadic CRT carry, equivalently the plaquette flux `(26)` |

The smallest odd squarefree composite is `15`, and it supplies the minimal
numerical hostile to multiplicativity:

```text
S(3)S(5)=0,                    S(15)=8.                (27)
```

Thus the proof extracts cancellation from a factorized kernel while retaining
the nonfactorizing parity sidecar.

## 5. Reproducibility and boundary

The primary audit uses THM-4059's inverse-parity formula and a Cartesian CRT
enumeration. It checks all odd squarefree `q<=5000`, rational-square
consequences of `(15)` and `(22)`, exact divisor stars, 388 selected complete
spectra, and the smallest-prime/cofactor canonical plaquette for every composite
`q` in that range. The independent
audit instead uses Euclidean continued-fraction depth, direct lower-star
enumeration, extended-Euclid half-box inversion, and iterated local histogram
convolution; it checks all odd squarefree `q<=3000` and every `(h,k)` on nine
hostile moduli.

Reproduce the two normal/optimized pairs from the repository root:

```text
python3 -B 04-computation/stern_squarefree_packet_apex_balance_thm4068.py
python3 -B -O 04-computation/stern_squarefree_packet_apex_balance_thm4068.py
python3 -B 04-computation/stern_squarefree_packet_apex_balance_thm4068_independent_audit.py
python3 -B -O 04-computation/stern_squarefree_packet_apex_balance_thm4068_independent_audit.py
```

The theorem proves nothing for a modulus divisible by `p^2`: extending `(8)`
requires a separately pinned and audited prime-power complete-sum input.
Odd nonsquarefree moduli, and hence a theorem over the general unrestricted
odd-composite family, therefore remain **OPEN** here.
There is also no Euler product or multiplicativity theorem, zero-packet
classification, cancellation beyond pointwise Weil, improved logarithmic
loss, Duffin--Schaeffer correlation estimate, Khinchin constant claim,
all-vertex tournament theorem, or LRC consequence. **QED.**
