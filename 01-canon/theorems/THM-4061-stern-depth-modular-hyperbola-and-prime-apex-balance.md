---
id: THM-4061
title: "Stern depth modular hyperbola and prime apex balance"
status: >
  PROVED + CITED WEIL INPUT + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED. For every odd denominator, the THM-4059 Stern depth
  packet is exactly four times a half-box modular-hyperbola count minus the
  totient. At prime denominators, Fourier completion and the classical Weil
  bound give |S(p)|=O(sqrt(p) log^2 p), so the normalized lower-star bias of
  the prime apex tends to zero. This proves no composite asymptotic, zero-
  packet classification, Euler product, Khinchin claim, or LRC consequence.
source: codex-frontier-synthesis-breakthrough-20260825 / depth-packet cancellation lane
audit: >
  PASS. The primary companion independently compares canonical continued-
  fraction depth, THM-4059 inverse parity, and the modular hyperbola for all
  2,499 odd denominators through 5,000, including 668 primes. A second
  implementation generates 2,736,187 Stern--Brocot fractions directly and
  checks all 1,499 odd denominators through 3,000. A separate hostile
  referee reconstructed every Fourier coefficient and Kloosterman completion
  for all 21 odd primes through 79. Every normal/optimized pair byte-matches.
depends_on:
  - THM-4059-stern-brocot-depth-packet-character-and-divisor-star-convolution
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
related:
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
  - THM-873-ramanujan-fourier-expansion-of-interval-core-good-sets
script: 04-computation/stern_depth_modular_hyperbola_prime_balance_thm4061.py
output: 05-knowledge/results/stern_depth_modular_hyperbola_prime_balance_thm4061.out
independent_audit_script: 04-computation/stern_depth_modular_hyperbola_prime_balance_thm4061_independent_audit.py
independent_audit_output: 05-knowledge/results/stern_depth_modular_hyperbola_prime_balance_thm4061_independent_audit.out
hostile_audit_script: 04-computation/stern_depth_modular_hyperbola_prime_balance_thm4061_hostile_audit.py
hostile_audit_output: 05-knowledge/results/stern_depth_modular_hyperbola_prime_balance_thm4061_hostile_audit.out
script_sha256: ef6af342a844ab219dd4ec0db93898540510c99cdf5b9d5848343ddc1bedceab
output_sha256: a358382749869609b7f48f29090e2f380de2cea16f105097ec98e8e0c3e333ec
independent_audit_script_sha256: 7cdb75d9b4ade320a037dd0a8b46c09a9cc7640f5bbe6c3eb977958930d0f98c
independent_audit_output_sha256: c0076d6b623cf34cb9937f2037fffa74fdb3694822c30973051fe3112aca17f6
hostile_audit_script_sha256: a0b62b2fa5889e58530e1d76245d203510446ae26582dc4bb3bc7c88583f3ca0
hostile_audit_output_sha256: 1fcb7d7a8b2335c07fd6bb4780369d64c70884280092b9e261fde4f8548ea752
hash_basis: raw LF bytes
---

# THM-4061 -- odd Stern packets are modular-hyperbola discrepancies

**PROVED + CITED WEIL INPUT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.** THM-4059 converts Stern--Brocot depth parity inside an odd
exact-denominator packet into parity of a residue and its modular inverse.
The first result below identifies the resulting sum with a modular-hyperbola
discrepancy. The second imports only the classical prime Kloosterman estimate;
every other step is proved here.

## 1. The exact half-box identity

Let `q>1` be odd. Use representatives `1<=a<q` for `U_q`, and let
`a^(-1)` denote the representative modular inverse. Retain THM-4059's

```text
epsilon(a,q)=(-1)^(a+a^(-1)),
S(q)=sum_(a in U_q) epsilon(a,q).                     (1)
```

Put `m=(q-1)/2` and define

```text
N(q)=#{(r,s):1<=r,s<=m and 4rs=1 mod q}.              (2)
```

Then

```text
boxed: S(q)=4N(q)-phi(q).                             (3)
```

### Proof

Partition `U_q` into its even and odd representatives. The mirror
`a -> q-a` toggles parity, so each class has size `phi(q)/2`. Write
`n_EE,n_EO,n_OE,n_OO` for the four counts according to the parity of
`(a,a^(-1))`. Inversion is a bijection, so `n_EO=n_OE`. Equality of the
two row sums now gives `n_EE=n_OO`. Therefore

```text
S(q)=n_EE+n_OO-n_EO-n_OE
    =4n_EE-phi(q).                                    (4)
```

Writing an even unit and its even inverse as `a=2r`, `a^(-1)=2s` identifies
`n_EE` exactly with `(2)`. This proves `(3)`.

The identity is exact for every odd composite as well as every odd prime.
It does not assert that `S` is multiplicative; THM-4059 already proves that
it is not.

## 2. Fourier completion at a prime

Let `p` be an odd prime, `zeta=exp(2*pi*i/p)`, and define the canonical
parity function on `F_p` by

```text
eta(x)=(-1)^x,             0<=x<p.                    (5)
```

Its normalized Fourier coefficients are

```text
c_h=(1/p)sum_(x=0)^(p-1) eta(x)zeta^(-hx)
    =2/[p(1+zeta^(-h))].                              (6)
```

Indeed `(6)` is a finite geometric sum, since `p` is odd. In particular,

```text
c_0=1/p,
|c_h|=1/[p|cos(pi h/p)|].                             (7)
```

Define the complete prime Kloosterman sum

```text
Kl_p(h,k)=sum_(a in F_p^*) zeta^(ha+k a^(-1)).        (8)
```

Fourier inversion in both parity factors gives the exact completion

```text
S(p)=sum_(h,k mod p)c_h c_k Kl_p(h,k).                (9)
```

There is no omitted conjugate or factor of `p`: the sign in `(6)` is chosen
so that `eta(a)=sum_h c_h zeta^(ha)`.

The degenerate sums are

```text
Kl_p(0,0)=p-1,
Kl_p(h,0)=Kl_p(0,h)=-1,        h!=0.                 (10)
```

For `hk!=0`, import the classical Weil estimate

```text
|Kl_p(h,k)|<=2sqrt(p).                                (11)
```

This is the sole cited analytic input. Put

```text
C_p=sum_(h=1)^(p-1)|c_h|.                             (12)
```

Separating the origin, the two axes, and the nondegenerate square in `(9)`
gives

```text
|S(p)| <= (p-1)/p^2+2C_p/p+2sqrt(p)C_p^2.            (13)
```

## 3. The explicit harmonic bound

For `n>=1`, write `H_n=sum_(j=1)^n 1/j`, and put

```text
H*(p)=2H_(p-1)-H_((p-1)/2).                           (14)
```

Write `p=2m+1`. Pairing `h` with `p-h` in `(12)` gives

```text
C_p=(2/p)sum_(h=1)^m sec(pi h/p)
   =(2/p)sum_(j=0)^(m-1)csc(pi(j+1/2)/p).             (15)
```

Since `sin x>=2x/pi` on `[0,pi/2]`,

```text
C_p <= sum_(j=0)^(m-1)1/(j+1/2)
    =2H_(p-1)-H_m
    =H*(p).                                           (16)
```

Substitution into `(13)` proves the explicit estimate

```text
boxed: |S(p)| <= 2sqrt(p)H*(p)^2
                  +2H*(p)/p+(p-1)/p^2.               (17)
```

The elementary bound `H_n<=1+log n` now yields

```text
S(p)=O(sqrt(p)log^2 p).                               (18)
```

No cancellation between distinct Kloosterman sums was used. Sharpening the
two logarithms would require additional input beyond the pointwise Weil
bound.

## 4. Prime-apex tournament balance

THM-4059 defines the lower-star imbalance `B(q)` of the maximal vertex in
the depth tournament on `{1,...,q}` and proves

```text
B(q)=sum_(d|q,d>1)S(d).                               (19)
```

For a prime `p`, this reduces to `B(p)=S(p)`. Hence

```text
indeg(p) =(p-1+S(p))/2,
outdeg(p)=(p-1-S(p))/2,                               (20)
```

and `(18)` gives

```text
indeg(p),outdeg(p)=(p-1)/2+O(sqrt(p)log^2 p),          (21)
|indeg(p)-outdeg(p)|/(p-1)=O(log^2 p/sqrt(p))->0.     (22)
```

This is the full star only because `p` is the apex of the initial tournament
on `{1,...,p}`. It is not the degree of vertex `p` inside a larger
tournament, nor an estimate for every vertex.

## 5. Imported input, audits, and boundary

The imported estimate `(11)` is André Weil's prime Kloosterman bound; its
primary citation and exact consumer boundary are pinned in
[the reference sidecar](../../05-knowledge/reference/stern-depth-kloosterman-weil-pin.md).
The modular-hyperbola identity, parity Fourier coefficients, harmonic loss,
and tournament consequence are in-repo proofs.

The primary exact companion checks `(3)` through `q=5000` by three routes:
canonical continued-fraction depth, inverse parity, and modular inversion in
the half box. The independent companion generates the Farey tree without
calling either packet implementation. The hostile companion numerically
reconstructs every Fourier coefficient and `(9)` through prime `79`.

Reproduce all three normal/optimized pairs from the repository root with

```text
python3 -B 04-computation/stern_depth_modular_hyperbola_prime_balance_thm4061.py
python3 -B -O 04-computation/stern_depth_modular_hyperbola_prime_balance_thm4061.py
python3 -B 04-computation/stern_depth_modular_hyperbola_prime_balance_thm4061_independent_audit.py
python3 -B -O 04-computation/stern_depth_modular_hyperbola_prime_balance_thm4061_independent_audit.py
python3 -B 04-computation/stern_depth_modular_hyperbola_prime_balance_thm4061_hostile_audit.py
python3 -B -O 04-computation/stern_depth_modular_hyperbola_prime_balance_thm4061_hostile_audit.py
```

The theorem does not classify `S(q)=0`, treat prime powers or composite
moduli asymptotically, produce an Euler product, improve a Duffin--Schaeffer
correlation estimate, identify Khinchin's constant at a named number, or
imply LRC(14). **QED.**
