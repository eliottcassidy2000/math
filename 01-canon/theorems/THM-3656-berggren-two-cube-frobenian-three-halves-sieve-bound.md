---
id: THM-3656
title: "Berggren two-cube Frobenian three-halves sieve bound"
status: >
  PROVED + ANALYTIC-SIEVE + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED;
  CITED for the arbitrary-dimensional upper fundamental lemma and the prime
  number theorem in fixed progressions.  The three correlated bad lines have
  sieve dimension 3/2, and the number of local-support survivors in
  0<V<n<=N is O(N^2/(log N)^(3/2)).  No lower bound, asymptotic constant,
  primitive-tail estimate, Pell/compiler survivor, or two-cube existence
  result is claimed.
source: kps-s191 / THM-3645 support-density continuation, 2026-08-21
audit: >
  PASS -- an independent hostile audit verified the local union sizes, CRT
  multiplicativity, uniform triangle residue remainder, summed level of
  distribution, character average, fractional-dimension sieve hypotheses,
  small-prime treatment and set inclusion.  It made the fixed parameters
  D=N^(1/3), s=16 and z=N^(1/48) explicit; no lower sieve is used.
depends_on:
  - THM-3645-berggren-two-cube-local-prime-support-and-mod24-cones
related:
  - THM-3375-berggren-positive-two-cube-pell-ray
script: 04-computation/berggren_two_cube_frobenian_sieve_dimension_thm3656.py
output: 05-knowledge/results/berggren_two_cube_frobenian_sieve_dimension_thm3656.out
script_sha256: b2d7b544639fe7591e9b170f4b014ad7f97c065d176485d72d37c19a8a1876b6
output_sha256: b5535603962fa86d108ddb1403a0c7ee0d16b602f5f5113ca1f3a5c30ebfe903
semantic_sha256: 2ea01388cebc907078d58844ff4799e8433d509f342195c6ef0240829476ae82
hash_basis: raw LF bytes
external:
  - https://doi.org/10.1016/0022-314X(88)90046-7
  - https://terrytao.wordpress.com/2015/01/21/254a-notes-4-some-sieve-theory/
---

# THM-3656 -- Berggren two-cube Frobenian three-halves sieve bound

**PROVED + ANALYTIC-SIEVE + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED;
CITED for the upper fundamental lemma and PNT in fixed progressions.**

For `N>=3`, let

```text
T_N={(n,V) in Z^2:0<V<n<=N}.                           (1)
```

Write `chi_24` and `chi_-8` for the quadratic Dirichlet characters attached
to discriminants `24` and `-8`.  Let `B(N)` count pairs in `(1)` such that,
for every prime `p>=5`,

```text
p|nV       implies chi_24(p)=+1,
p|(2n+V)   implies chi_-8(p)=+1.                       (2)
```

Then

```text
B(N) << N^2/(log N)^(3/2).                            (3)
```

The implied constant is absolute.  In particular, the primitive,
parity-correct, congruence-restricted slopes satisfying all mod-prime support
conditions of THM-3645 obey `(3)`.

## 1. Three correlated bad lines

For a prime `p>=5`, let `Omega_p` be the union of

```text
n=0 and V=0     if chi_24(p)=-1,
2n+V=0          if chi_-8(p)=-1.                      (4)
```

The three possible lines are distinct over `F_p`.  Put

```text
e_D(p)=(1-chi_D(p))/2,            rho(p)=|Omega_p|.   (5)
```

In the four Frobenius classes one obtains exactly

| `(e_24,e_-8)` | `(0,0)` | `(0,1)` | `(1,0)` | `(1,1)` |
|---|---:|---:|---:|---:|
| `rho(p)` | `0` | `p` | `2p-1` | `3p-2` |

Equivalently,

```text
rho(p)=e_24(2p-1)+e_-8 p-e_24 e_-8.                  (6)
```

With `g(p)=rho(p)/p^2`, equations `(5)--(6)` give

```text
g(p)=[3/2-chi_24(p)-chi_-8(p)/2]/p+O(1/p^2).          (7)
```

This is a single coupled density.  Treating the three lines as independent
half-dimensional sieves would lose their intersections and is not used.

## 2. The sieve dimension is `3/2`

Both characters in `(7)` are fixed, distinct and nonprincipal.  The prime
number theorem in the eight reduced residue classes modulo `24` gives,
uniformly for `5<=w<z`,

```text
sum_(w<=p<z) g(p) log p=(3/2)log(z/w)+O(1).            (8)
```

Partial summation yields the standard dimension condition

```text
product_(w<=p<z)(1-g(p))^-1
  <=K (log z/log w)^(3/2),                             (9)
```

and, in particular,

```text
mathcal V(z)=product_(5<=p<z)(1-g(p))
  asymp (log z)^(-3/2).                               (10)
```

Thus the coupled system has real sieve dimension `kappa=3/2`.  The
Diamond--Halberstam--Richert upper sieve applies to arbitrary real dimension
greater than one; the cited fundamental-lemma formulation supplies the
corresponding truncated upper weights.

## 3. Exact CRT density and lattice remainder

For squarefree `d` coprime to `6`, let `Omega_d` be the simultaneous inverse
image of the `Omega_p`, `p|d`.  The Chinese remainder theorem gives

```text
rho(d)=|Omega_d|=product_(p|d)rho(p),
g(d)=rho(d)/d^2=product_(p|d)g(p),
rho(d)<=3^omega(d)d.                                  (11)
```

Let

```text
A_d=#{(n,V) in T_N:(n,V) mod d in Omega_d},
X=|T_N|=N(N-1)/2.                                     (12)
```

For one residue pair modulo `d`, slicing the triangle by the `n` coordinate
gives uniformly

```text
#{(n,V) in T_N:(n,V)=(a,b) mod d}
  =X/d^2+O(N/d+1).                                    (13)
```

Summing `(13)` over `Omega_d` gives

```text
A_d=Xg(d)+r_d,
|r_d| << 3^omega(d)(N+d).                             (14)
```

Since `3^omega(d)<=tau_3(d)`, the standard divisor-sum estimate yields

```text
sum_(d<=D,(d,6)=1)|r_d|
  << ND(log D)^2+D^2(log D)^2.                        (15)
```

## 4. Explicit upper-sieve parameters

Take

```text
D=N^(1/3),               s=16,
z=D^(1/16)=N^(1/48).                                  (16)
```

The quantitative upper fundamental lemma is admissible here because
`s>=9kappa+1`.  Applied to `(9),(12)--(15)`, it gives

```text
S(N,z)
 << X mathcal V(z)+sum_(d<D,d|product_(5<=p<z)p)|r_d|. (17)
```

Equations `(10),(15)--(16)` imply

```text
X mathcal V(z) << N^2/(log N)^(3/2),
sum_(d<D)|r_d|
 << N^(4/3)(log N)^2+N^(2/3)(log N)^2
 =o(N^2/(log N)^(3/2)).                               (18)
```

Every pair counted by `B(N)` avoids `(4)` for all `p>=5`, hence survives the
truncated sieve through `z`.  Therefore `B(N)<=S(N,z)`, proving `(3)`.

## 5. Exact companion and boundary

Reproduce the local and CRT ledger with

```bash
python3 04-computation/berggren_two_cube_frobenian_sieve_dimension_thm3656.py
python3 -O 04-computation/berggren_two_cube_frobenian_sieve_dimension_thm3656.py
```

The assertion-free companion raw-pins THM-3645; verifies the complete joint
character table modulo `24`, its exact mean `3/2`, every bad-line union for
all 52 primes `5<=p<=251`, two direct CRT products, and the multiplicative
remainder majorant controls.  Normal and optimized streams match the stored
transcript.  The analytic application is the proof in Sections 2--4, not a
claim inferred from the finite checks.

Dropping the conditions at `2,3`, primitivity, parity and the two mod-`24`
cones only enlarges the counted set, so the THM-3645 consequence has no
hidden local gap.  This is an upper bound only.  It proves no matching lower
bound, asymptotic, Euler-product constant, primitive Moebius-tail estimate,
complete Pell class, compiler orbit, positive two-cube ray, or density among
actual two-cube representations.  **QED.**
