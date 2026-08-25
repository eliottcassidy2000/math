---
id: THM-4076
title: "Stern inverse-orbit mod-eight law and apex-degree congruence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. For every odd
  q>1, the Stern depth packet satisfies
  S(q)=-phi(q)+4*1_{omega(q)=1} (mod 8). Its divisor-star imbalance obeys
  B(q)=1-q+4*Omega(q) (mod 8), so the maximal-vertex indegree and outdegree
  are respectively 2*Omega(q) and q-1-2*Omega(q) modulo 4. Here omega and
  Omega count distinct and repeated prime factors. The proof is an inversion-
  orbit count; it gives exact necessary gates for packet and apex balance but
  neither a zero classification nor the mod-sixteen residue from fixed-point
  data alone.
source: codex-frontier-synthesis-creative-20260825d / Stern 2-adic orbit lane
audit: >
  PASS. The primary inverse-parity path checks every odd q through 4999,
  including 2,499 packets, 6,247,500 root-scan candidates, 630,947 nonfixed
  even-even inversion-orbit pairs, and 10,071 divisor-star terms. An
  independent Euclidean-depth path checks every odd q through 2999, performs
  1,823,912 primitive depth evaluations and 4,497,000 direct lower-star depth
  evaluations, reconstructs the half-box modular hyperbola without importing
  the primary path, and counts both degree classes directly. Normal and
  optimized outputs byte-match; both scripts have zero assert nodes and zero
  floating literals. The q=13/q=29 hostile proves that the orbit parity does
  not determine the next mod-sixteen bit.
depends_on:
  - THM-4059-stern-brocot-depth-packet-character-and-divisor-star-convolution
related:
  - THM-4061-stern-depth-modular-hyperbola-and-prime-apex-balance
  - THM-4071-stern-prime-power-stationary-phase-and-all-odd-apex-balance
script: 04-computation/stern_inverse_orbit_mod8_apex_congruence_thm4076.py
output: 05-knowledge/results/stern_inverse_orbit_mod8_apex_congruence_thm4076.out
script_sha256: 173a9d9dd5b346231327708e38eb485b325c6b9b5a44f0fe04a7de265669c33b
output_sha256: ab1aaf164b419e87b0bf1891824cd5cd63c13da009d6604a5dcc6ef87031a14e
independent_audit_script: 04-computation/stern_inverse_orbit_mod8_apex_congruence_thm4076_independent_audit.py
independent_audit_output: 05-knowledge/results/stern_inverse_orbit_mod8_apex_congruence_thm4076_independent_audit.out
independent_audit_script_sha256: 9e4236a284482217fea72218352d408e0c744dae44db3b7d85bd2b27769e92e0
independent_audit_output_sha256: cb2c269d88c4eed90fc874a08414b5adb6fb768823ec6f63966d8e5e1e03f4ef
hash_basis: raw LF bytes
---

# THM-4076 -- one more exact bit from Stern packet inversion orbits

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** THM-4059
identifies the odd-denominator Stern sign with the parity of a unit and its
modular inverse. Retaining the fixed points of that inversion, rather than
only pairing numerator parities, determines one further binary digit of every
packet, every divisor-star imbalance, and both apex degrees.

## 1. Statement

For an odd integer `q>1`, retain THM-4059's packet and lower-star notation

```text
S(q)=sum_(a in U_q) (-1)^(a+a^(-1)),
B(q)=sum_(d|q,d>1) S(d),                                      (1)

indeg(q) =(q-1+B(q))/2,
outdeg(q)=(q-1-B(q))/2.                                      (2)
```

Write

```text
omega(q)=#{p:p|q},             Omega(q)=sum_(p^e || q) e.     (3)
```

Then the complete mod-eight laws are

```text
boxed: S(q) = -phi(q)+4*1_(omega(q)=1)       (mod 8),         (4)
boxed: B(q) = 1-q+4*Omega(q)                 (mod 8).         (5)
```

Consequently the maximal vertex of the initial THM-4057 tournament has

```text
boxed: indeg(q)  = 2*Omega(q)                 (mod 4),
boxed: outdeg(q) = q-1-2*Omega(q)             (mod 4).        (6)
```

These formulas include every odd prime power and every odd composite; no
squarefree hypothesis is present.

## 2. The inversion-orbit mechanism

Let `N(q)` be the number of units `a`, represented in `1<=a<q`, for which
both `a` and its canonical inverse `u=a^(-1) mod q` are even. Negation sends

```text
(a,u) -> (q-a,q-u).                                           (7)
```

Because `q` is odd, `(7)` interchanges even-even with odd-odd and even-odd
with odd-even. Thus the four parity boxes have sizes `N,M,M,N` for some
`M`, whence

```text
phi(q)=2N(q)+2M,
S(q)  =2N(q)-2M=4N(q)-phi(q).                                (8)
```

This reproves the THM-4061 half-box identity at exactly the resolution needed
here: writing `a=2r,u=2s`, the even-even box is precisely

```text
N(q)=#{1<=r,s<=(q-1)/2 : 4rs=1 mod q}.                       (9)
```

Now let modular inversion act on the finite even-even set. Every nonfixed
orbit has two elements. Its fixed elements are exactly the even canonical
solutions of

```text
x^2=1 mod q.                                                  (10)
```

For each odd prime power `p^e`, the only roots of `(10)` are `+1` and `-1`.
The Chinese remainder theorem therefore gives exactly `2^omega(q)` roots
modulo `q`. Negation pairs them, has no fixed point for odd `q`, and flips
their canonical parity. Exactly half, namely

```text
2^(omega(q)-1),                                               (11)
```

are even. Hence

```text
N(q) congruent 2^(omega(q)-1)
     congruent 1_(omega(q)=1)                 (mod 2).         (12)
```

Substituting `(12)` into `(8)` proves `(4)`. In particular `(4)` reduces to
THM-4059's complete mod-four law on odd denominators, but its extra bit comes
from the fixed-point sidecar `(11)`.

## 3. Divisor stars and apex degrees

Sum `(4)` over all nontrivial divisors of `q`. The standard totient identity
and a direct divisor count give

```text
sum_(d|q,d>1) phi(d)=q-1,                                    (13)

#{d|q:d>1 and omega(d)=1}
  =sum_(p^e || q) e=Omega(q).                                (14)
```

Equations `(13)` and `(14)` yield `(5)`. Substitution in `(2)` is legitimate
modulo four because the numerators are known modulo eight, and gives `(6)`.

The distinction between `omega` and `Omega` is load-bearing: inversion fixed
points see the support of one packet, while the divisor-star sum counts every
prime-power divisor and therefore sees multiplicity.

## 4. Exact zero gates and their failure boundary

The congruences give the following necessary conditions.

1. If `S(q)=0` and `omega(q)=1`, then `q=p^e` with `p=5 mod 8`.
2. If `S(q)=0` and `omega(q)=2`, at least one of the two distinct prime
   divisors is `1 mod 4`.
3. If `B(q)=0`, then

   ```text
   q=1+4*Omega(q) mod 8.                                    (15)
   ```

For `omega(q)>=3`, Euler's totient is already divisible by eight, so `(4)`
adds no zero obstruction. These gates are not classifications. For example,
`q=25` passes the prime-power gate but `S(25)=-8`; `q=15=3*5` passes the
two-prime gate but `S(15)=8`; and the same `q=25` passes `(15)` but has
`B(25)=-8`.

Oddness is load-bearing. At `q=2`, THM-4059 has `S(2)=-1`, whereas formally
substituting the odd-modulus right side of `(4)` would give `3 mod 8`.
Negation no longer flips canonical parity, so the parity-box proof itself
identifies the failed implication.

Nor does the fixed-point parity quotient determine `S(q) mod 16`. Put

```text
F(q)=2^(omega(q)-1),       T(q)=(N(q)-F(q))/2,                 (16)
```

so `T(q)` is the number of nonfixed even-even inversion-orbit pairs. In fact,

```text
S(q)=-phi(q)+4F(q)+8T(q).                                     (17)
```

The prime controls `q=13` and `q=29` have the same `q mod 16`, `phi(q) mod
16`, `omega`, `Omega`, and fixed-root count, yet

```text
N(13)=3 mod 4,      S(13)=0 mod 16,
N(29)=1 mod 4,      S(29)=8 mod 16.                           (18)
```

Thus the next bit requires `T(q) mod 2`, equivalently `N(q) mod 4`: genuine
modular-hyperbola data that the fixed-point parity quotient discards. This is
a sharp boundary for the present quotient, not a claim that no different
mod-sixteen law can exist.

## 5. Connection and loss ledger

The source is THM-4059's inverse-parity packet. The map decomposes the
even-even set into inversion orbits and then takes its cardinality modulo two:
every two-orbit contributes zero and every fixed point contributes one. The
preserved predicate is `N(q) mod 2`, and its fixed-point sidecar is the CRT
root count `(11)`. Passing from a packet to its divisor star retains one
further coordinate, the number `Omega(q)` of prime-power divisors. The parity
projection destroys numerator locations and the nonfixed-orbit count `T(q)`;
`(17)`--`(18)` prove that `T(q) mod 2` is already necessary modulo sixteen.

The result supplies no estimate for `|S(q)|` or `|B(q)|`, no sufficient zero
criterion, and no new Mahler, even-graph, or LRC consequence.

## 6. Exact audits

The primary path evaluates `(1)` from the THM-4059 inverse-parity formula,
enumerates the roots `(10)`, checks `(8)` and every congruence through
`q=4999`, and applies all zero gates. The independent path reconstructs every
sign from Euclidean continued-fraction depth, counts `(9)` in half-box
coordinates, computes `B(q)` directly over all lower vertices before comparing
the divisor sum, and counts the two degree classes directly through `q=2999`.
It imports no primary implementation.

Replay from the repository root:

```bash
python3 -B 04-computation/stern_inverse_orbit_mod8_apex_congruence_thm4076.py
python3 -B -O 04-computation/stern_inverse_orbit_mod8_apex_congruence_thm4076.py
python3 -B 04-computation/stern_inverse_orbit_mod8_apex_congruence_thm4076_independent_audit.py
python3 -B -O 04-computation/stern_inverse_orbit_mod8_apex_congruence_thm4076_independent_audit.py
```

The normal and optimized outputs byte-match their frozen companions. Both
scripts replace optimization-sensitive assertions by explicit runtime gates
and contain no floating-point literal.
