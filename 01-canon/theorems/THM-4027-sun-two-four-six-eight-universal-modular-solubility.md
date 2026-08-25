---
id: THM-4027
title: "Sun 2-4-6-8 universal modular solubility"
status: >
  PROVED + FINITE-EXACT SMALL-PRIME SIDECAR + INDEPENDENTLY AUDITED. Every
  residue modulo every positive modulus is represented by the admissible
  2-4-6-8 binomial sum. Thus THM-4026 has no congruence obstruction. Exact
  local densities are positive but nonuniform; N lies in the unique worst
  class 20 mod 33, with relative factor 16/33.
source: root + independent local/modular and Pascal-period audits, 2026-08-24
depends_on: []
related:
  - THM-2412-delta-exponential-and-central-newton-layer-split
  - THM-2642-cyclic-difference-relation-saturation-and-thick-holotopy-no-go
  - THM-4026-sun-two-four-six-eight-binomial-counterexample
  - THM-4028-sun-two-four-six-eight-average-order-criticality
  - MISTAKE-363
script: 04-computation/sun_2468_universal_modular_solubility_thm4027.py
independent_audit_script: 04-computation/sun_2468_universal_modular_solubility_independent_audit_thm4027.py
independent_audit_output: 05-knowledge/results/sun_2468_universal_modular_solubility_independent_audit_thm4027.out
local_density_script: 04-computation/sun_2468_local_density_thm4027.py
local_density_independent_audit_script: 04-computation/sun_2468_local_density_independent_audit_thm4027.py
local_density_independent_audit_output: 05-knowledge/results/sun_2468_local_density_independent_audit_thm4027.out
direct_crosscheck_script: 04-computation/sun_2468_local_density_direct_crosscheck_thm4027.py
direct_crosscheck_output: 05-knowledge/results/sun_2468_local_density_direct_crosscheck_thm4027.out
---

# THM-4027 -- every congruence class is soluble

**PROVED + FINITE-EXACT SMALL-PRIME SIDECAR + INDEPENDENTLY AUDITED.** For
every positive integer `m` and every residue `r mod m`, there are integers

\[
w\ge2,\qquad x\ge3,\qquad y\ge5,\qquad z\ge7
\]

such that

\[
\binom w2+\binom x4+\binom y6+\binom z8\equiv r\pmod m.       \tag{1}
\]

In particular, the exact counterexample in THM-4026 is not explained by any
congruence obstruction, even one combining arbitrarily many fixed moduli.

## 1. Inheritance and the typed connection

THM-2412 supplies the lawful integer-valued Pascal coordinate

```text
Delta C(n,k)=C(n,k-1).                                  (2)
```

THM-2642 is the closest proved use of Cauchy--Davenport saturation. Here the
vertices are not runners or atom values: the intrinsic objects are the four
**role-labelled value sets** of the functions `C(n,k)` over a finite cyclic
input period. MISTAKE-363 is load-bearing: for `p<=k`, reducing
`n(n-1).../k!` by a modular inverse is undefined. Every finite sidecar below
uses exact integer binomials or Pascal addition.

The source-to-target map is

```text
(w,x,y,z) |-> sum_k C(t_k,k) mod m,  k in {2,4,6,8}.    (3)
```

It preserves congruence existence and exact periodic multiplicity. It
destroys magnitude, equality to a fixed integer, the bounded search box, and
the root in the triangular discriminant. THM-4026 is the canonical hostile:
`(3)` is onto for every `m`, yet one positive integer has no exact preimage.

## 2. A denominator-safe period theorem

Let `p` be prime, `a>=1`, and `k>=1`. Put

```text
e=floor(log_p k),            P=p^(a+e).                 (4)
```

Then `P` is a period of `n -> C(n,k) mod p^a`. Indeed, Vandermonde gives

\[
\binom{n+P}{k}-\binom nk
=\sum_{j=1}^k\binom Pj\binom n{k-j}.                   \tag{5}
\]

For `1<=j<=k<P`,

\[
v_p\!\binom Pj=v_p(P)-v_p(j)\ge a,                    \tag{6}
\]

so every term in `(5)` vanishes modulo `p^a`. The independent period audit
also proves minimality: if an arbitrary positive shift `h` is a period, take
finite differences of `(5)` in the binomial basis to force
`p^a|C(h,j)` for all `j<=k`; choosing `j=p^min(v_p(h),e)` forces
`p^(a+e)|h`. Thus `(4)` is the exact least period.

For a composite modulus, the lcm of the prime-power periods is a period. For
the four roles, the extra exponents `e` are

| `p` | `k=2` | `k=4` | `k=6` | `k=8` |
|---:|---:|---:|---:|---:|
| 2 | 1 | 2 | 2 | 3 |
| 3 | 0 | 1 | 1 | 1 |
| 5 | 0 | 0 | 1 | 1 |
| 7 | 0 | 0 | 0 | 1 |
| `p>=11` | 0 | 0 | 0 | 0 |

This table is exactly where a naïve `p^a`-sample would corrupt the small
primes.

## 3. Odd-prime base coverage

For an odd prime `p`, let

```text
A_2^reg={C(w,2): w in F_p, 2w-1 != 0},                 (7)
```

and let `A_k` be the image of `C(t,k)` modulo `p`. The involution
`w -> 1-w` and

```text
C(u,2)=C(v,2)  <=>  (u-v)(u+v-1)=0                    (8)
```

give

```text
|A_2^reg|=(p-1)/2.                                      (9)
```

If `p>8`, then `C(t,k)` is a nonconstant degree-`k` polynomial over `F_p`.
Every fibre has at most `k` elements, hence

```text
|A_k|>=ceil(p/k).                                       (10)
```

Iterated Cauchy--Davenport now yields

\[
|A_2^{reg}+A_4+A_6+A_8|
\ge\min\!\left(p,{p-1\over2}+\left\lceil{p\over4}\right\rceil
+\left\lceil{p\over6}\right\rceil+\left\lceil{p\over8}\right\rceil-3\right).
                                                               \tag{11}
\]

The second argument is at least

```text
25p/24-7/2,                                             (12)
```

so `(11)` is `p` for every prime `p>=89`.

The remaining odd primes are a finite exact sidecar. Using the complete
periods from `(4)`, the companion constructs a representation with regular
triangular coordinate for every residue at

```text
3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,
61,67,71,73,79,83.                                      (13)
```

It freezes one witness per target residue and checks the same finite gate by
an independent Pascal-state implementation. Thus every residue modulo every
odd prime has a solution with `2w-1` nonzero.

The reciprocal-degree sum is not decorative. In the general typed packet
`{2=d_1,d_2,...,d_s}` the same argument has lower slope

```text
sum_i 1/d_i.                                            (14)
```

When this exceeds one, large-prime saturation follows after finitely many
small-prime gates. For `2,4,6,8`, the excess is only `1/24`; THM-4028 shows
that this is also the mean-growth exponent.

## 4. Lifting every odd prime power

Fix a target integer `R` and a regular solution modulo `p`. Hold `x,y,z`
fixed and write

```text
G(W)=C(W,2)+C(x,4)+C(y,6)+C(z,8)-R.                    (15)
```

For `h=t p^a`, direct integer expansion gives

\[
G(w+h)-G(w)\equiv t p^a{2w-1\over2}\pmod {p^{a+1}}.    \tag{16}
\]

The coefficient is a unit modulo the odd prime. Therefore exactly one
`t mod p` corrects the next target digit. Induction lifts the representation
to every `p^a`; no factorial denominator from the higher roles is inverted.

## 5. The prime two

For `q=2^a`,

\[
\binom{w+2q}{2}-\binom w2=q(2w+2q-1)\equiv q\pmod{2q}. \tag{17}
\]

Thus pairing `w` with `w+2q` preserves its residue modulo `q` and toggles
exactly the new binary digit. Since triangular numbers hit both residues
modulo two, induction proves that `C(w,2)` alone is onto modulo every `2^a`.
In fact it is exactly two-to-one over its period `2^(a+1)`, so convolving the
other roles leaves the complete two-adic distribution exactly uniform.

An independent connection supplies a second route: the prime-power Newton
coordinate `n -> C(n,p^e) mod p^a` is balanced on its exact period. Taking
`p=2,e=3` makes the octic role `C(z,8)` itself uniform. The triangular-toggle
and balanced-octic proofs are independent controls on the delicate
denominator prime.

## 6. Coordinatewise CRT and the lower bounds

Write `m=prod_i p_i^(a_i)` and choose one local solution at every prime
power. For each role `k`, impose its chosen index modulo

```text
p_i^(a_i+floor(log_(p_i) k)).                           (18)
```

For distinct `p_i` these periods are coprime. Coordinatewise CRT therefore
produces a single quadruple whose binomial values agree with every local
quadruple. Adding a common positive multiple of the composite period to each
index preserves all congruences and enforces `w>=2,x>=3,y>=5,z>=7`. This
proves `(1)`.

## 7. Exact local density: rarity without obstruction

For `q=p^a`, let `P_(q,k)` be the exact periods and define

\[
\sigma_q(r)=q\,{\#\{(t_k)\in\prod_k\mathbf Z/P_{q,k}:
\sum_k\binom{t_k}k\equiv r\pmod q\}\over\prod_kP_{q,k}}. \tag{19}
\]

This is the factor relative to a uniform target residue; its mean over `r`
is one. The universal theorem says it is always positive, not that it equals
one.

For the counterexample `N`, exact complete-period convolution gives

| modulus | `N mod q` | `sigma_q(N)` | minimum |
|---:|---:|---:|---:|
| 3 | 2 | `22/27` | `22/27` at 2 |
| 9 | 8 | `68/81` | `62/81` at 5 |
| 11 | 9 | `72/121` | `72/121` at 9 |
| 33 | 20 | `16/33` | `16/33` at 20 |
| 99 | 53 | `544/1089` | `496/1089` at 86 |

Thus `N` belongs to the unique worst class modulo `33`, exactly the class
identified as sparse in the 2019 MathOverflow discussion. The deeper
modulo-99 class `86` is thinner than the target's actual class `53`; prime-
power refinement need not monotonically improve a search prior. Every entry
is still positive.

For the pairwise-coprime packet

```text
Q=9*11*25*49*13*17*19*23=11712375675,                  (20)
```

the target's exact factor is

```text
5359658991616/22179335403225 = 0.241651018580... .      (21)
```

This is an exact periodic density and a useful search ranking, not a
probability of nonrepresentation.

## 8. Smooth-factor boundary

At a prime `p>8`, if the hypersurface fibre

```text
C(w,2)+C(x,4)+C(y,6)+C(z,8)=r                          (22)
```

has no point where all four first derivatives vanish, then each solution has
exactly `p^3` lifts for each prescribed next target digit. Consequently

```text
sigma_(p^a)(r)=sigma_p(r) for every a>=1.               (23)
```

The exact audit verifies smoothness for the target at
`p=11,13,17,19,23,29,43`. The hypothesis is essential. At `p=31` there are
four singular target points,

```text
(16,17,8,11), (16,17,8,27),
(16,17,28,11),(16,17,28,27),                            (24)
```

and the factor changes from `942/961` at `31` to `29198/29791` at `31^2`.
This is the required Hasse-jet hostile against blind lift stabilization.

## 9. Reproduction and scope

Run the primary and independent modular proofs with

```text
python -B 04-computation/sun_2468_universal_modular_solubility_thm4027.py
python -B -O 04-computation/sun_2468_universal_modular_solubility_thm4027.py
python -B 04-computation/sun_2468_universal_modular_solubility_independent_audit_thm4027.py
python -B -O 04-computation/sun_2468_universal_modular_solubility_independent_audit_thm4027.py
```

and the local-density companions named in the frontmatter. The direct
Cartesian audit reproduces the selected histogram convolutions without using
convolution.

This theorem settles all fixed congruence quotients. It proves no pointwise
asymptotic, no exact integer representation, and no finite bound for a global
preimage. The missing coordinate is archimedean alignment; THM-4026 proves
that it can fail even when every local quotient succeeds. **QED.**
