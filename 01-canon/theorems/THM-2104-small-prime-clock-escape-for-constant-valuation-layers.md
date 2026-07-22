---
id: THM-2104
title: "Every constant small-prime quotient layer has a universal clock escape"
status: >
  PROVED. Let p,q be independent integer characters, take the guard g=p,
  and write every terminal as c_i=a_i p+n_i q with n_i nonzero. If, for
  any ell in {2,3,5}, all n_i have one common ell-adic valuation v, then
  p.X=1/2 and q.X=1/(2*ell^(v+1)) gives guard distance 1/2 and terminal
  distances at least 1/(2ell)>1/14. For ell=7 the same construction is
  sharp: it gives at least 1/14, and gives at least 1/7 if every normalized
  quotient unit avoids +/-1 modulo 7. Consequently a cover by the guard and
  these transverse bands must have quotient-valuation diversity at each of
  2,3,5. This statement does not control any additional guard-proportional
  bands and does not prove LRC(14).
source: codex-2026-07-22-LRC-small-prime-clock-valuation-layers
depends_on:
  - THM-2103
related:
  - THM-1065
  - THM-2069
  - THM-2073
  - THM-2095
  - THM-2099
  - THM-2120
---

# THM-2104 -- small-prime clocks for constant valuation layers

Let `p,q:T^2->T` be independent integer characters, so the torus map

```text
X |-> (p.X,q.X)                                         (1)
```

is surjective. Let

```text
g=p,
c_i=a_i p+n_i q,              a_i in Z, n_i in Z\{0}.  (2)
```

Fix a prime `ell` and suppose there is a common integer `v>=0` such that

```text
nu_ell(|n_i|)=v                         for every i.     (3)
```

Then there is an `X in T^2` satisfying

```text
||g.X||=1/2,
||c_i.X|| >= 1/(2ell)                    for every i.    (4)
```

For `ell in {2,3,5}`, (4) is strictly larger than `1/14`. Thus `X` is
strictly safe for a radius-`1/7` guard and all the radius-`1/14` terminals
in (2).

## Proof

Put `b=ell^v`. By surjectivity in (1), choose `X` with

```text
p.X=1/2,                       q.X=1/(2ell b).          (5)
```

Write `m_i=n_i/b`. Hypothesis (3) says that `m_i` is an integer not divisible
by `ell`. Equation (2) gives

```text
c_i.X=(a_i ell+m_i)/(2ell)                  in T.       (6)
```

The numerator in (6) is nonzero modulo `ell`. Its residue modulo `2ell` is
therefore neither `0` nor `ell`; in particular its circular distance from
`0` modulo `2ell` is at least one. Dividing by `2ell` proves (4). Finally,

```text
1/4 > 1/14,             1/6 > 1/14,             1/10 > 1/14,
```

for `ell=2,3,5`, respectively. QED.

## The sharp prime-seven boundary

The same proof with `ell=7` gives only

```text
||c_i.X|| >= 1/14.                                      (7)
```

This is exactly the terminal radius, so it does not by itself produce a
strict escape. The equality cases are transparent. If

```text
m_i=n_i/7^v is not congruent to +/-1 mod 7               (8)
```

for every `i`, then the numerator in (6) has circular residue at least two
modulo fourteen, independently of the parity of `a_i`. Hence

```text
||c_i.X|| >= 2/14=1/7.                                  (9)
```

Thus a strict cover whose quotient coefficients occupy one `7`-adic layer
must contain a normalized unit congruent to `+1` or `-1` modulo seven. The
prime seven is the exact endpoint of this one-clock argument; for larger
primes the universal lower bound `1/(2ell)` falls below `1/14`.

## Dyadic and multi-orbit corollaries

The earlier dyadic tower has

```text
n_i=b_0 2^(e_i).                                        (10)
```

All coefficients in (10) have the same `3`-adic and `5`-adic valuations.
Either the third-clock (`ell=3`) or fifth-clock (`ell=5`) construction gives
a strict escape. More generally, the theorem permits arbitrary unit factors,
signs, repetitions, and any number of multiplicative orbits. The clock point
depends only on one constant small-prime valuation layer.

The six THM-2103 exceptions were therefore neither sporadic geometries nor a
phenomenon specific to powers of two. They lie in full `3`-adic and `5`-adic
valuation layers on which a threshold-labelled residue evaluation is
uniformly safe.

## Coefficient-plane consequence and exact scope

For a torus cover whose complete band list consists of the guard `p` and the
transverse terminals in (2), the theorem forces all three inequalities

```text
#{nu_ell(|n_i|):i} >= 2,              ell=2,3,5.        (11)
```

This is a basis-aware statement. Replacing `q` by another transverse
character changes the quotient coordinates `n_i`; one may use (11) only in
a named integral two-character presentation. The invariant formulation is
that the terminal images in the rank-one quotient lattice by `Zp` cannot all
lie in one `ell`-adic valuation shell for any `ell in {2,3,5}`.

There are two scope guards.

1. An additional terminal proportional to `p` has `n_i=0` and is not covered
   by the proof: at `p.X=1/2` an even multiple of `p` may be unsafe. Therefore
   (11) is not, without another argument, a certificate against a larger row
   containing extra guard-proportional bands.
2. When `Zp` is not saturated in the ambient character lattice, the quotient
   has torsion. The theorem still applies whenever the decompositions (2)
   exist; passing to a saturated quotient requires retaining the torsion class
   as a sidecar rather than silently identifying it with an integer coordinate.

## Frontier effect

The carrier ladder sharpens to

```text
pair weights
 -> signed affine rank
 -> small-prime clocks
 -> the joint (nu_2,nu_3,nu_5) quotient profile
 -> unit residues at the prime-seven boundary.          (12)
```

THM-2103 found the sixth-clock residue empirically in a finite cube. The
present theorem makes it all-height and reveals a three-prime obstruction:
any all-transverse cover must cross a valuation wall at each of `2,3,5`.
THM-2095's p-adic deck warning remains essential because layer
multiplicities and unit residues, not only prime support, must be retained.

The next lawful target is to combine the three simultaneous valuation-wall
requirements with THM-2069's deletion-code/cogirth constraints. Merely
counting valuation signatures cannot suffice: two coefficients, for example
`1` and `30`, already vary at all three primes. A successful deletion argument
must use which terminals witness each wall and how those witnesses persist
after peeling, not just the number of occupied shells.

## Assumption challenge and Tournament Analysis

The challenged assumption is that powers of two cause the sixth-clock escape.
They only make the `3`- and `5`-adic layers constant. The actual invariant is
constant quotient valuation at any small prime.

Candidate tournament vertices were terminals, valuation signatures, valuation
walls, unit residues, clock states, and deletion obligations. Orienting
terminals by one valuation gives a complete tie on the branch closed here;
any tie-break tournament then has artificial scores, cycles/SCCs, edge flips,
and Hamiltonian paths. A more faithful future tournament has **wall witnesses**
as vertices and directs one witness to another when deleting the first erases
the second's valuation contrast. Its quotient preserves persistence of (11)
under peeling but destroys the actual clock residues, which must remain a
sidecar.

This closes only constant small-prime layers in the all-transverse model. It
does not bound the number or separation of layers, control added
guard-proportional bands, treat quotient torsion without its sidecar, discharge
the vertical branch, or prove LRC(14).
