---
id: THM-2839
title: "Prime-power unit-mass full spectrum and q11 response provenance boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For N=p^d, every integral cyclic mask whose total mass is prime to p
  is a unit in Q[C_N] and F_p[C_N].  Hence all N Fourier characters are
  nonzero, its N translates span the full regular module, and any linear
  equivariant realization observing that mask has dimension at least N.
  More generally, augmentation detects units in the modular group algebra
  of every finite p-group; nonabelian Fourier blocks are matrix-valued.
  A coset-supported cyclic mask has the equivalent cyclotomic-at-one proof.
  Applied to THM-2835's 449-sheet q11 horn in one residue mod169, this
  gives full ancestry spectrum and rank 13^5=371293.  The exact
  20-state response, 44-state one-step-memory, and 260-state
  residue-response graphs are nondeterministic quotients; H->C loses
  provenance, and only the unrotated H response is physical.  THM-2835's
  E3/endpoint-gauge clutch remains empty.  This is a signed algebraic
  reconstruction theorem, not a positive physical multiplier, current
  action, row exclusion, or LRC(14) proof.
source: >
  root/q11-outer-word-scan-2026-07-28;
  independent root/audit-2809-augmentation-unit-2026-07-28
depends_on:
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
related:
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - THM-2845-local-residue-versus-split-trace-unit-observability
script: 04-computation/lrc14_prime_power_unit_mass_q11_response_thm2839.py
output: 05-knowledge/results/lrc14_prime_power_unit_mass_q11_response_thm2839.out
script_sha256: 68ae72b62b7974e4f2c2bf7d570615c8c524746978c57cf120f6372a7250ece4
output_sha256: 495829603ea0c3944f83d7ae269cbbc5cbdec9fdc452395e96c78de8b2e7e24b
hash_basis: LF-normalized bytes
---

# THM-2839 -- prime-power unit mass gives full cyclic spectrum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2835 proves the exact 449-sheet q11 semantic horn and its E3/gauge
boundary.  The new theorem here is not another horn theorem.  It extracts
a reusable algebraic principle from the horn's mass:

```text
integer mass prime to p on a p-power orbit
  => group-algebra unit
  => every Fourier character survives
  => full translate/Hankel rank.                          (1)
```

For the q11 horn this rank is `13^5=371293`, not the `13` or `26` of
THM-2835's aggregate q-response modules.  The price is sharp: the inverse
is signed, and canon does not supply the positive physical multiplier or
typed ancestry action needed to use it.

## 1. Prime-power unit-mass lemma

Let

```text
p be prime,              d>=1,              N=p^d,

A(Z)=sum_(j in Z/NZ) a_j Z^j in Z[C_N],

s=A(1)=sum_j a_j.                                         (2)
```

Assume

```text
p does not divide s.                                     (3)
```

Then:

1. `A` is a unit in `F_p[C_N]`;
2. `A` is a unit in `Q[C_N]`;
3. `A(xi) != 0` for every complex character `xi^N=1`;
4. the `N` translates `Z^j A` form a basis of `Q[C_N]`;
5. the cyclic convolution/circulant matrix of `A` has rank `N`;
6. every equivariant scalar realization

   ```text
   a(q)=O rho(q) B,             q in C_N,                 (4)
   ```

   over a characteristic-zero field has `dim rho>=N`.

No coset-support assumption is needed for this form.

### Proof

In characteristic `p`,

```text
F_p[C_N]
 =F_p[Z]/(Z^N-1)
 =F_p[e]/(e^N),                    e=Z-1,                 (5)
```

because `Z^(p^d)-1=(Z-1)^(p^d)`.  The augmentation ideal `(e)` is
nilpotent.  Write

```text
A=s+n,                    n in(e).                        (6)
```

Since `s` is a nonzero element of `F_p`,

```text
A^(-1)
 =s^(-1) sum_(k=0)^(N-1)(-s^(-1)n)^k
                                                   in F_p[C_N].  (7)
```

Thus the integral multiplication matrix `M_A` has determinant nonzero
modulo `p`; in particular its determinant is a nonzero integer.  Hence
`M_A` is invertible over `Q`, proving (1)--(2) and the translate-basis
claim.

Over `C`, Fourier diagonalization of `M_A` has eigenvalues

```text
A(xi),                    xi^N=1.                         (8)
```

Invertibility makes every eigenvalue nonzero and proves (3)--(5).
For `(4)`, the cyclic Hankel matrix factors as

```text
H_(u,v)=a(u+v)
       =(O rho(u))(rho(v)B),                              (9)
```

so `rank H<=dim rho`.  Up to column reversal, `H` is `M_A`; it has rank
`N`.  Therefore `dim rho>=N`.

### Finite p-group extension

The unit statement is not restricted to cyclic groups.  Let `G` be any
finite `p`-group and let

```text
epsilon:F_p[G]->F_p
```

be augmentation.  Its kernel `I` is nilpotent.  One elementary proof is by
induction on `|G|`: choose a central element `z` of order `p`, put
`u=z-1`, and note

```text
u^p=0,                  F_p[G]/(u)=F_p[G/<z>].          (10)
```

If the augmentation ideal of the quotient has nilpotence exponent `m`,
then

```text
I^m subset (u),                  I^(mp)=0.              (11)
```

Thus `I` is the Jacobson radical and

```text
A in F_p[G] is a unit  iff  epsilon(A)!=0.              (12)
```

For an integral, or more generally `p`-integral, lift with `p`-unit
augmentation, regular multiplication is consequently invertible over
`Q`, with the sharper congruence

```text
det L_A
 =epsilon(A)^|G|
 =epsilon(A)                         mod p.              (13)
```

For finite abelian `p`-groups every scalar complex character value of `A`
is nonzero.  For nonabelian `G`, the correct statement is that every
matrix Fourier block is nonsingular; it is not a scalar-character claim.
In either case the group translates span the full regular module.

For a rational mask, the safe projective formulation is: clear
denominators to its primitive integral coefficient vector and require
that vector's augmentation to be prime to `p`.  Merely asking that the
original rational augmentation be a `p`-unit is false.  For example,

```text
(1/p) sum_(r=0)^(p-1) Z^(r p^(d-1))                    (14)
```

has augmentation one but the same cyclotomic zeros as the subgroup sum.
Primitive-integral scaling preserves the Fourier zero set and regular
rank.

### Coset/cyclotomic form

Suppose additionally that the support lies in one coset

```text
c+p^r Z/NZ.                                               (15)
```

Then

```text
A(Z)=Z^c G(Z^(p^r)),               G(1)=s.                (16)
```

At an `N`th root, the argument of `G` has order `p^t` for some
`0<=t<=d-r`.  If `t>=1` and `G` vanished there, its integral minimal
polynomial `Phi_(p^t)` would divide `G`, forcing

```text
p=Phi_(p^t)(1) divides G(1)=s,                            (17)
```

contrary to `(3)`.  At order one the value is `s`.  This proves full
spectrum directly and records the effective root order on each
valuation stratum.

## 2. Sharp boundaries

The prime-to-`p` condition is sharp as a universal hypothesis.  The
one-`p`-cycle mask

```text
1+Z^(p^(d-1))+...+Z^((p-1)p^(d-1))                       (18)
```

has mass `p` and vanishes at every character which is nontrivial on that
cycle.  Divisibility of the mass by `p` permits zeros; it does not force
them.

The `p`-power group hypothesis is equally load-bearing.  If `ell!=p`
divides a cyclic order `N`, the `ell`-subgroup sum has augmentation `ell`,
a `p`-unit, but vanishes on characters nontrivial on that subgroup.  The
smallest control is

```text
N=6, p=2,             A=1+Z^2+Z^4:

A(1)=3 is odd,        A(primitive cube root)=0.
```

Thus the theorem is a finite-`p`-group theorem, not an arbitrary finite
group theorem parametrized by one chosen prime.

The inverse supplied by `(7)` is algebraic and generally signed.  If a
nonnegative mask with at least two support points had a nonnegative
convolution inverse, one positive inverse support point would pair with
the two mask support points and create positive mass at two distinct
sums, contradicting a delta convolution.  Thus:

```text
nonnegative mask + nonnegative inverse
    => mask support is a singleton.                       (19)
```

This is the load-bearing physical boundary.  Faithful signed recovery is
not a lawful Boolean/positive multiplier.

The finite-group extension also resolves an apparent Heisenberg rank
paradox.  For the order-`p^3` Heisenberg group in THM-2779, a genuinely
group-indexed `p`-unit-mass mask would pay the full regular rank `p^3`.
The natural faithful transitive carrier has only `p^2` points because it
is `G/H` for an order-`p` stabilizer.  Its lifted stabilizer-coset mask
has augmentation `p`, lies exactly on the hostile boundary `(18)`, and
has regular convolution rank `p^2`.  The lost stabilizer coordinate,
not a rank contradiction, separates the two carriers.

## 3. THM-2835 horn application

Import from proved THM-2835 the fixed q11 horn label set

```text
L={a:S_QB(a)=T_QA(11,a)=1} subset Z/13^5Z,

|L|=449,                min L=59306,       max L=311961,

a=156 mod169             for every a in L.                (20)
```

Put

```text
h_L(Z)=sum_(a in L) Z^a.
```

Since

```text
449=7 mod13,                                                (21)
```

the lemma gives

```text
h_L is a unit in Q[C_(13^5)] and F_13[C_(13^5)],

h_L(xi)!=0 for all 371293 ancestry characters xi,

rank circulant(h_L)=371293.                               (22)
```

The coset form is especially explicit:

```text
h_L(Z)=Z^156 G(Z^169),

min supp G=350,          max supp G=1845,

G(1)=449.                                                  (23)
```

The exact frequency strata are

| ancestry frequency | count | effective order in `G` |
|---|---:|---:|
| 13-unit | 342,732 | 2,197 |
| valuation 1 | 26,364 | 169 |
| valuation 2 | 2,028 | 13 |
| valuation at least 3 | 169 | 1 |

Every coefficient in every row is nonzero.  In the modular proof,

```text
F_13[C_(13^5)]=F_13[e]/(e^371293),

augmentation(h_L)=7,       augmentation inverse=2,

det circulant(h_L)=7 mod13.                               (24)
```

Thus the horn loses no ancestry coordinate in signed group-algebra
terms.  Its unique rational convolution inverse recovers `delta_0`.
The absent datum is multiplier access to `h_L` inside the typed physical
current, exactly the carrier-versus-observable distinction also exposed
by THM-2815 in the finite Laguerre setting.

## 4. Exact response provenance

Encode each depth-five ancestry label by its three source-word bits and
all `13*3` target-word bits from the THM-2835 three-word atlas.  The exact
response has 20 codes.  The horn is precisely the complete state

```text
H=628342390802,                 frequency 449.             (25)
```

Every `H` has successor

```text
C=628342690404.                                           (26)
```

But `C` has global frequency 10,328 and ten predecessor states.
The marked edge `H->C` has multiplicity 449.

Set-theoretically, a predecessor bit recognizes the horn arrivals.
Equivalently,

```text
{a:response(a)=C and a=157 mod169}=L+1,       size449.     (27)
```

The low residue `a=1 mod13` selects 898 `C` labels and is insufficient.
The response quotients have:

```text
response graph:             20 states, 44 edges,
one-step-memory graph:      44 states, 68 edges,
normal residue response:   260 states, 476 edges.          (28)
```

They are nondeterministic quotients, not smaller equivariant
realizations of `1_L`.  In particular, the least cyclic period of the
complete response, the horn indicator, and the marked `H->C` edge is
`13^5`.  For the horn this also follows immediately from `(22)`.

There is no physical cyclic 13-orbit of `H`.  Rotate its thirteen target
response rows.  All thirteen formal codes are distinct, but their exact
frequencies in the physical ancestry bank are

```text
(449,0,0,0,0,0,0,0,0,0,0,0,0).                          (29)
```

Only the unrotated profile exists.  A formal regular `C13` capacity does
not supply the missing physical basepoint.

## 5. Why this does not duplicate THM-2835 or close LRC

THM-2835's rational ranks `13` and `26` concern the aggregate
thirteen-residue q-response after the depth-five ancestry phase has been
forgotten.  Equation `(22)` concerns the indicator of one exact
depth-five horn subset.  There is no contradiction:

```text
aggregate q-response module:          rank 13 or 26;
exact ancestry horn observable:       rank 371293.         (30)
```

THM-2835 also proves that every one of the 63 repaired q11 cells loses
macro factor `E3` at q7, no one of all 567 cells repairs the full
`I/J0/J7` packet, and the outer-word gauge and endpoint-positive gauge
have empty fibre product.  The present companion rechecks those fixture
facts to prevent a false spectral-to-physical promotion.

The theorem supplies:

1. a reusable prime-power unit-mass/full-spectrum lemma;
2. a sharp p-cycle hostile and positive-inverse boundary;
3. exact full ancestry spectrum and rank for the THM-2835 horn;
4. an exact provenance invoice for the 20/44/260 quotients; and
5. a precise missing operation: lawful positive access to the horn
   multiplier while preserving semantic word, E3, endpoint origin, and
   current phase.

It supplies no positive inverse, physical ancestry-shift action,
same-current word transition, E3 transport, row exclusion, ledger
decrement, or LRC(14) proof.

## 6. Exact evidence

Run

```text
python 04-computation/lrc14_prime_power_unit_mass_q11_response_thm2839.py
python -O 04-computation/lrc14_prime_power_unit_mass_q11_response_thm2839.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_prime_power_unit_mass_q11_response_thm2839.out.
```

The companion reconstructs the inherited THM-2835 fixture, all 371,293
depth-five response labels, exact response/pair/normal-residue graphs,
all possible 13-power period shifts, the mod-169 provenance selector,
the full-spectrum valuation strata, the characteristic-13 unit
certificate, all thirteen cyclic q-rotations of `H`, and the inherited
E3/gauge-clutch hostiles.

## 7. Independent hostile audit

The independent algebra audit proved the cyclic lemma from the local
group-ring model, rederived the determinant congruence and Hankel-rank
factorization, and extended the argument to every finite `p`-group by the
central-order-`p` induction above.  It separately checked cyclic groups
of orders `2,4,8,3,9,5`, the nonabelian group of order eight, and
`C_3 x C_3`.  Subgroup-sum controls have exact rank `|G|/p`; `p delta_0`
has full rank and confirms that prime-to-`p` mass is sufficient but not
necessary.

The same audit found and repaired the only load-bearing scope hazard:
arbitrary rational masks cannot be classified by their displayed
augmentation without `p`-integrality or primitive denominator clearing.
It also proved the bi-positive inverse boundary `(19)`.  A separate full
replay reconstructed all `371,293` ancestry labels; normal, optimized,
and stored transcripts agree, and the LF-normalized hashes in the
frontmatter match.

**QED.**
