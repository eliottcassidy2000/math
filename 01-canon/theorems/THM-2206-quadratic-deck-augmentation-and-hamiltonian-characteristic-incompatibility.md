---
id: THM-2206
title: "Quadratic-deck augmentation and Hamiltonian characteristic incompatibility"
status: >
  PROVED algebraic obstruction plus an explicitly OPEN continuation target.
  If Delta=sigma-1 in k[C_2], then Delta^j=(-2)^(j-1)Delta. In
  characteristic different from two the deck algebra is semisimple and every
  augmentation power has the same image and kernel, so higher deck/Hasse
  differences contain only parity. In characteristic two the augmentation
  ideal is square-zero, but polynomial squares are Hamiltonian-central and
  the quadratic root algebra is not an etale C_2-deck; the quartic square
  prefix and Faber flux package disappear. Thus no field characteristic
  simultaneously supplies a nontrivial higher quadratic-deck jet and the
  square-prefix Hamiltonian information used in the quartic Jacobian
  reductions. Over Z_2, by contrast, I^j=2^(j-1)I and the filtration is
  separated. Turning that identity into a Keller proof requires an unproved
  integral order-raising lemma; none is asserted here.
source: codex-2026-07-24-quadratic-deck-characteristic-audit
depends_on: []
related:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2189-nonsplit-quartic-deck-forces-the-remaining-pole-congruence
  - THM-2194-uniform-degree-six-quartic-pole-closure
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2202-uniform-all-degree-quartic-pole-closure
---

# THM-2206 -- quadratic-deck characteristic incompatibility

The triangular Hasse carrier in THM-2201 works because the cyclic group has
prime order equal to the coefficient characteristic. It is tempting to copy
that mechanism to the quadratic deck in the quartic Jacobian reductions. The
following exact calculation shows why no field supports the desired copy.

Let `k` be a field, let

```text
C_2=<sigma>,             sigma^2=1,
A=k[C_2],                Delta=sigma-1.                (1)
```

## 1. Exact augmentation law

In every characteristic,

```text
Delta^2=(sigma-1)^2=-2Delta,                           (2)
Delta^j=(-2)^(j-1)Delta,              j>=1.            (3)
```

Equation (3) follows immediately from (2) by induction.

### Characteristic different from two

Assume `char(k)!=2` and put

```text
e_+=(1+sigma)/2,             e_-=(1-sigma)/2.          (4)
```

These are orthogonal idempotents with sum one, and

```text
A=ke_+ direct_sum ke_-,
Delta e_+=0,                 Delta e_-=-2e_-.          (5)
```

Consequently, on every `A`-module `M`,

```text
M=e_+M direct_sum e_-M,
im(Delta^j)=e_-M,            ker(Delta^j)=e_+M         (6)
```

for every `j>=1`. Equivalently, the augmentation ideal `I=(Delta)` satisfies

```text
I^j=I,                       j>=1.                     (7)
```

Thus iterated deck differences are scalar rescalings of the first parity
difference. They cannot form a higher jet or expose information hidden from
the invariant/anti-invariant split.

### Characteristic two

Assume `char(k)=2`. Then

```text
A=k[Delta]/(Delta^2),             I^2=0.               (8)
```

The abstract deck algebra now has a nilpotent augmentation coordinate, but
the Hamiltonian carrier degenerates. In the polynomial Poisson algebra
`k[x,z]`, for all `H,G`,

```text
d(H^2)=0,                 {H^2,G}=0.                   (9)
```

Hence squares lie in the Poisson center, and

```text
{H^2+B,-}={B,-}.                                        (10)
```

The exact square prefix is invisible to the Hamiltonian derivation. Moreover,
the algebra `K[U]/(U^2-V)` is never an etale quadratic deck in characteristic
two: its defining polynomial is inseparable, it becomes nonreduced after
splitting, and the proposed sign involution is the identity because `-U=U`.
Thus (8) is not the augmentation algebra of a nontrivial geometric quadratic
deck carrying the original square-prefix data.

Equations (5)--(10) prove:

> **Field-characteristic obstruction.** No field characteristic provides
> both a nontrivial higher augmentation-depth coordinate for a quadratic deck
> and faithful Hamiltonian access to a polynomial square prefix.

This is an incompatibility theorem, not a statement that mixed-characteristic
methods are impossible.

## 2. Exact interface with the quartic Jacobian canon

THM-2129 works over a characteristic-zero differential field. For

```text
P_0=z^4+pz^2+qz+r,                 alpha=n/4,          (11)
```

it expands `P_0^alpha` and defines

```text
Phi_n=4c_1,               Psi_n=4c_2,
R_n=4c_3+pc_1,                                            (12)
```

with the exact Hamiltonian identity

```text
L(E_n)=(z^2+p/4)Phi_n'+zPsi_n'+R_n'.                   (13)
```

In characteristic two, `n/4` and `p/4` are undefined, the displayed
fourfold flux normalization vanishes, and (9)--(10) erase the square carrier.
Thus this package does not descend to the nilpotent algebra (8).

THM-2189 instead works over the genuine characteristic-zero quadratic field

```text
L=K(U),                    U^2=V,
sigma(U)=-U.                                               (14)
```

It uses

```text
sigma(E_m)=(-1)^m E_m                                    (15)
```

to delete odd Faber seeds. On the nonsplit field deck, anti-invariant
constants vanish because the constant field is fixed. On the split etale
deck `K x K`, however, the constants contain the nonzero anti-invariant line

```text
{(c,-c):c in C}.                                         (16)
```

Equations (3), (6), and (16) show that higher powers of `Delta` cannot repair
that split constant-flux obstruction: they only rescale it.

THM-2194 closes reduced degree six by retaining the ordered boundary,
first-flux, and second-flux valuation phases rather than relying on parity;
this distinction is load-bearing on the split deck. THM-2202 extends that
phase-bank closure to every twice-odd reduced degree, so no quartic finite-pole
survivor remains. Neither theorem is a hidden higher augmentation-jet
argument. The remaining quartic target is the terminal nonmonic square-prefix
descent.

## 3. The mixed-characteristic survivor

Let

```text
R=Z_2,               A_R=R[C_2],             I=(Delta). (17)
```

Equation (2) now gives the genuine descending carry filtration

```text
I^j=2^(j-1)I,                     j>=1,
intersection_(j>=1) I^j=0.                               (18)
```

The second equality uses `2`-adic separatedness of the finite free
`R`-module `A_R`. More generally it holds after acting on any finite free
`A_R`-lattice. As an abstract positive control, on `R x R` with the swap
action,

```text
I(R x R)={(a,-a):a in R},
I^j(R x R)=2^(j-1)I(R x R).                            (19)
```

This is the only one of the three coefficient regimes that retains both
characteristic-zero differentiation and a separated notion of quadratic-deck
depth. It is not an integral splitting theorem for `R[U]/(U^2-V)`: at the
prime two the factors `U-u,U+u` need not be comaximal, so the actual root
algebra can remain ramified/non-etale. Equation (18) by itself proves no
Keller constraint.

The honest open target is:

> **Integral order-raising target (OPEN).** Construct an integral,
> `2`-torsion-free, `I`-adically separated lattice of denominator-cleared
> square-prefix/Faber observables such that every applicable boundary and
> flux constraint sends an anti-invariant defect
>
> ```text
> xi in I^j M       to       xi in I^(j+1) M.           (20)
> ```

If (20) held uniformly, then `xi` would lie in every `I^jM` and (18) would
force `xi=0`. Tensoring with `Q_2` destroys this route because `2` becomes a
unit and `I^j=I`; reducing modulo two destroys it because (9) erases the
Hamiltonian square prefix.

## 4. Cheapest decisive test and boundary

The least expensive hostile test is the exact degree-six bank of THM-2194.
Retain

```text
E_6,E_5,E_3,E_2,E_1
```

together with its boundary, `Phi`, and `Psi` rows; clear the rowwise powers of
two; and compute the induced elimination map channelwise on

```text
I^j A_R/I^(j+1) A_R isomorphic_to F_2.                 (21)
```

If eliminating any odd seed fails to raise the surviving anti-invariant flux
by one `I`-order, the proposed route fails already in degree six. Success is
only a finite positive control. A proof would still need:

1. compatible clearing of all `n/4` Faber denominators;
2. `2`-torsion-freeness of the observable lattice;
3. an all-degree integral identity stable under arbitrary complex
   coefficient specialization;
4. control of ramification in the actual integral quadratic root algebra; and
5. the uniform order-raising implication (20).

The exact connection ledger is

```text
source:       quadratic deck parity in the quartic square-prefix reduction;
map:          powers of the augmentation ideal;
preserved:    invariant/anti-invariant split over characteristic zero;
destroyed:    higher depth over fields of odd characteristic,
              Hamiltonian square data in characteristic two;
needed sidecar:
              an integral 2-adic observable lattice and Keller order raising;
decisive test:
              the degree-six map on I^j/I^(j+1).                       (22)
```

Therefore ordinary field-valued Hasse jets cannot repair the terminal
square-prefix descent. The `Z_2[C_2]` filtration is a precise, nonvacuous
research target, but it remains OPEN. QED for the obstruction.
