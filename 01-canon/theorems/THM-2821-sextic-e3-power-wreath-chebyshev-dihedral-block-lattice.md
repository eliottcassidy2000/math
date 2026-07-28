---
id: THM-2821
title: "Sextic e=3 power-wreath and Chebyshev-dihedral block lattice"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  THM-2817's power carrier has monodromy
  (C3 x C3) semidirect C2 of order 18 and a unique nontrivial block
  system, of block size three.  Its only nontrivial rational
  decomposition type is therefore degree three followed by degree two.
  The Chebyshev carrier has dihedral monodromy C6 semidirect C2 of order
  12 and exactly two nontrivial block systems, of sizes two and three.
  Besides the common degree-three-then-two decomposition, its evenness
  gives an exact degree-two-then-three decomposition.  This classifies one
  rational response's intermediate fields, not polynomial Keller-map
  decomposition, JC(2), or DC(2).
source: root/sextic-e3-monodromy-block-lattice-2026-07-28
depends_on:
  - THM-2817-sextic-e3-maximal-pole-power-chebyshev-accessory-classification
related:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2816-maximal-pole-clean-nielsen-ribbon-tree-prufer-vandermonde-count
script: 04-computation/jc_sextic_e3_monodromy_block_lattice_thm2821.py
output: 05-knowledge/results/jc_sextic_e3_monodromy_block_lattice_thm2821.out
script_sha256: 4e5cfd3c4286e2a5c409184ee4e977a04c85e2ea22aa7933d708efae82554f46
output_sha256: f3bc6b1384c528d8de8bad8c0eb7f11ae7c5b52cff38d02bdf38cb72551f6318
hash_basis: LF-normalized bytes
---

# THM-2821 -- the two sextic carriers have different 2/3 block lattices

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2817 proves that the minimal `e=3` maximal-pole response layer has only
two unmarked rational maps.  Both visibly have degree three followed by
degree two.  Their full intermediate-field lattices are nevertheless
different: the power map has only that decomposition, whereas the
Chebyshev map also factors in the reverse degree order.

## 1. Exact branch-cycle representatives

Use the right-to-left permutation convention on six sheets and fix

```text
rho=(0 1 2 3 4 5).                                      (1)
```

The two noncrossing zero-inertia representatives are

```text
tau_P=(0 1)(2 3)(4 5),
tau_C=(0 1)(2 5)(3 4).                                  (2)
```

Their pole permutations have cycle types

```text
type(tau_P rho)=(3,1,1,1),
type(tau_C rho)=(2,2,1,1).                              (3)
```

These are precisely THM-2817's power and Chebyshev pole multisets.  Hence
their geometric monodromy groups are

```text
G_P=<rho,tau_P>,                  G_C=<rho,tau_C>.       (4)
```

The five noncrossing perfect matchings of the hexagon form exactly two
rotation orbits.  Directly computing `(3)--(4)` on every matching gives the
same two group/block signatures proved below, so the conclusion is
independent of the chosen representatives.

## 2. The power carrier is the imprimitive wreath group of order 18

Put

```text
a=rho tau_P=(0 2 4),
b=tau_P rho=(1 3 5).                                   (5)
```

The two three-cycles have disjoint support, commute, and generate

```text
K=<a,b> isomorphic to C3 x C3.                          (6)
```

Moreover

```text
tau_P a tau_P=b,                 rho=a tau_P.           (7)
```

Thus

```text
G_P=K semidirect <tau_P>
   isomorphic to (C3 x C3) semidirect C2
   =C3 wreath C2,                         |G_P|=18.      (8)
```

The involution exchanges the two `C3` factors.

There is exactly one nontrivial block system:

```text
{0,2,4} | {1,3,5}.                                    (9)
```

Here is a proof that does not rely on the finite enumeration.  A proper
block in a transitive degree-six action has size two or three.  If a block
containing `0` also contains an odd sheet, the independent cycle `b`,
which fixes `0` and rotates all three odd sheets, forces that block to
contain every odd sheet; this already has size at least four.  Hence the
block lies in `{0,2,4}`.  The cycle `a` then forces the unique possibility
`{0,2,4}`; in particular no size-two block exists.

## 3. The Chebyshev carrier is dihedral of order 12

The second involution satisfies

```text
tau_C rho tau_C=rho^(-1).                               (10)
```

Therefore

```text
G_C=<rho> semidirect <tau_C>
   isomorphic to C6 semidirect C2,
|G_C|=12,                                               (11)
```

with the nontrivial `C2` acting by inversion.  This is the dihedral group
of the hexagon (order twelve).

The regular cyclic subgroup `<rho>` makes the complete block calculation
transparent.  A block containing `0` is a subgroup of the regular
`Z/6Z` sheet coordinate and must be stable under inversion.  The only
proper nonzero choices are

```text
{0,3},                    {0,2,4}.                     (12)
```

Thus the Chebyshev carrier has exactly two nontrivial block systems:
three blocks of size two and two blocks of size three.

## 4. Complete rational decomposition classification

Let

```text
phi(z)=z^2/(z^2-1),
P_3(x)=2x^3-1,
T_3(y)=4y^3-3y.                                       (13)
```

THM-2817 gives

```text
F_P=phi(P_3(x)),                 F_C=phi(T_3(y)).       (14)
```

These are the block-size-three decompositions in `(9)` and `(12)`:
degree-three inner map followed by degree-two outer map.

For the Chebyshev response, put `u=y^2`.  Since

```text
T_3(y)^2=u(4u-3)^2,
```

one also has the exact reverse-order factorization

```text
F_C=psi(y^2),

psi(u)=u(4u-3)^2/[u(4u-3)^2-1],                      (15)
```

where `deg(y^2)=2` and `deg(psi)=3`.  This realizes the block of size two
in `(12)`.

For a separable rational map over characteristic zero, intermediate fields
between `C(F)` and `C(x)` are equivalent, by Lüroth's theorem, to
monodromy block systems containing a chosen sheet.  Equations `(9)` and
`(12)` are the complete block lists.  Consequently:

```text
F_P: exactly the nontrivial degree pattern 3 then 2;
F_C: exactly the nontrivial degree patterns 3 then 2 and 2 then 3,       (16)
```

up to Möbius changes of the intermediate coordinate.

## 5. What the two/three co-occurrence does and does not mean

This is a literal binary/ternary statement on one rational map:

- in both carriers, the ternary inner coordinate is retained before the
  common binary outer response;
- only the Chebyshev carrier permits the order to be reversed, because the
  full response is even in its centered source coordinate;
- the difference is detected by the complete monodromy block lattice, not
  by a slogan about degrees.

It is not an action of

```text
C2*C3=PSL2(Z).
```

The actual monodromy groups are the finite groups `(8)` and `(11)`, and the
intermediate coordinate is a load-bearing sidecar.

More importantly, a block system of the one-variable rational response does
not furnish a decomposition of a two-variable polynomial Keller map.  Such
a transfer would have to preserve the Faber multiplier, the source
polynomial, and the constant-Jacobian one-form.  None of those data is
encoded by `(9)` or `(12)`.  Thus `(16)` proves no instance of `JC(2)` or
`DC(2)` and gives no counterexample.

## 6. Exact companion

The companion uses exact permutation and symbolic rational arithmetic to:

1. verify the two passports `(3)` and enumerate both generated groups;
2. exhibit the `C3 x C3` kernel and semidirect law `(5)--(8)`;
3. verify the dihedral relation and normal form `(10)--(11)`;
4. exhaust every subset containing sheet zero and recover exactly the
   blocks `(9)` and `(12)`;
5. enumerate all five noncrossing matchings and their two rotation orbits;
6. verify both factorizations `(14)` and the reverse factorization `(15)`;
   and
7. verify normal, optimized, and stored transcript identity.

It contains no Python `assert` node.  Run

```text
python 04-computation/jc_sextic_e3_monodromy_block_lattice_thm2821.py
python -O 04-computation/jc_sextic_e3_monodromy_block_lattice_thm2821.py
```

The finite computation is exhaustive in the declared degree-six
permutation universe.  Higher degree, other `e=3` passports, positive
characteristic, polynomial Keller-map decomposition, `JC(2)`, and `DC(2)`
remain outside the theorem.
