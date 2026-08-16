# Three coordinate cubics, one sign class

**Status: proved elementary deduction from THM-2546 for the fixed sporadic
map; classification synthesis only beyond that fixed map.**  This note does
not extend the discriminant identity to the outside infinite family.

## 1. The three discriminants collapse to one square class

Let

```text
K = Q(a,b,c)
```

be the generic target function field of the fixed sporadic Keller map.  The
THM-2546 referee computation proves, for its three coordinate cubics,

```text
Delta_x = -4 A_x^2 L,
Delta_r = -4 A_r^2 L,
Delta_z = -4 A_z^2 L,                                (1)
```

with the three explicit nonzero polynomial square factors recorded there.
Since `4` and every `A_i^2` are squares in `K`, equation (1) gives

```text
[Delta_x]=[Delta_r]=[Delta_z]=[-L] in K*/K*2.          (2)
```

Equivalently, on the generic locus,

```text
K(sqrt(Delta_x)) = K(sqrt(Delta_r))
                 = K(sqrt(Delta_z))
                 = K(sqrt(-L)).                       (3)
```

Thus the three coordinate eliminants have one common quadratic sign
resolvent.  Their finite affine odd-valuation carrier is the same Jelonek
divisor `V(L)`.  A statement about the complete branch divisor would also
have to audit valuations at infinity.

The equality `-4=-(det JF)^2` explains why the constant disappears in square
class.  It does not turn (2) into a degree-four statement: `4` is a square
coefficient, while the surviving information is the binary class `[-L]`.

## 2. Exactly what the common class remembers

For a separable cubic, its discriminant square class is the sign character of
the permutation action on its three roots.  Consequently (2) synchronizes
the three sign quotients.  It does **not** identify:

- the three cubic root sets;
- their full splitting fields or labelled `S3` actions;
- the rational maps between coordinate roots;
- boundary-effective polynomial sections; or
- the multiplicities of the discriminant divisor.

Different cubic extensions can share the same quadratic resolvent.  The
lawful picture is therefore a common-source cospan,

```text
                         x-coordinate cubic
                       /
[-L] sign class  ----  r-coordinate cubic
                       \
                         z-coordinate cubic,          (4)
```

not three pairwise identifications.

This also identifies the correct tournament boundary.  There is no intrinsic
pairwise orientation among the three coordinate cubics: all three expose the
same binary observable.  A three-vertex tournament would manufacture order
that (2) does not contain.  The useful finite carrier is one binary hub with
three views, plus the sidecars lost by the sign quotient.

## 3. Why the `V4`/four-object analogy stops

A generic `V4` quartic-resolvent package has three nonzero binary characters.
The exact effectivity hostile in
`jc-torsor-killing-versus-boundary-effectivity-exact-hostiles-codex-20260814.md`
shows that even all three Kummer characters and their class relations can be
present while natural rational views fail to be polynomial.

Equation (2) is a different phenomenon: the three displayed cubic views land
on the **diagonal**

```text
([-L],[-L],[-L]) in (K*/K*2)^3.                       (5)
```

So neither the three nonzero characters of `V4` nor a tournament on four
vertices may be read off from the three discriminants.  Recovering such a
four-object carrier would require an explicit resolvent map and an oriented
boundary-valuation sidecar.  The numerical coincidence with `-4` supplies
neither.

## 4. Composition proves that the class is dynamical

THM-2576 proves for the fixed map

```text
S_(F o F)=V(LH),                                      (6)
```

with two distinct irreducible Jelonek components.  THM-2582 then proves that
the degree-nine **x-coordinate** eliminant of `F o F` has square class

```text
[Disc(E_2)]=[H],                                      (7)
```

not `[LH]`: the old `L` contribution cancels in parity between the norm
denominator and the odd-block alternating factor.  Hence the atom-level
class `[-L]` is not a static label that can simply be multiplied along the
Keller monoid.  At the first composite rung the sign quotient selects the
new image component instead.

The all-degree norm-product lemma explains the parity mechanism, but only
the fixed-map, level-two x-coordinate application is proved.  In particular,
(7) does not prove a three-coordinate common-class law for `F o F`, nor a
newest-component law at every iterate.

## 5. Relation to the outside infinite family

THM-1605 verifies that the outside family `E_m` is strictly broader than the
repo's constructed examples: `m=2` is the fixed sporadic map, while the
higher members realize fibre cardinalities

```text
2m-1 = 1+2(m-1).                                     (8)
```

THM-1350 explains the necessity of the shape in (8) for the relevant
equivariant fibre: one fixed point plus free involution pairs.  This is an
orbit-count theorem, not a discriminant theorem.  Nothing currently proved
in the repo shows that the higher `E_m` have three coordinate eliminants with
one common class, the same Jelonek geometry, or the same composition tower.

Therefore (2) contributes one exact invariant at `m=2`; it does not classify
the infinite family.  A serious classification record needs at least

```text
degree/fibre grade
+ fixed/free involution orbit profile
+ full monodromy group and block systems
+ quadratic sign-resolvent class
+ Jelonek divisor and its image tower
+ oriented boundary valuations/effective cone.        (9)
```

The first and second entries can agree while the remaining four differ.
Conversely, equality of the sign class alone is far too coarse to identify a
map: it retains only one `C2` quotient of the root permutation data.

## 6. Connection contract and decisive next tests

| field | exact content |
|---|---|
| source | three coordinate cubic discriminants of the fixed sporadic map |
| target | the single class `[-L]` in `K*/K*2` |
| map | take discriminant and quotient by squares |
| preserved | root-permutation sign and finite affine odd divisor `V(L)` |
| destroyed | labelled roots, full monodromy, multiplicities, poles at infinity, effectivity |
| required sidecar | splitting-field/block data plus oriented boundary valuations |
| cheapest hostile | two nonisomorphic cubic extensions with the same discriminant square class |

The next exact tests are correspondingly typed:

1. compute the other two coordinate square classes for `F o F` and ask
   whether the atom-level diagonal (5) survives at grade nine;
2. compute the level-three image divisor and norm law, then test whether only
   the newest component has odd valuation;
3. for one explicit `E_m` with `m>2`, compute a generic coordinate eliminant,
   its sign class, and the first Jelonek divisor before proposing any family
   law; and
4. keep the sign class separate from the effective boundary cone throughout.

The strongest current conclusion is therefore precise and modest: the fixed
sporadic map has three cubic views of one quadratic branch character, while
its exact classification—and still more the classification of the outside
infinite family—requires the information discarded by that character.
