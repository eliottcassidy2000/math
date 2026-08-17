# Keller inertia and the Berggren wall share one signed `C4` class

**Status: VERIFIED-EXACT cohomological cospan over proved THM-3536 and
THM-3537; no LRC-current or D5-flux identification.**

THM-3537's old-`L` inertia has cycle type

```text
(4)(2)(1)^3.                                               (1)
```

The quartic orbit is a genuine four-sheet cyclic carrier.  Its determinant
line has monodromy `-1`; in mod-two notation it is the nonzero class in

```text
H^1(D*,F2)=Hom(pi_1(D*),F2)=F2.                           (2)
```

This gives a lawful Jacobian-side `H^1` object.  It is the sign local system
of a normalized inverse cover, not a physical flux.

## 1. One class, three gauges

Put a bit `1` on a negative edge sign.  A local system on a four-edge cycle
is determined up to vertex switching by the XOR of its four bits.  In a
one-cut gauge the Keller quartic determinant class is

```text
k=(1,0,0,0).                                               (3)
```

THM-3536's Berggren wall class is

```text
b=(1,0,1,1),                                               (4)
```

and the raw Fibonacci projective cycle has

```text
f=(0,1,0,0).                                               (5)
```

All three have odd XOR.  The vertex gauge

```text
v=(0,0,0,1)                                               (6)
```

satisfies `k+delta(v)=b`.  Two complementary gauges transport `f` to `b`:

```text
(1,0,1,0),             (0,1,0,1).                        (7)
```

Thus the cospan is exact at the level of the unique nonzero signed-`C4`
class:

```text
Keller quartic determinant  -> [1] <- Berggren wall <- raw Fibonacci.
                                                                  (8)
```

The arrows forget different data.  On the Keller side they forget the
quadratic orbit, three fixed sheets, ramification indices, and coordinate
orders.  On the Berggren side they forget ancestry labels, angle chambers,
and density.  Equality in `(8)` is not an identification of the underlying
trees.

## 2. XOR is exactly the old-prime parity cancellation

The determinant bits of the five inertia orbits in `(1)` are

```text
quartic:  4-1 = 3 = 1 mod 2,
quadratic:2-1 = 1 = 1 mod 2,
fixed:                         0,0,0.                    (9)
```

Therefore

```text
1 xor 1 xor 0 xor 0 xor 0 = 0.                           (10)
```

Equation `(10)` is the cohomological form of THM-3537's even old-`L`
discriminant exponent.  The class is not absent locally: two nonzero orbit
classes cancel in the determinant of their direct sum.  At the newest prime
there is only one transposition orbit, so its bit is `1`; this is why the
newest divisor survives in THM-3531's square class.

This also clarifies the three cubic discriminants `-4 A_i^2 L`.  Each
primitive cubic observation sees the same one-step determinant line.  Their
common square class records that line, not a classification of the maps or
coordinates producing it.

## 3. Four vertices, six pairs, two missing edges

The quartic orbit has successor edges

```text
0->1, 1->2, 2->3, 3->0.                                  (11)
```

Among the six unordered vertex pairs, the two antipodal pairs

```text
{0,2}, {1,3}                                             (12)
```

are unobserved by the successor relation.  Orienting each missing pair gives
exactly four tournament completions.  None is invariant under the phase
rotation `i->i+1`, because two rotations reverse either antipodal edge.

This is the same six-pair carrier shape found for the Berggren signed `C4`:
four observed cyclic transitions plus two missing antipodals.  It is not a
canonical tournament.  The cover supplies `(11)`, not orientations for
`(12)`.

## 4. What this contributes to the D5 program

The cospan `(8)` supplies an explicit target class for any proposed
word-current/JC-flux map:

```text
source word cocycle -> determinant-line w_1 of inverse monodromy. (13)
```

Any lawful map must pass two hostile gates.

1. It must respect direct sums: the quartic and quadratic nonzero classes in
   `(9)` must cancel globally as in `(10)`.
2. It must carry the missing-edge sidecar: a `C4` class alone cannot choose a
   tournament completion or recover the lost quadratic orbit.

THM-3534's provisional LRC cospan does not yet provide such a source map.  It
has a static middle duality, while `(13)` lives on the determinant line of a
punctured-divisor monodromy representation.  No chain map between those
complexes has been constructed.  Consequently this reflection is not the D5
map, but it sharpens the target and supplies a decisive cancellation test.

## Connection contract

| field | exact answer |
|---|---|
| Keller source | quartic inertia orbit at old `L` |
| Keller class | `w_1(det Q[O_4])=1` in `H^1(D*,F2)` |
| Berggren source | four-phase wall sign word |
| common quotient | unique nonzero signed-`C4` switching class |
| explicit gauge | `(0,0,0,1)` from Keller one-cut to wall word |
| global Keller XOR | quartic `1` plus quadratic `1` equals `0` |
| pair carrier | four successor pairs plus two missing antipodals |
| tournament status | four completions, none rotation-invariant |
| destroyed | orbit multiplicity, ancestry, amplitudes, coordinate index |
| still missing | chain map from LRC word-current complex to determinant line |

Reproduce the finite cochain, switching, XOR, and completion census with

```text
python -B 04-computation/keller_inertia_signed_c4_h1_xor_bridge_20260816.py
python -B -O 04-computation/keller_inertia_signed_c4_h1_xor_bridge_20260816.py
```

Both transcripts match the stored output exactly.
