# Independent audit of the `U_full` owner-boundary `K4 x F_13` endpoint factorization

**Verdict: SOUND, with a repaired packaging defect and a sharper endpoint
boundary statement.**  The mathematical claims in candidate commit
`e92508714abec0f054e033ef933e04b61ae3f1e1` survive an independent central
tensor reconstruction and are promoted as
`THM-3514-ufull-owner-boundary-k4xf13-endpoint-factorization-and-walsh-spectrum.md`.
The result remains an endpoint factorization in the coboundary/`B^1` lane.
It is not a common-ancestry realization, physical current, nonzero absolute
`H^1` class, scalar-row exclusion, or LRC(14) theorem.

## 1. Inheritance and audit board

The closest proved mechanism is THM-3479's fixed `U_full` refined endpoint
bank.  The canonical hostile is early common-sheet marginalization, which
kills every primitive translated-sheet mode.  The corrected near miss is the
minimum-joint-address checkerboard: equal endpoint marginals do not determine
their pairing.  The least-used live sidecar is the speed-13 owner-in boundary
already present in `PATTERN_E`.

The audit kept five live objects distinct:

| object | faithful representation | decisive predicate | lost coordinate / cheapest test |
|---|---|---|---|
| owner factor | strict circle inequality plus half-open interval representative | middle chamber empty before Fourier summation | singleton boundary convention / compare owner-only intervals |
| guard factor | 39 translated sheet/chamber atoms | exact `g_C(a+tau)` covariance | common sheet / restore before summing |
| endpoint bank | 169 materialized left/right atom tables | reconstruct all `13^3` frozen values | pair support / direct guarded controls |
| drift tensor | 117 `(C,D,d)` types | exact 52-support and nonzero bucket values | absolute sheet / retain common-sheet phase |
| THM-2471 stalk | Boolean `(y,u,v,a,b,e')` fibre | common support before endpoint marginalization | absent map / reject Cartesian pairing as ancestry |

This used the META-PATTERNS cards “refine and saturate before transporting a
shadow” and “join compatibility fibers before comparing marginal sizes.”  The
first succeeds at the endpoint level; the second remains the stopping
boundary.

## 2. Source-level owner law and the boundary repair

In current order,

```text
w=(1,183,27,131,53,313,13,2197,742586).
```

Write `s=91t=7a+r`, with `0<=r<7`.  The owner factor in `PATTERN_E` is

```text
distance(13t,Z)<1/14.
```

Because `13t=a+r/7` modulo one, literal strict point support is

```text
[0,1/2) union (13/2,7).
```

The endpoint engine uses the measure-equivalent representative

```text
[0,1/2) union [13/2,7).
```

The distinction is thirteen singleton left endpoints `r=13/2`, one in each
sheet.  The audit did not infer this from output.  It built the owner-only
`PATTERN_E` primitive and compared it exactly with the thirteen cyclic
intervals

```text
[14a-1,14a+1) in units T_DEN/182.
```

Thus the candidate's half-open statement is the correct endpoint-engine
representation, while the strict statement must retain the open left boundary
in ordinary point-set language.  Endpoint integrals and their exponential
boundary differences are unchanged by those singletons.

For the guard partition

```text
[0,1), [1,6), [6,7),
```

the engine assigns `r=1` to the middle chamber and `r=6` to the right chamber.
Direct midpoint and cut calculations recover the danger arcs

```text
left   {12,0,1},
middle {11,12,0,1},
right  {11,12,0}.
```

The owner factor kills all thirteen middle atoms before any transform.  The
actual endpoint support therefore uses 26 atoms.

## 3. Independent aggregation route

The audit imports only the promoted THM-3479 endpoint primitives, pinned at

```text
ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b.
```

It imports no candidate guard-sheet or drift-bucket companion.  Its flow is
different from the candidate worker accumulation:

1. remove only the H guard for each `(alpha,beta)`;
2. use a separately written two-pointer intersection with the contiguous
   39-atom partition;
3. materialize every `AX` and `BY` atom table;
4. restore `g_C(a+tau)` centrally for every `tau`;
5. compare 78 reconstructed endpoint pairs against direct guarded endpoint
   evaluation;
6. reconstruct the frozen gamma bank and inverse role values;
7. build the 1,521 guard-pair kernels directly;
8. contract all endpoint atom products centrally into the 117 drift buckets;
9. compute the four corner rows, Walsh transform, row rank, and drift DFT.

The 169 unguarded sets contain `7,107,008` intervals before atom splitting and
`7,108,460` after splitting.  The complete reconstructed bank has digest

```text
1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682
```

and exact split-field values

```text
q_H    320618948602619577408,
q_q5   503604956476841920373,
bridge 389266878372286537904.
```

These agree with the frozen THM-3479 bank.

## 4. Exact bucket and spectral verdict

All `1,521/1,521` abstract pair kernels are nonzero for each role frequency.
After actual endpoint support is imposed, exactly the 65 types containing a
middle chamber vanish.  The remaining

```text
52=(left/right)^2 x F_13
```

types are nonzero separately for `q_H`, `q_q5`, and their difference.  The
independent bucket digests are exactly

```text
H      c5bb013836492665e766a0423a71c45983f5226b99b2ea74ad7f53e33c2a8629
q5     933c6f0b768f22ee993dc02af171d3b731e0c1ab3acd7120bed157a491a29ec9
bridge 553cfe7289b0556a19a8bcd1a0382dc1372545358feac62e5229adca315f8a26
```

The four corner rows have rank four.  The constant, left, right, and mixed
Walsh rows are nonzero at all thirteen drifts, and every row has all thirteen
nonzero drift-Fourier coefficients.  Their digests are

```text
corner c651f1deb71258001213c32b812005fb1f07b3a834fbe0526feb5b9c1c54435a
Walsh  6f8ce43a58c5c05552fa4a815c955e10873ed3535ba8cb436284fe4ad18d996c
DFT    ef4c512d3e0444cce1cacee42ebbbb2d2edef7bd866517814484d1e8e921eaa9
```

Same sheet, same chamber, and same guard atom give bridge values

```text
324498447313453607031,
341590019311254368788,
 50212089322152546627,
```

none equal to the full bridge.  These are useful geometric hostiles, not
ancestry substitutes.

## 5. `K4` carrier, root-label gauge, and exact loss ledger

The four chamber-pair states are the four vertices of `F_2^2`, hence of an
undirected complete graph `K4`.  The Walsh rows are its constant, two
coordinate, and checkerboard channels.  No intrinsic binary orientation is
present, so this is not a tournament.

Incoming commit `4e6021ca3` sharpened the root-label debt.  Both the guard
sheet and the THM-2471 root label carry regular `+tau` translation.  Every
equivariant label bijection is

```text
u=a+c, c in F_13.
```

There are thirteen gauges.  A common gauge on both endpoint legs preserves
`u_R-u_L=a_R-a_L`; independent gauges shift drift by `c_R-c_L`.  Changing a
common gauge only phases primitive Fourier coefficients and cannot change
nonvanishing.  This removes absolute root-label origin as a spectral
obstruction, but it does not identify a guard atom with a physical root or
supply a common base.

The exact connection contract is therefore:

| field | audited content |
|---|---|
| source | actual pre-merge `U_full` endpoint atom tables |
| target | `F_2^2 x F_13` chamber-pair/drift table |
| map | `s=7a+r`, owner boundary, `d=b-a`, and optionally one common torsor gauge `u=a+c` |
| preserved | frozen endpoint bank, both role values, bridge, guard characters, bucket support, Walsh and drift spectra |
| destroyed | physical root node, common base, owner/word/source sheets, horizon, chronology, lawful left/right support relation |
| required sidecar | a character-independent THM-2471 support predicate before either endpoint sum, carrying one common root gauge |
| cheapest test | restrict the 676 active Cartesian atom pairs by that predicate and require recovery of the frozen inverse value before testing graph factors |

The edge differences remain in the exact coboundary image `B^1`.  Full bucket
support and rank do not create a nonzero absolute `H^1` class.

## 6. Candidate packaging issue and repair

The mathematical candidate survived, but commit `e92508714` was not a fully
immutable package:

- its source still had `EXPECTED_SEMANTIC_SHA256="TO_BE_PINNED"`; and
- its committed output omitted the `guard_equality_hostiles` line printed by
  that source.

Thus that exact source/output pair was not line-for-line reproducible even
though its stored numerical claims were correct.  Incoming commit
`4e6021ca3` repaired the pin, output, equality hostiles, root-label gauge, and
scope prose.  The independent audit does not depend on that repair and has its
own pinned package.

## 7. Reproduction and hashes

```text
python -B 04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py
python -B -O 04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py
```

Normal and optimized runs reproduce the stored transcript exactly.  The
script/output LF SHA-256 hashes are

```text
f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc
7684cdb6bb1641780977d9e2def3753d802bf026e126dbc5351bd8f8ddebd906
```

and the semantic ledger is

```text
d52c9f0a56c14a83e1e6b175c7b725314c99f09d44509bc8582847a5857f7da6.
```

The next stronger endpoint-only question is the announced all-role extension:
whether all five refined role classes retain 52 support and whether every one
of the 72 weighted-tree charts has complete chamber/drift spectral support.
Even a positive answer would still live in `B^1`; it would not supply the
missing ancestry predicate.
