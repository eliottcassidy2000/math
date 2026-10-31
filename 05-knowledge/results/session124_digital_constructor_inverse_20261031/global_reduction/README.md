# Exact reduction of the global odd spectrum to targets not divisible by three

**Proved logical reduction, not a global attainment proof.** This addendum
connects the digital non-3 inverse to the unrestricted tournament conjecture.
It does not change any of the sealed digital packet's 44 files.

## Ingredients and scope

Write H(T) for the number of directed Hamiltonian paths of a finite tournament.
We use:

1. `H(A▷B)=H(A)H(B)`, where every cross-edge in the ordinal sum points from A
   to B;
2. a directed triangle C3, with H(C3)=3;
3. the explicit nine-vertex tournament with H=63 retained in `witness63.json`.

Items 1–3 are justified below. The already established **global exclusions of
7 and 21** are treated as external known results, not re-proved here. They are
needed to identify the exact forbidden set, rather than merely prove the
attainment implication. All coverage statements below concern positive odd
integers; this is not a new proof of the general parity theorem.

## Equivalence of the attainment obligations

The following two statements are equivalent:

**A.** Every positive odd integer N other than 7 and 21 is attained by a finite
unrestricted tournament.

**B.** Every positive odd integer m with `3∤m` and `m≠7` is attained by a finite
unrestricted tournament.

A implies B immediately. Conversely assume B. Write uniquely

```
N = 3^k m,   k>=0,   3∤m.
```

If `m≠7`, choose a tournament with count m using B, then take its ordinal sum
with k directed triangles. Multiplicativity gives count `m*3^k=N`. For k=0,
use the original tournament with no appended triangles.

If `m=7`, the cases k=0 and k=1 are precisely N=7 and N=21, the two exceptions.
Every remaining case has k>=2 and is obtained from the retained H=63 tournament
by appending k−2 directed triangles:

```
63 * 3^(k-2) = 7 * 3^k = N.
```

This proves A. Thus, given the known exclusions, solving the unrestricted
non-3 attainment problem B resolves the claimed global forbidden set `{7,21}`.
The equivalence of A and B itself uses the explicit 63 witness and product
identity; it does not assume the exclusions to derive new attainments.

## What would suffice from the restricted digital inverse

The digital inverse is complete as a **decision/reconstruction algorithm in
its stated finite-actual-seed singleton/triangle language**. Completeness does
not say that every eligible target receives a witness.

If that restricted inverse can be proved to succeed for **every** positive odd
m with `3∤m` and `m≠7`, its actual witnesses establish B, hence A by the proof
above. This is a sufficient route to the unrestricted conjecture.

Do not infer the converse: the unrestricted conjecture does not require its
witnesses to belong to this particular seed bank or grammar. This reduction
provides no theorem that unrestricted existence forces universal success of
the restricted inverse. The missing obligation is still a uniform
positive-answer proof, not another finite coverage computation.

## Elementary product proof

A directed path in A▷B cannot visit B and subsequently return to A because
there is no edge from B to A. A Hamiltonian path is therefore uniquely a full
Hamiltonian path of A followed by one of B. Every such pair concatenates,
since the connecting edge points from A to B. This is a bijection, proving
`H(A▷B)=H(A)H(B)` and its iterated version.

The triangle `0->1->2->0` has exactly the three cyclic rotations of `(0,1,2)`
as directed Hamiltonian paths, so its count is 3.

## Explicit H=63 constructor and independent verification

The [million-prefix constructor](../../session124_million_path_spectrum_20261031/construct_witness.py)
returns the following actual recipe for target 63:

```
seed (order, edge-code, apex-cut) = (1,0,0)
singleton blocks, chronologically = 1,0,0,1,1,0,0
final order = 9
```

The initial parent is a single vertex a with edge a->x. Append seven parent
sinks ordinally in the displayed order; block bit 1 means x beats the new
vertex and bit 0 means the new vertex beats x. Put x last in the adjacency
matrix. `witness63.json` retains that entire matrix and the constructor/input
source hashes. No strongness claim is made or needed.

`verify.py` independently reconstructs every edge, checks the tournament
axioms and exact agreement with the retained matrix, and counts full paths
in two different ways: backward subset endpoint DP and literal enumeration
of all 9! = 362,880 vertex permutations. Both give **63**. The same script
checks H(C3)=3. It does not use the digital formula or transfer coordinates.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 verify.py
PYTHONDONTWRITEBYTECODE=1 python3 -O verify.py
```

The normal and optimized receipts agree byte-for-byte. This is a fixed-witness
control for the general product proof, not a new coverage census. The parent
packet's existing whitelist/hashes are unchanged; this addendum has its own
`SHA256SUMS` and verification record.
