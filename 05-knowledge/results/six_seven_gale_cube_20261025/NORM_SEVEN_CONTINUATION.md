# Two norm-seven ideals select two complementary Campedelli pairs

**Status: PROVED composition of inherited maps, with a FINITE-EXACT marking
audit over the fixed branch labels.** This is a narrow continuation of
[Session61's torus/Fano construction](../session61_hadwiger_six_seven_20261025/torus-fano/TORUS_FANO_CAMPEDELLI.md)
and [our normalized fiber-product theorem, section 6](COMBINATORIAL_BRIDGES.md).
The lattice/Fano construction and its four free codes are inherited;
the additional deduction identifies its two pairs inside the existing
32-pair census and applies the proved reconstruction theorem.

## The same labels throughout

Use the branch order

    (E, E+x, E-x, E+z, E-z, E+x+z, E-x-z)

and the corresponding dual coefficient vectors

    (0,0), (1,0), (-1,0), (0,1), (0,-1), (1,1), (-1,-1).

The actual concurrent triples are 012,034,056. With omega^2+omega+1=0,
the conjugate ideals (2-omega) and (3+omega) both have norm seven.
Their residue maps are c_k(a,b)=a+kb mod7 for k=2,4. On this fixed branch
order the residue lists are respectively

    c_2: (0,1,6,2,5,3,4),
    c_4: (0,1,6,4,3,5,2).

For each k, the two face systems on residues are

    U_(k,t)={t,t+1,t+1+k},
    D_(k,t)={t,t+k,t+1+k},       t in F7.

Pull both systems back through c_k to the seven branch indices. For k=4 the actual face formulas use 1+k=5. As collections of residue
triples, U4=D2 and D4=U2 (translation reindexes t). Reusing the k=2
formulas would therefore swap the U/D face marking while retaining the
same unordered pair; we use the actual formulas to preserve that marking.

The standalone audit recovers the following binary markings. These are
representatives up to GL3(F2), not additive maps from F7 to F2^3. In each
row the Fano lines are exactly the triples whose binary labels xor to zero.

| Ideal residue root | Face type | Binary labels in the fixed branch order |
|---:|---|---|
| 2 | U | (1,2,4,5,6,3,7) |
| 2 | D | (1,2,4,6,3,7,5) |
| 4 | U | (1,2,4,5,7,3,6) |
| 4 | D | (1,2,4,7,3,6,5) |

Each row defines a map G=F2^7/<1> -> F2^3 by sending the i-th branch
inertia generator to its binary label. Denote its order-eight kernel
by H_(k,U) or H_(k,D). The receipt checks both the syndrome kernel and
the independent decoder whose nonzero even words are complements of
the Fano triples. It also directly checks the single-, pair- and actual
triple-inertia conditions for freeness.

## The exact additional conclusion

For each k, the U and D Fano systems are disjoint. Nonzero elements of
their kernel groups correspond bijectively to Fano triples modulo
complementation. Therefore

    H_(k,U) intersect H_(k,D) = {0},
    G = H_(k,U) direct-sum H_(k,D).

All four systems avoid the three actual branch triples, so all four kernels
act freely on the hourglass resolution. They are four distinct embedded
subgroups on the same branch labels. The two unordered pairs are distinct
members of the **32 complementary free pairs** from our earlier census.
They select exactly two members of that census, not all thirty-two.
Session61 identifies all four individual codes with its orbit of size four;
no new claim of distinct abstract Campedelli isomorphism classes is made.

Let X be the normal seven-root cover of P2 and let
Z_(k,U)=X/H_(k,U), Z_(k,D)=X/H_(k,D). These are the finite normal quotient
models whose resolutions are Campedelli surfaces. Applying the prior
fiber-product theorem to each complementary pair gives

    X = normalization of (Z_(2,U) x_(P2) Z_(2,D)),
    X = normalization of (Z_(4,U) x_(P2) Z_(4,D)).

Indeed, the two invariant fields intersect in C(P2) and their composite
is C(X). Their finite flat tensor algebra is integral, and X is its finite
normalization. Normalization remains essential at the branch locus:
after separating unramified local branches, the raw product has u^2=v^2.
This statement concerns the finite normal cover models; it does not replace
them by smooth resolutions inside a finite fiber product.

The retained data are the ideal, the face type, the branch marking and
the kernel. No tournament orientation or chromatic assertion is transported,
and no rational or integral point-lifting statement follows.

## Reproduction

From the repository root:

```sh
python3 -B 04-computation/six_seven_gale_cube_20261025/check_norm_seven_pairs.py > /tmp/norm_seven_pairs.normal.json
python3 -B -O 04-computation/six_seven_gale_cube_20261025/check_norm_seven_pairs.py > /tmp/norm_seven_pairs.optimized.json
cmp /tmp/norm_seven_pairs.normal.json /tmp/norm_seven_pairs.optimized.json
cmp /tmp/norm_seven_pairs.normal.json 05-knowledge/results/six_seven_gale_cube_20261025/norm_seven_pairs.json
```

The script prints full JSON without requiring an output argument; optional
`--output PATH` additionally saves the same receipt. It imports no producer
code and uses Python's standard library. Its universe is all 5,040 markings,
all 30 Fano systems, all 16 free kernels and their 120 unordered candidate
pairs, plus the two specified ideals and four pulled-back face systems.
The exact result is **219 explicit gates**, with byte-identical normal and
optimized receipts. Receipt SHA-256:
`4fa307e77b655f640b74f23848f5827f3aed0a7843bc90d5b6acf65aa7565b5e`.
