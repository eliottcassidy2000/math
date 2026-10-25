# Dyadic re-entry and the failure of divisor monotonicity

25 October 2026. A structural explanation of the frozen quotient-ring classification, not a novelty or priority claim. Zero is allowed as a square. “Parker” means that no 3×3 magic array of nine distinct square residues exists over Z/nZ. Modulus 1 is included as trivially Parker.

## 1. Exact dyadic profile

Let

    O = {1,3,5,7,9,11,13,15,17,19,21,23,25,31,33,
         39,43,47,49,51,57,67,69,75,93,129,141,147,201},
    E = {25,49,75,147}.

For n=2^a m with m odd, n is Parker exactly in the following cases:

- a=0 or 2 and m∈O;
- a=1,3,4 and m∈O\E;
- a=5,6 and m∈{1,3}.

There are no cases with a≥7. The seven exponent-layer counts are

    29,25,29,25,25,2,2,

giving 137 moduli including 1, or 136 for n≥2. The maximum is 16·201=3216.

The odd set O and the elimination of other odd parts use the frozen CLASSIFICATION.md and its explicitly stated prime-field theorem dependency. The new explanation here is why the dyadic layers repeat and then collapse. The finite binary compatibility calculation is independently re-enumerated in verify.py; no source checker is imported.

## 2. The exact small-dyadic array identities

Fix the row-major cell labels 0,…,8. For an array T, C(T) is the set of unordered pairs of labels carrying equal entries. It is a relation on fixed labels, not a histogram. Let M_q denote ALL square-residue magic arrays modulo q, including repetitions, represented by entries in [0,q−1]. Let B=M_2. There are eight binary arrays; their collision relations consist of the full relation and three minimal relations.

**Lemma.** The following are identities of actual arrays, not merely of collision antichains:

    M_4  = {constant 0, constant 1},
    M_8  = {4T : T∈B} ∪ {constant 1},
    M_16 = {1+8T : T∈B} ∪ {constant 0, constant 4}.

Proof. The square residues modulo 4 are 0,1. The integer line sums of an array using only 0,1 lie in [0,3]. Equality modulo 4 therefore forces integer equality. For any integer 3×3 magic array, each opposite pair sums to twice the center. If the center is 0 every opposite entry is 0; if it is 1 every opposite entry is 1. Hence the array is constant.

Modulo 8 the square residues are 0,1,4. Reduction modulo 4 is constant. A reduction of 1 forces every entry to be 1. A reduction of 0 forces entries in {0,4}; divide all entries and line-sum differences by 4 to obtain exactly a binary magic array modulo 2. Conversely every such scaled binary array works.

Modulo 16 the square residues are 0,1,4,9. Again the reduction modulo 4 is constant. A reduction of 0 gives entries in {0,4}; dividing by 4 gives a magic array with entries 0,1 modulo 4, which is constant by the preceding argument. A reduction of 1 gives entries 1+8t with t∈{0,1}; subtracting 1 from all nine entries changes every line sum by the same 3, and dividing line-sum differences by 8 gives precisely the binary magic equations. Conversely all these arrays work. ∎

Consequently the FULL realizable collision-relation sets for q=2,8,16 agree. Their three minimal masks, using lexicographically ordered pairs of cell labels, are

    0x495742ea2, 0x2e1cebb88, 0x98b456ad5.

The only mask for q=4 is 0xfffffffff. The independent enumeration gives:

    q       arrays    different masks    minimal masks
    2          8             4                 3
    4          2             1                 1
    8          9             4                 3
    16        10             4                 3
    32        73            10                 6
    64       522            20                 8

For coprime factors CRT intersects collision relations on the same labels. Thus for EVERY odd m, not just m∈O,

    m is Parker  ⇔  4m is Parker,
    2m is Parker ⇔ 8m is Parker ⇔ 16m is Parker.

This is the mechanism of re-entry: passing from modulus 2 to modulus 4 removes all nonconstant binary local arrays. Passing to modulus 8 restores precisely their collision patterns via the values 0,4; modulus 16 realizes the same patterns via the values 1,9.

## 3. The small odd compatibility check

For each m∈O, compute the minimal odd collision antichain and intersect it with the three binary masks above. A zero intersection occurs exactly for m∈E. The certificate includes every intersection, all local masks, realizing arrays, square-residue sets, and local array counts. The following compact independent diagnostic gives, for each odd part, the minimum NUMBER of surviving collision pairs across all such intersections:

    m : minimum surviving pairs
    1:16, 3:16, 5:3, 7:3, 9:3, 11:3, 13:3, 15:3,
    17:1, 19:3, 21:3, 23:1, 25:0, 31:3, 33:3, 39:3,
    43:1, 47:1, 49:0, 51:1, 57:3, 67:1, 69:1, 75:0,
    93:3, 129:1, 141:1, 147:0, 201:1.

This diagnostic is derived from the exact masks, not used instead of them. Modulo 3 every square magic array is constant, so multiplication by a coprime factor 3 cannot change the collision intersection. This explains the paired exceptions 25,75 and 49,147.

## 4. Why the layers collapse at 32

For every a≥5, the three residues 1,9,17 are distinct squares modulo 2^a in arithmetic progression. They are all 1 modulo 8, and each such unit has a square root modulo every power of 2: an odd root r modulo 2^k can have the next square bit toggled by replacing r with r+2^(k−1), for k≥3. This gives induction from modulus 8.

For every odd m≥5, the integer squares 1,25,49 remain a distinct arithmetic progression modulo m: their differences are 24 and 48, neither divisible by an odd m≥5. A progression x,y,z with x+z=2y produces the magic array

    y z x
    x y z
    z x y.

Its three collision blocks are {0,4,8}, {2,3,7}, {1,5,6}. Reflect a second such array vertically; its blocks are {2,4,6}, {0,5,7}, {1,3,8}. Every block of the first meets every block of the second in exactly one cell. Thus CRT combining the dyadic first array with the odd reflected array has nine distinct square entries. This is legitimate CRT between the COPRIME factors 2^a and m; it is not a lift between 2^a and 2^(a+1).

This excludes all m≥5 for a≥5. Complete small tables show that 32 and 64 themselves are Parker; the modulo-3 constant-array fact includes 96 and 192. Finally for a≥7, take 1+8L where L is the ordinary Lo Shu array with entries 1,…,9. All entries are square units, the array is magic, and its differences 8d, 1≤|d|≤8, are nonzero modulo 2^a. Hence these dyadic rings are non-Parker. CRT with constant zero modulo an odd m then excludes every remaining case. If m is already non-Parker, the same argument in the opposite factor excludes 2^a m for every a. Together with Sections 2–3 this proves the profile.

## 5. Explicit loss of a witness: 50 → 100 → 200

A non-Parker witness modulo 50 is

    41   6  31
    16  26  36
    21  46  11.

Every line sums to 28 modulo 50. Corresponding square roots modulo 50 are

    29  34  41
    46  26  44
    39  36  31.

The entries are distinct. The array reduces modulo 2 to

    1 0 1
    0 0 0
    1 0 1,

and modulo 25 to

    16  6  6
    16  1 11
    21 21 11.

The two masks are 0x495742ea2 and 0x300000104; their intersection is zero. Modulo 25 alone the repeated entries are unavoidable for distinctness, but the binary labels split every repeated pair.

There is NO lift of this array to a magic array of square residues modulo 100. Such a lift would reduce modulo 4 to a constant array by Section 2, so its parity would be constant. Reduction modulo 50 preserves parity, contradicting the displayed nonconstant binary array. More strongly, no alternative nine-distinct array exists modulo 100: its mod-4 mask is full and its mod-25 component is Parker. Nevertheless 200 is non-Parker because modulo 8 the binary mask returns.

The same phenomenon occurs for 98→196→392, 150→300→600 and 294→588→1176. Every non-Parker divisor of a Parker modulus occurs in exactly the following six ordered pairs (Parker multiple, non-Parker proper divisor):

    (100,50), (196,98), (300,50), (300,150), (588,98), (588,294).

The finite complete Parker list makes this exhaustive: enumerate every proper divisor of every listed Parker modulus and test membership in the same complete list. Thus 50→100 is the smallest loss, both in its starting modulus and its Parker target. This disproves general upward preservation of non-Parkerness under divisibility, and equivalently downward preservation of Parkerness. It does NOT contradict persistence under adjoining a coprime factor by CRT. Nor does it assert arbitrary liftability or imply an integer magic square of distinct squares.

## 6. Reproducibility and dependencies

Run `python verify.py --verify` from this directory (or pass its path from the workspace root). The checker enumerates every small prime-power local table independently via top-row/center line equations, checks the exact array identities, recomputes all binary compatibility intersections and 203 candidate flags, validates the explicit witness and roots, and enumerates all six divisor counterpairs. It agrees with the frozen candidate certificate without modifying or importing its checker. `python -O verify.py --verify` also passes; validation uses explicit exceptions, not assertions.

SOURCE_PINS.json records SHA-256 hashes of exact dependencies. The unbounded classification depends on the frozen reduction proof and separately audited prime-field theorem, not merely on the finite test. README.md describes the file layout. No remotes were used and no literature-priority assertion is made.
