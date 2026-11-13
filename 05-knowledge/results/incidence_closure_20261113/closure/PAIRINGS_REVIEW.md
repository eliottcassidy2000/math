# Independent review of signed converse transport

**Verdict: PASS after the recorded random-walk basis clarification.** The all-order claims are the precisely scoped finite-matrix statements S1–S4. The orders3–7 conclusions remain FINITE-EXACT, conditional on the complete matrices and their stated provenance. This review does not promote the candidate or assert a preferred presentation action.

Reviewed sources:

- `pairings/SIGNED_LIFT.md`, corrected SHA-256 `1710c4730bb2512e8574984ddb4a5a398fe46c699824571bd309b40210407ac7`.
- `pairings/check.py`, SHA-256 `e4df424debed4c4d4f2903405fbb7e8f90c63083d3b402ff0eef855abc300cf8`.
- `pairings/normal.json` and `pairings/optimized.json`, each SHA-256 `afd9c477812121bd7df67b970a12a330313e04f093bc99df7c91dc0fcabe982f`.

## Correction lineage and source weights

The initial proof snapshot, SHA-256 `1c10b5771df3ddb7650f89d9aaaed6d863c7bbbc73dd0c978328c157ec37920c`, discussed the orthonormal plus basis for C and then called the merged stochastic P its analogous random-walk block without specifying a basis change. In the **natural fibre-constant basis**, the plus block is indeed P_Y. In the preceding orthonormal basis it is M^(1/2)P_YM^(-1/2). The author inserted this distinction; the normalized minus block remains diag(d(x_u))^(-1)D. No matrix, computation or finite claim changed.

The raw C plus block is M^(-1/2)AM^(-1/2) in the orthonormal sum basis. Equivalently, C acting on unnormalized fibre-constant coordinates has matrix M^(-1)A. Source ordinary degree is d(x)=H/Aut, while merged degree is m(u)d(x). The corrected normalization is consistent with all of these conventions.

## Mathematical audit

**S1: reconstruction.** Converse commutation and symmetry force a double–double block [[a,b],[b,a]], with A=2(a+b) and D=a−b. Solving gives a=A/4+D/2 and b=A/4−D/2. Singleton–double entries are A/2, and singleton–singleton entries are A. The chosen sign gauge changes D by diagonal sign conjugation. Integer/parity, nonnegativity and loop conventions are appropriately stated as reconstruction guards, not as sufficient conditions for realizing an arbitrary formal pair by tournaments.

**S2: kernel.** A fibrewise exchange multiplies signed row and column u by the same epsilon_u. Therefore a nonzero off-diagonal D(u,v) forces epsilon_u=epsilon_v, while diagonal entries impose no condition. Components, including isolated vertices, are exactly the independent choices. The C2^c kernel, its fixed polynomial x^s product(1+x^(2m_i)), and the unique minimizing global converse follow. This is conditional on typed merged vertex rigidity when promoted from a kernel to the entire group. The note keeps that hypothesis explicit.

**S3: liftable base group.** A candidate g first preserves fibre type and the required merged matrices/labels. It must also preserve zero positions, magnitudes and diagonal values of D. Subject to those necessary guards, the ratios on nonzero off-diagonal support edges are signs. Their products on support cycles equal1 iff they are a vertex-sign coboundary. Spanning-tree propagation proves sufficiency and counts exactly2^c solutions. The note correctly allows the liftable base group to be smaller and does not claim its group extension necessarily splits. The ratios, rather than the integer D entries, are the invertible transports.

**S4: representation and cycle index.** The plus/minus spaces have dimensions s+t and t. Their displayed bases are independent and together span the ordinary-class space, so the verified intertwining identities imply the determinant factorization. The kernel fixes the plus space and acts by one sign character per signed-support component on the minus space. The cycle-index factors correctly distinguish2m_i fixed elements from m_i transpositions. Burnside substitution counts colorings of **ordinary classes**, not tournaments on a labelled vertex set. For nontrivial base cycles the sign product determines two cycles of length l or one cycle of length2l; summing only valid sign solutions and normalizing by the actual group size is correct. The extra involution condition epsilon_u epsilon_(g(u))=1 is needed and is present.

These are commuting sign characters, not Clifford generators or a periodicity theorem. The text claims only a class-matrix action for arbitrary compatible lifts. It retains the separate inherited physical reflection lift for global converse and explicitly declines a preferred fixed-spine lift for arbitrary matrix symmetries.

## Computational audit

The checker's source was reviewed in full. The complete order3–6 orbit construction uses all n! vertex permutations, independently of the closure lane's adjacent-transposition orbit compiler. Every physical Hamilton path is tested; native complement is c XOR full_arc_mask XOR path_mask. Dividing by source Aut is justified by the free action on full path orders. All four ordinary/merged and total/blue matrices agree with the closure lane's separate complete construction.

The automorphism refinement preserves all graph automorphisms: it records original labels, weighted degree, loop weight and a multiset of weighted neighbors in each current color. This is stronger than row-sum equitability and must not be conflated with the closure lane's coarsest lumpable partition. Exhaustive backtracking within its surviving cells then verifies all matrix entries at completion. Discrete refinement therefore certifies the stated vertex rigidity.

The signed reconstruction check visits every ordinary matrix entry. Support components include isolates. The code checks every component-flip generator as an ordinary C symmetry, and each nonzero signed edge as a one-fibre-flip hostile. The lift solver validates diagonal and absolute values before propagating signs, and each returned solution is checked directly. Its synthetic four-double-fibre control compares all sign solutions and then independently enumerates all8! ordinary permutations; both the untwisted and changed representative gauges have24 base candidates,4 liftable base permutations and8 ordinary symmetries.

At7, the ordinary matrices remain the shared complete closure census, with two internal construction paths. The pairing checker independently recounts six full low-H rows by literal Hamilton permutations and unrestricted canonicalization. These six rows are correctly described as spot checks rather than a second complete order-seven census.

## Additional independent consequences replay

[check_signed_review.py](check_signed_review.py) imports no producer module. From the frozen matrix receipt it checks:

1. C intertwines both fibre-constant and fibre-difference bases with the stated operators, and the stochastic blocks have the correct source degrees.
2. Converse commutation and confinement of blue incidence to singleton fibres hold entrywise.
3. A separate weighted-vertex refinement reproduces the rigidity certificates.
4. Signed components computed by repeated set union agree with the producer BFS; every component action preserves C.
5. Literal cycle decompositions of all small resulting kernel elements give the fixed polynomials and the Burnside numerators at2 and3 colors.

These are independent consequence checks on the same frozen matrices, not another order-seven tournament census. Normal and `-O` runs pass **1260321 always-active checks** and are byte-identical. Source SHA-256 is `ae8247cb233bc15d600c6881b25e4d7e8344bd5983bd5961e51c04004d6e6734`; both receipt files have SHA-256 `d55f89bc9677a650f31a9a0519840662aa1924e6eff9fa19c697bf8b66b75f56`. The receipts pin both input and source hashes.

Reproduce:

```bash
python3 05-knowledge/results/incidence_closure_20261113/closure/check_signed_review.py
python3 -O 05-knowledge/results/incidence_closure_20261113/closure/check_signed_review.py
```

The resulting kernel ranks are0,1,2,1,1 at n=3..7; the connected components have sizes22 and184 at6,7. Combined with the audited vertex-rigidity certificates, this supports the complete ordinary-matrix groups and polynomials asserted in S5. No all-order connectivity or rigidity claim follows.
