# Independent audit of physical boundary-return collapse

**PASS in the stated scope.** Every nonempty triangle-ray return maps the full raw vector space into the balanced plane and maps that plane into the triangle ray. Two return blocks therefore collapse the complete linear map to rank one and freeze the normalized auxiliary readout. The resulting prime obstruction applies to the repeated-return scheme, not to arbitrary tournaments or the full synchronized inverse. The one-return prime139 control is valid and makes the threshold sharp.

## Frozen input hashes

- `RETURN-COLLAPSE.md`: `d224a40d84e8b772a5ec9d539798627a199ac78b274a54b1557c482f4cb4733a`.
- `verify_return_collapse.py`: `cf1f58429c578e98e662717ddba25109f6f25905de8e1f397634de4e0b7e8cee`.
- Both author receipts: `700564074549bc981957681d86e3ea82e32f0bd768c4d85afc2c13732194fd9a`.

No correction to the mathematical argument was required. Any later status-only promotion needs a recorded final diff or addendum; these are the exact reviewed proof bytes.

## All-order argument checked

1. **Global suffix reduction is valid.** For a physical return beginning at b, a final O cannot be balanced. If the word ends in a nonempty Z-run, its last preceding non-Z letter cannot be O, since O creates strictly positive f and Z preserves f. Hence it is J. The identity JZ=OJ is an identity of full four-coordinate linear maps, so the terminal JZ-run can be replaced by an O-run followed by J. If the word already ends in J, balance requires predecessor E=0; from the physical seed b only a final O can make E zero. A pure Z-run never returns to a scalar multiple of b. Thus every nonempty return is globally equivalent to P OJ.

2. **The positivity step is sound.** OJ maps every raw vector to `(D+2U,D,D+2U,D)`, so its image is the balanced plane L. If `b.P.OJ=lambda b`, the actual predecessor b.P has D=U=lambda. In nonnegative source coordinates, b has a=e=1 and f=p=0. Since P is a product of nonnegative matrices, vanishing f and p outputs at b forces both the a and e coefficients of those rows to vanish individually. These rows therefore vanish on the entire plane f=p=0, including formal signed vectors. It follows that P(L) lies in L. OJ then maps `(x,y,x,y)` to x b. No unjustified extension from one point to a plane is made.

3. **Rank collapse is exact.** The first return maps all of R^4 into L, and the second maps L into span(b). Their composite has rank at most one and cannot be zero, since it sends b to the positive multiple lambda1*lambda2*b. Thus its rank is exactly one. The auxiliary coefficient mu is an integer because the E-coordinate of b is1; nonnegativity follows from the extended positive cone. The physical synchronized joint cone makes the ratio strictly between1/2 and1. Later return choices multiply both boundaries by the same scalar and do not alter this ratio.

4. **Prime exclusion is properly scoped.** For at least two returns, the target factors as `G(v) * product(later multipliers) * [lambda1*lambda2 + L_k*mu]`. The first factor is at least5, and the bracket at least9. Every factor is a positive integer. This excludes primes from this construction alone, including its balanced readout. It does not exclude primes from C1, general joint pairs, or one-return constructions. The finite-menu statement concerns normalized sidecars and is consistent with differing later scale factors.

5. **The one-return positive control is genuine.** C=ZOJ sends b to5b but sends z to `(11,3,11,3)`, which is not on the triangle ray. With v=t this is `(m,h)=(25,19)`, and k=2 gives139. The adjoint reconstruction is the source word JZOOZZ, followed by the closing J. The resulting nine-vertex tournament has exactly139 Hamilton paths and one Hamilton cycle. Thus the two-return threshold cannot be replaced by one return.

## Replays and independent checks

The author checker was replayed under normal Python and Python -O. Both outputs match each other and both frozen receipts byte-for-byte. It exhausts all88,572 nonempty O/Z/J words through length10 and finds5,476 returns. Its short return menu has26 words;676 pair and17,576 triple checks pass.

The separate `independent.py` imports no author code. It constructs exact integer generator matrices and performs rational row reduction. For every return through length10 it verifies the global normalized OJ factorization, preservation of L by the prefix P, whole-space image in L, and plane image in span(b). In this finite universe5,458 return matrices have rank1 and18 have rank2. Every tested two-return composite has rank exactly1. The arbitrary-return triple checks retain only the predicted common scalar.

The A/B collision family has exactly three auxiliary states at every checked depth2 through12, agreeing with the all-order scalar law. For the139 control, the independent implementation builds a dense tournament adjacency matrix and enumerates every vertex permutation, rather than using subset dynamic programming. It independently obtains139 Hamilton paths and one cycle. C has rank2, providing an additional sharpness check.

The independent normal and optimized receipts are byte-identical. From this directory:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

Finite checks validate the matrix implementation, rank computations and graph example. The all-order conclusions rest on the explicit factorization and positivity proof above. Universal target coverage remains OPEN.
