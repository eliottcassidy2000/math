# Independent audit of the common-access two-cusp braid proof

**Status: independent analytic/source audit PASS; complete normal, optimized,
and frozen-output replay PASS. Audit frozen.** The exact computation is
FINITE-EXACT; the argument that its actual paths force the complement group
and the stated whole-support exclusion is proved separately below.

This audit concerns the exact four-loop argument in
[the primary note](planar_jc48_sep08_two_cusp_braid.md), for the literal curve

    U=t^4-2t^2,  V=t^6-3t^4/2+t^3/3-t.

Its audited conclusion is the actual affine complement group Z and hence
exclusion as a whole irreducible Keller nonproperness curve. It is not a
classification of other two-cusp curves. The independently audited
[geometry supplier](planar_jc48_sep08_two_cusps_audit.md) is retained as a
dependency; its former OPEN global-group obligation is what these new
common access paths address.

## 1. Literal polynomial, complete paths, and frozen universe

I read the entire default verifier, its two separate heuristic/production
modes, all of the primary proof, and the inherited marked-braid convention.
The stored coefficients are exactly the actual resultant of U-u,V-v;
the default verifier reconstructs that resultant and F(U,V)=0 over the
rationals. Its monic degree four and minimality are also supplied by the
geometry proof. The discriminant contains all six critical projection
values, retaining the two cusps together over u=-1.

The accepted universe is precisely four closed paths based at 4+3i:
centres 0,-1,3,-401031/2^20, each using the common radius 1/32 and the
six ordered edges stated in the primary. I independently decompressed all
four witnesses, checked their exact names, equal u/root-row counts, exactly
four two-coordinate root centres on every row, the common labelled base
fibre, all byte lengths, and all raw/compressed hashes. No extra roots or
hidden fifth coordinate occurs in the data.

A separate literal Fraction check, without importing the producer,
checked every one of the 77,174 base segments by a cross-product
collinearity test and strictly increasing dot-product edge parameter.
It recovers the six declared directed edges in order, ending at the row
indices

    smooth:    4648,4712,4776,4840,4904,9553;
    cusps:     5983,7869,9703,11536,13423,19412;
    node3:     10345,14603,18855,23107,27364,37713;
    node_left: 4771,5008,5247,5486,5723,10496.

These indices include the complete outbound and return stems. In particular the
witnesses are not just arbitrary regular loops assigned the desired names.
The two complex-node scouts are excluded from the proof universe. No
statement that each rational diamond encloses exactly one critical value
is necessary: the actual braid of any certified regular loop supplies a
necessary relation.

## 2. Exact disk certificate and endpoint homotopy

For a rational complex number z, the verifier's lower bound
max(|Re z|,|Im z|) is at most its Euclidean norm, and its upper bound
|Re z|+|Im z| is at least that norm. These inequalities have the directions
required by every disk and Taylor estimate.

The radius attached to a centre is one thirty-second of the smallest
pairwise lower distance. For any two centres the sum of their radii is
at most one sixteenth of that pair's lower distance, hence strictly less
than their Euclidean separation. Their closed disks are disjoint.
The exact Taylor expansion retains every term through degree six in u
and four in v, the full bidegree of the actual polynomial. The strict
Rouche inequality compares the polynomial to B01 w on a fixed circle.
All omitted terms have their full upper bound on the right, including
base displacement terms. Thus it gives exactly one root, counted with
multiplicity, in each of the four disks throughout the entire straight
base segment. All roots are simple, and degree four plus monicity proves
that there are no additional roots. This also certifies that every base
segment stays in the regular base.

At a joint, let R be its incoming disk radius, R' the radius around the
new endpoint centre, and epsilon=min(R,R')/64. The centre moves by
upper distance less than R/2. Therefore the separately isolated endpoint
root in the epsilon disk lies in both full adjacent disks, since
R/2+epsilon<R and epsilon<R'. Its identity agrees with the unique root
in each disk. This is the exact label-continuation certificate; numerical
nearest-neighbor choices made during production have no logical role.

On each segment the straight centre strand lies in its incoming convex
disk, as does the actual root strand. Straight interpolation gives a
configuration-space homotopy through four distinct points. The two
homotopies at a common vertex agree because both use the same endpoint
centre and the just-certified same actual root. Initial roots are
independently isolated in disks of radius R/64. All four paths use the
same ordered initial centres; their final unordered centre sets are
exactly the same initial set. Uniqueness in the small base disks makes
the endpoint identification agree with the initial one, even when the
individual path permutes the four labels. Thus these homotopies give one
common identification of the actual base fibre with the rational base
configuration. There is no uncontrolled path-dependent conjugation.

The generic projection is multiplication by 1+i/4, an orientation-preserving
complex linear map. All endpoint projected differences are nonzero. For
each pair a real crossing occurs precisely when its endpoint differences
have opposite signs; its exact time is the rational number a/(a-b).
The source rejects simultaneous crossings, checks adjacency in the
current order, and checks nonzero imaginary separation at the crossing.
These exhaust crossings of the polygonal strands. Adjacent inverse
cancellation is a valid reduction and no other unsupported braid
simplification is used in the numerical-to-exact interface.

Default replay uses rational arithmetic for the witness tests. NumPy
root proposals are imported only when an explicit scout or producer mode
is called. The word is reconstructed from the checked strands and is
then compared to both the stored and declared words; it is not accepted
from certificate metadata alone.

## 3. Marked direction, section, and omitted relations

The crossing sign agrees with a positive counterclockwise half twist:
an initially left point moves below the right point. With below-stem
meridians, its geometric pushforward on the two free generators is the
inverse Artin automorphism

    (a,b) -> (b,b^(-1)ab).

Contravariant transport of monodromy values uses its inverse, so the
chronological tuple rule is

    H+(a,b)=(aba^(-1),a),
    H-(a,b)=(b,b^(-1)ab).

This is the same coherent convention as the previously audited braid
suppliers. It fixes the relation between the literal crossing direction,
the chronological word, and the evaluated group tuple. It does not appeal
to an unselected sign/reversal convention or merely compare local cycle
conjugacy classes.

A continuous section above a Cauchy root bound exists over the whole
u-plane. On any compact filling disk one may use a single sufficiently
large vertical basepoint. Consequently the lifted base loop is null in
the full complement. Over the regular base, fibre loops and section lifts
generate the bundle group. Restoring exceptional fibres kills the latter;
loops in the full complement can be perturbed to avoid the finitely many
exceptional vertical fibres. The positive fibre meridians therefore
surject onto the full affine group, and the checked braid loops impose
exact fixedness on their tuple.

Only necessary relations are used. The four selected base loops need not
generate the punctured base group. Omitted node loops or other relations
can only further quotient the group already forced by the selected
relations; they cannot invalidate the surjection or weaken a proved
cyclic upper bound. The final winding argument below prevents such an
extra relation from turning the actual group into a finite cyclic group.

## 4. The arbitrary-group calculation gives an infinite cyclic group

I independently checked the following exact reduction in an arbitrary
group, without a permutation census. Let the four meridians be (a,b,c,d).

- The smooth word [3,1,-3] equals [1], since the two disjoint generators
  commute. Fixedness gives a=b.
- The node3 word [3,3] gives [c,d]=1.
- The node_left word [3,2,2,-3] is a square after the prefix [3]. Its
  exposed pair is (b,cdc^(-1)), so it gives [b,cdc^(-1)]=1. The preceding
  relation reduces this to [b,d]=1.
- The cusp word has prefix [3,2], middle [3,1,3,1,1,3], and the inverse
  prefix suffix. The middle uses only commuting generators 1 and 3 and
  is exactly [1,1,1,3,3,3]. The prefix sends the tuple to
  (a,bcdc^(-1)b^(-1),b,c), which becomes (a,d,b,c) **after** the two
  commutation relations. Fixedness of its two disjoint cubes gives
  ada=dad and bcb=cbc.

Now a=b and [b,d]=1 make a,d a commuting braid pair, hence a=d. Since
[c,d]=1, b,c commute as well; their braid relation forces b=c. Thus all
four original positive meridians agree. The source independently checks
the free-group prefix identities and the cusp-word rearrangement.
Every constant tuple fixes all four words, so no hidden finite-order
condition has been added in this algebraic reduction.

The surjection from the four fibre meridians makes the actual group
cyclic. The reduced defining equation F induces a map to C*. Around a
smooth curve point F is a unit times a local normal coordinate; a positive
meridian has winding number one. Therefore that cyclic generator has
infinite order, and the affine complement group is exactly Z.
This lower bound is actual and global, rather than a conclusion of the
selected presentation by itself.

## 5. The whole-support Keller consumer is correctly typed

Assume this C is the whole irreducible nonproperness curve of a
nonautomorphic polynomial Keller map Phi. Off C the map is a finite
etale cover. Its source is the complement of an algebraic curve in A2,
hence connected, so the monodromy is transitive. The primary now also pays d>1 directly: a degree-one etale map is a
birational open immersion by Zariski Main. An omitted divisor would pull
back to a nonconstant unit, impossible on A2; if the omitted set is finite,
normality extends both inverse coordinates and gives a polynomial inverse.
I checked this added argument, including the use of dominance to preserve
nonconstancy of the pulled-back equation.

There is an actual retained sheet over a general smooth point of C:
F composed with Phi is nonconstant by dominance. Its zero set has a
divisor, and quasi-finiteness prevents that divisor from mapping to a
point. Its image is therefore dense in the irreducible C. Etaleness
gives a genuine local inverse and a fixed label of a transverse positive
meridian. This is an actual-sheet statement, stronger than merely having
an abstract inertia-fixed label. Equivalently the earlier finite-envelope
supplier gives 1<=a<=d-1; its positive deleted defect also rules out d=1.

A transitive cyclic permutation action on d labels is one d-cycle. For
d>1 its generator has no fixed label, contradicting that retained local
sheet. This proves exclusion for this literal **whole support**; the
complete finite replay below has now been independently accepted. It does not assert a global Keller source
from an abstract passport or transfer the one-cusp Euler ledger. The
three-generator/local-A6 control in the companion passport note is used
only as a hostile to forgetting the actual common access relations; none
of its provisional classification is a dependency of this exclusion.

## 6. Frozen pins and completed independent replay

The source, read in full, has 10,866 bytes and SHA256
`d5f0d1ce16b439c13d936f35ce38f558b1c5210f892a11218ab7f4f4695198e1`.
Each raw and compressed witness was independently hashed:

| Path | Segments | Gzip bytes | Raw bytes | Gzip SHA256 | Raw SHA256 |
|---|---:|---:|---:|---|---|
| smooth | 9553 | 651158 | 2522124 | `3ae1dfc0ada5d2fe148cb4da1b4643fcfd10280e66ac8a73ad783b6a007b22b9` | `06d050671046d9145ce23c82fa7fbbbeaa2d8a73ba62fd664b74cdc8046820be` |
| cusps | 19412 | 1253330 | 5051033 | `fa2c749fd1cd3cf065973ad8fe5bfc8c0aa725055564f3d25f404c3d8405e018` | `e77a3b3d313d8a5f9d0eb87ba241318fef531004e13255d0c3245a28180caa76` |
| node3 | 37713 | 2395456 | 10053504 | `f86f82e2bbcd6f701a50c0a70961ef9b2e1828b3f840e3b5a9d266289f2068bc` | `7638249409e0db9cd66c58f802e54e5285b607883a102057474e76ea1e11d49a` |
| node_left | 10496 | 736165 | 2861579 | `2994c91af1c1de6d4f3aa8555f2899dfa8c24465070314352ec00d38421cf662` | `794de455b150396b31c6e96154d4403a0d7dbfb5e2c28d69f2f107e532d6e45a` |

Independent reproduction commands:

    python3 -B 04-computation/planar_jc48_sep08_two_cusp_braid.py
    python3 -B -O 04-computation/planar_jc48_sep08_two_cusp_braid.py

Both independent runs finished with **694,669 always-active gates**,
covering all **77,174 segments**. Their complete outputs are byte-identical
to each other and to the frozen primary output: **1,925 bytes**, SHA256
`40e8efddb55d83269a51d190c7dbd0d85c61355cca0a76123892a52b0e48cd68`.
The per-path gate counts are 85,994, 174,740, 339,433, and 94,484; the
remaining 18 gates verify the literal polynomial and group identities.
No source or witness correction was found in the complete proof/source
audit above. The final primary prose, including its self-contained
birational-Keller argument, has also been accepted. The scope is the
stated actual curve and its whole-support exclusion, not the broader
family. No prior frozen artifact was edited for this audit.
