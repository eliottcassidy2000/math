# One ternary action, two quadratic shells: Pythagorean seeds and square-sum absorbers

**2026-10-19. Exact structural bridge, with no historical novelty claim for Berggren trees or almost-Pythagorean triples.** The same three classical integer matrices organize both the primitive Pythagorean seeds and the complete positive-integer coefficient family of the prescribed-edge absorbers developed in this session. The connection is a branch-preserving bijection, rather than a coincidence of individual numbers.

## 1. The common matrices

Use ordered column vectors and put

    A = [[ 1,-2, 2], [ 2,-1, 2], [ 2,-2, 3]],
    B = [[ 1, 2, 2], [ 2, 1, 2], [ 2, 2, 3]],
    C = [[-1, 2, 2], [-2, 1, 2], [-2, 2, 3]].

Each preserves Q(a,b,c)=a²+b²−c² and is invertible over the integers. The classical Berggren theorem says every primitive Pythagorean triple with odd first leg and even second leg occurs exactly once by applying a finite branch word to(3,4,5). See Keith Conrad's [Pythagorean descent notes](https://kconrad.math.uconn.edu/blurbs/linmultialg/descentPythag.pdf) for the classical matrix/descent construction. The source uses a different ordering of the matrix letters; the displayed matrices fix the convention here.

We now prove the analogous statement on Q=1 directly.

## 2. Complete odd norm-one ternary tree

Let

    S₁={(a,b,c): a,b,c positive odd integers, a²+b²=c²+1}.

**Theorem 1.** Every element of S₁ has a unique finite address in A,B,C starting from(1,1,1). Every node has precisely three children, and the third coordinate strictly increases along every edge.

### Positive children

For a norm-one triple, c≥a and c≥b, since c²−a²=b²−1≥0 and similarly for b. Thus every displayed matrix sends a positive triple to a positive triple. For example A's first coordinate is a+2(c−b)>0, its second is2a+b+2(c−b)>0, and its third exceeds c by2(a+c−b)>0. The other two cases follow by symmetry or direct positivity. Oddness and the norm are preserved.

### A unique smaller parent

Set

    r=a+2b−2c, s=2a+b−2c, h=3c−2a−2b.

The three inverse images are

    A⁻¹(a,b,c)=(r,−s,h),
    B⁻¹(a,b,c)=(r,s,h),
    C⁻¹(a,b,c)=(−r,s,h).

The numbers r,s are odd and therefore nonzero. They cannot both be negative. For example if a≥b, then

    c²=a²+b²−1 < (a+b/2)²,

because ab−3b²/4+1>0. Hence2c<2a+b and s>0. The other order similarly gives r>0.

Except at(1,1,1), we also have0<h<c. To prove positivity, if one leg is1 then c is the other leg, and h=c−2>0 unless both legs are1. Otherwise a,b≥3; with S=a+b≥6,

    c²≥S²/2−1>4S²/9,

so3c>2S. Finally c<S since S²−c²=2ab+1>0, and therefore h=3c−2S<c.

Exactly one inverse image is positive: choose B when r,s>0, A when s<0, and C when r<0. Its norm and oddness are preserved. Repeating the unique descent strictly lowers the positive integer c, and it can stop only at(1,1,1). This proves existence and uniqueness of every branch address. It also excludes branch collisions. ∎

## 2a. The full positive norm-one forest has exactly three roots

The parity restriction can be removed completely. The equation admits exactly three parity sectors: odd/odd/odd, even/odd/even, and odd/even/even. Reduction modulo4 excludes even/even/odd, and the other patterns fail modulo2.

**Theorem 2.** All ordered positive integer solutions of a²+b²=c²+1 form three disjoint rooted ternary trees under A,B,C, with roots

    (1,1,1), (2,1,2), (1,2,2).

The same unique-parent proof applies outside these roots. Indeed every nonroot has c≥3, and

    9c²−4(a+b)² ≥ c²−8 > 0,

so h=3c−2a−2b>0; again h<c. At least one of r,s is positive by the earlier comparison. The only possible zero r comes from a+2b=2c. Then a is even, say a=2k, and the norm equation gives

    k(3k−2b)=1.

Positive integrality forces k=b=1, hence the root(2,1,2). The case s=0 gives(1,2,2). Thus every nonroot has exactly one positive inverse image. Its parity sector is preserved by all matrices. The only positive triples with c≤2 are the three displayed roots, where h is−1 or0. Descent therefore terminates uniquely at the root in the appropriate sector.

All three matrices produce positive children of increasing third coordinate in every sector. Consequently the forest has no collisions within a sector, and parity separates the sectors. This proves the full classification.

## 3. An actual tree bijection with the primitive Pythagorean triples

Map a primitive Pythagorean triple M_w(3,4,5) to M_w(1,1,1), using the same unique branch word w. Classical uniqueness on Q=0 and Theorem1 on Q=1 make this a bijection. It commutes with every A,B,C child operation.

This is an isomorphism of the two labeled rooted ternary trees. More generally, the full norm-one forest is in branch-preserving bijection with three labeled copies of the primitive-Pythagorean tree, one copy for each parity sector. It is not a single linear transformation between their quadratic surfaces: it is defined by their unique finite addresses. Nor does it automatically transport square-sum Hamiltonicity or any other property not encoded by those addresses.

The ordered convention matters. On Q=0 the first leg is odd and the second even. In the first Q=1 tree all three coordinates are odd; the other two trees have the parity sectors stated in Theorem2. Swapping unequal legs corresponds to a different address in the odd norm-one tree; diagonal nodes are fixed by that swap. In the other parity sectors, leg exchange swaps the root sector. Thus the domains have been specified so no quotient by leg exchange is silently taken.

## 4. Every non-boundary node is a square-sum absorber

Take any positive norm-one triple(a,b,c) with a,b>1, and take distinct positive integers x,y with x+y=t². Define

    u=a²(x+y)−x, v=b²(x+y)−y.

Then

    u+x=(at)², x+y=t², y+v=(bt)², u+v=(ct)².

Thus an existing host edge u–v can be replaced by u–x–y–v. Both hosts are larger than x and y and are distinct, exactly as proved in [the prescribed-edge absorber note](square_sum_prescribed_edge_absorbers_20261019.md).

In the odd sector, the only norm-one nodes with a leg equal to1 are the two boundary rays

    A^k(1,1,1)=(1,2k+1,2k+1),
    C^k(1,1,1)=(2k+1,1,2k+1).

This follows both from the norm equation and the unique parent. These degenerate coefficients reproduce one omitted vertex as a host and are excluded from the insertion theorem. Every other node gives an admissible gadget. Once both legs exceed1, all children retain that property.

The two additional boundary rays are

    C^k(2,1,2)=(2k+2,1,2k+2),
    A^k(1,2,2)=(1,2k+2,2k+2).

Together the four rays exhaust all positive norm-one triples with a leg equal to1. Every other node in the three-tree forest is a nondegenerate absorber. Thus the forest enumerates every ordered positive-integer coefficient gadget of this linear form, without duplication. For example A(2,1,2)=(4,7,8), recovering the asymmetric gadget as well.

## 5. The Pell lane is exactly the central branch

On the symmetric subspace a=b, the central matrix B acts as

    (a,c) -> (3a+2c,4a+3c), c²−2a²=−1.

Starting from(1,1), this gives

    (a,b,c)=(1,1,1),(5,5,7),(29,29,41),(169,169,239),...

Every diagonal norm-one node has a pure B address: at any nonroot diagonal node r=s>0, so its unique positive inverse is again diagonal under B. Thus this lane exhausts the positive symmetric coefficients. The first nondegenerate node produces u=24x+25y and v=25x+24y. Thus the previously found Pell absorber is literally the central Berggren branch acting on the norm-one shell. The same matrix acts on the norm-zero Pythagorean shell, where its first image is(21,20,29).

## 6. The requested seed has an exact address in this dictionary

The primitive triple underlying482,3362,6242 is(9,40,41), and

    A³(3,4,5)=(9,40,41),
    A³(1,1,1)=(1,7,7).

Its corresponding odd-sector norm-one node is therefore on a degenerate boundary ray; this is a limitation of transferring that particular address within the odd sector. One central step enters the useful absorber region:

    B(9,40,41)=(171,140,221),
    B(1,7,7)=(29,23,37),
    29²+23²=37²+1.

The same original seed address in a different parity sector already gives

    A³(2,1,2)=(8,31,32), 8²+31²=32²+1,

which is nondegenerate without the extra central step. These yield explicit host insertions through any square-sum pair. The AP seed, the absorber family, and the Pell lane share an exact matrix action, while their admissibility conditions remain distinct.

## 7. Scope

The theorem supplies the complete positive-coefficient gadget forest and precise branch-preserving bijections with three copies of the primitive-Pythagorean tree. It does not ensure a chosen host edge lies in an interval-spanning backbone. That separate requirement is the subject of the prescribed-edge gluing work. It proves no Collatz or graceful-tree conjecture and makes no priority claim for the norm-one tree construction.

The companion exact verifier checks matrix identities, inverse choices, unique tree generation, a bounded exhaustive norm-one comparison, the absorber equations, and the displayed addresses. The infinite classification is proved by the descent above, not inferred from those checks.

## 8. Reproduction and independent checks

The standard-library-only [exact verifier](../experiments/verify_berggren_absorber_bijection_20261019.py) and [saved JSON receipt](berggren_absorber_bijection_controls.json) accompany this note. From the repository root, run:

    python 05-knowledge/experiments/verify_berggren_absorber_bijection_20261019.py --height 2001 --depth 8

To reproduce the receipt:

    python 05-knowledge/experiments/verify_berggren_absorber_bijection_20261019.py --height 2001 --depth 8 > /tmp/berggren_absorber_bijection_controls.json

The run uses exact integers and checks all 10,987 positive norm-one triples with c≤2,001 against independent exhaustive enumeration: 5,395 in the odd sector and 2,796 in each other sector. Exactly 4,001 are boundary nodes, leaving 6,986 nondegenerate coefficient triples. It also checks all 319 parity-oriented primitive Pythagorean triples under that height bound, 9,841 addresses through depth8 in each of the four trees, and 2,584,820 absorber instances. The zero inverse-coordinate exceptions, diagonal Pell classification, and displayed seed addresses are checked explicitly. The receipt includes the verifier's SHA-256.

The [independent proof audit](square_sum_berggren_absorber_bijection_20261019.audit.md) records the proof checks and scope. These finite computations do not replace the infinite descent proof.
