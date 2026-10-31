# Constructive integral reflection normal form and coherent numeral growth

**PROVED for the stated graph involutions; FINITE-EXACT order-five basis.** This extends Session130's index calculation to a complete integral module decomposition. A basis is constructed explicitly; no classification of abstract involution modules is imported as a black box. No literature-priority claim is made.

## 1. Objects and normal form

Let G be a finite multigraph with involution R fixing every vertex. Each fixed edge is fixed pointwise, not reversed; every other edge belongs to an exchanged pair. Loops and parallel edges are retained. Choose coherent orientations on each exchanged pair. Let Q=G/R and let B be its spanning blue (fixed-edge) subgraph. Let P be the number of black edge pairs, and c(-) count connected components including isolated vertices.

The integral cycle lattice is `Z_1(G)=ker(partial:Z^E->Z^V)`. These are signed integral chains with zero boundary, not necessarily simple or Hamiltonian cycles.

**Theorem 1 (explicit integral C2-module).** There is an integral basis in which R has

- `beta_1(B)` fixed basis vectors;
- `c(B)-c(Q)` sign basis vectors, each sent to its negative;
- `P+c(Q)-c(B)` exchanged pairs of basis vectors.

Equivalently,

`Z_1(G) ≅ Z_triv^beta_1(B) ⊕ Z_sign^(c(B)-c(Q)) ⊕ Z[C2]^(P+c(Q)-c(B))`.

The multiplicities are intrinsic. The basis depends on spanning forests and orientations.

**Proof by forests.** Choose a spanning forest T_B inside every component of B. Extend it to a spanning forest T_Q of Q using black edges; exactly c(B)-c(Q) black edges are needed. Lift each black forest edge to one of its twins, obtaining a spanning forest T_G of G.

A blue edge outside T_B gives a fundamental cycle contained entirely in B: its endpoints already lie in one blue component and their unique T_Q path is the blue-tree path. This is a fixed vector. Each black forest edge has a second twin outside T_G; its difference from the chosen twin is a sign cycle.

For each black quotient edge outside T_Q, let z be the fundamental cycle using one twin and the T_G path. Use z and Rz as an exchanged pair. In the usual G fundamental-cycle basis, z has coefficient1 on its own chord, and Rz has coefficient1 on the other twin chord. Replacing the usual second twin cycle by Rz changes it only by a combination of black-forest twin differences. Hence the change of basis is integral triangular with diagonal units, so unimodular. The displayed vectors exhaust all fundamental chords and form a basis. QED.

A black bridge of Q contributes a sign block, not an exchange block. A black loop of Q contributes an exchange block. Their effects on integral gluing are therefore different even though both came from black twins.

## 2. Complete normal form at n5

Use the Session130 literal merged metagraph, with its owner display-node atlas: ten vertices, eight blue edges, twelve exchanged black pairs. The blue subgraph has three components and one cycle; Q is connected. The formula gives

**`Z_1(G_5) ≅ 1 trivial ⊕ 2 sign ⊕ 10 exchange blocks`.**

Thus the ranks of the even and odd eigenspaces are11 and12. Each exchange block has diagonal coordinates

`s=x+y`, `d=x-y`,

whose integral inverse exists exactly when `s=d (mod2)`. Its even-plus-odd eigensublattice has index2. The fixed and sign blocks have no such defect. Hence the full index is2^10, recovering the previous result while identifying the source of all ten compatibility bits.

The checker prints all23 basis vectors in the32 oriented-edge coordinates, the block types, and an exact determinant of±1 against a separate fundamental-cycle basis. It tests boundary and reflection identities directly. The independent auditor additionally exhausts4,165 small colored quotient multigraphs and360 seeded larger orientation/multiedge controls.

## 3. Interwoven numeral strands: the lawful operation

Let A_j=F_(2j+2), so A=(1,3,8,21,...). The inherited canonical Chung–Graham representation uses digits0,1,2 and forbids a block2 1* 2. Its shift is

`P(n)=ceil(phi^2 n)` for n>0, with P(0)=0.

Extend by sign-magnitude to all integers:

`K(n)=P(n)` for n>=0 and `K(n)=-P(-n)` for n<0.

This is a strictly increasing odd injection. No floating arithmetic is needed:

`P(n)=(3n+floor(sqrt(5n^2)))//2+1` for n>0.

The underlying numeral theorem is inherited from Session125. Ordinary Zeckendorf binary digits, which forbid adjacent1s, and this even-Fibonacci digit system are related but different canonical languages; the companion [carry packet](../zeckendorf/PROOFS.md) gives the explicit conversion.

**Theorem 2 (coherent equivariant growth).** Fix the integral normal-form basis from Theorem1. Apply K separately to every coefficient. The resulting map on Z_1(G) is injective, preserves zero boundary, and commutes with R.

**Proof.** Fixed coefficients remain in their fixed lines. On sign blocks, K(-n)=-K(n). On exchange blocks, applying the same function to the two coefficients commutes with swapping them. Every image is still an integral combination of cycle basis vectors and hence has zero boundary. Strict injectivity of K gives injectivity of the coefficient map. QED.

The map is not a group homomorphism: K(1)+K(2)=9 but K(3)=8. It depends on the selected basis and does not preserve arbitrary simple cycles, nonnegative edge chains, cycle lengths, or Hamiltonian paths.

**Why independent eigenstrand shifts fail.** For a single exchange pair, take `(s,d)=(3,1)`, which corresponds to `(x,y)=(2,1)`. Shifting s,d separately gives `(8,3)`, whose parities differ, so the inverse half-sum is not integral. The coherent operation shifts x,y first, producing `(6,3)` and hence `(s',d')=(9,3)`.

In general the correct operation is

`s'=K((s+d)/2)+K((s-d)/2)`,
`d'=K((s+d)/2)-K((s-d)/2)`.

The companion packet proves that its correction relative to independent K(s),K(d) is a pair of carries in{-1,0,1}; in the positive normalization case a strategy-independent bottom carry explains the defect. Thus the two strands are coupled by an explicit parity and carry law, not merely described as interwoven.

This numeral map is not a square-sum graph embedding:1+3=4 is a square, butP(1)+P(3)=3+8=11 is not. The [square-sum packet](../square_sum/README.md) retains the actual edge-reflection and interval/injectivity constraints and classifies the uniform affine maps that preserve every square-sum edge.

## 4. Frobenius: the same gluing mechanism, a different module

For a finite Frobenius-stable set of supersingular points over F_p^2, Frobenius has f fixed points and g exchanged pairs. Its integral permutation lattice is exactly

`Z^f ⊕ Z[C2]^g`.

Its odd rational sector has dimensiong, and its even-plus-odd integral splitting has index2^g. This is an elementary permutation-module statement, independent of the difficult arithmetic determining f,g. The [moonshine packet](../moonshine/README.md) supplies the version-specific cited genus/rationality theorem and exact arithmetic verification.

In the classical supersingular-j setting, the cited relation identifies g with the genus of X_0^+(p). The Monster/Ogg primes are those with g=0. At the level-two quotient, the four Monster-but-not-Baby primes each have one new Frobenius pair. These are genuine arithmetic instances of an exchange-block defect, not a claim that the metagraph has acquired a Monster action.

The difference is provable: the n5 metagraph lattice has odd rank12 and even rank11. A permutation lattice always has even rankf+g at least odd rankg. Therefore **no Frobenius permutation lattice can be isomorphic to this metagraph module, even over Q**, respecting the involution. The two extra sign blocks are essential.

Moreover index1 alone need not mean no odd sector for a general involution lattice. At n4 the metagraph has one trivial and one sign block, hence index1 but odd rank1. Ogg's vanishing criterion uses the permutation-module structure; it cannot be exported to arbitrary reflection lattices merely from matching index values.

## 5. Evidence and scope

`python3 check_normal_form.py` and `python3 -O check_normal_form.py` produce identical receipts. The graph atlas is pinned to the prior Session130 file, not downloaded at runtime. The checker constructs all basis vectors, verifies the unimodular determinant and all reflection blocks, and checks101 deterministic signed coefficient vectors under coherent growth. Explicit nonadditivity, diagonal-parity and square-sum hostiles are included; `--test-failing-gate` must exit1.

Theorems1–2 are general symbolic proofs with finite controls, not extrapolations from these101 vectors. The square-sum and moonshine transfers retain their own arithmetic hypotheses and cited inputs. Neither an all-order tournament classifier nor a numeral-based Hamiltonicity proof follows.
