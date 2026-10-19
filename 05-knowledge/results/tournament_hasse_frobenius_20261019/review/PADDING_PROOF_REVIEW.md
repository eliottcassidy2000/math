# Independent review: Hasse deletion, p-complexes, and integral reconstruction

2026-10-19. **PASS: proved identities and explicitly bounded exact controls.**

Reviewed producer artifact: `session6-automata/PADDING_FROBENIUS_PROOF.md`, SHA-256 `7efc364d5bcf48ee1d3cb4cb10943119db4c0fccf5c99b36a0e850d46720bedd`. The producer changed only the status header after the mathematical review; this hash includes that status-only update. This review covers the all-prime and all-digit contraction, carry and borrow identities, pointwise embedding observability, integral embedding-lattice warning, source-free reconstruction basis, integral cokernels, and rank/string consequences. No priority claim is assessed here.

## 1. Type and algebra audit

The function module includes the empty tournament and uses zero negative degrees. The subset-deletion maps increase degree; source evaluation h decreases degree. These conventions make

    D_a D_b=binom(a+b,a)D_(a+b),
    hD_m-D_mh=D_(m-1)

valid even at the bottom degree. The integral proofs preserve all vertex-subset multiplicities.

For d=D_1, d^p=p!D_p=0 over F_p. Normal ordering of d^a h^(p-1) and the hockey-stick identity leave only the coefficient

    (-1)^(p-1)(p-1)!=-1 modulo p.

Therefore the minus sign in the displayed p-contraction is correct, including p=2. If d^k f=0, precisely the terms indexed j>=k disappear. Factoring d^(p-k) from the remaining terms gives the producer's primitive of degree n-p+k. This proves every kernel/image equality, including the negative-degree injectivity boundary, without a dimension count.

For s=p^e, the coefficient of D_(ps) in D_s^p has p-adic valuation one. The commutator [h^s,D_s]=I follows because every internal binomial coefficient binom(s,j) vanishes modulo p. The same normal-order proof gives every all-digit contraction and correctly preserves degrees modulo s.

The reconstruction of D_m from the base-p digits of m has local factorials m_e!, each invertible modulo p. It does not divide by m! when m>=p. The carry criterion concerns the composition scalar, or the global composition operator: an individual vector may be annihilated even without a carry.

The cross-digit commutator is exactly zero for e>f, identity for e=f, and D_(p^f-p^e) for e<f. This follows by setting r=p^e,m=p^f in the integral higher commutator and retaining only j=0 and, when present, j=p^e. Independent digitwise source coordinates would therefore be a false strengthening.

## 2. Embedding carrier audit

For an n-vertex host T, its automorphism group acts freely by postcomposition on spanning embeddings of every oriented graph S. Thus every such count is divisible by |Aut(T)|. Taking S=T realizes precisely |Aut(T)| times the class indicator. Both conclusions follow:

- Over F_p, the visible span consists exactly of functions vanishing at hosts with p-divisible automorphism order.
- Integrally, the embedding lattice is the direct sum of |Aut(T)| times each class-indicator lattice.

The ordinary d and h preserve the embedding span by adjoining a labelled isolated vertex and summing deletions of source-eligible pattern vertices. Their contraction therefore restricts to the pointwise-reduced visible span and descends to its quotient.

The odd-prime higher-Hasse escape is genuine: D_p of the empty indicator is the constant function at degree p, which is nonzero on a cyclic p-vertex tournament, whereas every spanning embedding count there vanishes modulo p. The formula for adding m labelled isolates has the factor m!, not its inverse.

The added integral-lattice warning is essential and correct. The pointwise image of E tensor F_p is not generally an injective copy of that tensor product. For the cyclic triangle, d^3 of the empty indicator has lattice coordinate 2 in the basis vector 3 times the cyclic-triangle indicator. It survives in E tensor F_3 although the same function is zero pointwise modulo three. Thus d^p=0 on ambient pointwise functions does not imply d^p=0 on the abstract reduced embedding lattice.

## 3. Integral reconstruction and cokernels

The Taylor projector

    P=sum_j (-1)^j D_j h^j

is finite on each homogeneous component. The commutator makes hP telescope to zero, and P fixes ker h. Expanding sum_j D_j P h^j leaves the coefficient sum_k (-1)^k binom(ell,k), which vanishes for every positive ell and equals one for ell=0. This proves reconstruction.

For a in ker h, repeated commutation gives h^i D_j a=D_(j-i)a when i<=j and zero otherwise. Applying P then gives the Kronecker-delta extraction identity. Hence the decomposition is a direct integral sum and its change of basis is unimodular, with the producer's explicit inverse.

The unique universal source gives a bijection between source-bearing classes at degree n and all classes at degree n-1. Therefore ker h has precisely the source-free class-indicator basis claimed. Applying D_(n-r) turns these indicators into induced source-free subtournament counts. The inclusion-exclusion coefficient formula has exactly n vertices in every augmented host, so it is correctly typed. The warning that same-order source-free cores remain present prevents a false tournament-reconstruction claim.

D_m acts on each primitive string by multiplication by binom(n-r+m,m). Therefore its integral diagonal presentation and full cokernel direct sum are correct. The entries are all positive, proving integral injectivity; unused target coordinates give the stated free rank. Prime-adic sorting of the complete diagonal list gives the actual Smith invariant factors. The producer correctly distinguishes this from the unsorted binomial diagonal.

The modular rank formula follows by reducing the unimodular bases. The first injectivity failure at n=p-m is correct only for 1<=m<p, exactly as qualified. Ordinary p-string starts and higher-digit p-string starts follow either directly from the digit coefficients or from the graded dimension identity T(z)=B(z)(1+z^s+...+z^((p-1)s)). The displayed generating functions and coefficient formulas match both derivations.

## 4. Independent exact controls

The independent verifier is `verify_padding_independent.py`; run it normally or with Python optimization. It constructs tournament classes independently through order five by canonicalizing vertex-extension graphs, obtaining counts 1,1,1,2,4,12. It checks every integral divided-power identity within the bounded degree range, the source commutator, every source-free basis determinant, hP=0 and projector idempotence, and the integral diagonalization of every D_m with source and target order at most five. All six independent basis determinants are one, and source-free ranks are 1,0,0,1,2,8.

It also checks seventeen full-function p-contractions within the independently constructed degree range, the no-carry scalar atlas for p=2,3,5,7 and 0<=a,b<=64, 1,152 exact-binomial cross-digit borrow identities, and embedding blindness for all 27 oriented graphs on the cyclic three-vertex host. Normal and optimized JSON receipts byte-match. These are controls, not replacements for the universal proofs.

No mathematical defect remains in the reviewed freeze. The two important distinctions to retain in any summary are: pointwise reduction versus abstract embedding-lattice reduction, and diagonal presentation versus actual Smith invariant factors.
