# Independent audit: Niemeier genus-5 restrictions

Audit date: 2026-10-22. No repository changes or Git writes made.

## Verdict

The attached manuscript's proposed hyperelliptic separation is wrong. All five equal-Coxeter pairs coincide on the hyperelliptic locus in every genus. There are exactly **19**, not 23, hyperelliptic restriction classes, including in genus 5. Its desired Jacobian result **gJac(N24)=5** nevertheless survives, unconditionally, by an existing theorem, not by its proposed hyperelliptic argument.

## Verified primary source and exact applicability

Giulio Codogni, *Hyperelliptic Schottky problem and stable modular forms*, Documenta Mathematica 21 (2016), 445–466. [Publisher](https://ems.press/journals/dm/articles/8965467); [full text](https://emis.de/ft/9558). A downloaded copy and text extraction are saved beside this review.

Theorem 6.1, printed p461 (PDF page17), applies to rank-24 forms with equal norm-two-vector counts. Its conclusion is stable hyperelliptic vanishing. Definition1.1, p446, means vanishing in every genus. The standing assumptions, pp452–453, are even, positive definite, unimodular; Definition3.6 uses exp(πi∑Q(xi,xj)τij). This is precisely manuscript Eq1.1, not a different theta normalization. The theorem explicitly lists all five pairs below. Theorem6.4, p464 (PDF page20), states genus≤4 Jacobian vanishing and a slope-12 divisor in genus5. For distinct pairs, this entails nonzero genus5 Jacobian restriction. As usual the theorem's latter clause requires distinct isometry classes; literally taking P=Q gives zero, not a divisor.

The printed p461 table lists 72,120,144,216,360 under “# roots”; these are half the full signed counts. This apparent convention/typographical issue does not affect equal-count applicability. The manuscript uses full counts24h.

| h | Root systems of the two Niemeier lattices | Full root count of each |
|---|---|---:|
|6|D4^6; A5^4 D4|144|
|10|D6^4; A9^2 D6|240|
|12|E6^4; A11 D7 E6|288|
|18|D10 E7^2; A17 E7|432|
|30|E8^3; D16 E8|720|

Here root-system notation labels the unimodular Niemeier overlattice, not an unglued root lattice. In particular D16E8 denotes D16+⊕E8. Reversing the source's pair order only changes a sign.

## Stronger diagnosis: Hypothesis4.6 is actually false

This is an independent elementary coefficient computation, not an inference that merely one of the manuscript's arguments must fail.

Let a,b,c be an ordered A3 simple-root chain, and let d be a root orthogonal to all three. Its Gram matrix is Cartan(A3)⊕[2]. Put T0 equal to half that Gram matrix. Fourier coefficients count ordered tuples with this Gram matrix, without dividing by Weyl groups or automorphisms.

In a D4 component there are 576 ordered A3 chains and no root inside that component orthogonal to a whole chain. In an A5 component there are 720 such chains and exactly two roots inside that component orthogonal to each chain. Both facts can be checked directly with roots ±ei±ej for D4 and ei−ej for A5; the accompanying script exhausts exactly these finite root sets using integer arithmetic.

Both h=6 Niemeier root systems have 144 roots. Every connected A3 chain lies in one root component. Therefore:

- D4^6: 6×576×(144−24)=414720 tuples.
- A5^4 D4: 4×720×(144−30+2)+576×(144−24)=403200 tuples.
- Thus a(T0,F6^(4))=11520, with the manuscript's ordering F6=Θ(D4^6)−Θ(A5^4D4).

This frame gives a six-root zero-sum configuration

u=(d,−d,a,b,c,−a−b−c).

The last entry is a root of A3. The six roots span rank4 and contain a zero-sum pair. In the manuscript's cumulative variables Eq4.6 they give x=(d,0,a,a+b,a+b+c). Its rank4 nonzero basis has Gram matrix integrally equivalent to Cartan(A3)⊕[2]. Hence the configuration meets Hyp4.6's stated rank4/zero-pair condition, while its Fourier coefficient is nonzero. Manuscript Eq1.4 identifies F6^(4) with a scalar multiple of E4^(4)J8^(4), so the corresponding E4J8 coefficient cannot vanish either. This explicitly refutes Hyp4.6 as written.

The conceptual mistake behind Remark4.7 is that a modular form vanishing upon restriction to Jacobians does NOT force selected ambient Fourier coefficients to vanish. The pullback involves cancellations among terms and period corrections.

## Page-bound consequences for the attachment

- p11, Hyp4.6: false as stated, by the explicit rank4 counterexample above.
- p12, Thm4.8: its elementary component-partition argument survives. Six zero-sum roots distributed among D4 components cannot span rank5.
- pp12–13, Thm4.9: formally states an implication under Hyp4.6. A conditional with a false premise is not a counterexample to Codogni, but yields no geometric separation theorem. Its conclusion about actual H5 is false.
- p14, Thm4.11 and Eq4.15: hyperelliptic nonvanishing is false. Ambient nonzero root-count coefficients, even if correct, do not establish a nonzero pullback to H5.
- p14, Table5.1: H5 row must say that none of the five pairs separate, and the count is19. Jacobian row can say all five separate unconditionally, citing Codogni6.4. The trigonal row's assertion that the other four separate is not established by these arguments or by the two Codogni theorems checked here; do not preserve its 23-class count without a separate proof. This review does not independently decide that trigonal count.
- pp14–15, Thm5.1: numerical conclusion gJac(N24)=5 survives, but the proof through H5 must be replaced. Calling that conclusion new would also require addressing the prior source.
- p1 abstract: proposed nonzero H5 restrictions and conditional proof mechanism require correction, while ambient coefficient results are not disproved by hyperelliptic vanishing.

## Exact strongest corrected result

For every genuine equal-h pair h∈{6,10,12,18,30}, Fh^(g)|H_g=0 for every positive genus g, Fh^(g)|M_g=0 for g≤4, and Fh^(5)|M5 is not identically zero. Consequently the 24 Niemeier theta series have exactly19 hyperelliptic restriction classes in every positive genus, exactly19 Jacobian restriction classes for genera1–4, and24 Jacobian restriction classes in genus5. Thus gJac(N24)=5 unconditionally.

The exactness of19, rather than merely an upper bound, follows from distinct genus-one theta series for distinct h and degeneration: smooth hyperelliptic curves can approach stable hyperelliptic boundary curves whose Satake limits contain an elliptic Jacobian, so equality of restrictions would force equality of genus-one theta series. Equivalently use a compact-type hyperelliptic chain of elliptic curves, attaching at involution-fixed points, and theta factorization, then degenerate all but one elliptic factor to the cusp. The elliptic q coefficient is24h, distinguishing all19 values (including h=0 for Leech). The same boundary argument distinguishes unequal-h pairs on M_g. This derivation is separate from the theorem's same-h statements.

No claim is made here that genus5 pair differences are linearly independent: pairwise distinct restriction classes are a different, weaker assertion.
