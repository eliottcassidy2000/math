# Session131: reflection, nonconsecutive numerals, square-sum seams and moonshine

**PROVED scoped extensions; exact finite certificates; cited arithmetic classifications.** The session makes the even/odd analogy precise through integral exchange blocks and carries. It distinguishes four carriers: metagraph cycle chains, canonical numerals, actual square-sum paths, and Frobenius permutations. No global tournament-spectrum or new moonshine theorem is claimed.

## 1. Reflection has a complete integral normal form

[The constructive module theorem](integral_normal_form/THEORY.md) extends the previous index1024 result. For any vertexwise-fixing graph involution with fixed edges and exchanged edge pairs, choose a blue spanning forest and extend it with black edges. Its cycle lattice has an explicit unimodular basis of fixed, sign and exchanged blocks.

At order5 the exact decomposition is

**`10 exchange blocks + 1 trivial block + 2 sign blocks`.**

Each exchange pair(x,y) has eigen-coordinates s=x+y,d=x-y. Recovering x,y integrally requires s=d mod2. The ten exchange blocks explain all ten parity constraints; the two sign blocks require no such gluing. This is more information than a rational eigenspace decomposition.

Applying the same odd sign-magnitude CG shift to each normal-form coefficient gives an injective, reflection-equivariant map of integral cycle chains. It preserves boundary zero but is basis-dependent and nonadditive; it need not preserve simple or Hamiltonian cycles. This is the lawful representation-level bridge.

## 2. Nonconsecutive representations need a carry law

[The numeral packet](zeckendorf/PROOFS.md) retains both canonical systems:

- ordinary Zeckendorf binary digits have no adjacent1s;
- Chung–Graham digits on1,3,8,21,... use0,1,2 and forbid2 1* 2.

For equal-parity eigen-coordinates write s=2u+epsilon,d=2v+epsilon. Then

`x=u+v+epsilon`, `y=u-v`.

Adding two such triples has the exact carry law

`(u,v,epsilon)+(u',v',eta)`
`=(u+u'+epsilon*eta, v+v'+epsilon*eta, epsilon XOR eta)`.

Reflection sends v to−v−epsilon. The ten-bit atlas code and the graph basis preserve the constraints that independent strands would lose.

Both paired canonical languages admit a **minimal17-state automaton** recognizing equal numerical parity. Full reachable-state and distinguishing-suffix certificates establish minimality. Splitting ordinary Zeckendorf positions into rows e_j=z_(2j+2),o_j=z_(2j+3) is legal exactly when the rung and diagonal exclusions hold. A reversible two-state converter sends these to signed even-Fibonacci digits `c_i=e_i+o_(i-1)-o_i`; this is not an assertion that simply pairing bits gives canonical CG digits.

The crucial operation theorem is sharper than a parity observation. For two nonnegative canonical CG numerals, the shift defect

`P(a)+P(b)-P(a+b)`, with `P(n)=ceil(phi² n)`,

counts exactly the normalization carries touching the lowest weight, independently of rewrite strategy. It is0 or1, so there is at most one such boundary carry. Signed correction terms lie in{-1,0,1}.

Independent one-place eigenstrand shifts can fail: `(s,d)=(3,1)` becomes `(8,3)`, while shifting the actual twin coefficients first yields `(9,3)`. There is also an exact repair without changing coordinates: **three-place shifts are the smallest universal shifts preserving parity gluing**, in both systems. An initially compatible pair either remains compatible at every shift or precisely at multiples of3. This follows from the order-three Fibonacci parity recurrence, not an observed coincidence.

## 3. Square-sum cycles have a sharp residuewise seam budget

[The square-sum theorem](square_sum/PROOFS.md) concerns actual Hamiltonian cycles on[1,N], N>=8. Call an edge an even-root edge when its square sum has even positive square root. The number of these edges is at least

**`(4,7,8,7,8,7,8,7)[N mod8]`.**

Explicit cycles atN=80,...,87 attain all eight entries. Thus the bounds are sharp uniformly by residue class; attainment at everyN is not claimed. In particular the inherited uniqueG32 cycle has16 even-root edges, not4.

The proof removes the even-root edges and counts every path endpoint incidence in the four odd-square channels

`{0,1}, {2,7}, {3,6}, {4,5} mod8`.

Matching the remaining endpoints gives exact flow and parity equations. AtN≡4, the apparent four-seam solution necessarily closes one block separately, so connectedness forces eight. The marked-path version retains its two exposed ports and the actual endpoint matching. This avoids the earlier false step of treating local sewn cycles as every global Hamiltonian cycle.

The odd-root subgraph has exactly these four connected components for everyN>=44. This is eventual stabilization, not the first isolated occurrence: N=8 already gives four matching components, then arrivals9..12 create new components.

A root word s_1,...,s_m acts by partial reflections `x -> s_j²-x`. The packet gives exact closure, interval and collision tests, retaining the initial point. Uniform affine maps preserving every square-sum edge are exactly `x -> r²x`. By contrast, the CG shift sends the square edge1--3 to3--8, whose sum11 is nonsquare. Reflection/numeral representations therefore do not automatically transport Hamiltonicity.

## 4. The1–2–4 decomposition must include carries

[The binary refinement](binary_carries/PROOFS.md) recovers the audited1–2–4 antecedent and THM-466's existing warning that H does not determine its odd-cycle packing coefficients.

For the previous reflection chargesq andi, the actual low binary digits are

`b0=1`, `b1=q mod2`, `b2=(floor(q/2)+i) mod2`.

Thus H≡7 mod8 exactly when either q≡1 mod4 andi is odd, or q≡3 mod4 andi is even. The coefficientsq,i are counts, not independent bits. H9 withq2,i1 andH15 withq3,i2 are concrete carry controls.

The positive lift is exact:

`H=1+2q+4i+8R`,
`R=j+p+sum_(k>=3) 2^(k-3) alpha_k`.

Herej counts disjoint mirror-cycle pairs, p counts nonfixed reflection orbits of disjoint two-cycle packings, andalpha_k counts k-cycle packings. If t counts fixed three-packings, then

`H=1+2q+4i+8(j+p+t) mod16`.

The next bit includes the carry `floor((q+2i)/4)`. At order6, a transitive tournament hasH1 andC3->C3 hasH9, although both have(q,i)=(0,0); q,i cannot determine modulo16. At order9 the three-triangle stack forces the fixed-three-packing term. These are refinements, not a realization theorem for arbitrary coefficient lists or a proof of the forbidden7/21 spectrum.

## 5. The prime groups: exact inherited meanings and one conditional correction

[The primary-source audit](moonshine/README.md) retains the supplied sets

`A={2,3,5,7,11,13}`,
`Btail={17,19,23,31,47}`,
`C={29,41,59,71}`,
`D={37,43,53,61,71}`.

Established group orders give prime supportsFi22=A, Baby=A∪Btail, Monster=A∪Btail∪C. The prior2A-centralizer theorem gives the same three shells: forH=2.B andI=|M|/|H|, A is the common prime support, Btail occurs only inH, andC only inI. These results were already in the repo and are credited, not rediscovered as new classifications.

The first six primesA are precisely those whose squares divide the Monster order. The four PRIMEGAME data labels{2,3,5,7} are precisely those whose squares divide the Baby order. The remaining six program labels are{11,13,17,19,23,29}. These numerical identities are true, but prime-renaming conjugacy shows the algorithm does not require these arithmetic labels.

D intersects the Monster set at71. If the intended object was its complement among primes<=71, the set is{37,43,53,61,67}. No typo is assumed: the supplied D instead contains four genus-one primes and genus-zero71, while the missing67 has genus two.

The absolute-supersingular convention used is Ogg's: every supersingular j-value in characteristicp lies inF_p. It is equivalent to genus zero for the **Fricke quotient**X_0^+(p), not ordinaryX_0(p). The exact fifteen-prime classification is cited. In particular11 belongs even thoughX_0(11) has genus1.

## 6. A real Frobenius bridge, with a provable obstruction

On the integral permutation lattice of supersingular points overF_(p²), Frobenius hasf fixed points andg exchanged pairs:

`L=Z_trivial^f + Z[C2]^g`.

Its odd rank isg and its even-plus-odd index is2^g. The cited arithmetic theorem identifiesg with the appropriate quotient-genus defect. Forp=37,43,53,61,67,71 the classical genera are1,1,1,1,2,0. At leveltwo, each of29,41,59,71 contributes one nonrational Frobenius pair; the packet gives the four exact quadratic polynomials and independently verified finite-field computations.

This is the same **exchange-block gluing mechanism** used above. It is not a whole-module identification: a permutation lattice has even rank at least odd rank, while the n5 metagraph has ranks11 and12. Its two extra sign blocks obstruct such an isomorphism even overQ. At n4, one trivial plus one sign block has index1 but a nonzero odd sector, showing why the Ogg criterion cannot be exported to arbitrary reflection lattices from an index alone.

## 7. Verification and next boundary

[Verification](VERIFICATION.md), [independent audit](independent/AUDIT.md), and `verify_session.py` bind final proofs, data and checkers. The square-sum witnesses are checked directly with the standard library; their search software is not trusted for validity or optimality. The universal lower-bound proof supplies optimality. Automaton minimality has explicit certificates. Primary arithmetic inputs remain cited, while finite polynomial checks give independent controls rather than claiming a new all-prime proof.

The synthesis supplies precise representation and flow laws. It does not identify numeral forests with Hamiltonian forests, binary packing coefficients with bits, prime-register labels with intrinsic moonshine data, or rational eigenspace splitting with integral independence. Those distinctions are the information that makes the connections useful.
