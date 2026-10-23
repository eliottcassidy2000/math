# A three-channel cycle gas from one reversed inter-block arc

Status: PROVED by explicit cycle decomposition; small cases independently checked. The old n9 real-root failure is credited in `CYCLE_GAS_REAL_ROOT_BOUNDARY.md`. This note supplies a structured infinite strong family and exact generating functions, not a new threshold claim.

## 1. Construction and a support constraint

For m>=2 take the ordered join of m directed triangles. Reverse one arc a->c, where a lies in the first block and c in the last. Call the result U_m. It is strong: each block is internally strong, the forward joins pass along the blocks, and c->a closes the block-level circuit.

All directed cycles other than the original m triangles contain c->a. Therefore at most one such new cycle can occur in a vertex-disjoint family. This fact, not a free-particle approximation, makes the exact cycle gas tractable.

Put I_m(z)=sum_F z^|F| over vertex-disjoint directed odd-cycle families. The inherited OCF gives H(U_m)=I_m(2).

## 2. Exact partition polynomial

THEOREM. For every m>=2,

I_m(z)=(1+z)^m + z/2 [9(z+10)^(m−2) − (z−2)^(m−2)].

Proof. Families without a new cycle select arbitrary original triangles, giving(1+z)^m. A new cycle starts at a after using the reversed closing arc and follows increasing blocks until it reaches c. In each visited block it follows a nonempty directed simple path.

Let t mark the number of vertices on that new cycle. In an endpoint triangle, exactly one path of each length1,2,3 starts at a or ends at c, giving P(t)=t+t²+t³. In an intermediate triangle, a path of each length has three choices, giving3P(t). An unvisited intermediate triangle may independently be absent or included as an original odd cycle, giving1+z. Hence the new-cycle contribution is z times the odd-t-degree part of

P(t)² (1+z+3P(t))^(m−2).

Evaluate its odd part by(f(1)−f(−1))/2. Since P(1)=3 and P(−1)=−1, the displayed formula follows. It counts each directed cycle once because the distinguished reversed arc fixes its cyclic starting point. QED.

For m2: I_2=1+6z+z². For m3: I_3=1+49z+7z²+z³. Although the parity-filter expression uses subtraction and a half, its coefficients are nonnegative integers by the counting proof.

## 3. Every member from m3 onward fails the free positive-mode test

The leading coefficients are

a_m=1,
a_(m−1)=m+4,
a_(m−2)=binom(m,2)+46(m−2).

If the monic polynomial had all real roots, write it as product_i(z+r_i). Comparing the top coefficients and applying Cauchy–Schwarz gives

a_(m−1)² >= [2m/(m−1)] a_(m−2).

For our coefficients the right side minus the left side equals

4(21m²−48m+4)/(m−1)>0 for m>=3.

Contradiction. Thus every U_m with m>=3 has a nonreal zero pair in its odd-cycle fugacity polynomial. In particular none admits I_m(z)=det(I+zK) with finite Hermitian positive semidefinite K. Interacting, signed or non-Hermitian representations remain possible and are not excluded.

## 4. Hamiltonian specialization and the disappearing parity channel

At z=2 the negative parity-filter channel z−2 vanishes. For m>=3,

H(U_m)=3^m+9·12^(m−2)=3^m[1+4^(m−2)].

The exceptional m2 value is17, not18: the exponent is zero and its channel contributes−1. Values starting at m2 are17,135,1377,15795,187353.

For m>=5 the Hamiltonian values obey

H_m=15H_(m−1)−36H_(m−2).

The starting index matters: the order-two recurrence fails at m4 if the exceptional H_2 is used. This is a precise parity cancellation, not a claim of a physical supersymmetry.

## 5. Generating function and generic recurrence

As a formal series in t,

sum_(r>=0) I_(r+2)(z)t^r
 = (1+z)²/[1−(z+1)t]
   + (z/2){9/[1−(z+10)t]−1/[1−(z−2)t]}.

Thus the generic three characteristic channels are z+1,z+10,z−2. Their elementary symmetric functions give

I_m=(3z+9)I_(m−1)−(3z²+18z−12)I_(m−2)
     +(z³+9z²−12z−20)I_(m−3),    m>=5.

The coefficients may have exceptional specializations, such as z2 where one channel vanishes. The finite matrix/recurrence interpretation does not erase the cycle-overlap interaction; it works because all new cycles share one explicitly retained arc.

## 6. Relation to the session's main question

The gas has two distinct mechanisms: arbitrary compatible original cycles give a product of binary occupations, while a mutually exclusive family of new cycles couples all blocks. Odd/even length is selected by t=±1; component fugacity z is a different variable. Conflating those two parities would destroy the formula.

This gives a concrete way to use the fermionic/bosonic analogy without assuming all cycle occupations are independent: declare the support exclusion, retain the common-arc interaction, and only then perform the parity filter. The formal transfer channels are exact. No continuum field theory, number-theoretic classification or metagraph-wiring theorem is asserted.
