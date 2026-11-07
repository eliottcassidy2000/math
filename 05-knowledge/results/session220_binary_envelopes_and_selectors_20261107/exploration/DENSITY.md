# Density lane: no promoted advance

Status: OPEN; no stronger support lower bound or density theorem proved.

Read-only input snapshot: acbad6c8636cf8cbb675b4c947fa19e33ac0df23.

Inherited and checked sources:
- Session216 overnight constructions density/RESULT.md and its copies of THM-4887-one-apex-energy-at-scale-two-to-the-n-over-three.md and THM-4890-one-apex-collision-balance.md.
- Session217 exact inverse and gate limits run-entropy/PROOF.md.

The valid unconditional inherited conclusion remains
A(X) >= X^(2/3)/(60 (ln X)^(2/3)) for X>=4096.

## Explored directions and exact stopping points

1. Reversal-product approximation. Writing g(v)=v(1-rev_n(v)/2^n)+e(v)+delta(v) gives an O(n) error. After multiplication by 2^n, the product uncertainty is O(n 2^n), not O(n). Consequently an ordinary divisor bound for a single integer product cannot bound a fiber of g. The reversal constraint and the exact defect must both be retained.

2. Symmetric two-block products. The exact inherited expression is
   g(h)+g(c)+hc+[h odd][c odd]-1.
   Treating either g-coordinate as an arbitrary bounded perturbation returns to the already optimized strip estimate. No additional bounded-correlation estimate was proved. The inherited random-micro model is not evidence sufficient to replace true micro values.

3. Multiplicative closure. Joining tournaments permits multiplying counts, but a power-law lower bound for the initial set does not on its own prove a larger exponent for its multiplicative closure. No extra separation or factorization argument was established.

4. Affine maps g(v)+cv. Large slopes do force injectivity, but the increase in output scale pays for this and does not improve the inherited exponent. No physical rooted family with enough independently adjustable small slopes was constructed. Arbitrary arithmetic coefficients cannot be silently substituted for realizable rooted states.

## Small exact side calculation, not a density advance

Let v>=0, t=nu_2(v+1), and when t>=1 let q=(v+1)/2^t (odd). Direct recursion gives

- If v is even: g(v+1)-g(v)=-v/2-((v/2) mod 2).
- If v is odd: g(v+1)-g(v)=(2^t-1)q+1-(q-1)/2-(((q-1)/2) mod 2)>0.

For the second identity, g(v)=g(q-1) after deleting t trailing ones;
g(2^t q)=g(q)+(2^t-1)q+1; and
 g(q-1)-g(q)=(q-1)/2+(((q-1)/2) mod 2).
These identities explain the sawtooth behavior but provide no nontrivial
fiber control. No novelty is claimed for this elementary consequence.

No finite scan was used as an all-length density argument. No repository
files were changed, no theorem ID was reserved, and no publication was made.
