# No aliquot three-cycle with a dyadic core and one odd prime at every node

19 October 2026. **PROVED, independently reviewed.** This resolves the scoped hypotheses A3 and D3 in the earlier dyadic frontier, not the general aliquot three-cycle problem. Historical priority has not been established.

## Result and inheritance

The earlier [dyadic frontier](../aliquot_dyadic_frontier_20261014e.md) proves that repeated dyadic cores cannot support a primitive aliquot three-cycle, classifies their sole nonconstant positive-odd formal exception, and leaves two explicit hypotheses:

- A3: three pairwise distinct powers of two admit no integral cyclic prime-variable solution, even before primality is imposed.
- D3: no primitive aliquot three-cycle lies entirely in the class 2^a p, a>=1 and p an odd prime.

Both follow from the theorem below together with that proved repeated-core classification. No larger finite exponent census is needed.

**Theorem.** Let c,d,e be three pairwise distinct powers of two, each at least two. The cyclic affine system

    d q = (c-1)p + 2c-1,
    e r = (d-1)q + 2d-1,
    c p = (e-1)r + 2e-1

has no solution with p,q,r all integers. Indeed the variable attached to the largest core is already nonintegral.

**Corollary.** There is no primitive aliquot cycle of length three all of whose members have the form 2^a p with a>=1 and p an odd prime. The prime may vary between nodes, as may the exponent; there is no height or exponent bound.

The source identity is s(2^a p)=(2^a-1)p+2^(a+1)-1. It requires an actual odd prime p; the corollary retains that guard. The theorem itself is stronger and purely arithmetic.

## 1. Sort the cores without discarding cyclic orientation

Write c<d<e for the sorted core values. Then

    c>=2,   d>=2c,   e>=2d.

After cyclic rotation, the two possible directed orders are (c,d,e) and (d,c,e). In either order, call r the variable attached to e. The common denominator is

    D=cd+ce+de-c-d-e+1>0.

It is odd and D=1 mod c. The forward orientation (c,d,e) has r=N_f/D, with

    N_f=2c²d-2c²+2cd²+2cde-3cd-2ce+2c-2de+d+2e-1.

The reverse orientation has r=N_r/D, where

    N_r=N_f+(c-d)(2(c+d)-1).

These are exact Cramer/elimination formulas for the three displayed equations. In particular N_r<N_f. All three cyclic variables have denominators dividing the same odd integer D. Therefore integrality of any one variable propagates around the cycle: its edge equation makes the next variable times a power of two integral, and a denominator dividing odd D cannot be canceled by that power of two. This is the predecessor's integrality-propagation lemma. Consequently proving the largest-core variable nonintegral also proves A3's original first-variable nondivisibility statement. Since c divides both d and e, both numerators are -1 mod c. Therefore any integral r obeys

    r=-1 mod c.                                             (1)

## 2. A two-sided size trap leaves exactly one integer candidate

We prove, simultaneously for both orientations,

    c-1 < r < 3c-1.                                        (2)

It suffices to prove N_r>(c-1)D and N_f<(3c-1)D.

The coefficient of e in N_r-(c-1)D is

    (c-1)(d-c-1)>0.

At e=2d, the same expression is

    d[4(c-1)d-c(c+1)+3]+c(c-1).

Because d>=2c, the bracket is at least 7c²-9c+3>0 for c>=2. Thus the lower inequality is strict for every e>=2d.

The coefficient of e in (3c-1)D-N_f is

    3c²+cd-2c+d-1>0.

At e=2d, this expression is

    c[(7c-5)d-c+2]+2d(d-1)>0.

For example d>=2c makes the first bracket at least 14c²-11c+2>0. This proves the strict upper inequality and hence (2).

Combining (1) and (2), the only possible integral value is

    r=2c-1.                                                (3)

The bounds are important: merely approximating r by a harmonic mean would not give this exact candidate. Both directed cyclic orders have been retained.

## 3. Forward orientation: its 2-adic valuation is impossible

Substitute (3) into the forward equation N_f=(2c-1)D. Rearrangement gives

    e(2c²-c+d-1)=c(2d²-1).                                (4)

The parenthesis on the left is odd; so is 2d²-1. Consequently v_2(e)=v_2(c). Since both are powers of two, e=c, contradicting c<d<e. There is no forward integral candidate.

## 4. Reverse orientation: one exceptional valuation case also fails

The reverse substitution N_r=(2c-1)D gives

    e(2c²-c+d-1)=2(c-1)(c+d²)+d.                         (5)

Write c=2^a and d=2^b, with b>a>=1. The left parenthesis is odd.

If b>a+1, the three terms on the right of (5), written as

    2c(c-1) + 2d²(c-1) + d,

have valuations a+1, 2b+1, and b. The first is uniquely smallest. Thus v_2(e)=a+1<b, contradicting e>d.

It remains that b=a+1, so d=2c. Equation (5) becomes

    e(2c²+c-1)=2c²(4c-3).                                (6)

Both 2c²+c-1 and 4c-3 are odd. If e is a power of two, valuation in (6) forces e=2c². Cancellation then requires

    2c²+c-1=4c-3,
    2c²-3c+2=0.

The last expression is positive for every real c (its discriminant is -7), and in particular for c>=2. This is impossible. The reverse orientation has no integral candidate either, proving the theorem.

## 5. Close the native aliquot statement

In a hypothetical primitive cycle, mark each source as 2^a p with p an odd prime. The aliquot formula gives the cyclic affine system. If its three cores are pairwise distinct, the theorem excludes even integrality.

If two cores coincide, the earlier dyadic frontier's Section 3 proves that the only positive-odd formal possibilities are a fixed point or the nonconstant core pattern (2,2,16), up to rotation, with variables (23,13,1). The latter is not native because 1 is not prime; its displayed last affine edge is not an aliquot edge. The repeated-core theorem is independent of any finite census or unproved part of A3. Thus every possible core pattern is excluded, proving D3.

This does not exclude three-cycles containing a node with two odd prime factors, an odd prime power of exponent greater than one, or an odd node outside the specified class. It does not classify other cycle lengths or decide trajectory termination.

## 6. Mechanism and verification

The useful transfer is a combination of an ordered size interval and a retained congruence. The size interval reduces an unbounded prime variable to one residue-compatible candidate; its valuation then contradicts the ordering of the cores. It resembles the earlier repeated-core trap, but works for both orientations with three distinct exponents. Neither a numerical prime search nor a statistical abundance argument is used.

The companion verifier checks the exact elimination formulas, the four coefficient/boundary identities in Section 2, the two forced equations, and finite controls. The proof above supplies the all-exponent quantifiers. Prior repeated-core results and their native-prime exception are explicitly inherited rather than reclassified by a larger census.
