# Independent audit of block-uniform modulo-eight blindness

**PASS.** The theorem covers all 64 block-uniform physical queries for C3[A1,A2,A3] versus quotient converse, with arbitrary and potentially unequal even-order children of order at least four. It asserts zero residue modulo eight, not an exact response collision or equality on nonuniform named masks. No mathematical corrections were required.

## Universal profile reduction

I checked the dependency's even-order profile derivation. The classical input is the modulo-four Hamilton-path/odd-cycle congruence, verified directly in Grinberg and Stanley's [detailed paper](https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge-long.pdf), Theorem 7.1, printed page 119, June 13, 2026 version.

Applying it partwise to the exact unordered-partition formula gives the asserted ordered-cover congruences. For even child order every odd cycle is proper, and its contracted atom set has even order; the required Stirling parities therefore produce q2=2H-4 and q3=4 modulo eight when the order is at least four. Ordered slot permutations give divisibility by k!, hence by eight for every k>=4. H is odd. Thus the reduced profile (0,1+2s,-2+4s,4) is valid for every eligible child, independently of all other children.

The m>=4 condition is substantive. It excludes the order-two profile, whose q3 is zero and whose q2 has a different reduced value. The checker retains its nonblind modulo-eight control.

## Full query functional and all parameters

The run-word functional includes both port orders. The fixed x->y edge makes the middle segment optional only when x comes first; the other order requires a nonempty old segment. Empty prefix/suffix segments are permitted. The functional takes one ordered-cover entry per child, indexed by its total number of runs. Consequently each independent parameter appears affinely and the complete contrast is genuinely multilinear in s1,s2,s3.

After reduction, indices exceeding three vanish modulo eight, so truncating the formal run series at the stated finite indices is exact for this modular calculation. The source's eight Boolean parameter evaluations determine all eight coefficients of each multilinear polynomial. Divisibility of all 512 coefficients proves the identity at every integer parameter triple, not just the tested Boolean values. The parameters are independent: no equality of child profiles or orders has entered.

## Independent certificate and physical checks

The independent checker imports no code from coarse_mod8_boundary or covered_response. Instead it enumerates bounded child-run words directly, allowing each port exactly once and each child one to three runs. It expands the affine profile factors symbolically, without using the producer's rational matrix calculation or Möbius inversion. Every one of its 512 coefficients agrees exactly with the retained producer certificate and is divisible by eight.

It separately checks the profile reduction on all 64 labelled four-name tournaments and six declared six-name controls. Nine physical query pairs are recomputed by a subset DP over the actual old vertices, including children of unequal orders (4,4,6), unequal profiles, and nontransitive interiors. All residues vanish. The independent physical engine also recovers the sharp excluded TT2 query values 509 and 545, whose difference is four modulo eight.

Normal and optimized audit receipts are byte-identical. Integer bounds are safe: the literal controls use at most sixteen augmented names, below the engine's declared limit of seventeen, and 17! fits in unsigned 64-bit arithmetic. All mathematical checks use explicit exceptions.

## Ambient and exact-equality boundaries

The total eligible old core order is even, so adjoining the two ports leaves even order. Giving both ports the core module's exterior incidences makes the augmented core a common module. The separately proved even-order transfer theorem multiplies its Hamilton difference by an odd unit modulo eight, so zero residue remains zero in every such canonical ambient lift.

This argument does not cover exterior attachments that destroy the augmented module. It also does not constrain nonuniform masks inside a child. The explicit nonzero integer gap for three TT4 children is retained and is compatible with every modular statement. The exact named-response problem has not been replaced by this coarse observer.

The receipt pins both reviewed source files, the retained coefficient certificate, the even-transfer proof dependency, and independent checker/helper sources.
