# Smooth lift trees and exact 2-power/3-power quotient-ring thresholds

Status: PROVED elementary statements; finite controls in certificate.json. No integer magic square of distinct squares is constructed. No external novelty claim is made.

## 1. Two different infinite families

Let L=(8,1,6;3,5,7;4,9,2) be the ordinary Lo Shu square. Its common sum is 15. Entrywise affine changes preserve its magic equations and distinctness over the integers.

A=3L-2=(22,1,16;7,13,19;10,25,4), common sum39.
B=8L-7=(57,1,41;17,33,49;25,65,9), common sum99.

**Theorem 1.** The quotient ring Z/3^k Z is non-Parker exactly when k>=3. The quotient ring Z/2^k Z is non-Parker exactly when k>=7. Here k is a positive integer and square entries must be nine pairwise distinct residues; zero is permitted.

Proof of existence. Every integer congruent to1 modulo3 is a square modulo every3^k. Begin with root1 modulo3, and given a root r modulo3^k choose the unique t modulo3 making (r+t3^k)^2 equal to the chosen integer modulo3^(k+1). The coefficient2r is invertible. Every entry of A meets this condition. The entries are distinct modulo27 (all lie between1 and25), hence remain distinct modulo every larger3-power. Its eight line equations hold over the integers already.

Likewise every integer congruent to1 modulo8 is a square modulo every2^k. For induction from modulus2^j to2^(j+1), j>=3, replace an odd root r by r+t2^(j-1). Its square changes by t2^j r modulo2^(j+1), so t in{0,1} chooses the next bit. This establishes existence at each level (it is not a claim that every residue root has a lift preserving all its bits). The entries of B are1 modulo8 and are distinct modulo128, giving the second existence claim.

Proof of nonexistence. For3-powers, the square residues modulo3 are0 and1. Each opposite pair sums to twice the center. If the center is0 then all entries are0 modulo3; all square roots are divisible by3, so all entries are divisible by9. Division by9 converts a distinct square magic array modulo3^k into one modulo3^(k-2), preserving distinctness, whenever k>=2. If the center is1, all opposite entries are1 modulo3. There are exactly3^(k-1) unit-square residues modulo3^k, so nine distinct entries require k>=3. The modulus3 base case has too few elements. Thus descent excludes every k<3.

For2-powers with k>=2, square residues modulo4 are0 and1. Opposite pairs force the entire array to be0 modulo4 or1 modulo4. In the first case divide by4 to obtain a distinct square magic array modulo2^(k-2). In the second, all roots are odd, and all entries are1 modulo8 for k>=3. There are at most2^(k-3) such residues, so k>=7 is necessary. The smaller base rings have fewer than nine square residues. Descent proves the claim.

Division by p² above is performed on residue classes known to be multiples of p²; it is a well-defined bijection from that ideal to Z/p^(k-2)Z. The roots divide by p. Both magic equations and entry distinctness survive.

**Scope.** These rings are not the extension fields F_(2^k) or F_(3^k). In particular every characteristic-two field is Parker, whereas its same-cardinality quotient ring becomes non-Parker at128.

## 2. Infinite no-wrap residue towers still do not solve the integer problem

A is a canonical no-wrap residue-square magic grid for every3^k>=81, since39<81. B has that property for every2^k>=128, since99<128. Their entries remain fixed integers. Many of those fixed integers are not integer squares. Therefore even infinitely many compatible p-adic square-root constraints, with fixed positive magic entries and arbitrarily large no-wrap moduli, do not imply integer squareness.

For A, odd-prime root lifting gives an actual compatible3-adic root sequence for every entry. This is not a coherent sequence of ordinary bounded integer roots. The archimedean root-height condition is exactly the missing coordinate. This statement does not assert an analogous obstruction when all primes and real places are jointly imposed on a single fixed rational number.

## 3. Every nonsingular odd-prime square grid has a regular p³ lift tree

Use the integral magic parametrization

M(c,a,b)=(a,3c-a-b,b,c-a+b,c,c+a-b,2c-b,a+b-c,2c-a).

Let p be odd. Suppose a magic square modulo p has all nine entries nonzero squares. Fix one choice of each root modulo p.

**Theorem 2.** For each k>=1, exactly p^(3(k-1)) root-vector solutions modulo p^k reduce to that fixed root vector and have square entries satisfying the magic equations. Every such solution has exactly p³ children at the next level. If the nine initial entries are distinct modulo p, they remain distinct at every level.

Proof. Choose the three pivot entries c,a,b arbitrarily lifting their fixed residues: there are p^(3(k-1)) choices modulo p^k. The other six entries are uniquely M(c,a,b). Every entry remains a unit square modulo p. Every unit entry has exactly one square root modulo p^k reducing to its fixed selected root modulo p: at each step the derivative2r is a unit and determines one correction digit. This yields existence and uniqueness. Conversely any solution yields precisely its three pivot entries and these unique roots. Each parameter triple has p³ lifts. Distinctness already detected modulo p cannot disappear at a higher modulus.

This is also the Jacobian calculation: the six primitive linear magic residuals have full rank6, and differentiating their square-root equations multiplies columns by the invertible2r_i. Thus the modular root variety is smooth of relative dimension3 at such a solution. The elementary parametrization proof above does not require algebraic geometry.

For p=59, the fixed integer grid (28,3,17;5,16,27;15,29,4) from the independently checked weak no-wrap example has nonzero square residues. The certificate verifies101 independent parameter-lift samples and explicit compatible roots through59^9. This is a counterexample to interpreting an endlessly extensible local root tree as a bounded integer-root witness, not to the correctness of Hensel lifting.

## 4. Search implications and connection contract

Source: a modular magic square with unit roots at an odd prime.
Target: its inverse-limit root-lift tree.
Map: lift three pivot entries and uniquely lift nine chosen roots.
Preserved: magic equations, unit-square status, initial labels, and any distinctness already seen modulo p.
Lost if only the tree is retained: integer height, global rational/integer root status, and zero-carry equality of the actual least-root squares.
Restoring sidecar: an archimedean height bound plus the exact carry syndrome from the sibling carry-lattice packet.
Prediction tested: p³ descendants per layer, not a pruning phenomenon. Therefore merely deepening a nonsingular single-prime Hensel search cannot by itself eliminate these local branches.

The p³ branching is a genuine dimension/lift law, not a correspondence with a ternary/Fibonacci fractal or tournament count. Additional maps would be needed for those claims.

## 5. Verification

Run `python3 -B verify.py` and `python3 -O -B verify.py`. The standard-library checker independently enumerates all magic square-residue arrays by their three pivot entries at moduli2,3,4,8,9,16,27,32,64,128; it checks the exact witnesses at five levels of each prime-power family, their eight ordinary line sums, and odd-prime unique root lifting. All load-bearing checks use explicit exceptions, not optimization-removable assertions. Complete finite census is only a control for the all-k proofs.

Context: the user supplied Cain's paper https://arxiv.org/html/1908.03236v2. The particular ring thresholds and elementary lifts here are rederived with full proofs; no assertion of first discovery. The usual linear magic-square parametrization is classical. The integer-square existence problem remains open in this packet.
