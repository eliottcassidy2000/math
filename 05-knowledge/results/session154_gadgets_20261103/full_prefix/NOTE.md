# Mixed-prefix orbit removes the 7 × 13² phase obstruction

Status: exact finite-state computational theorem with executable closure certificate. This concerns residues of actual strong-tournament families, not exact integer coverage or primality. Matrices, source words, seed and parity identity are inherited from `golden_inverse_defects_20261103/TRANSPORT.md` at commit `f6d5304849e50b29df6e368239e92fffa92e57e3`.

Let m=1183=7·13². Write q for q2 and let b be the triangle. On the integer coordinates (Y,X,g), the literal U seed is (25,34,6), with h=Y+g. The two prefix actions, applied one at a time, are

    q: (Y,X,g) -> (11Y+4X,18Y+7X,g),
    b: (Y,X,g) -> (5Y+2X,8Y+3X,-g).

We retain parity p of the number of b letters. In each parity layer the reached pairs (Y,X) modulo m are exactly

    S = {(Y,X): (Y,X) != (0,0) mod 7,
                  (Y,X) != (0,0) mod 13}.

There are |S|=(7²−1)(13⁴−13²)=1,362,816 states per parity, 2,725,632 total. The q and b matrices have determinants 5 and −1, both units modulo 7 and 13, hence preserve S. The seed belongs to S. The executable breadth-first search in `orbit.c` starts at the exact seed and stores each visited triple (parity,Y,X) once. After the queue is exhausted it independently scans all 2·1183² possible triples, asserting that membership in the computed orbit is exactly the displayed predicate, and asserts the resulting cardinality. The loops enumerate all possibilities with exact integer arithmetic; the maximum intermediate multiplication is below 32-bit limits. Run with assertions enabled:

    cc -O3 -std=c11 -Wall -Wextra -o orbit orbit.c
    ./orbit

The source and recorded transcript are in this directory. This is a finite exhaustive certificate, not a claimed conceptual proof of transitivity. A closed-form group-theoretic proof of this orbit remains a separate question.

For even parity g=+6. The two literal positive companions have H=3(Y+6)+2 and H=3(Y+6)+4 modulo m. For every prescribed Y modulo m, choose X=1; (Y,1) belongs to S. Since 3 is invertible modulo m, each companion separately attains every residue modulo 1183. In particular, each separately attains every residue 7j modulo 1183, j mod169. For any corresponding actual count H divisible by 7, this is exactly the statement H/7 ≡ j (mod169). It does not assert that H/7 is prime.

The q-only orbit has exact period 4368. It yields 934 distinct residues for each offset. For 3h+2 the omitted j modulo169 for H≡7j are

    1,11,24,27,37,50,53,63,66,76,79,89,102,105,115,118,128,141,144,154,157,167.

For 3h+4 they are

    2,5,15,18,28,31,41,54,57,67,70,80,93,96,106,109,119,122,132,145,148,158.

For a concrete recovered q-only omission, breadth-first search returns the forward transition sequence `qbqqbbbbbqq`, i.e. the literal prefix `qqbbbbbqqbq` when each action prepends a symbol. The seed count becomes h=11,639,406,485,921 and its +2 companion is H=34,918,219,457,765=7·4,988,317,065,395, with quotient ≡1 mod169. This is a genuine strong-tournament count from the inherited literal word construction, but the quotient ends in 5 and is composite. It is not a 7p witness.

An independent reverse-transition breadth-first computation reported the same 2,725,632-state primitive-pair orbit and full residue coverage. This certificate establishes the stated modulus 7·13² only: it makes no claim for higher 7- or 13-power moduli, arbitrary moduli, or an infinite profinite orbit.

SHA-256 (LF bytes):

    orbit.c   27b2c55a1082fbab37cc9043382a91b3b63d616989eb873572471c20e89d1c28
    orbit.out d57d21c6ee23a5a4ae46f148d5c940a33de4fbf20ebd677b7116bdfb9ef49625
