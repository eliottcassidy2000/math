# Independent review: exact marked path duality

Verdict: PASS. This review binds `../marked_path_duality/PROOF.md` at SHA256 `b92f9c0b7ddd6ba95b772af876e6ec3ef144b52e6765df97e1e05aca999c3e33`, verifier `fc97d44e0e2c8d46378ea360005796c27ca314f98e973f5cac14c8483ced7d85`, and paired producer receipts `47adec547640b26fe6a4f6288c1d08b927ef8b3ba5342df56e45cf7717df7e71`. No correction is required.

## All-order algebra

The matrix factorization uses the diagonal identity Z^2=0 correctly: `(I+Z)(I+ZB)-ZJ=I-ZB^t`. Its inverse applied to Z is N before the rank-one update. All denominators have constant term one, so they are units in the integer squarefree ring. The rank-one formula gives the stated transpose matrix and H=h/(1-h), hence exact reciprocity without dividing by an integer. Its row and column sums have the displayed signs. Antisymmetrizing the matrix formula gives `D/(1+H)=E(-z)-E(z)`, which is exactly minus twice the odd-support projection. No overlapping physical supports or independent endpoint clones enter these products.

The arbitrary-child quartic specialization is also exact. Independently inserting the two exterior singleton names u,v into the complete eight-term disjoint quartic gives the stated six-term source expression. The standalone `formal.py` reproduces this calculation with seven algebraically independent child-path symbols and squarefree u,v, then subtracts the transposed expression to obtain `(1-H)D+(2H+H^2)E` identically. Applying the proved marked identity gives equation (6). This is a universal symbolic check, not merely a test on small children.

The scope is correct: the marked duality holds for every tournament A, while the response specialization uses the particular two-singleton cyclic exterior. Neither non-annihilation for all A nor arbitrary-exterior positivity nor arbitrary-ambient sign transfer follows. The Paley example disproves a universal nonzero-mod8 claim for this chosen coefficient matrix only; it is not a full-response collision. The text makes no historical-priority claim for the inherited reciprocity method.

## Independent executable evidence

The auditor imports no producer algorithms. It enumerates all directed paths as vertex permutations, uses dictionary-based disjoint-support multiplication, and counts actual augmented graphs with a separate predecessor-removal Hamilton dynamic program.

- All 1,103 complete algebra controls agree: every labeled tournament through five names and four larger named controls. There are 33,423 independently verified coefficientwise reciprocity, row/column and marked identities, together with the full kernel-polynomial identities.
- All 820 literal fourth-difference controls and all 10,144 distinct actual graph counts agree.
- Every one of the 2,116 complete named Paley-bank coordinates and every retained full-mask reconstruction agrees. The contrast retains 420 nonzero cubic and 420 nonzero quartic coordinates, with all zeros checked.
- The Paley marked matrix has exactly the claimed off-diagonal values +8 and -8. The actual query (a,c)=(1,2) gives 16,185 versus 16,419, so the full response is visibly different even modulo four.
- Both audit programs pass in normal and optimized Python modes with byte-identical paired receipts. Checks are explicit exceptions, not removable assertions.

These controls support the independently reviewed all-order algebra; finite evidence is not used to assert an additional universal nonzero observer.
