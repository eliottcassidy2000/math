# One extra cycle has an exponential path tax in the binary comb

Status: PROVED elementary family identity, independently reviewed. This is an addendum to `PROOF.md`, using exactly its comb definitions and the source identities in `SOURCES.json`. It is not an arc-reversal monotonicity claim for arbitrary tournaments.

Let rev_m(q)=sum_i epsilon_i 2^(m−1−i), with m bits including leading zeros, e0=epsilon_0 and et=epsilon_(m−1).

**Theorem.** For every m>=1, 0<=q<2^m and tau in {0,1},

    H(B+)−H(B−)=(4−tau)2^m−2−q−rev_m(q)+e0−et.         (1)

Consequently the sharp bounds are

    (2−tau)2^m <= H(B+)−H(B−) <= (4−tau)2^m−2,         (2)

with lower equality exactly at q=2^m−1 and upper equality exactly at q=0. Meanwhile the inherited cycle difference is C(B+)−C(B−)=1.

**Proof.** Only Hamiltonian words containing the two apices consecutively can change under the reversal. Let N_ab count otherwise valid words containing adjacent a,b, and let N_ba count those containing adjacent b,a, independently of whether that apex arc is currently present. Then H(B+)−H(B−)=N_ab−N_ba.

A word counted by N_ab has form I,a,b,J, with I,J complementary increasing chain runs. Ignore the condition on the last vertex of I for the moment. If J is empty there is one word. If J starts at chain index j, all smaller indices lie in I, while all larger indices have two choices. Summing over the allowed outneighbors of b gives

    1+2^(m+1)+(1−tau)2^m+sum_i(1−epsilon_i)2^(m−1−i).

The forbidden last vertex of I is x1. There are exactly two possible such runs: I={x1}, always giving an otherwise valid word, and I={x0,x1}, giving one exactly when epsilon_0=0. Subtracting 2−e0 gives

    N_ab=(4−tau)2^m−2−rev_m(q)+e0.                    (3)

A word counted by N_ba has form I,b,a,J. If J is empty, I is the full increasing chain and its last vertex must beat b; this contributes et. If J is nonempty, it starts at x1, because a has no other chain exit. Thus x0 belongs to I, and I must end at some bit vertex x(i+2) with epsilon_i=1. Each earlier bit vertex may be in I or J, and every later bit vertex must lie in J. The count is sum_i epsilon_i2^i=q. Therefore

    N_ba=q+et.                                        (4)

Equations (3)-(4) prove (1). The coefficient of epsilon_i in q+rev_m(q)−e0+et is

    w_i=2^i+2^(m−1−i)−[i=0]+[i=m−1],

which is strictly positive, including m=1 where the indicators cancel. Its sum over all i is 2^(m+1)−2. Hence its unique minimum and maximum occur at the zero and all-one words. Substitution proves (2). QED.

**Interpretation and limits.** The cycle code's easiest unit move is not a bounded-cost move in total path count: it has a tax at least 2^m, and padding makes that tax arbitrarily large even with q,tau fixed. This complements, but does not replace, the density obstruction. It does not exclude a different path-preserving surgery, deleting the sole-exit vertex, or a different construction family.

**Verification.** `python3 verify_arc_flip.py` and the same command with `-O` produce `ARC_FLIP_RESULT.json`. The checker independently enumerates all complementary chain-run choices for N_ab,N_ba through m=8, cross-checks the total-path difference by subset recursion through m=6, and checks all Hamiltonian permutations for m<=2. It uses explicit exceptions.
