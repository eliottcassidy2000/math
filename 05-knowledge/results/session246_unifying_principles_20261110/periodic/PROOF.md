# A common flag for switched apex tapes and polynomial periodic laws

**Status: PROVED universal upper bounds and exact physical mixed-period examples; independently audited.** The bounds are physically sharp at every positive eigenvalue for k=1,2. All-k sharpness is not claimed. No all-target spectrum coverage follows.

## 1. Inheritance and exact setting

Fix any tournament on k labeled apexes and a transitive chain of tape vertices. A letter a is a k-bit mask: bit i=1 means apex i beats the new chain vertex. The physical graph is T_A(w). The positive permutation/gap automaton is inherited from Sessions215/216, restated and audited in Session221 `automata/PROOF.md`. Its dimension is k! 2·3^k. Session245 `forest/PROOF.md` proves that a single letter is a sum of commuting idempotents, giving square-free annihilator product_(j=0)^(k+1)(X-j). We extend its constant-letter statement to products of different letters.

Overlap checks: THM-4878, `01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md`, already supplies finite-state rationality and Hankel calculations for graded spine grammars. THM-3181, `01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md`, already gives polynomial-exponential repeated-*join* laws for path-colour jets. Neither is the common-flag statement for arbitrary fixed-apex tape letters below. We do not claim novelty for finite-state recurrences or polynomial-exponential sequences in general.

The root permutation p specifies k+1 chain gaps. A nonfinal gap has states e0 (empty), e1 (occupied, last vertex cannot enter its right apex), e2 (occupied, last vertex can enter its right apex). The final gap has e0,e2. The initial vector is all-empty. A terminal state is accepted precisely when every nonempty nonfinal gap has status e2 and every empty internal gap has the required arc between consecutive apexes. Summing all permutations counts each Hamilton path once, by its unique apex order and chain-gap assignment.

A letter operator is A_a=sum_g E_(a,g), where E assigns a vertex to gap g, subject only to the empty-gap left-apex eligibility. The sum records multiplicities when distinct assignments coincide. Composition here acts on column vectors, so a chronological word a1...al has matrix M=A_al...A_a1.

## 2. A flag shared by every letter

For each nonfinal gap use the ordered basis

    b0=e1-e2,  b1=e2,  b2=e0.

Let b be 1 if the left apex beats the letter (always 1 in the first gap), and let delta be 1 if the letter fails to beat the right apex. The assignment map satisfies

    E b0=0,
    E b1=b1+delta b0,
    E b2=b(b1+delta b0).

Its matrix is triangular with diagonal (0,1,0), independently of the letter. In the final gap use (b0=e2,b1=e0); its assignment matrix has diagonal (1,0) and sends b1 to b b0. This proves the parent session's candidate common flag exactly, including disallowed empty-gap assignments.

Use the tensor-product basis for all gaps. Assign rank r=sum of the basis indices. Every nonzero off-diagonal entry of every A_a strictly decreases r. The diagonal label is

    j = number of nonfinal indices equal to 1
        + indicator(final index equals 0).

It belongs to {0,...,k+1} and is independent of a. Therefore a nonempty length-l block M has diagonal j^l, with the same strict-rank triangularity. In particular its only eigenvalues are 0,1,2^l,...,(k+1)^l. The result holds separately for each apex permutation and hence for their direct sum. Actual apex orientations affect acceptance, not this flag.

**Important distinction:** the individual A_a are diagonalizable, but a product of different A_a need not be. Sharing a triangular flag does not mean sharing an eigenbasis.

## 3. Universal Jordan bounds

**Theorem.** Let l>=1 and M be any length-l physical tape block. Then

    M^(k+1) product_(j=1)^(k+1) (M-j^l I)^(k-j+2) = 0.       (1)

Thus for any fixed physical prefix u and suffix v,

    H(T_A(u W^m v)) = sum_(j=1)^(k+1) P_j(m) j^(lm),
    degree P_j <= k-j+1,                                    (2)

for every m>=k+1. There can be a zero-eigenvalue transient before that index. Coefficients are rational. The universal recurrence order is at most (k+1)(k+4)/2, independent of period length and far smaller than the k! 2·3^k state count; this is an upper bound, not a minimality assertion. For k=0 the assertion remains valid (the actual graph is transitive and its count is always 1).

**Proof.** The filtration by tensor rank is invariant, and its associated graded action is diagonal. For j>=1, ranks occupied by basis states with label j range from j-1 to 2k-j+1. For j=0 they range from 1 to 2k+1. Moreover r+j is always odd: each nonfinal index contributes an even number to rank plus label, while the final index contributes 1. Thus at fixed j the allowed ranks are separated by at least two.

To justify using these counts after eliminating other eigenvalues, take the generalized eigenspace G_j of M at lambda=j^l. Its intersection with each rank-filtered subspace is invariant. The associated graded of G_j injects into the associated graded of the ambient space. On that ambient grade M is diagonal, so only label-j basis components can support a generalized lambda-vector. Consequently the filtration of G_j has at most k-j+2 nonzero grades when j>=1, and at most k+1 when j=0. On each such grade M-lambda I acts as zero. Its indicated power is therefore zero on G_j. Distinct j give distinct eigenvalues because l>=1, proving (1).

The primary decomposition and the binomial formula for powers of lambda I plus a nilpotent matrix give (2). The zero primary component dies after k+1 powers. Prefixes and suffixes alter the scalar coefficients but cannot increase any degree. QED.

The numerical upper bound is governed by the number of same-label grades, not the huge state-space dimension. We have not proved that every bound is attained for every k. The proof provides a universal annihilator, not necessarily a minimal one. For a constant block Session245 gives the stronger square-free law and all polynomial degrees vanish.

## 4. The smallest mixed physical example

Take one apex and W=(1,0), meaning an apex-beaten vertex followed by an apex-beating vertex. The actual graph T(W^m) has 2m+1 vertices. Exactly, for all m>=0,

    H(T(W^m)) = (4·4^m+6m+5)/9.                           (3)

The nonzero linear term is unavoidable: distinct exponential-polynomial modes are linearly independent, or simply apply the shift operators to isolate the base-1 component. The counts begin 1,3,9,31,117,459,1825. The six-dimensional block has nullities 2,3,3 for successive powers of M-I, so it has a size-two Jordan block at 1.

Formula (3) follows from (1), whose recurrence has degree 5 at k=1,l=2, and five exact initial values. Both sides obey that recurrence; the five values are independently obtained by the direct positive gap count. The standalone verifier also checks these small graph counts by endpoint-subset dynamic programming. This is a finite certificate for an all-m identity, not extrapolation from a fitted curve.

Minimality is precise: k=0 gives only transitive tournaments, so no nonconstant polynomial factor; period length one has the inherited diagonalizable single-letter transfer for every k. Thus k=1 and period length two are both the smallest possible parameters for a mixed-period polynomial term. No claim is made about the fewest total vertices needed to distinguish arbitrary competing fitted formulas.

## 5. Sharp bounds at two apexes

Take two apexes with arc 0->1 and W=(3,2,1). Its masks mean respectively: both apexes beat the vertex; only apex 1 beats it; only apex 0 beats it. For all m>=0,

    H(T(W^m)) = P_1(m) + P_2(m) 8^m + P_3 27^m,

where

    P_1(m)=236833/811538 + (2650/4459)m + (13/49)m^2,
    P_2(m)=426044/866761 + (1684/6517)m,
    P_3=26433/122018.

In particular the maximal allowed degrees 2,1,0 all survive in the physical scalar readout. Counts begin 1,13,225,4917,121381,3166885,84462961. For one apex-permutation block the nullities of (M-I)^r are 5,7,8,8, certifying a size-three block at 1. For the complete two-permutation matrix they double.

The universal annihilator has degree 9. Nine initial counts prove the displayed identity, and the verifier checks further counts as additional controls. Its exact rational linear solve recovers the coefficients rather than injecting them. There is no assumption that an ambient Jordan block must survive acceptance: it genuinely does here. As a hostile, a different mixed block can lose its highest scalar coefficient; triangularity alone never guarantees noncancellation.

For k=1 and k=2 the examples prove physical sharpness of every positive-eigenvalue degree bound. General-k sharpness, and sharp zero-transient lengths, remain unclaimed.

## 6. Consequence, limits, and reproduction

The common flag makes every periodic fixed-apex tape exactly analyzable using a bounded collection of integer exponential bases and low-degree polynomials. It strengthens the constant-run calculus without pretending that mixed letters commute. It also provides a cheap universal recurrence certificate for symbolic physical families. It does not establish that the union of these families realizes every allowed odd integer, or that fixed-prefix/suffix parameters can be freely chosen to solve a target equation.

The connection contract is explicit: physical tape block -> tensor-flag matrix -> primary decomposition; preserved quantity is every fixed-context Hamilton path readout. The triangular diagonal forgets extension/Jordan data; the rank filtration and actual prefix/suffix are the missing sidecars. The minimal mixed example refutes the tempting stronger claim that all periodic tape families are pure exponential sums.

Run `python3 -I -B verify.py` and `python3 -I -B -O verify.py`; compare their output to RECEIPT.json. `--fail` must exit nonzero. The standalone verifier checks all length-1..3 blocks for k<=2 and all apex permutations (185 block checks), the common flag for all 11 letter/permutation cases, and independent subset-DP counts for all tape words of length <=4 at k<=2 (377 graphs). It verifies exact scalar identities with more initial values than the proved recurrence order, including all transient indices. No new numerical spectrum-coverage cutoff is computed.
