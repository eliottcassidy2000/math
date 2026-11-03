# Three explicit feedback-two clone families

## Statement

Let T be the 7-vertex tournament whose out-neighbor bit masks, in vertex order 0,...,6, are

    [62, 60, 120, 16, 96, 8, 43].

For r in {0,3,5} and k >= 1, replace vertex r with k transitive clones having all r's orientations to the six other vertices. Write H_r(k) for its number of Hamilton paths. Then

    H_0(k) = 18*2^k - 1,
    H_3(k) = 19*2^k - 3,
    H_5(k) = 12*2^k + 11.

Each tournament is strong and has feedback vertex number exactly 2. These formulas give both prime-value candidates and 7p candidates, but neither infinitely many primes nor completeness is asserted. Each H_r(k) is divisible by 7 if and only if k is congruent to 1 mod 3. Explicit certified instances include H_3(4)=301=7*43 and H_3(7)=2429=7*347. They are not Hamilton-path counts of any tournament with feedback vertex number at most one. Also H_0(2)=71, H_3(3)=149, H_3(6)=1213, H_3(8)=4861, and H_5(2)=59 are prime counts outside that one-apex class (primality follows by trial division). Other examples H_0(4)=287=7*41 and H_5(4)=203=7*29 happen to be in the one-apex spectrum.

## Proof

The given T is a tournament, has H(T)=35, and is strong. Its every one-vertex deletion contains a directed triangle. Deleting {5,6} makes it transitive; deleting {2,4} also makes it transitive. These assertions are checked explicitly by the finite control. The clone blow-up retains an induced copy of T, so it has feedback number at least 2 and is strong (the retained strong copy reaches and is reached by each clone). Deleting {5,6} for r=0 or 3, or {2,4} for r=5, leaves a transitive tournament: cloning a vertex of a transitive tournament preserves transitivity. Hence the feedback number is exactly 2.

For each permutation of the six non-clone vertices, call a slot valid if a nonempty clone run can fit at it. An internal slot is a defect if the successive non-clone vertices point backwards. A path is possible only if each defect is valid and occupied. If the permutation has t valid slots and q defects, its number of extensions to k clones is

    sum_{j=0}^q (-1)^j binom(q,j) (t-j)^k.

Indeed each labeled clone independently chooses a valid slot; the transitive orientation fixes the clone order within each slot. Inclusion-exclusion enforces the q compulsory slots. Summing over the 6! permutations gives these exact multiplicities of (t,q), omitting impossible permutations:

    root 0: (1,0):3; (2,0):14; (2,1):4
    root 3: (1,1):5; (2,0):11; (2,1):8
    root 5: (1,0):7; (1,1):4; (2,0):12.

Thus, for example, H_3(k) = 5 + 11*2^k + 8*(2^k-1) = 19*2^k-3. The remaining formulas follow identically. This is the clone-slot identity already recorded in the S155 notes applied to a specific 35-path witness, not a claim of a new general clone theorem.

For the 7-divisibility claim, reduce the coefficients mod 7. In all three cases the defining congruence is 2^k=2 (mod 7), equivalent to k=1 (mod 3). The examples are established by direct evaluation and trial division.

For each claimed exclusion from the *entire* one-apex spectrum, enumerate its strong-core binary words by the S155 recurrence. The bound H >= 2^(m-2)+1 limits their length m. Source and sink padding cannot alter H, so the finite enumeration is exhaustive, not merely a bounded-order search. The run with threshold 4861 certifies every exclusion above. This relies on S155's core-reduction and recurrence proof and does not show that the all-k families avoid the one-apex spectrum.

## Reproduce

    python3 certificate.py

The standard-library certificate independently verifies the tournament, deletion obstructions, slot tallies, and subset-DP path counts for k=1,...,5, and performs the one-apex enumeration. The parent strongness claim may also be checked by directly enumerating outgoing edges of all nonempty proper vertex sets of T. The exhaustive slot count gives the formulas for all k, not only the k <= 5 controls. No assertion of literature priority is made.
