# Exploratory positive-cancellation scout

**EXPLORATORY / NOT A NEW PUBLICATION PACKET.** Exact identity controls pass;
no broader constructive target subcase was obtained. The frozen factorization,
digital-graph, Hammersley and bounded-error proofs were left unchanged.

For k-bit A=sum_i epsilon_i2^i, B=2^k and R=rev_k(A), set

    Q(A)=sum_(i<j) epsilon_i epsilon_j
                    (2^(j-i-1)+1_(j=i+1)).

The actual index of its chronological source word is u=2B-R, and

    g(u)=2A+1-Q(A),
    Q(rev_k(A))=Q(A).

Thus the reversed word's actual index v=2B-A satisfies

    g(u)-2u=g(v)-2v.

For k>=1, complementation gives the exact further identity

    Q(B-1-A)=Q(A)+B-2-A-R+epsilon_0+epsilon_(k-1).

These follow from the weighted bit-pair formula and preservation of pair
separations under reversal. For complementation, the sum of weights incident
to position i is2^i+2^(k-1-i)-1_(i=0)-1_(i=k-1), and Q(B-1)=B-2.

The hoped-for cancellation fails in the displayed physical-form expansion:
reversal/complement pairings express both g coordinates as a linear term
minus the SAME Q(A). Every nonempty central form retains its positive odd
AB coefficient multiplying Q(A)^2. Empty retains -2Q(A). This observation
concerns these explicit substitutions only; it is not a general impossibility
theorem for more complicated digit constructors or other actual-source
identities.

The precise remaining positive condition for any fixed physical row is

    n=alpha(2A+1-Q(A))+beta(2B-rev_k(A))
         +gamma epsilon_(k-1)+delta,

with the last-bit term set to0 at k=0. A uniform constructor would need an
actual binary A hitting this INTEGER equality, across the allowable rows.
Slope/discrepancy control places nearby real values but does not supply that
exact weighted binary-quadratic hit.

A sparse backup was checked but deliberately not promoted: for the word
(ZO)^m, m>=0,

    v=(4^(m+1)+2)/3,
    g(v)=(4^(m+1)+6m+5)/9,
    Psi_J(2,v)=(28*4^m+6m+17)/9.

This supplies only an exponential-linear sequence of targets, not an interval,
a fixed congruence class, or a broader integer target subcase. It therefore
does not justify another publication packet.

Reproduce the bounded algebra controls with

    python3 session201-positive-cancellation-scout.py

The script reads the existing frozen twelve-entry menu. It does not scan
target coverage. Outcome: stop this branch here; preserve the positive
approximate-slope constructor and the exact carry-hitting obligation.
