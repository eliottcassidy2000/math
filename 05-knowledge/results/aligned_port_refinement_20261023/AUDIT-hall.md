# Independent review of aligned actual-port refinement

23 October 2026. Reviewer: hall_bott_recovery. **PASS** for the precise scope
in PROOF.md; no whole-metagraph completeness claim follows.

- The integer identity C=A+B-(J-I)=F-F^T is correct on every matrix entry.
  Counts into common retained cells therefore determine the signed path
  counts. The unique +1/-1 divergence identifies the two actual endpoints;
  singleton columns peel the next owners from both ends. The middle vertex
  at odd order needs no extra round. Thus floor(n/2) synchronous rounds suffice.
- The only common equitable partition is discrete: starting from the uniform
  partition, every common equitable partition refines every computed round.
  The proof does not replace actual aligned vertices by separately sorted
  source/target colors.
- For a nonnegative doubly stochastic X commuting with A,B, it commutes with
  C and fixes C1. Its unique maximum/minimum force the endpoint rows. Column
  sums and nonnegativity force the corresponding columns. Commuting with C
  then propagates fixed basis vectors along the path. Thus X=I.
- Under the separate indegree/outdegree guards, D is exactly one endpoint
  path plus disjoint directed cycles. Nilpotence, det(I-D)=1, full traversal,
  and the (n-1)-step s-to-t entry are equivalent. The smallest detached-cycle
  counterfeit has five vertices because the path has at least two and an
  oriented tournament cycle at least three.
- Converse preserves the equitable partitions since A^T=J-I-A, reverses
  the agreement path, and exchanges the actual endpoints. The F2 binomial
  identity, dyadic product, and binary-subset rule all hold with the actual
  row-vector orientation used in the code.

An additional clock guard is worth retaining: for n>=2 the operator I+F
over F2 has exact order 2^ceil(log2 n). Indeed its 2^k-th power is I+F^(2^k),
which equals I exactly once 2^k>=n; its order is a power of two since one
such power is the identity. Consequently a snapshot forgets the integer
clock modulo this period. The explicit integer t in the theorem is a real
sidecar, not something recoverable from the finite state.

`audit_hall.py` independently uses sparse signed-neighbor records rather than
the two tournament matrices. It checks orders2..256: the C-only observer has
exactly floor(n/2) rounds, and chronological Pascal updates have the exact
stated return period. This does not prove actual A/B sharpness at every order:
the two full tournament records can split faster. The all-size statements
are established by the proofs above, not by this finite control.

Reproduce with `python -B audit_hall.py` or `python -O -B audit_hall.py`.

## Power-of-two modulus clock extension

**PASS, PROVED** for the incoming valuation-based extension. Let n>=2, s>=1,
let F be the actual directed n-vertex path matrix, and work modulo 2^s.
Then U=I+F and its start-row orbit both have exact order

    2^(s-1+ceil(log2 n)).

For 0<j<2^r,

    v2(binom(2^r,j)) = r-v2(j).

Indeed binom(2^r,j)=(2^r/j) binom(2^r-1,j-1), and the second factor is
odd: the row indexed by 2^r-1 has every binary digit equal to one. If
2^r<n, the coefficient of F^(2^r) in U^(2^r) is one, so that power cannot
be the identity. If 2^r>=n, the largest valuation of an index 1<=j<n is
floor(log2(n-1)). Thus every surviving coefficient vanishes modulo 2^s
exactly when

    r >= s+floor(log2(n-1)) = s-1+ceil(log2 n).

The path powers occupy distinct superdiagonals, so these coefficients cannot
cancel each other. A power of two does give the identity; the operator's
order therefore divides that power and is itself a power of two. The smallest
one passing the test is the claimed order. The start row sees every path
position separately, so it has the same exact return period.

The extended consumer independently performs chronological integer-residue
updates for every n=2..32 and s=1..4. All 124 first returns equal the stated
period; both normal and optimized execution retain the checks. This is an
arithmetic clock result on the authenticated path. It does not identify the
path with the native electrical mod-eight pairing or establish Bott
periodicity. No code from the incoming five-tap valuation implementation is
used in these chronological controls.
