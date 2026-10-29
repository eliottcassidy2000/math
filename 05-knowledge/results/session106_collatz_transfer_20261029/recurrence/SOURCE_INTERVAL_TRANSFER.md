# Source-interval transfer: an exact integer stopping-time sieve

**Status:** PROVED elementary all-word identities and interval closure; FINITE-EXACT implementation checks. **New-to-this-session assembly:** a closed interval carrier on the inherited division-exact cylinder/endpoint gauge for first-descent filtering, and its original-source generating-function recursion. No literature-priority, finite-state, speedup, or Collatz-convergence claim.

## Findings first

The useful transport from Session105 is **extend on the genuine carrier, then project**, not a claim that tournament spectra control Collatz. A Collatz prefix admits an integer state whose source and endpoint are affine in the *same* nonnegative integer parameter. Exact valuation extension is a congruence split of that parameter. Imposing every previous first-descent survival condition leaves just an interval, not an arbitrary subset. Thus an exact source-preserving first-descent sieve closes on six integers (with an optional infinite endpoint), even though their sizes and the number of states are unbounded.

This is a genuine interval compression of the surviving sources inside a word, not a compression of the word itself. It supplies a practical, auditable target for source-weight inequalities. Such an inequality still has to be proved; none follows merely from the recursion.

## 1. Inherited carrier and its integer gauge

Fix odd integer b>=3 and T_b(n)=(bn+1)/2^v2(bn+1) on positive odd integers. For a chronological word u of length k and sum H, inherited Sessions100/101 give

    A=b^k, B=2^H,
    T_b^k(n)=(A n+C)/B,
    n=r mod 2B,   r=(B-C) A^(-1) mod 2B,

where 1<=r<2B is odd. The exact source cylinder is n=r+2Bq, q>=0. Define

    M=2B,       t=(A r+C)/B.

Then t is a positive odd integer and the entire cylinder has

    source n=r+Mq,       endpoint y=t+2Aq.                (1)

The empty word has (r,M,A,t)=(1,2,1,1). The simultaneously recovered `history/ATOMIC_INCOMING.md`, §3 equation(3), already retains least source r and endpoint h and gives the arbitrary-block lift; (2) below is its single-native-letter specialization, not a new lift theorem. That incoming result is pinned to commit `201fa96604d51171909318cc482220e00df54c38`. Its five-return minimum search and the seven-window continuation (`history/SEVEN_INCOMING.md`, commit `63ef3ca1172665f7100d06d405a3aa34c42efa57`) are distinct inherited finite searches, not universal descent results. These coordinates are integers. They retain the source mark and the intercept: C=Mt/2-Ar. Given b, A and M also determine k,H; thus the gauge is not a lossless-word compression beyond the already inherited carry encoding.

**Inherited exact extension, native-letter specialization.** To append a>=1 let P=2^a and

    s = [P/2-(bt+1)/2] (bA)^(-1) mod P,   0<=s<P.

The sources with next exact valuation a are precisely those with

    q=s+Pq', q'>=0.

Their new state is

    r'=r+Ms, M'=MP, A'=bA,
    t'=(bt+1+2bAs)/P.                                  (2)

All divisions in (2) are exact and t' is odd. Moreover 1<=r'<M', so no recanonicalization shift is hidden.

**Proof.** The next numerator is bt+1+2bAq. Exact valuation a means it is congruent to P modulo 2P. Divide that congruence by 2; since bA is odd, its unique solution is q=s modulo P as stated. Substituting q=s+Pq' gives an odd constant quotient t' and endpoint t'+2bAq'. Source substitution gives r',M'. Positivity follows from positive t and s>=0. The bounds on r' follow from r<M and s<P. QED.

For each positive source exactly one a occurs. Hence the countably many child cylinders partition their parent, even though they have unequal moduli. Unlike a fresh random valuation, this split is tied to the original n.

## 2. Interval closure of the entire survival history

Keep a permitted source-parameter interval

    L<=q<=U,  L>=0,  U integer or infinity.

The extension in (2) sends it to

    L'=max(0,ceil((L-s)/P)),
    U'=floor((U-s)/P),                                  (3)

with U'=infinity when U=infinity; discard if L'>U'. In the new coordinates, survival of first descent at this endpoint means y'>=n, namely

    (2A'-M')q' >= r'-t'.                                (4)

Intersect (3) with this single half-line. If d=2A'-M' is positive, raise L' to ceil((r'-t')/d); if negative, lower U' to floor((r'-t')/d); if d=0, keep the entire interval exactly when t'>=r'. Discard empty intervals.

**Interval-closure theorem.** Starting from any interval of original positive odd sources, every exact word's sources that have avoided descent below their own initial value at all its prefixes are represented by one state (r,M,A,t,L,U) of the above form. Different words at the same length represent disjoint source sets.

**Proof.** Initially sources in an interval correspond to an interval of q. Congruence splitting maps an interval to (3). Every previous survival condition is already encoded in that interval. The only new condition is (4), whose intersection with an interval is an interval. Induction gives closure. Unique native valuations give disjointness. QED.

The distinction between immutable source and endpoint is essential: the test is y'>=r'+M'q', not y'>=the previous endpoint. A failed descent at one prefix can later become paid; only sources that have already achieved their first descent are removed.

For a finite source cap N, start with L=1,U=floor((N-1)/2), excluding the terminal source1. Each state has finite U. The largest endpoint is y_max=t+2AU, so all nonempty children satisfy

    1<=a<=floor(log2(by_max+1)).

This is an effective *finite* branch bound, computed by bit length. Thus the sieve halts for every prescribed depth and N without assuming Collatz. Each retained cell contains at least one source, so at a level its number is at most floor((N-1)/2). The elementary implementation makes no complexity advantage over directly iterating these sources.

## 3. Source-preserving weighted recurrence

For 0<z<1 attach to a cell its generating mass

    F(z)=sum_(q=L)^U z^(r+Mq)
        =z^(r+ML)/(1-z^M)                  if U=infinity,
        =z^(r+ML)(1-z^(M(U-L+1)))/(1-z^M) otherwise.    (5)

This is a nonnegative real series (or a formal power series). Before descent filtering, the child F sum equals the parent F **coefficient by coefficient**, because (2)--(3) partition its actual sources. After filtering the removed coefficients are exactly the sources whose first descent happens at that step. No endpoint-dependent resampling occurs and no cancellation is involved. Finite source caps produce polynomials, even though the displayed expression uses a rational function.

Let G_K be the sum of retained masses after K steps starting with all odd n>=3. Its coefficients are exactly 1[n has no first descent by K], so G_(K+1)<=G_K. Normalize to the full-support probability

    w_z(n)=(1-z^2)z^(n-3),  n=3,5,7,... .

The exact survival mass is S_K=(1-z^2)z^(-3)G_K. Each surviving n contributes its unchanged positive atom w_z(n). The inherited Session104 atomic argument therefore gives

    S_K ->0  iff every odd n>1 eventually descends below itself. (6)

For b=3 this is equivalent to Collatz termination by strong induction; the converse uses eventual arrival at1. Equation (6) is an **open target**, not an estimate established here. A uniform contraction S_(K+1)<=rho S_K for rho<1, or a proved source-selective tail modulus, would suffice, but the interval theorem supplies neither. Evaluating more finite levels cannot justify the limit. The exact recurrence makes the missing estimate explicit rather than replacing atomic weights by Haar cylinder sizes.

## 4. Quotients and hostile controls

* The clocks A,M alone discard t, which enters both the next exact congruence s and survival intercept r-t. Discarding r loses the initial-source comparison as well. A proposed coarser state must preserve the complete chosen continuation response (or give a proved weighted disintegration), not just a denominator or slope. These statements identify necessary obligations, not a minimality theorem for all possible algorithms.
* For b=5, word (1,1,5) has affine fixed point13 and native cycle13->33->83->13. Starting at immutable source13, every endpoint is >=13. It survives every repeated-word sieve level, so S_K>=w_z(13)>0. The implementation explicitly keeps13 for four repetitions. A contraction argument based only on this parameter-independent interval formalism would therefore be false.
* At b=3, inherited words (1,1,5),(1,3,3) have fixed points19/101,31/101. Append1 and **reclose**:37/35 and221/175. Both the state and code preserve the distinct carries; a denominator-only attachment recurrence fails. Reclosure is not a step of either native trajectory.
* The source-interval machinery includes all positive integer sources. It neither treats all rational periodic words as integer sources nor discards any finite prefix as impossible. Every unfiltered exact-word cylinder is nonempty and infinite, as required by Session100.

## 5. Verified scope, provenance, and reproduction

Directly inspected inputs:

* Session105 `metagraph/WEIGHTED_NATIVE_RECURSION.md`: actual-source extension before quotient, normalization and nonlumpability warning.
* Session105 `enumeration/INTEGER_PROJECTIVE.md`: exact integer gauges versus information compression; divisions must be proved exact.
* Session100 `bijections/VALUATION_CONJUGACY_AND_INTEGER_SIDECAR.md`: exact source cylinders and the hostile5-cycle; integer-source boundedness is not a finite-pattern filter.
* Session101 `continuations/CONTINUATION_STATE_AND_SYMMETRY.md`: affine carry, lossless decoding, 101->35/175 obstruction.
* Incoming `history/ATOMIC_INCOMING.md`, §3: source/endpoint least-residue lift and finite growing-return minima. The single-letter lift here is inherited; the new assembly is the exact interval filter for all-prefix original-source non-descent.
* Session104 `descent/POSITIVE_ATOMS_AND_INTEGER_DESCENT.md`: immutable-source debt and first-descent induction; atomic coverage is an inherited theorem, not a new convergence proof.

`check_recurrence.py` independently compares affine formulas with actual native iteration for all340 words of lengths1..4 on letters1..4, for b=3,5,7 (1020 states); tests exact cylinder sources and interval cuts; checks geometric source-mass conservation; replays the named hostile controls; and compares the interval sieve against a literal source-to-current-value dictionary on every odd3<=n<=4095 for18 steps at each multiplier. It uses no removable assertions, external libraries, floating-point arithmetic, or unbounded stopping-time search.

    python3 check_recurrence.py > checks.json
    python3 -O check_recurrence.py > checks_optimized.json
    cmp checks.json checks_optimized.json

Both receipts agree byte-for-byte. Tests validate the implementation on their finite domains; the proofs above establish the all-word statements. No global tail estimate or universal first-descent bound has been supplied.
