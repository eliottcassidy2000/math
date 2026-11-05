# Independent audit: physical gap-one cores and their density obstruction

**PASS after one wording clarification.** The gap-one compiler physically realizes its targets, and its entire infinite family has counting coefficient strictly below 17/35. Therefore it misses at least N/70-O(1) odd targets through N. This is an obstruction to this particular family, not to the unrestricted C1/C2/C3 spectrum or to a union with other core families. Universal {7,21} coverage remains OPEN.

This is a separate audit of the gap-one continuation. It does not enlarge or replace the earlier audits of the canonical-deficit-core or synchronized-joint-inverse packets.

## Frozen inputs and required clarification

The original source directory is `session200-deficit-core-release/session200/gap-one-core-continuation/`. `INPUTS.json` records SHA-256 hashes of its six frozen files. In particular:

- README: `b12fbfa49f71b0417cfc89e4ee65eb8f64a324943d530cbacf4fd1d1f4f35e2f`.
- Density calculation: `c76fee7f48026a831a4016db0b3c5a8b5679379c5f767fe1f572d3c667de3c5a`.
- Author verifier: `91047e89190b76692ef81751d247398638cdd95da6acc36d783e56181d333328`.

README Section 1 says that the balanced-core slice is “exactly” the OJ family `(G,G,D)`. That can incorrectly suggest that it includes every balanced core from nontrivial sources. The triangle has core `(3,3,1)`, whereas G(v)>=5 for every nontrivial source v. Replace the sentence in the reviewed release with “The OJ construction therefore supplies the balanced-core subfamily ...”, retaining the displayed coordinates. The author accepted this clarification. No density formula or physical construction changes. This report approves the mathematical claims with that clarification; a release promotion should separately verify its exact changes and final hashes.

## Physical proof audit

For any source v=(D,E,U,d), direct operation gives `v.OJ=(D+2U,D,D+2U,D)`. The identity `JZ=OJ` is a chronological identity on all four coordinates. Thus the one-return construction indeed stays within its displayed balanced-core subfamily. This is not itself a completeness assertion about all balanced cores.

Every integer cover t>=2 has one zero-J word: even t must have last letter O and predecessor t/2; odd t>=3 must have last letter Z and predecessor (t+1)/2. Both predecessors are smaller, and the recursion terminates at cover 2. The output deletion is always one, and its endpoint is zero for even cover and one for odd cover. Updating D through Z proves the displayed recurrence for g, including the boundary g(1)=1 used in the even recursion. No fictitious cover-one physical source is required.

For odd t>=3, v=(g(t),1,t,1) is therefore physical. Its J-image is `(g+t-1,g-1,g+t,g)`. The positive endpoint is g-1>=2 and the retained triple is `(g+t-1,g+t,g)`, with gap exactly one. These are actual nonbalanced cores, not integer points accepted solely by a cone. The physical endpoint-height bound gives `t<=2g-3`, hence the claimed stronger lower bound `g>=(t+3)/2`.

Applying the canonical core compiler with R=2^(a+1), M=2^(k+1)-1 gives

    M[R(g+t)-g]+(g+t-1)+g
      =[M(R-1)+2]g+[MR+1]t-1.

The source recipe is the zero-J word at cover t, then J, then O^(a+1), then Z^k, and finally a closing J. It uses one J before the closing operation. This is the cost of this supplied recipe, not a claim that every target has minimum precursor cost one. At k=0 it simplifies to `(R+1)(g+t)-1`. Neither the additive offset nor the infinite set of available t proves arbitrary target landing.

## Density proof audit

Put A=M(R-1)+2 and B=MR+1. Both are positive, and A>=3. The bound g>=t/2+1 yields `n>=((A+2B)/2)t+A-1`. Therefore n<=N implies `t<=2N/(A+2B)`. The number of odd t>=3 in that interval is at most N/(A+2B), without a positive additive rounding error. Summing presentation counts over all nonnegative a,k bounds distinct targets, even when presentations collide.

The denominator is exactly `A+2B=(2^(k+1)-1)(3*2^(a+1)-1)+4`. Its lower bound `5*2^(a+k)` proves absolute convergence. Split the index set into the four a,k<2 terms, all a>=2, and the remaining a=0,1,k>=2. These are disjoint and exhaustive. The four denominators are 9,15,19,37. After discarding +4 in each tail:

- `sum_(k>=0)1/(2^(k+1)-1)<=5/3` and `sum_(a>=2)1/(3*2^(a+1)-1)<=2/23`, giving 10/69.
- `sum_(k>=2)1/(2^(k+1)-1)<=2/7`, giving `(2/7)(1/5+1/11)=32/385` for the other tail.

Each geometric comparison follows by subtracting the proposed bound: the three differences are respectively `2^(k-1)-1`, `2^(a-2)-1`, and `2^(k-2)-1`, all nonnegative in the stated ranges. Thus these are all-order inequalities, not inferred from a finite check.

The exact total is `27199346/56025585`; its positive gap below 17/35 is `13081/56025585`. There are N/2+O(1) odd integers through N, so this family misses at least N/70-O(1), equivalently asymptotically at least 1/35 of odd targets. The supplemental 6-by-6 partial sum plus tail 1/40 gives a valid stronger coefficient below 0.471, but is not needed for the stated conclusion.

## Independent implementation and literal graphs

`independent.py` imports no author modules. It uses matrix operations, a forward queue of O/Z covers rather than the author's recursive word compiler, and separately built adjacency matrices with subset dynamic programming for literal Hamilton paths and cycles.

- Every cover from 2 through 5,001 is generated exactly once: 5,000 zero-J profiles.
- All 122,500 odd-cover/exponent choices (odd cover 3..5,001 and a,k=0..6) satisfy the physical core and target formulas.
- Thirty-two literal C1 graphs (odd cover 3..17, a,k=0,1), of orders up to 11, have exactly the predicted Hamilton-path count and exactly one Hamilton cycle.
- The full family through target 20,001 is enumerated with a complete cover cutoff and exponent cutoffs. The declared receipt contains the distinct-target and presentation counts. This finite count is implementation evidence, not the asymptotic proof.
- Exact rational arithmetic independently confirms both density bounds and the strict gap below 17/35.
- Independent normal and optimized Python receipts are byte-identical. Both author replays are byte-identical to both frozen author receipts.

From this audit directory:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

The all-order conclusions depend on the displayed proof and the separately audited physical grammar/core compiler. Finite graph controls check the actual target rather than only profile equations. No unbounded source family beyond the stated gap-one construction is excluded.
