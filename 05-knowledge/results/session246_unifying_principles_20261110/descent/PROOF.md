# Bounded-deletion dynamics and the binary shadow

**Status: PROVED elementary deductions from the inherited physical grammar; independently audited.** This is an infinite-dynamics theorem and a precise failure boundary for projective arguments. It does not prove Chan–Pak, a new coverage interval, or a missing-target predecessor.

## Inheritance, scope, and concept board

Use Session199 `source-chambers/PROOF.md` for the actual-source inequalities and complete inverse, `hidden-structure/PROOF.md` for the binary block compiler and the positive-density deletion-one gaps. Session244 `dyadic/PROOF.md` rules out the naive continuation-closed two-sign menu. Session245 `descent/PROOF.md` gives synchronized physical rays but not arbitrary-target selection. The chamber inverse is already known and is not claimed new here.

The concept board is: physical source; deletion growth; binary projective shadow; exact scalar readout; source certificate. The anchor is a lawful inverse; the niche is infinite physical trajectories; the wildcard is local-surjectivity versus global-lifting failure. Searches of the current source/canon routes found no statement of the quantitative shadow theorem below; this is a bounded novelty check, not a priority claim.

All trajectories below start at an actual nontrivial source, and all letters are literal graph extensions. Let its coordinates be (D,E,U,d), with

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

The seed is (1,0,2,1). We use the inherited inequalities 0<=E<=d, D>=d, U>=D and U>=2d. In particular d>=1. One can alternatively check their preservation directly under the three maps, together with the inherited stronger positive-cone inequalities.

## 1. A quantitative deletion-growth dichotomy

**Theorem 1.** In an infinite physical trajectory, d stays bounded if and only if there are finitely many J letters. If d_i is the deletion coordinate immediately after the i-th J, then

    d_(i+2) >= d_(i+1)+d_i.

Consequently d_i>=F_i for i>=1, where F_1=F_2=1. In particular a bound d<=M allows at most 2 floor(log_2 M)+2 root changes. This deliberately loose estimate is sufficient for finiteness.

**Proof.** O preserves D and d; Z preserves d and increases D by U-E>=d; J has the same D increment and replaces d by the old D>=d. Thus both D and d never decrease. Immediately after J number i+1, D is at least d_(i+1)+d_i. Until the next J it cannot decrease, so that next J gives the displayed inequality. Since d_1,d_2>=1 the Fibonacci bound follows. The recurrence also gives d_(i+2)>=2d_i, yielding the logarithmic estimate. Conversely after the last J every remaining letter preserves d. QED.

This is a growth statement about one supplied word, not a minimum-J complexity theorem about an unmarked graph. It does not infer that every graph with bounded deletion has a unique word.

## 2. Quantitative binary shadowing

After the last J, index an O/Z tail by n>=0 and write its constant deletion d. Put a=U_0-d>0, x_n=D_n/U_n, and epsilon_n=0 or 1 according as letter n is O or Z. Define

    y_0=x_0,   y_(n+1)=(y_n+epsilon_n)/2.

**Theorem 2.** For every n>=0,

    U_n >= d+2^n a,
    0<=x_n<=1,
    |x_n-y_n| <= n d/(2^n a).                         (1)

Moreover

    y_n = 2^(-n)x_0 + sum_(j=0)^(n-1) epsilon_j 2^(j-n). (2)

Thus x and the reversed binary-prefix reading y have the same subsequential limits.

**Proof.** O gives U'-d=2(U-d)+d, whereas Z gives U'-d=2(U-d). This proves the first bound. The second is the physical inequality D<=U. O gives x'=x/2 exactly. For Z, direct subtraction gives

    x'-(x+1)/2 = [d(D+U)-2EU]/[2U(2U-d)].

Since 0<=D<=U and 0<=E<=d, the numerator has absolute value at most 2dU. Therefore the one-step error is at most d/(2U-d), and hence at most d/(2^(n+1)a). Iterating a contraction of factor 1/2 gives n terms each at most d/(2^n a). Formula (2) is immediate induction; the bound tends to zero. QED.

The estimate compares normalized states, not actual integer counts. It is uniform over every binary tail and requires no genericity, randomness, or finite-order extrapolation.

**Constructive finite consequence.** For any t in [0,1] and any r>=1, choose A=min(2^r-1,floor(2^r t)) and use its r binary digits in least-significant-first order as epsilon_0,...,epsilon_(r-1). The resulting physical suffix satisfies

    |D_r/U_r-t| <= (1+r d/a)/2^r.

Indeed its affine shadow is (A+x_0)/2^r, within 2^(-r) of t. This is an explicit, uniform normalized-state controller. It chooses a lawful word, not an arbitrary cone point, but controls a ratio rather than an integer Hamilton count.

## 3. Exact consequences for infinite trajectories

**Corollary 3 (convergence rigidity).** On an eventually O/Z physical trajectory, D_n/U_n has a limit if and only if the letters are eventually constant. The limit is 0 for eventual O and 1 for eventual Z.

If x_n converges, the one-step error tends to zero, so epsilon_n=2x_(n+1)-x_n+o(1) converges to the same limit. A convergent zero/one sequence is eventually constant. The converse follows from the affine recursion and (1). In particular infinitely many O and infinitely many Z forbid convergence even though d/U tends to zero.

**Corollary 4 (periodic tails).** For a repeated length-l binary block epsilon_0,...,epsilon_(l-1), set A=sum epsilon_j 2^j. At block endpoints the limit is A/(2^l-1). Intermediate phases follow y'=(y+epsilon)/2. A mixed block produces a nonconstant rational cycle; its period can be a proper divisor of l. For instance repeated OZ has endpoint limit 2/3 and the other phase 1/3. This follows because one block maps y to (y+A)/2^l.

**Corollary 5 (dense shadow with positive-density arithmetic omissions).** There is a single physical d=1 trajectory whose normalized counts have subsequential-limit set [0,1]. Nevertheless the entire physical d=1 family, and therefore this trajectory's balanced counts G=D+2U, misses a positive proportion of odd targets.

For construction concatenate, for m=1,2,..., all binary words of length at most m in lexicographic order (within each length), and interpret bits as O/Z. Each finite binary word occurs arbitrarily late. To approach a chosen t in [0,1], choose successively longer binary prefixes of t, arrange that their reversals end at successively later positions, and use (2). Boundary t=1 uses prefixes of all ones. Bound (1) transports the limit to physical x. Since x lies in [0,1], the limit set is exactly that interval.

The arithmetic assertion is inherited Session199 hidden-structure Theorem 3: with c=sum_(k>=0)1/(2^(k+2)+1)<9/20, the full d=1 balanced image has at most cN+O(log N) members through N. Thus at least N/20-O(log N) odd targets are omitted. This is not a new proof of that density bound. The new conjunction gives an explicit physical infinite path with maximal continuum limit set inside a rigorously arithmetic-incomplete family.

## 4. Procedural analogy and what transfers

After the known inverse failed to yield an all-target selector, a reproducible sample used Python random.Random(246031).sample(sorted(git ls-tree theorem paths),4) at HEAD 0e140250d07c. Its four files were THM-3635-russell-cylinder-retained-curve-jet-plane-actual-rank-witness.md, THM-3261-projected-k3-z216-unique-gcd18-terminal-descent.md, THM-3817-quintic-r-repairs-of-nodal-carriers-have-critical-points.md, and THM-643-blue-black-line-parity-structure.md.

THM-3635 was the useful analogy: actual-ring local jet realizability still retains a global algebraization obligation. Here the source is physical word trajectories, the target is normalized real states, the map is (D,E,U,d)->D/U, and the preserved predicate is asymptotic binary shadowing. The map destroys integer scale, deletion and endpoint information, and exact attainment of G. The required sidecar is the original word/full integer state. Corollary 5 is the decisive test: even the strongest possible real limit-set coverage coexists with a proved positive-density integer omission. This is an analogy with a precise target-side theorem, not a claimed map between the Russell ring and tournaments, and no algebraic theorem is imported as a tournament proof.

## 5. Remaining frontier and verification

No scalar failure-preserving predecessor is obtained. The inherited balanced inverse still must select a physically realizable deletion coordinate. The result says that a projective limit/compactness strategy must retain an arithmetic sidecar; real density alone cannot supply that selector. A future infinite-J strategy escapes the fixed-deletion obstruction only by paying unbounded deletion complexity, at least with the growth proved above.

The standalone verify.py uses exact fractions, all O/Z/J words through length 7, all binary tails through length 8 from selected physical roots, independent subset-DP graph counts on all words through length 4, and periodic/dense-word controls. No range-coverage expansion is performed. Normal and optimized runs must match; --fail must fail.
