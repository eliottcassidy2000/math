# Source-paid residue certificates and the finite-sidecar barrier

**PROVED, scoped:** the arithmetic-cylinder verifier, the finite-sidecar state lower bound, and the residue-polynomial rank obstruction below. **FINITE-EXACT:** the atlas, its coverage, and the bounded first-descent experiment. **OPEN:** source payment for every positive odd integer. No Collatz solution or literature-priority claim.

## Findings

A compact, replayable atlas contains **751 exact valuation words**, each proving genuine descent below the immutable original source on an infinite arithmetic progression. Together they cover **29,864 of 32,768 odd residue classes modulo 65,536**. The remaining 2,904 classes are *unresolved by this fixed precision*, not divergent or exceptional classes. Unlike “eventually a non-growing return,” each accepted certificate pays the original source, including its exact carry.

The atlas uses the inherited guarded affine compiler. Its contribution here is a reproducible exhaustive realization and independent symbolic check, not a new all-source theorem. The sharper negative conclusion concerns the state interface: a rank that retains an arbitrary immutable source but otherwise only S finite states cannot strictly decrease through all failed returns. For sources below N, such an interface needs at least floor(log2(N+1)) distinct states. This does **not** exclude unbounded current-integer, carry, credit, or clock coordinates.

## 1. Inheritance actually inspected

* Session104 `descent/POSITIVE_ATOMS_AND_INTEGER_DESCENT.md`: local integer rank; exact original-source payment; paid-reset compiler; one-return time bound.
* Session100 `controllers/RETURN_CONTROLLER_TRAP.md`: exact (3,3) return termination, queue trap, and 5x+1 hostile.
* Session101 `continuations/CONTINUATION_STATE_AND_SYMMETRY.md`: ordered affine carry and its lossless word encoding.
* Session105 `README.md`: retaining native alignment and guards is necessary; no tournament-to-Collatz convergence transfer is claimed.
* Older local `math/05-knowledge/results/collatz_adaptive_return_budget_20261013a.md`, §§2,4: original-source budget, exact return cylinders, native-lattice repair. Its compiler already proves infinite guarded families; this atlas is a simpler first-descent specialization.
* Older local `math/05-knowledge/results/collatz_composed_weight_floor_20261011.md`: expanding actual-word hostile and source/flow distinction. The all-one family below is inherited, not rediscovered.
* Incoming `history/ATOMIC_INCOMING.md` at commit `201fa96604d51171909318cc482220e00df54c38`, and `history/SEVEN_INCOMING.md` at `63ef3ca1172665f7100d06d405a3aa34c42efa57`: five growing complete returns, fixed original-source atoms, and the lower-74 exclusion. Their finite exceptional-language bound is not replaced or improved here. This packet does not rerun that much larger search.

## 2. Exact paid-cylinder theorem

Let w=(a_1,...,a_J) be a nonempty positive valuation word, and recursively set

    (A,B,C)=(1,1,0),
    (A,B,C) <- (3A, 2^a B, 3C+B).

Thus A=3^J, B=2^D, D=sum a_i, and the endpoint is (An+C)/B.
The exact word guard is

    n = r mod 2B,  r=(B-C)A^(-1) mod 2B.

The final oddness bit is indispensable. Here is a direct proof that it authenticates the whole word. At a proper prefix of length j with carry C_j and division clock D_j, concatenation gives An+C=3^(J-j)(3^j n+C_j)+2^D_j C_suffix. The nonempty suffix carry is odd. Since the left side is divisible by 2^(D_j+1), the prefix numerator has exact valuation D_j. At j=J the guard itself says the quotient is odd. Thus every intermediate quotient is an odd integer; consecutive prefix identities give precisely the prescribed accelerated steps. The symbolic replay below checks the same fact independently.

Let n_* be the smallest n>1 in the guarded positive odd progression (r itself unless r=1, when n_*=1+2B). If

    B>A and (B-A)n_*>C,                                  (1)

then every n=n_*+2Bt, t>=0, satisfies T_3^J(n)<n. Indeed its payment is

    (B-A)n-C = [(B-A)n_*-C] + 2B(B-A)t > 0.

There is no slope-only inference: the constant term in brackets must be strictly positive. Equality is rejected.

The independent symbolic checker propagates a source x=u+vt, initially u=n_* and v=2B. On each letter a it checks that 3u+1 and 3v are divisible by 2^a, then replaces them by their quotients. It requires the new u odd and v even, proving the exact valuation for every t. At the end u<n_* and v<2B establish strict source descent. Three literal native replays are additional controls, not the justification of the infinite family.

## 3. Exhaustive finite atlas construction

For each odd residue r modulo 2^M, choose its smallest source n>1. Follow its true accelerated trajectory while a candidate word has D+1<=M. All sources in that class then share the exact prefix. Stop as soon as (1) holds and deduplicate its exact word cylinder. Otherwise mark the class unresolved when the next exact valuation needs more bits. No assumption about the eventual trajectory is made.

The resulting first-payment cylinders are pairwise disjoint. A prefix with positive slope difference has already paid at the smallest member of its cylinder iff it pays throughout that cylinder; thus a later certificate cannot overlap it. The verifier additionally enumerates coverage and checks disjointness and completeness of the covered/unresolved partition.

| Modulus bits M | Certificates | Covered odd classes | Total odd classes | Unresolved | Maximum accepted odd clock |
|---:|---:|---:|---:|---:|---:|
|8|19|95|128|33|4|
|12|99|1,766|2,048|282|6|
|16|751|29,864|32,768|2,904|9|

Coverage is an exact natural-density fraction among positive odd integers for these finite unions. It is **not** fixed-source atomic probability. Unresolved classes include sources needing more exact valuation precision, not just long first descents. For example a very large final valuation may pay immediately but be unavailable to this exact-prefix precision rule.

Separate bounded orbit checking establishes first descent for all **50,000 odd sources 3<=n<=100,001**. Maximum observed first-descent odd clock is85, first attained at35,655. A 10,000-step cap is explicit; no source hits it. This is an intentionally modest diagnostic, not a new verification record. Strong induction using these individual descents proves convergence for those bounded inputs, but says nothing about larger inputs not covered by a grounded dependency proof.

## 4. A sharp barrier for finite source-only sidecars

For every k>=1 set n_k=2^(k+1)-1. Its first k exact valuations are all1, and

    T_3^j(n_k)=3^j 2^(k+1-j)-1 > n_k,  1<=j<=k.          (2)

Proof: induction gives the formula, and before the next step its numerator has exactly one factor of two; the ratio of shifted endpoints to n_k+1 is (3/2)^j. Every letter1 is a leaf of the (3,3) controller, so each is a completed singleton return. All first k resets remain unpaid.

**Finite-sidecar theorem.** Suppose a putative source-preserving certificate keeps the original integer n0 immutable and a state q in a finite set Q, with |Q|=S. Suppose its rank R(n0,q), in any strict order, must decrease on every completed return until original-source descent or genuine termination. Then it cannot work for every positive odd source.

Take k=S in (2). The S return transitions occur while no original-source payment or terminal1 is available, requiring S+1 distinct ranked states. Two q values repeat by pigeonhole; their immutable source is the same, so their rank is the same, contradicting strict decrease. The rank may depend arbitrarily on n0; no continuity, computability, or polynomial assumption is used.

For sources n0<=N the same proof requires S>=floor(log2(N+1)) whenever N>=3. This is a necessary state-count bound, not a sufficient construction. If a native-step rank is used instead, the identical proof applies. Arbitrary fixed residue sidecars, clocks modulo a fixed period, and finitely many guard modes remain finite and are subject to it.

**Scope matters:** retaining the full current integer, unbounded carry, exact division clock, or unbounded debt is not a finite sidecar. The theorem does not rule out finite *descriptions* of unbounded ranks. Nor does it rule out a variable-length macro that independently proves payment before its sole ranked edge. Such a macro still owes a totality proof. This is a precise interface lower bound, not a general impossibility result for finite-state algorithms.

## 5. A second scoped no-go: endpoint residue-polynomial ranks

Fix M>=1. There is no nonnegative rank on all positive odd integers of the form

    R(n)=P_(n mod 2^M)(n),

with each P_r a real polynomial, that strictly decreases on every active accelerated Collatz step n>1.

Consider arbitrarily large n=-1 mod 2^(M+1). These have exact a=1; the endpoint y=(3n+1)/2 is larger than n and has the same residue -1 modulo 2^M. Write P=P_(2^M-1). Nonnegativity on its unbounded positive progression implies either P is constant, or its highest nonzero coefficient is positive. A constant cannot decrease. In the latter case P is eventually strictly increasing, so P(y)>P(n) for sufficiently large members, again impossible.

The same argument covers any eventually nondecreasing function on that residue progression, not just polynomials. It does not cover arbitrary oscillatory ranks, source-dependent two-variable polynomials, or a mode that changes with unbounded history. The obstruction is mathematical; sample checks are not its proof.

## 6. Hostiles and reproduction

For5x+1, word(1,1,5) has (A,B,C)=(125,128,39) and guarded source13. Its payment is (128-125)13-39=0; the certificate rejects it although its slope is below1. The rotated controller word(5,1,1) at83 likewise returns to83 with zero payment. Both native exact guards are replayed.

The verifier checks all-one witnesses through k=256 and rejects empty words, zero letters, noncontracting words, the unpaid5-cycle, and corrupt carry data. All checks use explicit exceptions, not Python `assert`, and run unchanged with optimization:

    python3 check_descent_transfer.py --output checks.normal.json --atlas atlas16.json
    python3 -O check_descent_transfer.py --output checks.optimized.json
    cmp checks.normal.json checks.optimized.json

Both modes pass with identical receipts. `atlas16.json` is generated evidence, while the symbolic inequalities above are the infinite-domain proof. The open next step is not another unmarked local return bound: it is an unbounded but controlled source-payment mechanism, or a genuinely total source-paid macro family.
