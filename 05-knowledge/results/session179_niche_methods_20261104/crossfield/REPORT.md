# Crossfield survey: transfer the certificate, not the analogy

Status: self-contained elementary proofs plus FINITE-EXACT controls; independent audit completed with two wording clarifications applied (spectral operator antecedent and positive-integer domain). No global Collatz or odd-perfect-number claim. Repository novelty is not asserted for the general methods. The concrete valuation separator and runnable crossfield controls are this packet's deliverables.

## Useful result first

The most useful transfer found is **prime-exponent program semantics -> exact divisor-sum columns -> separating linear functional**. It imports the repository's actual flow/Farkas certificate mechanism into a multiplicative arithmetic problem, without relying on visual resemblance or changing quantifiers.

For the 30 blocks `p^a`, `p in {3,5,7,11,13}`, `1 <= a <= 6`, the integer weights

    w(2)=-4, w(3)=-9, w(5)=-13, w(7)=-14,
    w(11)=-17, w(13)=-18, w(q)=256 for all other primes

have nonnegative pairing with every column

    c(p,a)[q] = v_q(1+p+...+p^a) - a*[q=p],

but pair to `-4` with the perfect-number target `e_2`. Consequently **no nonnegative multiset of these blocks**, even allowing repeated blocks with the same base prime, satisfies

    product sigma(p^a) = 2 product p^a.

In particular no odd perfect integer supported on these five primes, with each exponent at most six, exists. That restricted perfect-number consequence is deliberately small and not presented as new number theory. The stronger multiplicative-semigroup separation is the exact certificate being demonstrated. It can be applied to other catalogues; its usefulness at unbounded exponents or expanding support remains OPEN.

The other lane recovers (rather than rediscovers) a saturated Collatz max-cycle mechanism and supplies an exact rational affine-state compiler and hostile controls. It demonstrates why an averaged spectral certificate cannot be silently promoted to pointwise descent.

## Survey and inheritance: what actually transfers

All source paths below are relative to the repository `math/`. These are targeted repository-source findings, not a fresh literature survey. External paper claims in historical notes were not re-verified or imported as theorems.

1. **Logic / order topology:** `01-canon/theorems/THM-4623-collatz-functional-graph-order-topology-endomorphisms-and-koopman-torsion.md`. Reachability sobriety reformulates no divergence; it does not supply it. The usable positive mechanism is finite quotient nilpotence and exact parity conjugacy, not importing a lattice sobriety theorem to a poset lacking its hypotheses. The huge endomorphism monoid also warns that graph symmetry is not a faithful complexity surrogate.
2. **Realizability / category models:** `05-knowledge/results/collatz_five_papers_synthesis_20261006.md`, section 4; `05-knowledge/results/logic_lenses_jc2_collatz_hn_20261012.md`. The former itself labels the countable-real and measurement-independence comparisons ANALOGY. No map preserving integer termination is supplied. A useful future question would concern effective extraction of arithmetic certificates, not claiming a choice principle decides a fixed integer dynamical system. No such transfer was established here.
3. **Prime-exponent computation / FRACTRAN:** `05-knowledge/results/collatz_procgen_20260922_barrier_atlas.md`, Conway row (line 43 in the inspected version). Its source is secondary and it is an undecidability boundary for a general class, not for fixed Collatz. The operational feature retained here is elementary: multiplication/division becomes addition/subtraction of prime-valuation vectors. Program priorities and execution order are discarded, so Petri-style incidence constraints are necessary shadows, not a simulated execution. This is the source operation behind the divisor-sum certificate below.
4. **Farkas duals:** `01-canon/theorems/THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary.md`. The transferable object is an explicit verified rational dual separating a target from necessary columns. The theorem itself distinguishes certificate validity from solver-dependent metadata; this packet therefore replays a fixed integer vector, not a solver verdict or hash as proof.
5. **Valuation-sensitive obstruction:** `01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md`. Arbitrary nonlinear finite polynomial-valuation corrections cannot yield global stepwise Collatz logarithmic descent. This blocks an attractive but already closed transfer from the odd-perfect valuation functional back to a universal Collatz rank. It does not block finite arithmetic infeasibility certificates, whose quantifiers and targets differ.
6. **Spectral transfer / de Bruijn:** `01-canon/theorems/THM-4631-level-factorisation-of-the-compressed-twisted-transfer-operator-of-the-px-plus-1-maps-fast-modes-are-the-doubling-cycles-on-the-units-mod-p-i.md`. Weighted scalar spectral collapse and the prime-depth factors are already present. The packet's weighted checks are independently implemented inheritance controls, NOT a new theorem claim. Arithmetic primes can change fast-mode multiplicities while leaving the unresolved integer closure untouched.
7. **Max-plus / ergodic optimization:** `05-knowledge/results/procgen_tension_20260926_ranks_christoffel_duality.md`, Theorems R and D. The cycle/potential alternative and optimal rank defect already exist, including Collatz's expanding all-odd obstruction. Hence another proposed 'new max-plus solution' is not underused research. The artifact here repackages the mechanism as rational affine inequalities with an explicit additive-error core and checks a different small language, not as a new principle.
8. **Stochastic realization:** `01-canon/theorems/THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary.md`. Every finite terminal law admits a prefix Markov realization. Therefore merely constructing a chain with desired terminal statistics cannot be an arithmetic bridge. The exact compatibility with actual arithmetic transitions must remain in the certificate.
9. **Monoid framing:** `01-canon/theorems/THM-1885-the-half-and-plus-one-monoid-is-baumslag-solitar.md`. Its explicit affine maps are useful; broad 'all dyadic phenomena are the same generator' prose does not define maps between mathematical predicates. Here the prime-valuation monoid map is actually written down and checked. No complexity or group-theoretic claims from that file are imported.
10. **Perfect-number / tournament numerology:** `07-reflections/perfect-numbers-are-the-arc-counts-of-forcade-tournaments.md`. The even-perfect/triangular identity alone says nothing about odd perfect numbers or tournament spectrum. The needed replacement is the complete prime-power divisor-sum column, not an arc count. No spectral consequence is claimed from it.

The clearly less-used cell in this targeted survey is the **divisor-sum valuation cone plus certified separators**. The survey did not establish a repository-wide absence result or external novelty. The topology, generic realizability, weighted de Bruijn, and max-plus Collatz cells were already substantially occupied.

## Exact bridge A: divisor-sum incidence and separation

Let B be any finite collection of prime-power blocks `(p,a)`. Let Q contain 2, all base primes, and all prime divisors of every sigma(p^a). Define the integer matrix C with columns c(p,a) above.

**Proposition A1.** For nonnegative integers m_b,

    C m = e_2  iff  product_b sigma(p_b^a_b)^m_b
                       = 2 product_b p_b^(a_b*m_b).

Proof: take each q-adic valuation of the equality. Conversely all prime valuations agree, so the positive integers agree. No execution, ordering, or heuristic probability is involved.

For a genuine integer N, use at most one block per base prime and set that block's exponent to v_p(N). Multiplicativity of sigma then makes `C m=e_2` exactly sigma(N)=2N. Dropping the one-block condition is an explicitly stated relaxation; repeated same-prime blocks no longer compute sigma of their product, but an impossibility for the relaxation remains sound.

**Proposition A2.** If a rational row w has `w C >= 0` coordinatewise and `w e_2 < 0`, no nonnegative real solution Cx=e_2 exists, hence no nonnegative integer solution exists. Proof: `w Cx >=0`, contradicting `w e_2<0`.

This is the exact source-to-target bridge:
- source: additive prime-exponent transition effects and dual separation;
- target: multiplicative divisor-sum identity;
- map: prime valuations of each block ratio sigma(p^a)/p^a;
- preserved predicate: equality of positive rational products, exactly;
- discarded information: unique exponent block per prime when relaxed; transition order if viewed as a program;
- retained sidecar: every output prime and its full multiplicity;
- decisive test: verify every integer dot product and the negative target pairing.

**Concrete certificate.** The displayed w verifies all 30 inequalities. There are 23 coordinates, with primes

    2,3,5,7,11,13,17,19,29,31,37,43,61,71,157,1093,
    2801,3221,4733,19531,30941,45319,5229043.

Exactly five columns have zero pairing: `(3,2),(5,1),(7,3),(11,1),(13,1)`. The others have strictly positive pairing. The `7^3` column is important: sigma(7^3)=400=2^4*5^2. The tempting weaker weights `(-1,-2,-3,-3,-4,-4)` on `(2,3,5,7,11,13)` produce pairing `-1`, so they are NOT a certificate. The corrected vector gives zero there. This is a cheap hostile probe rather than a post hoc claim that a solver succeeded.

Positive controls use the blocks of 6 and 28. They have exact target e_2; at least one newly introduced base-2 column has negative pairing, so the separator does not falsely reject a real perfect number. All factorization, column construction, dot products, and controls use integers. A floating-point LP was used only to discover a candidate rational vector; it is unnecessary for replay and supplies no proof evidence.

**Boundary.** An unbounded OPN argument needs a uniform column inequality for all allowed primes/exponents, or a logically covering branching scheme whose finite leaves carry such certificates. Enlarging the catalogue can destroy the separator. No claim of uniformity or eventual completeness is made. This certificate is stronger than enumeration in its permitted block multiplicities, not in its prime/exponent range.

## Exact bridge B: affine-state descent with a finite core

Let a finite directed graph have edges e:u->v labelled by affine maps n->a_e n+b_e, with positive rational a_e and nonnegative rational b_e. Actual arithmetic trajectories in this statement have positive integer states n. Fix rational 0<rho<1.

**Proposition B1.** Positive h on vertices satisfying

    a_e h(v) <= rho h(u)

exists iff every directed cycle gamma satisfies `product_(e in gamma) a_e <= rho^length(gamma)`.

Necessity follows by multiplication around the cycle. For sufficiency set h(u) to the maximum over walks starting at u of `product(a_e/rho)`, including the empty walk. Removing any repeated-vertex cycle cannot reduce this product; hence the maximum occurs on a simple path and is finite. Appending an edge gives the inequalities. For rational labels this h is rational and a max-product Bellman-Ford recurrence computes it exactly in at most |V| updates.

Put `B=max_e b_e h(v)` (B=0 if the edge set is empty) and `V(u,n)=h(u)n`. Then any accepted arithmetic step satisfies `V(next)<=rho V(now)+B`. Thus `V` strictly decreases whenever `V>B/(1-rho)`. An infinite accepted integer orbit must reach that finite weighted core: if it never did, the recurrence brings it arbitrarily close to the core, whereas finitely many h values and bounded integer n give a positive discrete gap above it. This proves entry into a core or an exit from the accepted language, not Collatz convergence without checking the core.

**Positive control.** Vertices are four-bit histories; append a bit only if the new five-bit window contains at most three ones. Edges act by n/2 or (3n+1)/2 according to the appended bit. At rho=99/100 the exact compiler produces 16 vertex weights and verifies all 26 edges, with B=1250/1089 and weighted core 125000/1089.

A direct independent block calculation is sharper: among all 26 five-bit words with at most three ones,

    T^5(n) <= (27 n+76)/32 < n for n>=16.

The simultaneous coefficient maxima are attained at word 00111. If this window condition holds forever, repeated five-step descent eventually reaches an integer at most15. The bound alone does not show every Collatz orbit satisfies that window condition.

**Hostile control.** Prohibiting only 111 is insufficient: the allowed periodic word (110)^k has affine multiplier (9/8)^k. The exact compiler rejects rho=1 for this graph. For k=1..8, independent binary lifting constructs positive integer witnesses, all above a fixed large starting floor; the final witness starts at1694498811 and ends at4347718816 after24 steps. These are changing starts with increasing k, not one divergent orbit.

## Averaged spectrum versus pointwise certificate

For states w of length L, let K send w to each shifted word with appended bit with weight a_(w_1)/2. There is exactly one L-step path from w to each endpoint, and its product is `2^-L product_i a_(w_i)`. Thus K^L has identical entries along each row and rank one. The matrix K has only one nonzero eigenvalue, `(a_0+a_1)/2`; K^L has its Lth power as its only nonzero eigenvalue; this is inherited from THM-4631, replayed here at L=1..5 for three rational weight pairs.

For true Collatz multiplier s-moments, `a_0=2^-s`, `a_1=(3/2)^s`. At s=1/2 the averaged eigenvalue is `(1+sqrt(3))/(2sqrt(2))<1`, yet the all-odd loop has multiplier sqrt(3/2)>1. Hence an averaged positive-operator contraction does not imply any pathwise inequality of B1 on the unrestricted graph.

Every finite shortcut parity word occurs in a positive residue class modulo 2^k, with arbitrarily large representatives. This follows inductively: adding 2^k to an initial representative preserves its first k bits and flips the next bit because its k-step affine slope numerator is odd. Therefore a finite-word language accepting prefixes of every positive integer orbit cannot reject any finite word at all. This is a safety-language obstruction, not an impossibility theorem for arithmetic-aware liveness/Buchi conditions or height-dependent observers.

## Reproduction and limits

    python3 verify.py > results.json
    python3 -O verify.py > results.optimized.json
    cmp results.json results.optimized.json
    python3 valuation_dual.py > valuation_results.json

Both scripts ran successfully normally and under -O, with byte-identical outputs. No installs, network requests, repository modifications, theorem-ID reservations, or commits were made by this worker. The scripts use the standard library only. The finite parity bijection is independently checked through depth12. These checks support the specified formulas and certificates; they do not extend the proved finite arithmetic catalogue or settle any prize problem.
