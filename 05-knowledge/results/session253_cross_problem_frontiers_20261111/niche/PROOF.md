# Session253 niche bridges: the first ramified Smith band and zero-carry fibers

Status: **PROVED IN THIS NOTE**, with small **FINITE-EXACT** new-claim controls. Not independently canon-audited; no theorem number or priority claim. Chan–Pak all-height coverage, general jet partitions, and Sun leastness/density remain **OPEN**.

## Findings first

1. An all-node/all-multiplicity jet theorem closes a genuine first-band frontier: for a uniform m-jet bank in one p-adic cluster of minimum depth e, the first positive Smith exponent is e+v_p(m). Its multiplicity is one, except at the shallow dyadic even-multiplicity boundary, where it is exactly two. This includes highly nested node sets and arbitrary ramified m, rather than assuming all pair distances equal.
2. Sun's universally full modular images admit an exact **minimum-integer-value** replacement. Each binomial atom has a finite first-period minimizer, and min-plus cyclic convolution gives the exact minimal integer in every residue fiber. The missing target coordinate is its nonnegative carry, not another support bit. This supplies a lossless, finite height test and a precise common carrier for the Chan–Pak costed-fiber program.
3. A universal forced-divisor invoice gives a collision-cap instrument for arbitrary two-cube targets, including nonprimitive representations. Unlike the restricted singleton decoder, it is stated directly from the target factorization.

## 1. Inheritance, sources, and concept board

Read from the clean supplied Session252 worktree. The missing September result files are tracked in Git and were recovered with `git show HEAD:path`, not treated as nonexistent because the checkout was sparse. Git was not mutated.

- Closest proved mechanism: `01-canon/theorems/THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall.md`, especially the first-positive-layer bound and residue-cluster CRT decomposition.
- Precision sidecar: `01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md`, cardinal Hermite denominators and dyadic unit dependence.
- Least-used companion: `05-knowledge/results/overnight9_20260906_jets_deuring.md`. Divided reflected solutions P,R and G=(P−R)/p satisfy a Wronskian identity that caps simultaneous top-jet cancellation at one digit. The exact largest loss is (3m−1)e−[H_p(a)=0], m=(p+1)/2.
- Hostile to retaining only metric and largest loss: `05-knowledge/results/overnight13_20260906_jets_p31_intermediate.md`: normalized a=3 and a=4 share metric, ordinary Deuring bit, and largest exponent47, but an intermediate determinantal ideal distinguishes them. This is inherited, not re-audited.
- Sun sources: THM-4027 `sun-two-four-six-eight-universal-modular-solubility`, THM-4036 `sun-2468-energy-and-support-exponent`, THM-4246 `sun-reflection-orbit-two-adic-tower-and-odd-cycle-firewall`.
- Two-cube sources: THM-463 `two-cube-representations-are-a-divisor-property-on-the-split-axis`; THM-3793 `inert-prime-sum-all-scale-two-cube-singleton`; THM-4018 `ramified-three-two-cube-singleton-extension`; THM-3825 `prime-colour-valuation-two-cube-decoder`.
- Chan–Pak source: `05-knowledge/results/session252_fivefold_bridges_20261111/SYNTHESIS.md` and `arithmetic/PROOF.md`: seven-coordinate compiler, K1 universal Fermat-product residue coverage, explicit missing small-output lift.
- Wildcard reconnaissance: `05-knowledge/results/overnight11_20260906_no3line_rowfreeze.md`. Its conditional row-label concentration is a productive sidecar principle, but no arithmetic map was found; no new no-three-line claim is made.

Concept board: complete Hasse banks / weighted binomial bands / divided companions / minimum-value residue fibers / square-discriminant divisor fibers / Chan–Pak bounded outputs. The operations used are clearing a complete jet block, rescaling, first associated-graded rank, replacing support by tropical minimum, and attaching valuation-plus-height invoices.

## 2. All-node uniform-jet first-band theorem

Fix a prime p, integers n>=2 and m>=1, and distinct integer nodes x_0,...,x_(n−1). Assume

    e=min_(i<j) v_p(x_i−x_j)>=1.

Let M be the nm by nm Hasse evaluation matrix on polynomials of degree<nm:

    M_((i,r),c)=binom(c,r)x_i^(c−r), 0<=r<m, 0<=c<nm.

Let its ordered p-Smith exponents be α_1,...,α_(nm).

**Theorem.** There are exactly m zero exponents. The smallest positive exponent is

    λ=e+v_p(m).

Its multiplicity is

    2, if p=2, e=1, and m is even;
    1, in every other case.

There is no bound on n, m, higher pair depths, or unit lifts. The exceptional multiplicity2 is possible already at two nodes and m=2. At e=0 this single-cluster statement does not apply; split into distinct residue clusters first.

### Proof

Translate x_0 to zero, an integral unimodular coefficient change, and write x_i=p^e h_i. All h_i are integral, and at least one h_i is a p-adic unit: otherwise every pair difference would have depth greater than e. The m rows at node0 are the first m coordinate rows. Clear these columns from all other rows. This leaves I_m and a residual matrix R with columns c=m+s and rows r=m−t, where s>=0 and 1<=t<=m:

    R_((i,t),s)=binom(m+s,m−t)(p^e h_i)^(s+t).

Put L=s+t. The binomial numerator is the product of the L consecutive integers from m−t+1 through m+s; it contains m. Hence

    v_p(R_((i,t),s)) >= eL+v_p(m)−v_p(L!).

For all L>=1, v_p(L!)<=L−1. Consequently every residual entry is divisible by p^λ. The entry t=1,s=0 at a unit h_i is m p^e h_i and has valuation exactly λ. This proves the smallest positive exponent. Since R vanishes mod p, the zero multiplicity is exactly m.

The number of residual Smith factors of exponent λ equals the rank over F_p of p^(−λ)R modulo p. If p is odd and L>=2, then L−1−v_p(L!)>=1. If e>=2 and L>=2, then e(L−1)−v_p(L!)>=1 for every p. Thus only L=1 survives in either regime. That band has one column s=0, one derivative level t=1, and entries equal to a nonzero constant times h_i. Its rank is one.

It remains to take p=2,e=1. If m is odd, λ=1 and any L>=2 entry has valuation at least L>=2, since the binomial is integral. Again only L=1 remains.

Suppose m is even. At L=2, the two possible binomials are binom(m,2) and binom(m+1,2). Both have valuation v_2(m)−1, so both survive after division by 2^λ. For L>=3, the consecutive numerator product contains, besides m, at least floor(L/2)−1 further even integers. Therefore

    v_2(R)−v_2(m) >= L−v_2(L!)+floor(L/2)−1
                              = s_2(L)+floor(L/2)−1 >=2,

where s_2(L) is the binary digit sum. Thus no higher band survives. The surviving two columns are s=0,1 and two derivative levels t=2,1. At an odd h_i the corresponding normalized block is

    [1 0]
    [1 1]  over F_2,

because h_i²=h_i and m/2^v_2(m) is odd. All other rows lie in these two columns. The normalized rank is exactly two. ∎

### Consequence: first measurable fiber defect

For the map M modulo p^a, every nonempty output fiber has cardinality

    p^(sum_j min(a,α_j)).

At precision a=λ and a=λ+1, the increase in the exponent of this cardinality is the number of Smith exponents at least λ+1. The theorem therefore measures a one-versus-two-channel release exactly at the first positive threshold. Neither determinant nor worst inverse denominator identifies this release in general.

For node banks in several residue classes, local polynomial CRT is unimodular over Z_p because the cluster node-polynomials have unit resultants. The full partition is their sorted multiset union. Apply the theorem inside each nonsingleton residue cluster, using that cluster's own minimum pair depth; singleton clusters contribute only zeros. This does not determine higher bands.

### New-claim controls

`probe.py` independently builds full integer Hasse matrices and obtains exact Smith forms with SymPy over ZZ. It tests 192 cases: p in {2,3,5,7}, e in {1,2}, m=1,...,6, and h in {(0,1),(0,1,2),(0,1,p),(0,2,3)}. The nested (0,1,p) family is included deliberately. Every case matches the first exponent and its multiplicity. Results are in `results.json`. These are tests of the new theorem, not repeat audits of inherited claims.

## 3. The minimum-value carry replacement for modular support

Let f_i(t)=binom(t,k_i), with fixed k_i>=1 and domains t>=b_i>=k_i−1. Thus each atom is nonnegative and nondecreasing on its domain. For modulus M>=2 put

    P_i(M)=product_(p|M) p^(v_p(M)+floor(log_p(k_i))).

The inherited Pascal/Vandermonde period theorem says f_i(t+P_i)≡f_i(t) mod M. Define the extended integer cost array

    c_(i,M)(a)=min{f_i(t):t>=b_i, f_i(t)≡a mod M},

with +infinity for an empty fiber, and residues represented by 0,...,M−1.

**Theorem (finite compression and exact carry).**

(a) Every finite c_(i,M)(a) is attained at some

    b_i <= t < b_i+P_i(M).

(b) The exact minimum represented nonnegative integer in residue r is the min-plus cyclic convolution

    C_M(r)=min_(a_1+...+a_s≡r mod M) sum_i c_(i,M)(a_i).

(c) Write c_(i,M)(a)=a+M h_i(a) at finite entries. Then

    C_M(r)=r+M κ_M(r),
    κ_M(r)=min_(sum a_i≡r) [sum h_i(a_i)+floor((sum a_i)/M)].

In particular, for 0<=N<M,

    N is represented exactly iff κ_M(N)=0 iff C_M(N)=N.

(d) If M divides M', then

    C_M(r)=min_(s mod M', s≡r mod M) C_(M')(s).

Thus these costs form an exact refinement system, not merely a list of modular feasibility bits.

**Proof.** Subtract multiples of P_i from any t>=b_i until it belongs to the displayed first period. Its residue is unchanged, and monotonicity cannot increase its value. This proves (a). Every tuple induces a residue tuple whose atom costs are no greater, and independently attaining each atom minimum gives an actual tuple, proving both inequalities in (b). Euclidean division of the sum gives (c). For N<M, a nonnegative integer congruent to N is at least N; hence minimum N is equivalent to exact membership. Finally the M-residue fiber is the disjoint union of its M'-residue fibers, giving (d). ∎

The ordinary modular image retains only whether C_M(r) is finite. Its exact lost coordinate is κ_M(r). Finite local multiplicities and reflection-orbit counts from THM-4246 are distinct from this coordinate: many high positive carries do not imply a zero carry.

**Sun application.** Use k=(2,4,6,8), b=(2,3,5,7). All C_M(r) are finite by inherited universal modular solubility. A hole N nevertheless has C_M(N)>N whenever M>N. This is an exact finite target test, not a proof of leastness and not a near-linear algorithm for the large known target. The first-period arrays can be large and their values far exceed M. The construction makes the height obstruction measurable without pretending it has disappeared.

**Cheap hostile.** For ranks(4,6) and starts(3,5), residue3 mod4 is represented by value7, but integer3 is not represented: the minimum cost is7 and the carry is1. Checking only residue existence would miss it. `cost_probe.py` compares costed convolution against independent bounded integer supports for both the four-role Sun packet and this two-role packet, at every modulus2,...,60 (118 cases), including every target0<=N<M. All comparisons pass. The four-role modular images are full throughout that finite control bank.

## 4. Back-and-forth with Chan–Pak

The actual transfer is the following carrier contract, not a claim that binomial atoms equal tournament core outputs.

- Source: a representation grammar with a nonnegative integer output H.
- Target: each residue fiber is labelled by its minimum integer output, and optionally the least construction order among its minimum-output witnesses.
- Preserved predicate: H=N, when modulus exceeds N, is exactly zero carry.
- Destroyed by ordinary residue support: all integer winding and witness order.
- Restoring sidecar: C_M(r), with a reconstructible minimizer; graph order is an additional optimization coordinate, not a replacement for integer output.
- Decisive test: compare minimum carries for the same target with and without each H9 core, rather than recomputing complete residue supports.

For Sun, monotonicity plus a proven input period gives the finite compression theorem above. For Chan–Pak, Session252 gives a finite modular state for the arithmetic compiler but no monotone input-period compression of minimum **integer** H. A shortest modular word need not minimize H. Leading-one extensions can preserve g while word length changes; the graph-order objective is genuinely separate.

Accordingly, an appropriate next experiment is the finite-window cost profile

    C_(M,L,S)(r)=min{H(w): grammar bank S, word cost<=L, H(w)≡r mod M},

retaining the actual witness and seven-coordinate core signature. Compare banks S, S\{S5a}, S\{S5b} at the known repair targets21901 and86641, with M>N. Any zero-carry repair witnesses an integer phenomenon invisible to Session252's universal K1 modular surjectivity. A uniform all-height theorem still requires a contraction or provable cost bound across L; a finite profile is not that theorem.

The jet comparison is equally precise: full Smith data gives finite congruence-fiber cardinalities and deformation thresholds, but still not a small integer representative. The Deuring largest-loss theorem controls only one end of that filtration. The p31 intermediate-ideal hostile explains why a residue bit, determinant, and maximum precision together can still be insufficient. Arithmetic carry is an additional archimedean coordinate, not another Smith factor.

## 5. Universal two-cube forced divisor and collision cap

Here is a modest native extension of the recovered prime-colour decoder. It is a general target-side invoice, not a new support asymptotic.

For n>0 put a_p=v_p(n). If a_3=1, there is no integer two-cube representation. Otherwise define

    D(n)=3^f(a_3) product_(p≡2 mod3) p^(a_p−2 floor(a_p/3)),

where f(3q)=q, f(3q+2)=q+1, and f(3q+1)=q+2 for q>=1; set f(0)=0.

**Proposition.** For every positive distinct representation n=x³+y³, its sum d=x+y is a multiple of D(n). Also

    n<d³<4n.

Consequently the number of positive distinct representations is at most

    #{j>=1:n<(jD(n))³<4n}.

In particular two representations require two multiples of D(n) in that open interval. This is an exact-integer collision-cap computation with no cube-root rounding ambiguity.

**Proof.** Write x=gX,y=gY with gcd(X,Y)=1. At an inert prime, the primitive norm cofactor X²−XY+Y² is a unit, so a_p=3b+u with b=v_p(g),u=v_p(X+Y)>=0. Thus v_p(d)=b+u=a_p−2b>=a_p−2floor(a_p/3). At3, either u=0 and a_3=3b with v_3(d)=b, or u>=1 and a_3=3b+u+1 with v_3(d)=b+u=a_3−2b−1. Minimizing these allowed expressions gives exactly f above; a_3=1 is impossible. Finally n=d³−3xyd<d³ and 4n−d³=3d(x−y)²>0. The inherited good-divisor bijection makes the sum d injective on unordered representations. ∎

This universal invoice is weaker than the full square-discriminant test, but cheaper. It generalizes the **domain** of the restricted inert/cube-free decoder, not its uniqueness strength. Candidate multiples must still divide n and pass the square test from THM-463. No claim is made that a nonempty invoice interval supplies a representation, or that this yields a singleton family outside the inherited ones.

## 6. Next closest frontiers and stopping reasons

- **Closest jet extension:** compute the next positive band after λ, retaining the rank-one versus dyadic-rank-two elimination and actual normalized residue configuration. The first band is now all-node; the next is where node-position polynomials can appear. Do not guess the complete partition from the first band.
- **Closest Sun extension:** implement branch-and-bound min-plus cost convolution at selected moduli with reconstruction, and use admissible partial-sum lower bounds. A zero-carry criterion is proved; efficient enormous-target execution and global leastness are not.
- **Closest two-cube extension:** combine D(n) with short-interval divisor statistics and the discriminant square condition. The invoice gives a rigorous restricted candidate pool, not independence of square events.
- **Chan–Pak:** costed integer fibers are the common actionable object. Neither p-adic ranks, Sun's modular universality, nor K1 Fermat coverage alone closes the all-height target problem.
- **No-three-line:** recovered the uniform frozen-row bound, but no honest arithmetic map; stop rather than force an analogy.

Reproduce from this directory with `python probe.py` and `python cost_probe.py`. Outputs are `results.json` and `cost_results.json`. Neither script writes into the repository or invokes Git. The full mathematical proofs above, rather than finite experiments, supply the unbounded quantifiers.

## 7. Final exchange with the compiler lattice lane

The parallel compiler lane reports that the first three cores already span the universal linear seven-response lattice, index16, whereas their affine difference lattice has index160; either H9 core drops the affine index to32. Both H9 charts therefore make the same affine torsion repair while having distinct exact integer-hole effects. These are **parallel-lane findings**, not recomputed here. They sharpen the shared next experiment: compare zero-carry fibers *after* the common affine obstruction has been discharged. Full affine lattice equality is not equality of the positive realizable semigroup or of bounded-output fibers.

Additional new-claim check: `cube_probe.py` verifies the forced-divisor and exact open-interval invoice for all7140 positive distinct pairs with 1<=x<y<=120, and checks the collision cap for the resulting target fibers. `cube_results.json` records the universe. That coordinate box is not claimed to exhaust every representation of its largest targets; the individual necessary conditions and upper bound are tested on the enumerated subset, which is the correct direction.
