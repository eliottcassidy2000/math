# Focused meta-ledger: triples, canonical codes, and ternary return

This is a targeted source→idea→lever→obstacle survey, not a claim of equivalence between the user's inspirations. All counts H below mean actual directed Hamilton paths, except when a marked/weighted object is explicitly identified.

## 1. Divisor balance {p, p³, p²qr}
Source: https://github.com/eliottcassidy2000/math/blob/99181f1707aff5f8323e21f1cfa86527a3413af1/05-knowledge/results/arithmetic_braids_20260917_divisors.md
Source: https://github.com/eliottcassidy2000/math/blob/99181f1707aff5f8323e21f1cfa86527a3413af1/05-knowledge/results/collatz_mod6_20260917_divisor_balance_family.md

For N>=2, F counts proper nontrivial divisors, S the squarefree ones, and U their prime elements. F=S+U iff N=p,p³,p²qr (last three primes distinct); (F,S,U)=(0,0,0),(2,1,1),(10,7,3). If N=1 is permitted it also satisfies equality. For N=∏p_i^a_i and r=ω(N), F=∏(a_i+1)-2, S=2^r-1-[N squarefree], U=r-[N prime]. Nonsquarefree equality is ∏(a_i+1)=2^r+r+1. The support-4 bound 3·2^(r-1)>2^r+r+1 eliminates r>=4. The later survey explicitly rejects interpreting “three shapes” as three polynomial roots without a map.

Idea: classify three states by retaining a full exponent box, a Boolean subcube and defect, rather than by an undifferentiated count.
Concrete lever: for return-search, carry the full path-cover/port state and a physical defect or guard. Classify low-defect outputs by integral constraints before searching constructors, with p² and p²q as failure controls.
Obstacle: no divisor-to-tournament port map is supplied. Numerical coincidence of 7 or a three-item classification cannot prove a ternary operation.

## 2. S151 and S152, history → compiler → tests
Sources: https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session151_spectrum_metagraph_20261103/spectrum_history/HISTORY.md
https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session151_spectrum_metagraph_20261103/strategy/TERNARY_COMPILER.md
https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session152_sequences_20261103/inductive_target/README.md

S151's state bank target has base odd counts 23..67 and three operations U_j with H(U_j(T))=3H(T)+2j and *closure of the decorated bank*. Arithmetic induction is proved; ordinal sum T→C3 realizes the first identity but not automatically closure, and this output generally is not strong. The +2 and +4 gadgets remain open. S152 exhibits an H=31 strong parent with one-ear outputs 93,95,97, but the chosen H=93 and H=95 children have no one-ear child at any of their own three targets. This is a concrete failure of that selected state continuation, not impossibility for all representatives.

Lever: do not stop at one-step triple availability. Enumerate source states, complete state vectors, legal successors and closure strongly connected components; inspect return words with verified actual witnesses. Hostile: 93 and 95 children; graph state and whether branch is strong must be tracked.

## 3. S156 physical return and S162 family obstruction
Sources: https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session156_spectrum_20261103/returns/NOTE.md
https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session162_mrb_square_20261103/square_transfer/MODULE_BARRIER.md
https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session162_mrb_square_20261103/square_transfer/PRIME_OBSTRUCTION.md

S156 retains three integer coordinates with physical cone h>=0, 0<=x<=y<=h+x, h<=x+y. Its quotient gives necessary slope tests, actual transported reset identities, and excludes 53 of 65 specified last letters by a ray test; 12 survivors are not thereby realizable. No closed ternary bank is shown. S162: square-sum-path transfer yields a block B whose one- and two-slot cover counts are divisible by 3; by binomial inversion and k! divisibility, *all* one-hole modular contexts give H divisible by 3. A named nonuniform attachment escapes and recovers existing 133=7·19 witness. This is a reusable causal contrast: preserving an interface as a module preserves a genuine obstruction; nonuniform wiring can break it.

Lever: perturb membership of a source module, retain endpoint-resolved cover profiles and test projected semigroup/cone closure. Obstacle: literal legal orientations and actual descendant state needed; a scalar or rational solution is not enough.

## 4. S158 Zeckendorf square-sum code
Source: https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session158_creative_20261103/square_overlap/ZECKENDORF.md

A canonical nonconsecutive Fibonacci support I maps to actual H=265∏_{i∈I}(1+2^(i+2)); the factors arise from chosen disjoint square-sum reversal intervals, and ν₂(product-1) recovers the smallest exponent successively. The code is injective, not onto, and is generally not strong.
Lever: seek an information-preserving normal form to make a constructive bank countable with explicit decoding.
Obstacle: a fixed multiplicative interval gadget qH cannot give 3H+2j on an unbounded bank if j≠0, since (q-3)H=2j has at most one solution. Integer encoding, mod-3 diversity, and exact target realization differ.

## 5. Exact Pythagorean endpoint-fiber bridge
Source: https://github.com/eliottcassidy2000/math/blob/f18780daea94f262eb1b7d6fd174015e2657ea86/05-knowledge/results/endpoint_pythagorean_fibers_20261019.md

For C>1 with r distinct prime divisors all 1 mod 4, the E_(r+2) endpoint-reversal tournament has H=1+2^r; its nonexceptional paths are binary Gaussian factor allocations, and complementary address/conjugation pairs correspond exactly to the 2^(r-1) primitive Pythagorean triples of hypotenuse C. Prime powers affect phase, not the number of coordinates. E.g. p and p³ carry r=1, and p²qr carries r=3, when p,q,r are distinct and all 1 mod 4. This is a conditional, marked application to the divisor shapes, not the definition of those shapes. Omitting arithmetic marks loses angles/labels; complement pairing is not generally a tournament automorphism.
Lever: build an exact paired-path or port quotient and retain all decorations that affect subsequent attachment; compare what remains before and after gluing.
Obstacle: no ternary operation or identification of all three mathematical objects follows. The 3 equal-looking target branches are not the 2^r Boolean allocations.

## 6. Collatz and 1,3,5 mod 6
Source: https://github.com/eliottcassidy2000/math/blob/dd7f1b5a1fe24353fac899d4af8c4947840b88ac/05-knowledge/results/collatz_mod6_20260921_berggren_edge_transport.md

For odd n, k=v₂(3n+1), y=(3n+1)/2^k: since 2^k y≡1 mod 3, y≡1 mod 6 for even k and y≡5 mod 6 for odd k; never 3 mod 6. For odd h, 3h,3h+2,3h+4 lie in classes 3,5,1 mod 6. Valuation parity therefore describes which of the latter two *arithmetic* branches an accelerated Collatz output occupies, but does not provide the zero-offset branch. The historical Berggren transport source independently warns that the three children of a primitive triple generally do not preserve its Collatz guard, so three tree children are not automatically three legal descendants.

Lever: label a proposed cellular/return rule by residue and valuation parity; test both forward legal transitions and a backwards preimage/port guard, with all three target classes required.
Obstacle: this is arithmetic inspiration; no Collatz-to-tournament dynamical conjugacy, no global convergence, and no gadget follows.

## 7. New all-size cyclic three-block obstruction
Source: https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session158_creative_20261103/cyclic_blocks/NOTE.md

Let C(a,b,c) be the strongly connected tournament with nonempty transitive classes A→B→C→A. For all positive a,b,c, H(C(a,b,c))≢1 mod 3. More precisely H≡0 mod 3 when the number t of even a,b,c is 0 or 3, and H≡2 when t=1 or 2.

Proof. In S158's run formula Q(n,k)=k!S(n,k), so Q(n,k)≡0 mod 3 for k>=3. For positive n, Q(n,1)=1 and Q(n,2)=2^n-2≡2 if n is even, 0 if odd. All k=0 terms vanish because all three classes are nonempty; all terms k>=2 vanish, as does the 3Q(a,1)Q(b,1)Q(c,1) term. Writing x=Q(a,2) mod 3 etc., H≡x+y+z+xy+yz+zx. If t sizes are even this is 2t+binom(t,2), which is 0,2,2,0 for t=0,1,2,3. ∎

Consequently none of these strong tournaments realizes the 3h+4 branch of S151. The theorem applies to arbitrary transitive class sizes, not arbitrary substituted tournaments; it does not claim novelty. `check_cyclic_mod3.py` independently counts all labelled cases a,b,c∈{1,2,3,4} via subset DP, and passes under python3 and python3 -O.

## 8. Pure-blue vocabulary and remaining type check
Source: https://github.com/eliottcassidy2000/math/blob/466113f5c067dd92b2cf2ece5a3122b461917c55/01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md

This source (pre-audit status) identifies the strong pure-blue alphabet K1,C3,T5, with H=1,3,15, where pure-blue means every Hamilton path is palindromic under positional reversal as an anti-automorphism. It does not identify the supplied numbers 7,17,55 as H counts. Their intended type and source need a separately grounded identification; neither label collision nor shared three-ness proves a correspondence.

