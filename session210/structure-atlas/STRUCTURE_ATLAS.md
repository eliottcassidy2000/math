# The 3/10 indices: native objects, inherited obstructions, and current spectrum routes

Session 210. Status: source-checked synthesis, an elementary corollary of the
known odd-cycle formula, and independent finite controls. No new completeness
theorem, new Hamiltonian-path bijection, or literature-priority claim.

The reference snapshot is main commit
`968dbe3af678bcfdb6c2a94fa5b8d5ce394f992c`. `SOURCES.json` records fixed-commit
URLs and Git blob IDs fetched through the GitHub connector. The shared dirty
checkout was not changed. Historical source claims are scoped below rather
than inherited wholesale.

## 1. A literal object counted by the odd index

Let `alpha_k(T)` count the sets of `k` pairwise vertex-disjoint directed odd
cycles of a tournament `T`; cycles have length at least three, and cyclic
rotation does not produce a different cycle. The inherited odd-cycle formula
(Grinberg-Stanley; repository THM-002/PROP-001) is

    H(T) = sum_(k>=0) 2^k alpha_k(T),   alpha_0(T)=1.

Define `Q(T)` as follows. Choose a nonempty such cycle packing and color each
cycle with one of two colors. Identify a coloring with its simultaneous
global color swap. The swap has no fixed points on a nonempty packing, so

    |Q(T)| = sum_(k>=1) 2^(k-1) alpha_k(T) = (H(T)-1)/2.

Thus the two proposed missing indices `3,10` are literal cardinalities of a
natural finite orbit set, not only a change of coordinates. Isomorphisms of
tournaments transport this object. It records actual cycle supports and
cycle words; merely recording their number forgets the intersections.

This proves a cardinality interpretation. It does not give an explicit
bijection between Hamiltonian paths with one distinguished path removed and
two copies of `Q(T)`. Such a bijection would need an actual path map, and
equivariance or a canonical distinguished path is a separate obligation.

For a one-apex tournament every directed cycle passes through its apex, so
only `alpha_1` is nonzero. In that special family the index is literally the
number of directed odd cycles. This is exactly the `r_core` coordinate of
HYP-1782 and HYP-2198, not a new interpretation of that family.

## 2. The operation law is not addition

In an order join `A -> B`, directed cycles stay inside their component.
Every colored packing either lies in one side or uses both. In the mixed
case the independent choices of the two global colorings leave a two-element
relative-color torsor after quotienting by simultaneous swap. Therefore

    q(A -> B) = q(A) + q(B) + 2 q(A) q(B),
    q(T) = |Q(T)|.

Equivalently this is the usual inherited multiplicativity of `H=2q+1`.
The relative torsor has two elements, but no preferred labeling by 0 and 1
is claimed without choosing normalized representatives.

An immediate screening rule follows. A class closed under disjoint union
with an additive size statistic cannot have size spectrum exactly
`N_0 minus {3,10}`: it would contain size 1 and size 2 and hence size 3.
In particular calling this an additive numerical semigroup of genus two is
incorrect. Its natural operation is `a+b+2ab`.

The control `C3 -> C3` has packing profile `(1,2,1)`, `H=9`, and `q=4`, not
`q=2`. Treating its two cycles as an additive pair loses the relative-color
choice and the mixed packings.

## 3. Why arbitrary independence complexes do not solve the problem

Let `Omega(T)` have the directed odd cycles as vertices, with two adjacent
exactly when they intersect. Then `H(T)=I(Omega(T),2)`. The orbit construction
above works for the nonempty independent sets of any graph `G`, but allowing
arbitrary graphs immediately loses the gaps:

    G=K_m gives I(G,2)=1+2m and |Q(G)|=m, for every m>=0.

Consequently the whole conjecture lives in the realizability constraint
`G=Omega(T)`, not in the formal independence polynomial or the color-swap
quotient alone. `K_3` and `K_10` are especially transparent forbidden formal
targets, but the `H=21` exclusion also has non-clique candidate profiles.
The graph/complex construction is a faithful counting carrier, not an
independent classification of its realizable values.

## 4. The single-core route was already completed, and is insufficient

For a binary word `s` describing an apex above a transitive chain,

    r(s) = sum_(i<j, s_i=1, s_j=0) w(j-i),
    w(1)=1, w(d)=2^(d-2) for d>=2,
    H(s)=1+2r(s).

Removing leading zeros and trailing ones preserves `r`. A nonzero normalized
word begins in 1 and ends in 0. If it has length `L>=3`, the extreme pair
already contributes `2^(L-3)`. Thus every index `r<=R` has a witness of length
at most `3+floor(log2 R)`. HYP-2198 used this to certify the all-length
single-core support through `R=131072`, including the absence of `3,10`.

The same construction also misses `6,14,17,20,21,24,...`. In particular
`H=13`, which is actually attained by a five-vertex tournament, has no
one-apex realizer. Conversely `r=31`, hence `H=63`, is realized by the words
`1001100` and `1100110`. Neither a seven-divisibility obstruction nor a
Mersenne-gap tower survives.

Our independent small replay certifies the one-core support through `R=128`
using the length bound, and compares the weighted formula with direct
Hamiltonian subset DP on all 255 words of lengths zero through seven.

The larger intercyclic route, in which no two directed cycles are disjoint,
is also already classified in the audited Session16 packet. The current
THM-4885 source note records its spectrum as the single-core spectrum plus
`3*2^t+1` (`t>=2`) and `9*2^t-3` (`t>=1`). Thus it is incomplete too. A
new complete-spectrum route must allow genuine compatible cycle packings.

## 5. The strongest inherited finite family and exact current scope

For a one-apex word with binary value `v`, the inherited rooted triple is

    (g(v), e(v), 1+v),
    g(0)=1,
    g(2v+1)=g(v),
    g(2v)=g(v)+v+(v mod 2),
    e(v)=1-(v mod 2).

If `(a,b,c)` is an actual rooted tournament triple, the one-arc join law gives
an actual tournament with

    Phi_(a,b,c)(v) = a g(v) - b e(v) + c(1+v).

The fixed bank of 2,732 distinct rooted triples from orders at most seven is
therefore a particularly concrete candidate: a finite union of integer
sequence ranges, every value backed by a literal tournament construction.
The source state, not an arbitrary triple of integers, is essential.

- THM-4875(D): the two join families together attain every odd value through
  `2^32`, except `7,21`; independently reproduced. Session203 recertified
  this with the redundant plain family disabled.
- HYP-9478: the fixed-bank family alone is independently reproduced through
  `2^30`.
- HYP-9478's November6 addendum: the fixed-bank family alone attains every
  odd value through `2^35=34359738368`, except `7,21`. This is explicitly one
  run at that largest scale, with lower-scale cross-checks, not an independent
  rerun at `2^35`.

The unbounded assertion that this fixed bank suffices remains open. The
two-one-apex family alone does not cover those finite intervals: at block
`[2^31,2^32)` it has 934,107 missing odd values. Do not conflate this family
with the union or with the fixed-bank family.

## 6. Other already-explored structural routes

The following are inherited, not fresh discoveries of this session.

1. **Strong-component arithmetic:** THM-4094 reduces global completeness,
   using actual carry witnesses 49,63,343, to strong realization of the two
   lanes `p` (odd prime other than 7) and `7p` (odd prime other than 3).
   Factorization alone cannot supply new prime realizers.
2. **One-vertex interpolation:** THM-4876 proves that an increment exactly 2
   forces the old count to be 1 or 3. A universal `H -> H+2` one-vertex
   induction is therefore impossible once `H>=5`.
3. **Finite-rank spine grammars:** THM-4878 organizes one-apex, two-apex,
   cutwidth-two and ribbon-like families as rational weighted automata.
   Hankel ranks four and seventeen recur. Their range questions are open;
   finite rank does not make the output range an interval.
4. **One stationary nontransitive core:** Session197's
   `ONE_CORE_BINARY_INVERSE.md` gives an exact two-integer normal form and
   finite target inverse. The order-at-most-five core bank first misses
   462661 in addition to 7,21; two actual cores repair it. This is a clean
   example of a restricted-language hole, not a global forbidden count.
5. **Analytic routes:** THM-4885 proves first-moment mass laws, including the
   fixed-bank mean 48.197612 and the two-apex `pi^2/12` logarithmic law.
   THM-4886 proves minor-arc bounds. THM-4887 proves
   `A(X)>=X^(2/3)/(13 ln X+102)` for `X>=4096`, where `A` counts distinct
   attained values. Positive density, density one and full completeness
   remain open in these sources. Mean mass and residue mixing do not imply
   every integer has a representation.

## 7. Corrections that matter when recovering older work

- THM-1370's corrected statement proves the omissions 7,21 and explicitly
  leaves every-other-odd attainment open. The omission and completeness
  discussion already appeared on MathOverflow in 2016; Chan-Pak records the
  conjecture later.
- HYP-2198's interpretation contains the stale phrase that only 7,21 are
  globally forbidden. Its all-length **single-core** theorem survives;
  the global phrase does not prove completeness.
- HYP-3617 claims that the global spectrum is known and invokes an invalid
  factor intuition for 35. A strong tournament can have `H=35` without any
  component having `H=7`; prime divisors of a path count need not themselves
  be component counts. Use THM-4094's actual strong-factor reduction.
- MISTAKE-055 and its repair replace false fitted strong-floor recurrences
  by Busch's residue law `f(n+3)=5f(n)` for `n>=3`, starting `3,5,9`.
- The THM-4878 source note's prior-art paragraph says the two-one-apex family
  alone covers through `2^32`. THM-4875 and HYP-9478 distinguish the families
  and disprove that reading with the explicit miss counts. This atlas uses
  the more specific family statements above.

## 8. Independent controls and remaining obligation

Run:

    python3 verify_orbit_carrier.py > RESULT.json
    python3 -O verify_orbit_carrier.py > RESULT.optimized.json
    cmp RESULT.json RESULT.optimized.json

The checker imports no repository computation. It counts Hamiltonian paths
by subset DP and odd cycles/compatible packings by literal enumeration on
all 1,099 labelled tournaments of orders one through five, plus 32 fixed-seed
samples at each order six, seven and eight. The output distinguishes these
universes. It also checks the two-cycle order join and the binary-word
controls. All checks use explicit exceptions and are identical under `-O`.

The most faithful response to the 3/10 idea is therefore twofold: the indices
already count a genuine natural orbit object, but its tournament-realizable
subclass is the hard part. For a new completeness proof, produce an
unbounded inverse into an actual physical family, such as the fixed-bank
one-apex readouts, or prove the two prime lanes. A matching scalar spectrum
in unrelated objects would be a lead, not a bijection or a proof.
