# Seven old bridges toward an operation-faithful theory of symmetry

**Session 88 archive lane.** These are reusable mechanisms, not a keyword
inventory and not claims that the source subjects are equivalent. The source
chain runs from March through August 2026. Every principal theorem below was
read at statement/proof/boundary level; its local bytes were checked against
commit `7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51`.

## The common question

A symmetry quotient should be judged relative to **an operation alphabet and
an observable**, not by whether its current fingerprint happens to be complete.
Three different objects must not be called the same state:

1. a present-time decision certificate;
2. enough information to propagate selected future observables;
3. enough information to reconstruct or compose the original objects.

The old sources prove actual separations between these. The new arithmetic
result below adds a particularly clean example: all future one-step response
values recover an entire permutation cycle type, although the current exact
edge-orbit count does not.

## 1. Burnside cycle types → divisor incidence algebra → response completeness

**Source:** `01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md`
(June 15). **Status:** inherited PROVED; new general theorem PROVED in
`poset_future_state_theorem.md`, with a separate FINITE-EXACT audit.

- **Source/target/map:** a unit-free odd permutation cycle multiset `mu`
  maps to divisor counts `M_d=#{p in mu:d|p}`. This is the zeta transform on
  the odd-divisibility poset, invertible by Möbius inversion.
- **Preserved:** adding a `p`-cycle increments every `M_d` with `d|p`; the
  edge-orbit exponent is `e=C(t,2)+1/2 sum_(d>=3)phi(d)M_d^2`.
- **Lost by the tempting quotient:** `(mass,t,e)` retains a current fixed-
  permutation count `2^e`, not the divisor-resolved derivative under adding
  a cycle. It also does not determine the Burnside denominator `z_mu`.
- **Cheapest hostile:** `(3,11)` and `(5,9)` both have `(mass,t,e)=(14,2,7)`.
  Adding 3 agrees; adding 5 gives increments `4,8`. The new fixed counts
  are `2048,32768`. This is the smallest such collision **among unit-free
  odd cores**; with units allowed the first is `(1,7)` versus `(3,5)` at
  mass 8, and adding 3 separates them.
- **Derived theorem:** for any finite positive weighted incidence Gram
  matrix `C=Z^T diag(w)Z`, the corrected add-atom response vector of
  `Q(m)=m^TCm/2+ell^Tm` is `Cm`, hence determines `m`. For unit-free odd
  parts `3..N`, all odd queries `p=3..N` suffice without a known length
  because the restricted GCD matrix is positive definite. Unit-containing
  scopes are separated explicitly in the theorem note. For a weight-preserving
  poset automorphism group, this response embedding is equivariant and
  preserves stabilizers exactly; response orbits recover state orbits.
- **Why it matters now:** it gives an exact meaning to *faithful future
  state* across arithmetic, species/Burnside symmetry, and incidence algebra.
  A short present certificate cannot be promoted to a closed recursion
  without checking these operation responses.

**Boundary:** cycle addition is a formal operation on permutation cycle
**types**, not one-vertex extension of a tournament. The theorem does not rule
out compressed aggregate recurrences for `B[m,t]`, analytic formulas, or
algorithms that avoid reconstructing individual cycle types. No significance
is attached to the hostile's mass coinciding with a different problem's 14.

## 2. Inverse-digit collisions → conditional expectation → a quantitative loss test

**Source:** `01-canon/theorems/THM-2519-last-digit-collision-drift-and-k13-dirichlet-boundary.md`
(July 27). **Status:** PROVED, exact-verified and independently audited.

- **Source/target/map:** a depth-`L` pair of inverse branches of `x->13x`,
  with last relative digit `u` and higher address `e`, maps to
  `h=P_(13^(L-1))F` by **averaging** the higher address. The future-owner
  weight is `W(y)=G(13y)`, with `G>=0`.
- **Preserved:** the collision averages become exactly
  `B_u=int W(y)h(y)h(y+u/13)dy`. The one-digit averaging operator `E_13`
  is an orthogonal conditional expectation because `W` is fibre-invariant.
- **Exact loss/equality law:**
  `B_0-(1/13)sum_u B_u=int W|h-E_13h|^2` is the weighted `K_13` Dirichlet
  energy. It vanishes iff the future owner sees equal mass at all thirteen
  predecessors. The rational-step case has an all-zero/all-positive
  nontrivial Fourier-colour dichotomy.
- **Lost:** high-address identity and antipodal orientation: `B_u=B_-u`.
  A positive square spectrum supplies neither a directed tournament edge
  nor an owner-oriented service certificate.
- **Corrected near miss/hostile:** a Boolean future-measurable
  delta-plus-replicas source can have complete bulk square spectrum but
  `h=E_13h`, hence **zero** last-digit drift. Bulk charge does not force a
  new collision colour.
- **Transfer to the new work:** on any finite quotient block, averaging is
  rank one, and `||r-Er||^2` is an exact response-loss certificate. The
  arithmetic hostile in §1 has loss `4` for the add-5 response and `0` for
  add-3. This gives a natural diagnostic for weighted forgetting maps or
  block-based interval-Gram certificates: ask which exact observable has
  zero conditional variance, rather than infer completeness from positivity.

**Boundary:** positivity needs the nonnegative, fibre-invariant weighting.
This theorem does not justify arbitrary reweightings or commuting two
forgetting maps. It provides the single-forgetting operator identity and
sharp nullspace that such a new commutation theorem must respect.

## 3. Marked tilings → groupoid spans → source-preserving gluing

**Sources:** `01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md`
(July 15) and `THM-830-b3-deletion-deck-mirror-current-calculus.md` (July 14).
**Status:** all-size structural laws PROVED; numerical non-lumpability and
classification ranges FINITE-EXACT.

- **Source/target/map:** `X_n`, tournaments with a fixed directed Hamiltonian
  observer path, has maps to merged unmarked nodes `M_n` and complement
  lines `L_n`: the exact object is the span `M_n <- X_n -> L_n`.
  A node fibre is `HP(T)/Aut(T)` (or its converse-pair union), not a unique
  marked representative.
- **Preserved gluing:** marked tilings satisfy
  `X_n ~= (X_(n-1) x_(X_(n-2)) X_(n-1)) x F_2(apex)`.
  The overlap is the actual common interior tiling. Lower **line** pairs
  instead leave a `C_2` coherent endpoint-phase torsor for `n>=5`.
- **Lost:** forgetting the path destroys a canonical parent map. Forgetting
  coherent phase destroys unique gluing. Independently sorting the upper
  and lower endpoint pairs destroys straight-versus-crossed coupling.
- **Cheapest hostiles:** already at `n=4`, simultaneous endpoint coupling
  has 4 support types while independently sorted marginals have 3. At
  `n=5`, composing normalized node-deletion kernels disagrees with actual
  two-step deletion in `16` entries across `8/10` high rows (maximum error
  `1/2`). The intermediate node does not induce a uniform marked fibre.
- **Further exact coordinate:** THM-830 identifies the tiling with a
  labelled tournament on the *gaps* of its observer path. Internal deletion
  must retain its newly adjacent shortcut/seam bit; the full A/B/C
  recursion is a three-axis address-simplex recurrence, not bare deletion.
- **Why it matters now:** one row per root-switch equivalence block is
  mathematically safe only for the declared equivalence/decision. For
  gluing or further operations, retain the source/target coupling and
  stabilizer/phase data, or prove the relevant operation descends. A span
  is the correct transfer object when a function is unavailable.

**Boundary:** equal weighted correspondences are not equal object maps;
complete static fingerprints need not be strongly lumpable. The July theorem
has exact ordinary/converse-sheet loop holonomy, but its source tree and
transport groupoid must not be confused with §5's thin ambient tree.

## 4. Even-path forests → labelled species → a minimal linear future-mean module

**Sources:** March `01-canon/theorems/THM-076-walsh-ocf-factorization.md`,
continued by July `THM-848-hamiltonian-path-walsh-needle-stalk.md`.
**Status:** all-size formulas/minimality PROVED; directional classification
through eight vertices FINITE-EXACT only.

- **Source/target/map:** Hamiltonian-path count `H` on the full arc cube
  decomposes into even Walsh layers `H_(2r)`. The March coefficient law is
  supported on disjoint even-length path forests, with a factor `2` per
  unrooted path component. This is an actual labelled assembly species,
  not just a resemblance between generating functions.
- **Preserved:** under uniform random arc flip `P`, each layer has eigenvalue
  `1-4r/binom(n,2)`. The vector `W_n=(H_2,...,H_(2q))` therefore propagates
  every conditional future mean of `H`; Vandermonde inversion gives its
  linear minimality. The invariant function space has dimension `q+1`
  including the constant.
- **What less data can do:** `(H,K)` determines one-step averaged drift,
  where `K=nH-2 sum_r rH_(2r)`. It does not determine every future mean.
  The unordered directional-derivative multiset retains the response
  spectrum but not which target node each response reaches; the
  target-coupled derivative orbit restores that incidence.
- **Corrected near miss:** MISTAKE-151 headed by the Walsh-level energy
  error retracts the extrapolated numerator `2(n-2r)^r`. Squaring the
  forest expansion weights components by `4`; the exact labelled EGF is
  `exp(x+2x^3 z/(1-x^2 z))`. The level numerator is
  `[w^r]((1+w)/(1-w))^(n-2r)`. At `n=8,r=3` it is `12`, not `16`.
  All tests through seven vertices had accidentally lain on equality
  boundaries.
- **Why it matters now:** a finite response module may be minimal for a
  *chosen linear observable* without being a full object state. This is
  the clean positive counterpart to the non-Markov quotient in §3 and the
  no-finite-state result in §6. It also warns sparse pattern/Gram
  reconstruction to retain support-component multiplicities, not only a
  low-order formula that happened to agree on boundary cases.

**Boundary:** the minimality is linear-dynamical for the uniform full-arc
flip operator, not arbitrary root switches, deletion, gluing, or Mahler
continuation. A finite-exact complete invariant through `n=8` is not an
all-size isomorphism theorem.

## 5. Gaussian factor toggles → Berggren LCA paths → flat source-dependent symmetry

**Source:** `01-canon/theorems/THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost.md`
(August 12). **Status:** PROVED with exact independently audited witnesses.

- **Source/target/map:** a Gaussian prime-allocation XOR torsor at fixed
  hypotenuse maps each external edge `(x,x+e_p)` to the unique undirected
  path in the Berggren ancestry tree. If addresses are `w_x,w_y` and their
  common-prefix length is `r`, the path is
  `reverse_inverse(w_x[r:]) w_y[r:]`; its length is
  `|w_x|+|w_y|-2r`.
- **Preserved:** source, target, inverse arrows, and coherent composition.
  In particular
  `a_p(x)a_q(x+e_p)=a_q(x)a_p(x+e_q)`; exact squares are flat.
- **Lost by prime-colour alone:** ancestry cost and actual path word. In
  the `1105=5*13*17` fibre, the two prime-5 edges have costs `8` and `27`,
  despite identical prime and folded factor weight. Other directions also
  have unequal costs. A rooted-tree automorphism cannot implement these
  matchings, because some paired vertices have unequal depths.
- **Corrected near miss:** freezing one path per prime at one basepoint
  produces apparent defect-word lengths `39,39,37`. Those are **not**
  intrinsic curvature: the second arrow has the wrong source. Restoring
  source-dependent arrows makes every closed ancestry path reduce to the
  identity. The external K4 has cycle rank three while its tree realization
  has zero ambient holonomy.
- **Why it matters now:** a root switch may be a genuine symmetry of an
  external fibre without being a uniform word in the source's native
  operations. LCA/address data are necessary for source-preserving
  composition even when a short external label decides equivalence. This
  is a useful hostile to any proposed block-to-native-operation dictionary.

**Boundary:** the theorem excludes fixed ancestry words and rooted-tree
symmetries, not every source-reading transducer. It proves an infinite
prime-5 family with unbounded dispersion and a rational transduction on one
unary branch, but leaves a uniform transducer across all Boolean fibres open.

## 6. Mahler carry safety × native termination → exact fibre product, no finite state

**Source:** `01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md`
(August 25). **Status:** PROVED; concrete minimal hostile FINITE-EXACT.

- **Source/target/map:** a carry word `c` maps compatibly, level by level,
  to its native binary word `e=theta(c)`, preserving the exact 2-adic point,
  depth and finite cylinders. The update requires an integer sidecar:
  `e_m=c_m xor (u_m mod 2)` and
  `u_(m+1)=(3u_m+c_m+3^(m+1)e_m)/2`.
- **Exact preserved predicate:** positive candidate Z-number codes are the
  synchronized fibre product of strict carry safety and positive ordinary
  stabilization. Safety is nonrejection plus infinitely many reset edges;
  stabilization is finitely many native ones, with at least one. The
  safe follower states are countably many and have distinct continuation
  languages.
- **Lost:** follower state alone forgets the native clock; native bits
  alone forget strict/equality phase; retaining clock zero still loses
  the parity and then unbounded value of `u_m`.
- **Cheapest hostile:** at terminal depth four, `A=8` has carry prefix
  `0001` and `A=13` has `1001`; both have follower state `q=1`, clock `h=0`,
  and positivity seen. At the next step, 8 rejects and 13 survives (until
  its first rejection at index ten). The common short certificate does
  not even determine the next decision.
- **Stronger obstruction:** every finite safe prefix has a native positive-
  terminal completion, so the naive finite terminal-prefix test prunes
  nothing. The rooted pair-prefix language is nonregular and `theta`
  has no finite synchronous transducer. The exact product is not thereby
  empty, nonempty, or a finite automaton.
- **Why it matters now:** finite per-instance certificate compression and
  a uniform finite future state are fundamentally different promises.
  The contrast with §4 is useful: specify the observable and operation
  before proclaiming a universal finite-state symmetry theory.

**Boundary:** Mahler's Z-number problem remains open. The unknown full
candidate intersection is not called nonsofic; nonregularity applies to the
specified rooted finite-prefix language. Inverse limits of finite relaxations
retain equality/nonterminal rays and need not equal the target language.

## 7. Sum of two cubes → Eisenstein norm factor → an exact divisor gate

**Source:** `01-canon/theorems/THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis.md`
(June 10). **Status:** elementary complete PROVED; large finite audits are
supporting evidence, not the proof.

- **Source/target/map:** unordered integer pairs `{x,y}` with `x^3+y^3=n>0`
  map bijectively to positive divisors `d=x+y` for which
  `s=(d^3-n)/(3d)` is integral and
  `Delta=(4n-d^3)/(3d)` is a square. Recovery is
  `{(d+sqrt(Delta))/2,(d-sqrt(Delta))/2}`; parity is automatic.
- **Preserved:** the full unordered solution, including negative members
  and the repeated-root boundary. This is a quotient by `S_2` with a
  complete elementary-symmetric coordinate, not a lossy support count.
- **Arithmetic bridge:** the complementary factor is
  `n/d=x^2-xy+y^2=N(x+y omega)` in the Eisenstein integers. For a primitive
  pair, it is odd, its primes other than 3 are `1 mod 3`, and its 3-adic
  valuation is at most one. This is the same norm/divisor-character axis
  as the older Moser-ladder transverse-vector count.
- **Corrected near misses:** split-prime support alone is not the square-
  discriminant gate. Doubly-primitive taxicab numbers need not be fully
  split: inert content must lie in the common sum-divisor, not the norm
  cofactor. The first recorded inert example is `4342914`. The number
  `1729` genuinely links the taxicab and Eisenstein/Moser mechanisms,
  but the separate tournament-1729 resemblance was refuted as coincidence.
- **Why it matters now:** this is a positive model of a *short complete
  certificate*: fixed `n` plus one good divisor determines an entire
  `S_2`-orbit. It contrasts with §1's current-energy certificate and §6's
  terminal-state certificate. One should seek such decoder proofs—not
  merely small fingerprints—when claiming reconstruction of root-switch
  patterns from sparse numerical data.

**Boundary:** the decoder solves this fixed two-variable cubic symmetry
problem, not general equal-sum-of-cubes problems or tournament enumeration.
The source/target arithmetic map, not the shared record integer, supplies the
actual connection.

## Recommended synthesis, not another census

The seven sources suggest a compact working language:

```
carrier + symmetry + operation alphabet + observable
    -> quotient + exact decoder or explicit loss operator
    -> continuation test / gluing test / equality boundary.
```

Two positive outcomes are different and both valuable:

- **certificate completeness:** a small witness decides or reconstructs the
  current target (the good divisor; potentially an interval-Gram decision);
- **operation completeness:** retained coordinates support every specified
  future update (divisor profile; Walsh/Krylov module for mean flow; marked
  fibre plus coherent phase for gluing).

Three inherited hostiles should accompany any new claimed bridge:

1. **wrong source:** the Gaussian frozen-arrow false holonomy;
2. **wrong coupling/weight:** independently sorted endpoint pairs and the
   false two-step Markov product;
3. **wrong future:** same current energy or terminal clock, different next
   response.

The new incidence-response theorem closes one concrete loop between these
threads. It does not claim that every mathematical symmetry is one group
acting on one state space; the archive repeatedly shows why typed groupoids,
spans, conditional expectations and observable-specific modules are needed.
