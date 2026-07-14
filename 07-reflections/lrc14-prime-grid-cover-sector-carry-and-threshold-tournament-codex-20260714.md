# LRC(14): persistent prime-grid covers, sector carries, and the threshold tournament

**Status:** exact reformulations and elementary lemmas proved below; no proof of LRC(14).

The corrected frontier (HYP-6780, THM-758, THM-760) says that raw speed size is not the
right compactness variable. This note supplies two complementary finite objects for the
remaining work:

1. the paper's initial sieve `I(13,p,1)` is exactly a cyclic translate-cover problem, with
   a forced supply of private blocker times and a Fourier pressure inequality;
2. the round tournament becomes threshold-faithful only after adding a 14-sector potential
   and one global microphase order. The binary tournament alone provably cannot decide even
   the instantaneous witness predicate.

The resulting hard object is a **projectively persistent blocker cover with scale/residue and
carry data**, not a raw runner tournament or a bounded speed box.

## 1. Exact prime-grid translate cover

Fix a prime `p>14`, put

`h=floor((p-1)/14)`,

and define

`G=F_p^x/{+1,-1}`, `n=|G|=(p-1)/2`,

`D={[1],...,[h]} subset G`.

Here `D` is the strict danger set: equality at distance `1/14` is safe. Since `p` is
not divisible by `14`, there is no prime-grid endpoint ambiguity.

For a tuple `v=(v_1,...,v_13)` of nonzero residues, let

`x_i=[v_i]^{-1} in G`.

### Proposition 1 (initial improper tuples are covers)

At lift level one,

`v in I(13,p,1)` if and only if `G = union_i D x_i`.

**Proof.** The omit-one gcd alternative in the definition of properness is impossible when
`l=1`. A nonzero grid time `a/p` is blocked by runner `i` exactly when

`||a v_i/p||<1/14`,

equivalently `[a v_i] in D`, equivalently `[a] in D[v_i]^{-1}=D x_i`. Thus every
nonzero grid time is blocked exactly when the thirteen translates cover `G`. Sign changes
vanish in `G`, permutation only reorders the blocks, and common scaling translates all
centers simultaneously. QED.

Although `G` is cyclic, `D` is the discrete-log image of a short **additive** interval, not
generally an interval in exponent coordinates. Ordinary cyclic interval-cover theorems do not
apply without further work.

### Proposition 2 (the tight cover is always present)

The tuple `(1,2,...,13)` gives a thirteen-translate cover for every prime `p>14`.

**Proof.** For any nonzero `a`, place `0,a,2a,...,13a` on the additive `p`-cycle. The
fourteen circular gaps sum to `p`, so one has integer length at most
`floor(p/14)=h`. Its endpoints differ by `d a` for some nonzero
`d in {-13,...,13}`. Therefore `d a` is strictly dangerous and `[a]` belongs to the
translate belonging to speed `|d|`. QED.

Consequently the `k=13` task cannot be to make `I(13,p,1)` empty. It is to classify these
covers and show that every projective lift tree eventually dies or enters the omit-one gcd
route. This sharpens the April 2026 paper's stated bottleneck.

## 2. Private blocker pressure

For a thirteen-center multiset `X`, define its covering multiplicity

`C_X(g)=sum_i 1_D(g x_i^{-1})`.

Let `A=13h` and `alpha=A/n`. For `p>14`, one has `1<alpha<2` whenever `X` is a cover.

### Proposition 3 (private-time lower bound)

Every thirteen-translate cover has at least

`2n-A`

points covered by exactly one translate. Hence some runner is the unique blocker of at least

`ceil((2n-A)/13)`

prime-grid classes.

If `p=14q+r`, with `r in {1,3,5,9,11,13}`, this becomes

`2n-A=q+r-1`.

**Proof.** Since `C_X>=1`, the total redundancy is

`sum_g(C_X(g)-1)=A-n`.

Every non-private point consumes at least one unit of this redundancy. There are therefore at
most `A-n` non-private points, and at least `n-(A-n)=2n-A` private points. Assign each
private point to its unique translate and apply pigeonhole. QED.

This is directly recursive information: when a cover is lifted, every private base fiber must
either remain with its unique owner or acquire a new blocker. A lift search should branch first
on the owner with the largest private fiber, rather than enumerate coordinate tuples uniformly.

### Proposition 4 (Fourier pressure)

Use the unnormalized transform

`hat f(chi)=sum_{g in G} f(g) conjugate(chi(g))`.

Let `mu_X=sum_i delta_{x_i}`. Every cover satisfies

`sum_{chi != 1} |hat 1_D(chi)|^2 |hat mu_X(chi)|^2`

`    >= n^2 (alpha-1)(2-alpha)`.

**Proof.** Since `C_X(g)` is an integer at least one,

`(C_X(g)-1)(C_X(g)-2)>=0`.

Averaging gives

`Var(C_X)>=3alpha-2-alpha^2=(alpha-1)(2-alpha)`.

Parseval and the translate formula for `C_X` identify the variance with

`n^{-2} sum_{chi != 1}|hat 1_D(chi)|^2|hat mu_X(chi)|^2`. QED.

Equality occurs only when every grid class has multiplicity one or two. This does not eliminate
the tight cover, but it gives an exact spectral pruning test and identifies the survivor as a
multiplicative-Fourier alignment between the short additive interval and the speed centers.

## 3. The 14-sector carry object

At a time `t`, let `p_i={v_i t}` (including the stationary observer `p_0=0`) and set

`c_i=floor(14 p_i)`, `rho_i={14 p_i}`,

`gamma_ij=floor(14 {p_j-p_i}) in {0,...,13}`.

### Proposition 5 (sector-carry factorization)

Exactly,

`gamma_ij = [c_j-c_i-1_{rho_j<rho_i}]_14`,

where brackets denote the least nonnegative residue modulo `14`.

**Proof.** Write `14p_i=c_i+rho_i`. Before reduction modulo `14`, the difference is

`c_j-c_i+rho_j-rho_i`.

Its floor loses one precisely when `rho_j<rho_i`; reduction of the integer part modulo `14`
then gives the formula. QED.

Thus the edge-labelled object is not an arbitrary gain tournament. It is a `Z_14` vertex
potential plus the inversion matrix of one global weak order of the microphases. The continuous
unfloored differences are a genuine circle-valued coboundary; flooring introduces one carry bit.

Immediate consistency laws include

`gamma_ij+gamma_ji=13` when `rho_i!=rho_j`,

with sum `14` on a noncollision sector wall and sum `0` at a collision, and

`gamma_ij+gamma_jk-gamma_ik in {-1,0,13,14}`.

This compresses the nominal directed edge table to the sector potential vector plus a weak
order. It is a natural dynamic state for a fixed-clock movie, but static sector realizability
alone cannot prove LRC: the constraint must come from how the state evolves with `t`.

### Closed-threshold observer readout

Away from sector walls, runner `i` is safe from the observer exactly when

`gamma_0i notin {0,13}`.

For the closed LRC inequality, the exact rule is

`gamma_0i in {1,...,12}`

or

`(gamma_0i,gamma_i0)=(13,1)`.

The exceptional pair is the safe endpoint `p_i=13/14`. A one-sided floor label merges this
endpoint with the dangerous interior `(13/14,1)`; the two directed labels or an explicit wall
flag are necessary.

## 4. Why the raw tournament is too small

Use the convention `i -> j` when `{p_j-p_i} in (0,1/2)`. The half-turn round tournament
records only whether a relative phase lies in a semicircle.
For a pair with speed difference `d`, its movie has `2|d|` collision/antipode walls per period.
The 14-sector edge has `14|d|` walls. The binary quotient therefore erases exactly `6/7` of
the wall events, including the LRC transitions at relative phases `1/14` and `13/14`.

This loss is not merely heuristic. At `t=1/100`, compare

`S_A={6,12,18,24,30,36,42,48,54,63,70,78,87}`,

`S_B={8,12,18,24,30,36,42,48,54,63,70,78,87}`.

Moving the first phase from `0.06` to `0.08` crosses no collision or antipodal wall with the
observer or any other listed phase, so the complete labelled half-turn tournament is identical.
Yet `t` is not a witness for `S_A` because `0.06<1/14`, while it is a witness for `S_B`
because every phase has circular distance at least `0.08>1/14`. No static statistic of the raw
tournament—scores, cycles, SCCs, Hamiltonian paths, or spectrum—can decide the witness predicate.

## 5. A faithful mixed threshold tournament

There is a useful tournament after changing what observer arcs mean. Off the finite wall set,
let

`B(t)={i: ||v_i t||<1/14}`, `F(t)=[13] minus B(t)`.

Orient the observer incidences by the LRC predicate,

`blocker -> 0 -> safe`,

and retain the half-turn orientation only between moving runners. This is the repository's
mixed observer-source construction; it must not be confused with the pure round tournament.

Unwrap blockers as `b in (-1/14,1/14)` and safes as `s in [1/14,13/14]`. The blocker
subtournament is transitive because its arc has length `1/7<1/2`, and across the cut

`s -> b` if and only if `s>b+1/2`.

The cut neighborhoods are therefore nested: its adjacency matrix is Ferrers. If
`tau_0^mix(t)` counts cyclic triangles through the observer in this **mixed** tournament, then
time reversal keeps the threshold incidences but reverses each moving-runner edge, so

`tau_0^mix(t)+tau_0^mix(-t)=B(t)(13-B(t))` almost everywhere.

Since each danger set has measure `1/7`, integration gives the exact bridge

`E tau_0^mix = 78/7 - sum_{i<j} mu(D_i intersect D_j)`,

equivalently

`E B^2 = 169/7 - 2 E tau_0^mix`.

This identifies the classical pair-correlation/second-moment calculation with a marked
tournament triangle observable. It does not by itself give the needed inequality: the Ferrers
area can vary from zero to the whole rectangle, and the missing information is the arithmetic
handoff/monodromy of blockers through time.

For the **pure** half-turn tournament, time reversal reverses every arc and preserves cyclic
triangle count; the displayed identity is false there.

## 6. Combined hard object and next tests

The useful state is not `speeds <= W`, nor a binary tournament. It is a packet containing

`(normalized core shape, scale residue, exceptional offsets,`

` prime-grid cover up to common multiplication, private-owner fibers,`

` sector potential c, microphase weak order, lift depth)`.

HYP-6780 supplies the scale action; THM-760 removes the one-exception common-factor lane;
HYP-6785 supplies the endogenous pair-sum obligation/blocker complex. The prime-grid cover
describes the published sieve's base layer, while the sector carry describes its continuous
wall movie. They meet at blocker ownership and handoff.

The most concrete next moves are:

1. canonicalize `I(13,p,1)` as cyclic set covers, then branch on the largest private owner;
2. record how private base fibers split under the `c=2` and `c=7` lifts forced by `14=2*7`;
3. search for an inverse theorem: low private mass or high overlap forces the AP/Goddyn-Wong
   cover orbit, while high private mass makes a lift die quickly;
4. implement the sector movie as `(c, weakorder(rho))`, not 156 edge labels, and measure the
   blocker-handoff monodromy of the exact survivors;
5. join this ledger to HYP-6785's pair-sum obligations and prove a scale-normal local-to-global
   blocker-pressure theorem.

These are proof-directed finite computations. None licenses promotion of a sampled raw cutoff
or a tournament correlation to a theorem about all speed sets.

## Source pointers

- `00-navigation/LRC14-FRONTIER-2026-07-15.md` (corrected status ledger)
- `05-knowledge/hypotheses/HYP-6780-lrc14-capped-band-scale-quotient.md`
- `01-canon/theorems/THM-760-coprime-exceptional-runner-sheet-dodge.md`
- `05-knowledge/hypotheses/HYP-6785-lrc14-endogenous-pair-sum-blocker-complex.md`
- `07-reflections/lrc-sector-vector-realizability-s540.md`
- `07-reflections/lrc-fixed-n-ferrers-interval-menu-s525.md`
- Sungkawichai--Trakulthongchai, arXiv:2604.23906
