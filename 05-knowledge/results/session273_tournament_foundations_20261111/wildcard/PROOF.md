# Arithmetic clock tournaments: local resonance, entropy, and a faithful lift

**Session273, 2026-11-11. Status: PROVED elementary statements below; FINITE-EXACT controls in `controls.json`. No LRC or Chan–Pak closure.** No literature novelty is claimed for the underlying square-wave or symbolic-dynamics tools. The contribution is their explicit composition into this tournament observer, including quantitative loss and recovery.

## Inheritance and the unexpected route

The starting point was not a recent Chan–Pak attack. A randomly selected old reflection, `07-reflections/three-imports-one-resonance-and-the-rigidity-is-diophantine.md`, places alternating projections, Christoffel kernels and circle dynamics around the same LRC feasibility problem. Its claims are inspiration, not dependencies. A second sampled reflection, `07-reflections/representation-lens-fundamentals-2026-05-30.md`, suggests sign assignments as common carriers. The precise usable mechanisms are:

* `01-canon/theorems/THM-374-half-turn-transitive-semicircle-criterion.md`: half-turn tournaments detect an empty semicircle, not an arbitrary loneliness margin.
* `01-canon/theorems/THM-4106-lrc-pair-owner-reciprocity-and-valuation-tree-decoder.md`: square-wave correlations retain primitive arithmetic ratios even when bare order is blind.
* `01-canon/theorems/THM-4602-tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips.md`: real rank-two determinant sign tournaments are exactly locally transitive tournaments (known result credited there). Its finite-field Legendre counterpart leaves the real chamber.
* Hostile: `01-canon/theorems/THM-383-lrc-half-turn-boundary-compactification-audit.md` shows that arbitrary tie completion changes the class menu; `THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density.md` warns that arithmetic type/order does not retain phase arrival.

The least-used sidecar is **the full dyadic time-address**, rather than the single tournament's H-value. This connects the 2-adic valuation controlling local correlations to ordinary base-two phase reconstruction, without identifying the two notions.

## 1. Intrinsic tournament, not a cosmetic completion

Take distinct integer speeds `v_0<...<v_(n-1)`. At generic `t in R/Z`, put

`i -> j iff sin(2 pi (v_j-v_i)t)>0`.

This is the half-turn orientation on phase points `{v_i t}`. Pair collisions and antipodal pairs are exactly the walls where this sine vanishes. They have Lebesgue measure zero. We leave them **unoriented**, rather than making an arbitrary tournament there. All probability statements below use uniform `t`. The displayed labelled orientation is the gauge; relabelling and converse remain actual symmetries.

## 2. Exact three-vertex law and sharp arithmetic bias

For three increasing speeds, let adjacent gaps be `a,b>0`. Set

`X=sgn sin(2 pi a t), Y=sgn sin(2 pi b t), Z=sgn sin(2 pi(a+b)t)`.

For positive integers `r,s`, define

`C(r,s)=gcd(r,s)^2/(rs)` if `r/gcd(r,s)` and `s/gcd(r,s)` are both odd, and zero otherwise.

Then for each `x,y,z in {-1,1}`,

`P(X=x,Y=y,Z=z) = (1 + xy C(a,b) + xz C(a,a+b) + yz C(b,a+b))/8`. (1)

**Proof.** The inherited square-wave Fourier calculation gives `E[sgn sin(2 pi r t) sgn sin(2 pi s t)]=C(r,s)`. The substitution `t -> 1-t` negates all three signs, so all odd moments vanish. Expand `(1+xX)(1+yY)(1+zZ)/8` and take expectation. This proves (1), not an independence approximation.

Write `a=dp,b=dq`, with `gcd(p,q)=1`. Exactly one pair-correlation survives:

* `p,q` odd: `C(a,b)=1/(pq)`;
* `p` even: `C(b,a+b)=1/[q(p+q)]`;
* `q` even: `C(a,a+b)=1/[p(p+q)]`.

Thus the exact labelled three-edge distribution differs from the uniform tournament cube in a **single quadratic Walsh character**. Its total-variation distance is half the surviving correlation. The unique pair is also the pair among `a,b,a+b` having equal 2-adic valuation. The law is never exactly uniform at finite integer gaps, but can approach uniform when that primitive denominator grows.

A directed triangle occurs precisely when `X=Y=-Z`. Consequently

```
P(cyclic) = (1+1/(pq))/4,                 p,q odd;
           (1-1/[q(p+q)])/4,             p even;
           (1-1/[p(p+q)])/4,             q even.
```

In particular

`1/6 <= P(cyclic) <= 1/2`,

with lower equality exactly for `a:b=2:1` or `1:2`, and upper equality exactly for `a=b`. Proof: the denominator in the negative cases is at least three, with equality only `(p,q)=(2,1)` or `(1,2)`; in the positive case `pq>=1`, equality only `(1,1)`.

For the induced three-vertex H-count, `H=1+2 1_cyclic`. Hence

`4/3 <= E H <= 2`,

and `E H` is never exactly `3/2`, the uniform labelled-tournament mean, at finite integer gaps. These are statements about a clock ensemble, not the full path-count spectrum. For any finite speed set the expected number of directed triangles is the sum of these explicit gap laws. No independence across triples is asserted.

## 3. Rank-two compression and an A000568 capacity bound

For arbitrary distinct nonantipodal phase points, choose a phase as origin and fold all points modulo `1/2`. Order their folded coordinates

`0=u_0<u_1<...<u_(n-1)<1/2`

and write original phases as `u_i+s_i/2 mod 1`, with `s_0=0` and `s_i in {0,1}`. For `i<j` the orientation is

`i -> j iff s_i=s_j`. (2)

This follows directly by subtracting phases. It is exactly a cut-switch of a transitive order. Conversely every such sign word has a realization by taking any strictly increasing folded coordinates. This recovers the inherited real rank-two/chirotope carrier from an explicit binary compiler; it is not a new rank-two classification theorem.

Two useful consequences:

1. The number `L_n` of isomorphism classes obtainable from half-turn circles obeys `L_n <= 2^(n-1)`.
2. Since `A000568(n) >= 2^(n(n-1)/2)/n!`,

`L_n/A000568(n) <= n! 2^(n-1-n(n-1)/2) -> 0`. (3)

The bound follows by assigning a binary word to each class, not by incorrectly dividing a labelled count by `n!`. Multiple words may yield the same class. This is an all-order capacity firewall: a real planar rank-two phase model is a vanishing fraction of tournament classes, even if every arithmetic speed row is allowed. It is useful precisely because the diagram with the finite-field Legendre construction does **not** commute with real positivity. A finite-field Paley pattern is not automatically a real planar determinant sign pattern.

This does not say that a sparse family cannot cover the odd H-spectrum. Class completeness and value completeness are different obligations. Nor does it say arbitrary planar JC or unit-distance graphs are rank-two determinant sign problems.

## 4. A fixed arithmetic row has exactly one bit per-scale entropy

Fix `V` with at least two speeds. Let `W_K(t)` be the labelled tournament word at scales `t,2t,...,2^(K-1)t`, at generic times. Write `N_K(V)` for the number of distinct such words. Set

`D={|v_i-v_j|:i<j}`, `G_d={j/d mod 1:0<=j<d}`,

`c(V)=| union_(d in D) G_d |`.

Then

`2^K <= N_K(V) <= c(V) 2^K`, and `lim_(K->infinity) log N_K(V)/K = log 2`. (4)

**Proof.** All walls through level `K-1` lie in `union_d G_(2^K d)`: finest-level grids contain all coarser grids. This union is the preimage of `union_d G_d` under multiplication by `2^K`, so it has exactly `c(V)2^K` points, hence that many complementary cells. The word is constant on each cell. For a lower bound, choose any edge of frequency `d`. Its signs at the K scales give the first K binary digits of `{dt}`. All `2^K` words occur on open intervals. This proves (4).

The finite alphabet of tournaments can be quite complicated, but for this fixed clock its allowable histories have only one bit of entropy per scale. Replacing them by independent A000568-valued symbols adds histories that no single clock can realize. This is a concrete symbolic-dynamics obstruction, not a general obstruction to automaton approaches to Chan–Pak.

### Exact finite-state dynamics, not an independent-symbol approximation

There is a computable finite stationary hidden Markov model for the **entire** scale-indexed tournament process. Take as states the `2c(V)` open cells of `P=union_d G_(2d)`, and let `ell_i` be the length of cell i. This endpoint set maps into itself under doubling. Each cell has length at most `1/2`, so doubling is injective on its interior and its image is a union of complete state interiors. Put `A_ij=1` exactly when that image contains state j. Label state i by its tournament.

Then

`sum_j A_ij ell_j=2 ell_i`, `sum_i A_ij=2`,

and the transition matrix and stationary distribution are

`P_ij=A_ij ell_j/(2 ell_i)`, `pi_i=ell_i`. (5)

The first identity is image length; the second counts the two generic preimages of a circle point. Hence P is stochastic and pi is stationary. A specified admissible state path of length K is an interval of length `ell_last/2^(K-1)`: each inverse branch divides lengths by two. Summing over paths with the desired tournament labels gives the **exact joint probability of every finite tournament word**. This also proves directly that every finite accepted path has a genuine open interval of clocks realizing it. It is not a relaxed transition graph.

Thus exact multitime enumeration and correlation prediction for a fixed arithmetic row reduce to rational matrix products. Collapsing states with the same instantaneous tournament or H-value is allowed as a hidden observation, but need not produce a Markov chain on those observations. This is the precise interface lesson shared with the repository's cellular-automaton and Hamiltonian-boundary work: a repeated scalar does not certify the same future.

A sharp physical hostile is `V=(0,1,3,7)`. On the entire cell `(1/14,1/12)`, `H(T_t)=5` and `H(T_(2t))=1`. On `(1/12,1/8)`, `H(T_t)=5` and `H(T_(2t))=5`. Both present tournaments are in the unique strong four-vertex isomorphism class. Thus even the present isomorphism class loses the next count. These rows prove failure of strong lumpability by H (and by class) of the physical-cell chain. They do not by themselves prove that its stationary H-process is not weakly Markov.

The state bound depends on the speed differences, potentially enormously. It is not a size-only finite controller for all LRC rows. The result is a standard Markov-partition mechanism applied to this specific arithmetic tournament observation; no general new symbolic-dynamics theorem is claimed.

## 5. Restore phase instead of asking the scalar H-observer to see it

Include speed zero, and let moving speeds be nonzero. The sign of edge `0 -> i` at scale `2^k` gives the `(k+1)`st binary digit of `{v_i t}`: negative sign means bit one. Therefore the full infinite labelled, scale-indexed root-edge word determines every phase. The kernel of `t -> ({v_i t})_i` has exactly `g=gcd(v_i)` elements, by Bezout; if `g=1`, the infinite observation determines `t` itself modulo one.

At finite depth K, decode `m_i=floor(2^K {v_i t})` and put `D=2^K`. The phase lies in the dyadic interval `I_i=[m_i/D,(m_i+1)/D]`. Since distance to the nearest integer is 1-Lipschitz,

```
L_K = min_i min_{x in I_i} ||x||,
U_K = min_i max_{x in I_i} ||x||,
L_K <= min_i ||v_i t|| <= U_K,
0 <= U_K-L_K <= 2^(-K).
```

So every strict LRC witness is certified at a finite scale. An exact equality witness need not be: for the single speed one at threshold `1/2`, the true maximum is `1/2`, but no dyadic interval lower bound reaches it. At wall times use exact rational phases or explicit boundary symbols, not arbitrary signs. This sharp equality hostile is the missing compactification sidecar.

Thus the bridge is operational: **labelled multiscale tournament observations reconstruct the actual margin**, while H, an unlabelled class, or a single time slice generally does not. It gives a faithful search/certificate representation, not a termination theorem, not an effective uniform depth, and not LRC(14).

## 6. Connection contract and useful next questions

| Source -> target | Map / preserved data | Loss / required repair | Next decisive obligation |
|---|---|---|---|
| Integer gap triple -> tournament probability | (1), full labelled law | scalar H loses which valuation pair resonates | use marked correlations, not H alone |
| Real planar frames -> tournament classes | determinant signs / (2) | lengths and most of A000568 disappear | add higher-rank or finite-field carriers only with an actual predicate transport |
| Fixed speed row -> symbolic tournament histories | scale doubling / (4) | independent-symbol model invents unattainable histories | characterize the allowed language, including equality-wall symbols |
| Multiscale rooted tournaments -> LRC phase | binary digit decoder | forgetting labels, scales or root destroys reconstruction | control the depth or the equality stratum uniformly |
| Chan–Pak finite-state constructions -> arithmetic clock models | ANALOGY: both are constrained positive languages | no map yet preserves exact H target and clock arrival | build such a map or retain the analogy label |

Promising **OPEN** question: can a small boundary-augmented language capture every equality witness of an LRC extremal row while retaining a finite-depth certificate for all nonextremal rows? This requires a genuine uniform theorem, not just one-bit entropy. For Chan–Pak, the corresponding design principle is to retain the entire continuation address before collapsing to a count; it does not supply an unrestricted H-selector.

## Reproduction

`python3 05-knowledge/results/session273_tournament_foundations_20261111/wildcard/clock_lift.py`

The script integrates over exact rational event cells, independently of Fourier formulas; checks all 900 adjacent gap pairs `1<=a,b<=30`, equality conditions, 24 exact finite-state language replays and growth sandwiches on four speed sets, Markov length/stationarity identities, and 2400 phase/margin controls. It uses no floating-point trigonometry and rejects ties. These controls support, but do not replace, the all-order proofs above.
