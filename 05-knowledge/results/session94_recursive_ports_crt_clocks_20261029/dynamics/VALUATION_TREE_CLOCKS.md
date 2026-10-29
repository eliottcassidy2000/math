# Exact arithmetic-source coding of the critical ternary tree

Session94, 2026-10-29. **PROVED (elementary, scoped)**; independent finite controls below. No priority claim for geometric valuation laws, Galton–Watson trees, or total-progeny formulas. The contribution here is their explicit source/clock coupling, recovery sidecars, and hostile boundaries in this repository's carriers. No Collatz convergence, Antihydra halting, or ROOT-forcing theorem follows.

## 0. Inheritance, source pin, and portfolio

Repository source pin: `88c969273f2b6bc4ea2a8c5d29bb36c558eb59d2`.

- Closest arithmetic mechanism: `05-knowledge/results/collatz_renewal_source_20261011.md`, R1, exact native cylinders with mass `2^-A`. This is inherited, not newly discovered here.
- Closest group mechanism: `05-knowledge/results/lucas_prime_reflection_towers_20261028/PORT_MOMENTS.md`, section3: uniform `Aut(W_k)` fixed points have the critical offspring law `(2+s^3)/3`. Its warning against an unqualified Collatz transfer remains correct.
- Antihydra rule and returned counter: `01-canon/theorems/THM-4771-odd-orbit-and-countdown-forms-of-bb6-cryptids-antihydra-a4-a11-a2-a5-a10.md`, part(a), and `05-knowledge/results/fixed_twists_conductors_reverse_potential_20261025/antihydra/ANTIHYDRA.md`. The latter's stale *all odd Z2* wording was already repaired by `session64_lca_reverse_valuation_20261026/session64-reverse-valuation/REPORT.md`, section2, with independent audit in `session64-audit-dynamics/REPORT.md`, and by `lca_ports_reverse_20261025/arithmetic/README.md`, section3. Section1 inherits those repairs; it does not claim to discover them.
- Ordered packet and first-hit guard: `05-knowledge/results/eleven_squares_routing_bridges_20261007.md`, sections3–5, and `collatz_carry_interfaces_20261004.md`.
- Golden period-five sources: `05-knowledge/results/session74_eleven_mechanisms_20261027/session74-clocks/RESULT.md`, especially “Elevenfold hostile”; its independent audit is the sibling `session74-golden-audit/AUDIT.md`.
- Actual ROOT caution: `05-knowledge/results/session75_overnight_black_root_spectral_20261027/session75-root-transport/REPORT.md`: native source measure and operator domain cannot be replaced by ambient Haar.

Anchor: an exact bridge from arithmetic sources to the critical ternary law. Niche: residual clock marks and named-root sensitivity. Wildcard: optimal variable-precision sampling rather than a nonexistent bounded-bit bridge. Concept board: **native cylinder / breadth-first root / residual gap / carry order / stopping guard / source measure**. Canonical hostile: the positive fixed source1. Corrected near miss: an almost-sure Haar statement is not a theorem on a prescribed positive integer. Least-used sidecar: the unused infinite valuation tail.

## 1. Three maps, one valuation law; a necessary puncture

On odd 2-adic sources define, whenever the numerator is nonzero,

```
S_e(x) = (3x+e)/2^a,       a = v2(3x+e), e in {+1,-1};
A(x)   = 3^(a-1)(3x-1)/2^a, a = v2(3x-1).
```

`A` is Antihydra's unrestricted odd-return map, before its counter guard is imposed. Negation conjugates `S_+` to `S_-`, preserving the complete valuation word. It does not preserve the positive half-line. Negation does **not** conjugate `A` directly to `S_+`: the multiplier `3^(a-1)` remains.

**Inherited domain repair (Session64 and LCA/ports/reverse).** The first map is undefined as a finite-valuation odd return at `x=-e/3`; the Antihydra map is undefined at `x=1/3`. These are odd 2-adic integers. Remove these points and all their finite backward iterates; call the remaining invariant domain `Omega`. The removed set is countable and Haar-null. Valuation coding is a homeomorphism from this *punctured*, noncompact space to positive-integer sequences. It cannot be a homeomorphism from all odd `Z2` to that noncompact sequence space, as asserted literally in the older Antihydra note's A3. The finite-cylinder calculation and all Haar laws survive the repair.

Proof without an appeal to compactness: the inverse branches, for odd `y`, are

```
I_e,a(y) = (2^a y-e)/3;
I_A,a(y) = (2^a y+3^(a-1))/3^a.
```

Each maps onto precisely one valuation-a cylinder, contracts 2-adic distances by `2^-a`, and makes the numerator's valuation exactly a because y is odd. Nested inverse branches determine one source for every infinite sequence of finite positive a's. Conversely every source in Omega has that sequence. Cylinders in either direction supply continuity. Each depth-r cylinder with `A_r=sum a_i` is a single odd residue class modulo `2^(A_r+1)`, hence has conditional odd-Haar mass `2^-A_r`. Thus the digits are iid with `P(a=j)=2^-j` on Omega. This is an exact source-measure assertion, not a random-orbit assumption about an integer.

A source packet is `(P,B,A,r)` with endpoint `(Px+B)/2^A` and native guard

```
x = r mod 2^(A+1),     r = (2^A-B) P^(-1) mod 2^(A+1).
```

Initialize `P=1,B=0,A=0`. Appending a digit a replaces

```
(P,B,A) by (pP, pB+b*2^A, A+a),
```

where `(p,b)=(3,e)` for `S_e`, and `(3^a,-3^(a-1))` for `A`. Exact oddness at the terminal point reconstructs the whole word by the inverse branches; merely reduced-denominator integrality is weaker and is not substituted for the guard.

## 2. The parity-of-valuation factor is exactly critical ternary branching

Set

```
b = 1[a is even],       G = floor((a-1)/2),       a = 1+b+2G.
```

Then b and G are independent, with

```
P(b=1)=1/3, P(b=0)=2/3;
P(G=g)=(3/4)(1/4)^g, g>=0.
```

Indeed the product is `2^-(1+b+2g)`. This factorization is coordinatewise and therefore preserves independence across successive digits.

Explore an ordered rooted tree breadth-first using the successive valuations. The current vertex receives three ordered children when b=1 and no child when b=0. Equivalently, the queue starts at1 and changes by `3b-1`. Stop at its first zero. Under odd Haar this is exactly the critical ternary Galton–Watson tree, with offspring `xi=3b`, mean1, variance2. Its generation-k population has exactly the same law as the number of fixed vertices of a uniform element of `Aut(W_k)` in the inherited port theorem.

This is the promised **explicit measure-preserving factor onto the fixed-point-tree distribution**. It is not an isomorphism of native arithmetic dynamics, not a full uniform automorphism sampler, and not an identification of arithmetic residues with tournament ports. In particular the two nonidentity macro rotations have been merged into the zero-offspring case. Group-valued holonomy is not recovered from this statistic.

For a rooted ordered finite tree, retain all G-marks in breadth-first order and the entire unused valuation tail. This restores every valuation digit, hence the original arithmetic source in Omega. Without these marks the factor forgets arithmetic information. The chronological breadth-first order is part of the inverse, not an optional layout convention.

The joint offspring/valuation generating function is

```
Phi(s,z) = E[s^xi z^a]
         = (z/2 + s^3 z^2/4)/(1-z^2/4).
```

At z=1 this becomes `(2+s^3)/3`. But the clock is not independent of offspring:

```
E[a | xi=0]=5/3;    E[a | xi=3]=8/3;
Cov(xi,a)=2/3.
```

The independent quantity is the residual gap G, not the whole valuation a. This is the precise repair for an unweighted branching analogy.

## 3. Exact all-order tree-size/division-clock recurrence

Let N be the number of processed vertices before extinction, and let D be the total valuation `sum a_i` at that stopping time. Use D here to distinguish the total division clock from the map A. Put

```
F(u,z)=E[u^N z^D].
```

It is the coefficientwise least solution with zero constant term of

```
F = u [2z/(4-z^2) + z^2/(4-z^2) F^3].                 (1)
```

This follows by independent root/child decomposition and is equally a formal power-series recurrence; no analytic convergence assumption is necessary for coefficients. The equation counts the recursively ordered tree while its recovered arithmetic word is read breadth-first. Reordering numerical instructions is not being performed.

The only possible sizes are `N=3m+1`. There are

```
C_m = binom(3m+1,m)/(3m+1)
```

ordered full ternary trees with m internal vertices. There are `2m+1` leaves. Their minimum total division clock is `4m+1`. Each vertex independently contributes an additional even amount `2G`. Consequently, for m,h>=0,

```
P(N=3m+1, D=4m+1+2h)
 = C_m binom(3m+h,h) 2^-(4m+1+2h).                  (2)
```

All other joint values have probability zero. One proof: fixed shape and weak composition of h among its `3m+1` vertices specify exactly one valuation word, hence one native cylinder of that mass. Alternatively use coefficient extraction in (1), a specialization of classical total-progeny/Lagrange inversion.

Summing h gives the familiar unclocked law

```
P(N=3m+1)=C_m (1/3)^m (2/3)^(2m+1).
```

Conditional on N, the tree shape is uniform among the C_m ordered shapes, and the gap marks stay iid with the above geometric law, independent of that shape. Thus

```
D = 4m+1 + 2 Z,     Z ~ NegativeBinomial(3m+1, success probability 3/4);
E[D | N] = 2N-1/3;  Var(D | N) = 16N/9.              (3)
```

Here Z counts failures before N successes. For one arithmetic odd step, the **standard** `3n+e` map consumes `a+1` steps, whereas the Terras map with its odd division included consumes a steps. Antihydra's floor(3n/2) map consumes a steps from one odd visit to the next. Therefore a stopped tree consumes `D+N` standard Collatz steps, D Terras steps, or D Antihydra steps. These clocks must not be interchanged.

Critical extinction is almost sure: the least fixed point in [0,1] of `q=(2+q^3)/3` is1. Nevertheless E[N] is infinite: each generation has mean1 and total progeny is their nonnegative sum. Also E[D] is infinite since D>=N. This does not contradict finite mean2 for an individual valuation digit; the stopping population has an infinite mean.

### 3.1 Finite generation clock recurrence

For height k, let `Z_k` be the population at generation k and let `D_k` be the sum of digits read at generations0,...,k-1 (the frontier generation is not read). Then

```
H_0(s,z)=s;
H_(k+1)(s,z)=[z/2+z^2 H_k(s,z)^3/4]/[1-z^2/4].       (4)
```

The exact identities `E[Z_k]=1`, `Var(Z_k)=2k`, and `E[D_k]=2k` follow. The same joint recurrence gives the exact mixed laws

```
Var(D_k) = (8/3) k^2(k-1)+2k;
Cov(Z_k,D_k) = (2/3) k(3k-2).                        (5)
```

For a direct proof, conditional on the first digit, write `Z_(k+1)=b sum Z_k^(i)` and `D_(k+1)=a+b sum D_k^(i)`, with three independent child pairs. Use `E[a]=2`, `E[a^2]=6`, `E[b]=1/3`, and `E[ab]=8/9`. Expanding products gives (5) by induction. The clocks thus have cubic variance even though the fixed-point population variance is merely linear. The nonzero covariance rules out attaching an independent clock after the unclocked population has been sampled.

Formula(4), not an independent clock attached after sampling, is the source-compatible refinement of the inherited fixed-point recursion. Independent copies can be used for the independent forest in the centralizer's fixed-count formula, but this requires separate source coordinates or an explicitly ordered parsing of independent remaining digits. It supplies no named-port transporter.

## 4. No bounded-precision free bit; optimal variable precision

Every event determined by a fixed number n of conditional binary source bits has dyadic probability with denominator `2^n`. Hence no such event has probability1/3. A bounded-precision exact factor to this branch law is impossible.

The valuation code reads fair bits until the first zero and returns the parity of its position. Its length has law `P(a=j)=2^-j`, expected2; even stopping positions have total mass1/3. These are **valuation-code bits**, not a claim that raw binary source digits can be read without the carry transformation. Native cylinders give exactly the same adaptive precision cost after the known odd bit.

This expected cost is optimal among deterministic exact samplers of a Bernoulli(1/3) variable driven by independent fair bits. At every depth n, at least one length-n prefix must remain unresolved; otherwise the output probability would be dyadic. Thus `P(L>n)>=2^-n`, and `E[L]>=sum_(n>=0)2^-n=2`. The first-zero parity sampler attains equality. Its unbounded worst-case cost is essential. No finite worst-case bit or arithmetic-step budget is hidden in the construction.

## 5. Hostiles: sources, roots, counters, and order

### Positive source1

`S_+(1)=1` with valuation2. Its tree is the full infinite ternary tree. The same tree is finite almost surely under Haar. A fixed positive integer is a Haar-null source, so the almost-sure statement cannot settle its arithmetic orbit or any other prescribed orbit.

### Density-one auxiliary extinction is not arithmetic descent

**Corollary (natural density among positive odd sources).** The set of positive odd `S_+` sources whose auxiliary tree becomes extinct has natural density1. For each finite clock cutoff M, the stopped words with total valuation at most M form a finite prefix-free family. Their native source cylinders are disjoint residue classes, so their natural density among positive odds is exactly their odd-Haar mass. As M increases, these masses tend to1 by almost-sure critical extinction. The full integer extinction set contains each such finite union; hence its lower natural density is at least every one of these masses and therefore is1. The upper density is at most1. This proof does not interchange an infinite union with natural density without a lower-bound argument.

This is a property of the **auxiliary queue**, not Collatz convergence, arrival at1, or even descent at the stopping time. The word `(2,1,1,1)` stops the tree after four vertices but has exact packet

```
x = 41 mod64;       S_+^4(x) = (81x+103)/32.
41 -> 31 -> 47 -> 71 -> 107.
```

The endpoint exceeds x for every positive source in that cylinder since `(81x+103)/32-x=(49x+103)/32>0`. Moreover that word contains an intermediate decrease, illustrating that tree stopping and first arithmetic descent are different stopping rules. Positive source1 still supplies a nonextinct exceptional source, consistent with density1 rather than universal coverage.

### The genuine golden-eleven connection needs the named root

The prior standard period-five family contains the rational cycle through `1/13` and the integer cycle through `-5`.

- The accelerated odd source `1/13` is fixed with valuation4, so its tree is infinite.
- `S_+(-5)=-7` has valuation1; `S_+(-7)=-5` has valuation2. Rooting the same orbit at -5 yields `(1,2)^infinity`, whose tree dies immediately. Rooting at -7 yields `(2,1)^infinity`, whose queue never empties (each pair increases it by1).

Both -5 and -7 belong to one standard period-five cycle, whose five golden labels lie in one nonzero orbit in `F11`. Thus the new factor is not invariant under moving the cyclic root and does not descend to an unmarked orbit. The finite-field eleven is an existing period-five carrier, not a reason to identify its uniform measure with Haar. No artificial twelve-fold counterpart is inferred.

### Antihydra halting is a different stop

At odd source27 the actual initial counter is6. After a valuation prefix its returned counter is

```
B=6+sum(2a_i-3).
```

This must be checked at each prefix, not just at the endpoint. It is unrelated to the queue being zero:

- word(1) kills the tree, while B=5 remains positive;
- the word consisting of six2's followed by twelve1's reaches B=0 for the first time while the tree queue is still1.

Each finite word has a nonempty native source cylinder, including positive integer representatives. These are hostile *source families*, not claims that the root27 orbit follows either word. Conditioning on survival changes distributions; the unconditioned iid law cannot silently replace the guarded machine.

### Equal clocks and even equal tree shape do not erase carries

For `S_+`, `(1,2)` gives `(9x+5)/8` on `x=11 mod16`; `(2,1)` gives `(9x+7)/8` on `x=9 mod16`. Same slope and division clock, different source guards and carries.

Even after retaining the tree shape and total clock, individual residual gaps matter. For example the one-internal-vertex tree has word parities `(even,odd,odd,odd)`; the words `(4,1,1,1)` and `(2,3,1,1)` shift a residual +2 from one vertex to another and preserve that shape and total clock, but change both the chronological affine carry and its source guard. Keep the ordered G-marks, not just their sum.

### ROOT and holonomy boundaries

Actual ROOT transport in the earlier work uses its native coordinate `X=3s+1` and positive affine carries; its measures and boundary guards are not S_± Haar. A common branch-count distribution does not transfer forcing bounds, source norms, source debt, first-hit membership, or integer arrival. Likewise the group loop-error theorem supplies legal re-marking paths in a groupoid of markings; our arithmetic tree factor does not identify these with arithmetic steps. Group law, named shared ports, native source, and clock remain separate coordinates.

## 6. Reproduction and external provenance

Run in this directory:

```
python3 check_valuation_tree.py
python3 -O check_valuation_tree.py
```

Both receipts agree byte for byte: **19,831 checks**. They cover all1,023 positive valuation words of total clock at most10 for all three maps, signed affine/guard checks, all stopped words of total clock at most16, and independent direct arithmetic enumeration of all4,096 odd source residues modulo8192 for each map. The stopped-tree mass within division clock12 is `3172/4096` for each map. A separate truncated bivariate Taylor implementation checks the generation-clock mixed laws through depth8. They check signs, the singular rational point, the ordered-carry hostile, the named-root golden hostile, and the distinct counter barrier. These finite tests verify formulas on declared universes; the proofs above establish all-order claims.

Primary literature checked2026-10-29:

- D. J. Bernstein and J. C. Lagarias, *The 3x+1 Conjugacy Map*, [author-hosted paper](https://websites.umich.edu/~lagarias/doc/bernstein.pdf). Its full binary-shift map uses the **Terras** convention `(3x+1)/2` on odds. Our induced odd-return cylinder proof is given explicitly above; the whole-space binary coding does not remove the induced-map exceptional set.
- Meyer Dwass, *The total progeny in a branching process and a related random walk* (1969), [publisher paper](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/B9505DF87F411C873D3AA511F21EBF8C/S0021900200026711a.pdf/total_progeny_in_a_branching_process_and_a_related_random_walk.pdf). Classical total progeny is acknowledged; the displayed clock refinement is proved here by exact cylinder counting.

**Open continuation:** use this decorated-tree carrier to transport specified arithmetic observables without erasing ordered carries, or prove that a given tournament-port observable factors through it. Neither a shared distribution nor equality of moment recurrences proves such factorization.
