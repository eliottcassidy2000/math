# Erdős Problems Through the Doubling Lens: the Anti-Counterexample Operator

**Source:** kind-pasteur-2026-06-09-S2 (THM-455/456/457, MISTAKE-068/069, T769,
HYP-2356..2359 resolved/advanced)
**Seed:** the human's directive — explore the repo's Erdős-problem work, exploit the vast pool
of connections, think outside the box.

## The box, and the outside of it

The repo's Erdős threads (64, 592, 625, u(22), LRC) had been approached through additive
machinery — Sidon ladders, covering systems, graded relations. Yesterday's skew-Sylvester
doubling (THM-447..454) looked like pure tournament algebra. The out-of-the-box move was to
notice that BOTH live on the same 2-adic seam, and to run the doubling machinery AT two Erdős
problems directly. Both times the doubling turned out to be the **anti-counterexample
operator** — the structure whose presence makes the Erdős statement trivially true, so that
counterexamples are precisely the structures that exclude it.

## Erdős 64: counterexamples are the twin-free, near-Turán-extremal core

A single edge of G plants a twin C₄ in the blowup G[K₂]; a k-cycle spreads into the FULL
interval [k, 2k] (every dyadic window). So no blowup is ever an Erdős–Gyárfás counterexample
— and quantitatively, a counterexample must be C₄-free (no two vertices with two common
neighbors: twin-free in the strongest sense) while δ ≥ 3 pushes its edge count INTO the
C₄-Turán ceiling. The corridor ex(n;C₄) − ⌈3n/2⌉ is NEGATIVE through n=9 (E-G vacuously true)
and opens by only +1, +1, +3 at n=10,11,12 — counterexamples live within O(1) edges of
C₄-extremality (THM-456). The exhaustive census found all 71 C₄-free δ≥3 graphs at n=10–12
uniformly killed by a forced C₈.

Above that, the hunt (THM-457) found the **dyadic gate ladder**: suppressing C₈ inflates C₁₆
(min count 614→970 across n=28–40); achieving {4,8,16}-freedom (provably impossible below
n=54, first realized at Exoo's G78) immediately hits C₃₂ — we verified the reconstructed
G78-type graphs contain C₃₂, new data beyond Exoo 2014. Each gate you close pressurizes the
next: the conjecture defends in depth. New exact rungs: girth-5 cubic C₈-free starts at
EXACTLY n=28; and the dihedral 3-reflection Cayley family yields specimens with dyadic
spectrum exactly {4, 32} — an extreme gap family nobody seems to have catalogued.

The parity-mixer reading: the OCF taught us (THM-454) that twin insertion converts even
cycles to odd. Here the same insertion converts ANY cycle into a full length-interval. The
doubling doesn't just mix parity — it mixes all residues, which is why its absence is the
precondition for lacunary (dyadic-avoiding) spectra. The lacunary extreme — theta graphs,
spectrum = 3-element sumset — is excluded from E-G precisely by the degree-3 floor:
**the entire difficulty of Erdős 64 lives in the gap between "no degree-2 vertices" and
"contains a twin-like structure".**

## Erdős–Moser: the tower's extremality window, and the third repetition of one pattern

trans(D(T)) ≥ trans(T) + 1 always (append the source's twin); equality in 17/18 small classes;
ONE exception where an alternating primed/unprimed chain gains +2 — the interleaving
obstruction binds chains, not mixed alternations. The Mersenne tower hits trans = 3, 5, 7 at
T7, T15, T31: extremal at 7 (Paley T₇, THE largest TT₄-free tournament), pointwise-minimal at
15, and EQUAL to Paley_31 at 31 despite being a different tournament. Then T63 jumps to 11
and the window closes.

This is now the THIRD time the same pattern has appeared, in three unrelated lenses:
- isomorphism: tower = Paley at 7, leaves at 31 (THM-448)
- Hadamard equivalence: tower = Sylvester through order 8, leaves at 16 (THM-451)
- extremal combinatorics: tower f-extremal through 31, leaves by 63 (THM-455)

Each lens sees the tower track the classical extremal object for a few doublings and then
depart — while |Aut| stays frozen at F₂₁ throughout. A reflection-worthy conjecture: the
tower is the "generic skew object with Paley initial conditions"; every classical coincidence
it inherits has a finite half-life, and the half-lives are ordered (iso dies at 31, Hadamard
at 16... in matrix order: 16 then 31+1=32 then 64 — the half-lives are CONSECUTIVE POWERS
OF 2: equivalence dies at 16, isomorphism at 32, extremality at 64). That ladder deserves
a sequel session: predict which property dies at 128.

## Corrections as connective tissue

Two mistakes entered canon this session and both teach the same lesson from opposite sides:
MISTAKE-068 (cycle-anchored DP reused for paths — rotation symmetry mistaken for path
symmetry) and MISTAKE-069 (first-in-enumeration mistaken for smallest — McGee has 34
eight-cycles, correcting S710's premise and re-grading the girth ladder). Both are
order-of-discovery artifacts. In a repo where agents inherit each other's claims, recording
the FULL profile (spectrum, not first hit) is what keeps the ladder sound.

## What to do next (ranked)

1. **The half-life ladder:** formalize "tower property half-lives at 16/32/64" and test the
   prediction at 128 (which classical coincidence survives to order 128? SNF? GF(2) rank?).
2. **HYP-2356 closure:** literature-confirm st(5)=13, st(6)=27, st(7) bounds; if f(31)=6,
   find which 31-vertex tournament beats trans 7 (the TT₇-free sub of the composition).
3. **The {4,32} dihedral family** (THM-457): characterize WHY refl(0,j,m/2) kills 8 and 16
   but not 4 and 32 — likely a 2-adic character condition; connects to the Sidon/B_h ladder.
4. **The alternating-chain exception** (THM-455(2)): classify which T admit the +2 gain;
   the answer is a tournament substructure dual to the interleaving obstruction.
5. **Erdős 592 bridge:** hand mac-mini the trans data (their R(n,2) witnesses are
   "transitive-leaning" per their checkpoint — our +1 law and exception may pattern-match
   their dyadic witnesses).
