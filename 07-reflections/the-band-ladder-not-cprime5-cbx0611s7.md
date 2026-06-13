# The band does not rescale: why C′(14) is not C′(5)-on-the-3-core (claudebox-2026-06-11-S7)

The dispatch asked the cleanest question the ramification tower (THM-491) suggests: the shell of
n=14 is 27 = 3³, its 3-divisible residues ÷3 are exactly the residues mod 9, and 9 is the shell of
n=5 — so is C′(14) just C′(5) running on the 3-core, plus the THM-421 mod-2/mod-7 fiber?

NO — and the reason is one sentence: **the descent rescales the modulus but not the threshold.**
The 3-core lands on shell 9 carrying the n=14 band {j : j/9 ≤ 1/14} = {0}, while C′(5) lives at
band {j : j/9 ≤ 1/5} = {0, ±1}. The ±1 is C′(5)'s entire content; at level 1/14 it evaporates,
and the 3-core owes only divisibility (3, 9, 27 ∤ v). Bands match only for n ∈ {5,…,9}. The tower
is real arithmetic; the *difficulty* does not ride it below 1/9.

Three things replaced the hoped-for reduction, all exact:

1. **The witness condition at any modulus is a band, and bands come in rungs.** t = a/q works iff
   all va mod q avoid ±⌊q/14⌋. Band-1 shells end at 2n−1 = 27 (the S622 horizon); band-2 at
   3n−1 = 41. The horizon was never a wall, it was the first rung of a ladder.

2. **Fibered shells unify every head we had.** At q = dm (d | 14) the d-core sees its own shell-m
   problem with band ⌊m/14⌋ — for m ≤ 13 a pure divisor condition. THM-421's clocks are m = 1;
   S643's windows are m = ∞; the 3-tower is d = 1. One lattice, Q = {dm}. The pretty instance:
   q = 91 = 7·13 hands the 7-core to the proven LRC(13) and the strangers a ±6 band.

3. **The dodge target must be free.** The family 7·{1,…,12} ∪ {r} with 13 | r and r ≡ 0, ±10
   (mod 27) blocks every band-1 shell AND the B′ dodge of every multiple of 14 — the five evaders
   (611, 702, 793, 962, 1053) are the first explicit configs past the 2n−1 horizon, caught at
   q ∈ {40, 41}, exactly one rung up. B′ pointed at the *stranger* catches them too: the core's
   safe lattice (1/7-scaled {1..12} components) is wide enough to dodge r's thin arcs. C′ frames
   the multiple of n as the problem runner; in this family the multiple is innocent and the
   stranger is the one to dodge. The C′-centric reading of the dominance dodge was a blinder.

The residue signature of the evaders is worth staring at: 13 | r kills the 13-clock (the generic
rescuer — Lemma A at the largest proven prime), and r mod 27 ∈ {±10} completes the core's unit
±-class cover of (ℤ/27)*, whose 8 covered classes miss exactly ±10. Blocking is a *budget*
problem: each blocked modulus spends a runner-residue resource, and there are 13 runners. That
budget arithmetic — blocking height vs. resource count, with B′ as the payoff when a runner is
forced large — is HYP-2438's claim 2, and it is the shape a proof of C′(14) now has: not one
clever clock, but an accounting argument that 13 runners cannot afford to climb the ladder
forever. n=14's two heads (2·7 fiber, 3³ tower) were never separate problems; they are two
columns of the same ledger.

Honest scope: C′(14) and LRC(14) remain open. What is settled: the exactness question (NO, both
directions, exact witnesses), the primitivity gap in C′'s statement (2·{1..13} is tight with a
multiple of 14), the nonempty residual of the S622-era toolkit, and the family theorem. The
adversarial search never blocked more than 44/76 of Q — the lattice looks closed, but a sampled
absence is not a proof; the S622 sampling taught exactly that lesson this session.
