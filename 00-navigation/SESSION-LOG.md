## boxeph-2026-07-21-S213 -- chirality + toothpick A139250 + tournament evenness = one Lefschetz parity (HYP-8850)

**Owner:** connect previous chirality work + toothpick A139250; tournament values that are EVEN (iso classes) and what toothpick diagrams correspond to their structure.

**ONE LAW:** a count graded by a Z/2 INVOLUTION obeys count=(#fixed)+2*(#free pairs) => count == #fixed (mod 2) [Lefschetz/Euler]. Master involution = reversal/converse R:T->T^op (=complement=antipodal, THM-584). Three realizations, all verified (chirality_toothpick_tournament_parity_boxeph_S213.py):

1. TOURNAMENTS EVEN: A000568=SC+2*(chiral pairs). THM-587 (PROVED): P_n(1)=A000568, P_n(-1)=SC=self-converse=antipodal Euler/Lefschetz number; SC=2,2,8,12,88,176 (n=3..8) all EVEN => A000568 EVEN for n>=3. Enumeration n<=6: (total2/SC2/0pairs),(4/2/1),(12/8/2),(56/12/22), count=SC+2pairs holds. THM-1830 BLUE self-comp atom=1 iff n odd = local single-fixed-point face.
2. TOOTHPICK A139250: sim=OEIS exactly; D4 mirror sym; =(#axis)+2*(#pairs); #axis ODD (central seed) => A139250 ODD all n>=1; first differences have exactly ONE odd term (seed); A139250(2^k)=(2^(2k+1)+1)/3=1,3,11,43,171 Jacobsthal.
3. LRC (S212): chi(G_delta) == #iota-fixed lonely pts (mod 2).

**TOOTHPICK<->CHIRALITY DICTIONARY:** central seed toothpick (axis-fixed) <-> achiral/self-converse classes; mirror toothpick PAIR <-> chiral pair {T,T^op}; D4 mirror <-> R=converse=complement=antipode; A139250 parity=#axis(mod2) <-> A000568 parity=SC=P_n(-1)(mod2); doubling A139250(2^k) <-> level-graded signed-cycle-index recursion. The toothpick fractal = the PICTURE of 'one achiral seed generating chiral mirror-pairs'. Unifies with S212: (P(+1),P(-1)) = Burnside count / Lefschetz fixed set = the chirality Euler class (S212 even chi + odd Gauss sum i*sqrt7).

**Honest:** THM-587, the toothpick OEIS/oddness/Jacobsthal facts, and S212 are each proved/verified; the contribution is the UNIFICATION (one reversal, one Lefschetz parity count==#fixed) + the explicit toothpick<->chirality dictionary. Artifacts: reflection chirality-toothpicks-and-why-tournament-counts-are-even-one-lefschetz-parity-boxeph-S213.md, HYP-8850, script (+.out).

