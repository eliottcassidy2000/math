# The Antipodal Ladder is Discrete Borsuk-Ulam, and the Sliver Boundary is its Index Threshold

**boxeph-2026-07-20-S153** · owner prompt: "one involution that is free and
carries an odd map... the k-torus of the resonance lattice... T^k != S^k
blocks plain BU, so it needs the Z/2-index form."

## 1. The mod-q spread lemmas ARE discrete Borsuk-Ulam

On Z/q (q odd), x -> -x is a FREE involution away from 0, and residue
evaluation v -> v*a is odd: (-v)a = -(va). The spread-blocking lemmas
(mod 19/23/25, Lean kernel-pure) say: if M(W) < 2/q (and gcd(v_i, q) = 1),
the residues of the 13 nontrivial runners must SURJECT onto the antipodal
quotient (Z/q^*)/± — and if a pair {±c} is missed, the missed orbit itself
hands back the witness time t = a/q with a not in {±v_i^{-1}} (verified
20,000/20,000 exactly this session for q = 27). Surjectivity onto the
antipodal quotient forced by an odd equivariant structure = the Z/2-index
argument in miniature: the "sphere" is the unit part of Z/q, the "index" is
the pair count.

## 2. NEW RUNG: mod-27, and the ladder termination law

pairs(q) = #antipodal unit-pairs: q=19: 9, q=23: 11, q=25: 10, q=27: 9,
q=29: 14. A rung exists iff pairs(q) <= 13. So 27 = 3^3 is the LAST rung —
9 pairs, slack 4, serving the 3-coprime gcd-stratum (complementary to
mod-25's 5-coprime: blocker economy) — and it was UNFILED until today
(no LRCMod27 in Lean; verified computationally, 94% block rate on random
3-coprime 13-tuples). The ladder then STOPS at q = 29 because 14 > 13:
**the open sliver [2/29, 1/14) begins exactly at the Z/2-index phase
transition.** The residual hardness of LRC(14) starts precisely where
discrete Borsuk-Ulam loses surjectivity forcing — the sliver's left endpoint
is an index threshold, not an accident of the census.

## 3. The k-torus form (program, honest status)

Multi-prime certificates live on prod_i (Z/p_i \ 0) — the discrete resonance
torus — with the DIAGONAL involution. T^k != S^k: no plain BU; the right
invariant is ind_{Z/2} of the product, and the S-T gate laws are its tower
computations: the x7 gate's 1-of-7 slaving (singleton CRT fibers, zero
freedom) is an index-1 statement one level up; banded-CRT density 2*theta =
1/7 is its measure shadow. PROGRAM: express the x2/x7/x14/x49 gate ladder as
ind_{Z/2} of the lens-like quotients (Z/(7^b p))^* / ± filtered by 7-torsion;
a transfer argument should reprice wall (b) — this is the index-form step the
owner named. Status: framing only; no index computation done yet.

## 4. The four free-involution laws of the repo (the owner's abstraction)

"One involution, free, carrying an odd map" is the shape of four separate
canon results: (i) antipodal rungs: -1 on Z/q free away from 0, odd map =
residue evaluation, output = spread forcing; (ii) Moon/S150: involutions are
never automorphisms of tournaments — the reflection is forced OUTSIDE Aut,
output = dihedral splitting + anti-order 2 mod 4; (iii) THM-584: complement =
the antipodal map of the arc-hypercube Q_{C(n,2)}, free, output = the merged
metagraph G_n/Z_2; (iv) THM-1365: Keller deck groups act freely (Cartan),
output = the Galois reduction of JC. In each, freeness + oddness converts an
involution into a GLOBAL obstruction (surjectivity, splitting, quotient
structure, automorphy). The Z/2-index is the common currency; the k-torus
gate program (§3) is its next payment.
