# Independent hostile audit: arithmetic and smooth witness bank

Status: **ACCEPT arithmetic; ACCEPT physical seed verification; ACCEPT all-exponent smooth-sector descent and fixed-bank cost theorem.** Chan–Pak remains OPEN. This audit owns only this directory. Graft obstruction is independently reviewed by a separate audit lane.

## Arithmetic verdict

The proof that every common integer denominator d of a positive rational square AP of difference 5 satisfies 12|d is valid without any primitive/reduced hypothesis. Independent enumeration of square residue triples gives d≡0 mod3 and d≡0,4,8,12 mod16. Thus the displayed d=12 is minimal. At d=12, exhaust all equal-parity divisor pairs (b−a,b+a) of720; all nine positive pairs are tested and exactly (31,41,49) survives. This exhausts unbounded numerators, rather than truncating a search.

The general squarefree-N denominator assertion follows correctly from step divisibility by8 and3. Squarefreeness is essential to the stated even-N conclusion.

The signed map identity Q=−2P is correct. The independent chord-law computation verifies P=(25/4,75/8), Q=(1681/144,62279/1728) and the curve equation. The producer's general tangent slope and doubled-coordinate formulas are algebraically valid.

The two positivity charts must remain distinct. Direct square-AP clique labels are ((X−10)/2,X/2,(X+10)/2), so Q is positive because X>10. The inverse triangle map at Q followed by the triangle-to-clique construction fails its different inequality (X²+25)²>40X(X²−25). There is no contradiction. Independent computation additionally shows 2Q has x=11183412793921/2234116132416<10. A vague assertion that Q itself exits the positive chamber without identifying the chart would be false; the final producer text avoids it.

The halving obstruction is rigorous and stronger than merely rejecting one chosen half. If Q=4R, then 2R=−P+T for one of all four rational two-torsion points. Independent chord addition gives their x-coordinates 25/4,−4,45,−5/9 and verifies every point doubles to Q. None has x−5,x,x+5 all rational squares. Every finite rational double with nonzero y has this necessary property by the displayed duplication identities, so Q is in2E but not4E. No rank, generator, or converse two-descent theorem is assumed.

Berggren reversal leaves area class5 and the positive graph-label chamber at the first step. Elliptic doubling preserves the curve but increases denominator valuations in this example. Neither is an infinite descent on the same forbidden predicate. The note correctly states that deleting a tournament already supplied cannot descend a hypothetical missing value: no tournament exists at that missing target to delete.

## Physical finite seed verdict

The supplied graph masks for all ten counts 3,5,31,35,41,49,63,217,287,343 are independently verified to be loop-free complete tournaments. A separately written, memoized start-vertex path recursion computes their Hamiltonian-path counts, independently of the producer's forward endpoint DP and scalar word recurrence. In particular the order10 seeds really have counts217 and287. This is a physical witness certificate, not just arithmetic support.

The independent finite semigroup check tests every odd {3,5,7,31,41}-smooth target through100000 and obtains precisely the allowed targets, excluding7 and21. This is FINITE-EXACT control; the all-exponent claim needs the elementary valuation proof, not extrapolation.

## Novelty and quantifier boundary

The prime/7p equivalence is already proved in THM-4094, `01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`, §6. Do not promote it as new. The finite basis closure is a specialization with explicit suppliers. A minimum order among joins of the named bank is not a minimum over all tournament witnesses.

A valid missing-target descent requires admissible smaller targets plus a constructor whose required rooted states are supplied, or a strengthened induction predicate which supplies those states. A finite inverse terminating on given graphs only tests membership in its construction class. It does not prove nonempty inverse fibers for arbitrary desired targets. The final arithmetic note respects these distinctions.

## Reproduction

Run from the repository root:

    python3 -B 05-knowledge/results/session242_congruent_descent_20261110/audit/verify.py
    python3 -B -O 05-knowledge/results/session242_congruent_descent_20261110/audit/verify.py

Output is `audit/verify.json`. The audit reads only the seed masks from `descent/CHECK.json`; it imports no producer code. All checks are explicit exceptions, not removable assertions. The arithmetic and descent producer scripts were additionally replayed in normal and optimized modes and byte-matched their stored outputs.

## All-exponent bank proof verdict

The supplied `descent/PROOF.md` §§2–4 pass hostile audit. Its universal join lift needs only the smaller target count, not a special rooted realization, so it genuinely preserves the failure predicate on the declared smooth sector. Every residual in the five-case descent remains in that sector and strictly decreases. The terminal1 is realized. This is stronger logically than testing inverse branches of supplied witnesses, but only inside the declared sector.

The fixed-bank cost formula accounts for all generator multiplicities, including the ordinary-prime costs consumed by63,35,217,287. The pure-seven carry minimization is valid: minimizing8v+11t with2v+3t=c is equivalent to minimizing v, whose least feasible value is0,1,2 for residues0,2,1. The c=1 exception is explicitly retained. An independent exhaustive carry-cost check covers0≤c≤200. For c≥K+2 every residual is at least2, so every candidate shifts by exactly11 under c→c+3; taking minima preserves that equality. No global minimum tournament-order claim follows or is made.

One documentary correction was sent to the parent: §0 initially cited THM-4094 with a nonexistent trailing `-reduction` in the filename. The actual path is given above. No mathematical blocker remains. Graft exclusion/sharpness is outside this audit and is assigned separately.
