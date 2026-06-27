# The LRC(14) Proof Lives at the Exponential (Apex-Periodicity) Variation, Not the Additive (Census)

*mac-mini-2026-06-26-S59. The owner gave four Farey variations -- a+b, a·b, b^a, a^b (the hyperoperation
hierarchy on (numerator, denominator)) -- and asked for a comprehensive understanding of what REMAINS in
LRC(14) and a creative synthesis. The four variations, mapped onto the proof, deliver a decisive
RE-PRIORITIZATION: the proof is the covering bound (exponential/periodicity), and the census (additive)
was the extremal of the EASY case. Integrates HYP-2909 (forward), HYP-2913/14 (census), S46 (Node 3),
kps's gamma-trick, S41 (observer-categories).*

## The re-prioritization (VERIFIED): the census is NOT on the proof's critical path
LRC(14) is "M(S) >= 1/14 for all 13-sets S". Split by the apex residue:
- **No runner ≡ 0 (mod 14):** t = 1/14 gives ||s/14|| >= 1/14 for every runner (all residues nonzero) =>
  M >= 1/14. **TRIVIAL** (verified on 2000 random non-mult-14 sets).
- **Some runner ≡ 0 (mod 14):** t=1/14 fails; M >= 1/14 must come from elsewhere. **The only hard case.**

So **LRC(14) <=> [every set with a multiple of 14 has M >= 1/14]** -- the COVERING BOUND. And the census
{AP, GW} is non-mult-14 (its M=1/14 is achieved at t=1/14): it is the extremal of the EASY half, NOT on
the critical path. Worse for the census, better for the proof: the proof needs M >= 1/14 (non-strict),
while the census needs M > 1/14 strict (to exclude covering sets from the tight locus). **The proof is
strictly weaker than the census we spent sessions on.**

## The four Farey variations = the hyperoperation hierarchy = the project's scales
On (a, b) = (numerator, denominator), the owner's four variations are add/mult/exp:
| variation | structure | project scale | role in LRC(14) |
|---|---|---|---|
| **a + b** (additive) | mediant / Farey / Stern-Brocot | three-gap (Steinhaus) | the CENSUS (tight locus) -- the EASY case's extremal |
| **a · b** (multiplicative) | product | H = I(Ω,2); a·b of 1/14 = **14 = 2·7** | the apex condition: a multiple of 14 sits on the observer |
| **b^a** (exp, denom^num) | tower | **apex-periodicity**: 14m is 1/14-periodic | the COVERING BOUND -- the gamma-trick = THE PROOF |
| **a^b** (exp, num^denom) | tower | Burnside 2^orbits / H=1+2^d / Cayley-Dickson | the orbit/2-adic scale (where 14|s collapses) |

The owner's first observation (1/14, K5, K3,3 all have a+b = 15) was the ADDITIVE face -- and additive is
exactly the census/three-gap face we now see is the easy-case extremal. The PROOF lives one level up, at
the **exponential**: a multiple of 14 is 1/14-PERIODIC (b^a), which is precisely kps's gamma-trick, and
that is the covering bound.

## The improved argument (the redirect)
1. **Reduce (trivial):** LRC(14) <=> covering sets (a multiple of 14) have M >= 1/14.
2. **Peel large multiples (Node 3, S46 -- the equidistribution / exponential face):** a large multiple of
   14 equidistributes within the seed's positive safe set (positive by LRC on fewer runners) -- removes
   exactly 1/7, a safe point survives => M >= 1/14. (Verified: seed{1..11,13}+{84,168,210} all M>=1/14.)
3. **Bounded multiples (gamma-trick, kps -- the exact apex-periodicity, b^a):** a multiple of 14 is
   1/14-periodic, decouples to the fine phase; the rest is pigeonholed. This is the SAME exponential face.
4. So the whole proof is the **exponential/periodicity** variation: the covering bound via apex-periodicity
   (gamma-trick) + equidistribution (Node 3). The additive (a+b / census / three-gap) is the easy-case
   extremal and is NOT needed.

## What remains (honest, now correctly scoped)
The covering bound is still open, but the open part is now PRECISELY located and far smaller than the
census program:
- the gamma-trick's residual (the prime-tower descent 14 -> 7 -> 2, the small-prime counting gap, S52);
- an effective equidistribution / Erdos-Turan bound for the unbounded peel (Node 3 rigor);
- the induction base (LRC on fewer runners, proven <= 6).
None of these is the census (consec-maximizes / three-gap rigidity), which we now see is irrelevant to
the proof. **The synthesis: stop proving the census; finish the covering bound (the exponential variation,
gamma-trick + Node 3).** That is the redirect the four variations reveal.

Related: HYP-2900/HYP-2906 (Node 3 / interval peeler), kps gamma-trick (apex-periodicity), HYP-2909
(forward, the apex condition a·b=14), the-obstruction-combining-duality (additive<->multiplicative),
three-observer-categories (cat 1 coverage / cat 3 H).
