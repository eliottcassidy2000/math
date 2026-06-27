# LRC as Tournament Analysis: the Winding Realization, and Where the Metric Lives

*mac-mini-2026-06-22-S57. The owner's technique: realize the LRC data as tournaments and ask which iso
classes are achievable. Applied with exact definitions -- and the honest finding is precise about WHAT
the tournament frame captures (the order) and what it does not (the metric M=1/14, which lives at the
optimum). Builds on HYP-2913 (three-gap/residue census), S48 (H=7 forbidden = the apex iso class).*

## The exact realization (the winding tournament)
The n "things" = the observer (speed 0) and the n-1 runners. The relation at time t: the cutoff at 1/2
turns the continuous positions into a tournament,
    u -> v   iff   frac((s_u - s_v) t) in (0, 1/2)   (u is "ahead" of v by less than half the circle).
For generic t this is a genuine tournament on n vertices; as t crosses a coincidence t = k/(s_u-s_v) it
flips an arc, winding through tournament iso classes. The achievable set A(S) = { iso[T(t)] : t generic }.

## What the frame captures -- and what it loses (verified)
- **H = 7 NEVER appears.** The winding tournament is a real tournament, and H=7 is forbidden (THM-029);
  the apex iso class (Omega = K_3) is simply not in A(S) for any S. (= S48's apex-7 from the order side.)
- **A(S) over all t does NOT distinguish tight from non-tight.** Computed (n=5, observer included):
  AP [0,1,2,3,4], GW [0,1,3,4,7], and the NON-tight [0,1,2,4,8] ALL achieve the SAME 4 iso classes
  (H = 1,9,11,15). At n=6, AP achieves 4 classes (H=1,17,41,45), GW achieves 6 (adds H=23) -- they
  differ, but not as a clean tight-invariant. **The reason:** ranging over t keeps only the cyclic
  ORDER; the metric M = max_t min ||s_i t|| = 1/14 is a DISTANCE, which the order forgets.

## Where the metric lives: the iso class AT THE OPTIMUM
The LRC's metric content is the single iso class T(t*) at the optimum t* = a/n. There the runners sit on
the n-grid (residues mod n), and the winding tournament is the **circulant** C_n({1,...,floor((n-1)/2)})
(with antipodal ties at the apex difference n/2 when n is even). For the AP this is the full rotational
tournament; for GW it is the circulant with one vertex's residue doubled and another skipped. So:

> **The LRC tight-locus census = the set of achievable optimum iso classes**, and that set is exactly
> the <=3-gap residue configs of HYP-2913 (AP-circulant + its single-swaps). Tournament analysis recovers
> the three-gap/residue characterization -- the census is a tournament-iso-class question, with the
> metric encoded as "the optimum class is the n-grid circulant."

## Honest assessment
Tournament analysis is the right strategic frame -- it brings the project's machinery (H, OCF,
H=7-forbidden, the apex) to bear on the LRC, and it makes the census a question about which iso classes
are achievable AT THE OPTIMUM. But it is a REFRAMING: the achievable-over-t classes lose the metric, and
the at-optimum class reproduces the residue/three-gap characterization (HYP-2913) rather than transcending
it. The open core is unchanged: proving that the achievable OPTIMUM classes are exactly {AP-circulant,
GW-single-swap} for n=14 is the same Steinhaus/consec-maximizes rigidity. The value added: the census is
now visibly a tournament-iso-class realizability problem ("which iso classes arise at the optimum for
some integer speed set"), the precise object the project knows how to study -- with H=7 the forbidden
apex on both the H side and the LRC side.

Related: HYP-2913 (three-gap census), HYP-2909 (binding pair = the apex antipode), S48 (H=7 = apex),
HYP-2605 (winding tournament), THM-029/200 (H=7 forbidden).
