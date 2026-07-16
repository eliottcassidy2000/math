# The parallel-class circle: quotient by rotation, extremal at contiguity

*kind-pasteur-2026-07-16-S128 (cont.33). Owner prompt: consider deeply the parallel-class
circle and what it represents abstractly; apply similar concepts to other niche problems
from the repo's history.*

## The abstraction

Death-star's THM-913 (+ LEM-030, THM-922) puts the chords of K_n (n odd) into n classes by
the invariant a+b mod n; the classes form a circle Z_n; interactions (crossings) descend to
a CIRCULANT form on that circle; and the extremal object (Guy's optimal 2-page drawing) is
the CONTIGUOUS split. Strip the specifics and a four-part pattern remains:

1. **A pair-structure with a Z_n rotation.** The primitives are pairs/arcs/chords/gaps;
   the symmetry is cyclic.
2. **A class invariant** — sum (a+b) or difference (a−b) mod n — labeling the orbits.
   The orbit space is a small circle: the parallel-class circle.
3. **Descent**: pairwise interaction costs are rotation-invariant, so they become a
   circulant matrix on the class circle. All spectral theory is Fourier at roots of unity.
4. **Contiguity extremality**: the optimizers of the natural extremal problems are
   CONTIGUOUS arcs of class-sets, and optimality statements become max-cut / interval
   statements on the circle.

The circle is the same object the covering theory calls torsion: the classes are the
n-torsion points, and "clean structure lives at torsion" (klein-S296's witnesses) is the
same slogan as "costs are circulant on the class circle."

## Instances already in the repo (some recognized only now)

**(a) The Farey/measure profile IS a parallel-class computation.** THM-826 writes
m({1..k}; λ) = Σ over consecutive Farey pairs of max(0, (1−λ(i+j))/(ij)) — the gaps organize
by the MEDIANT LEVEL s = i+j, which is precisely the sum-class invariant. The deep-well
corridor law (THM-853 II) is the statement that on [1/14, 1/13] exactly ONE parallel class
survives (s = 13), and the whole corridor is the single-class sum (2H₁₂/13)(1−13λ). THM-819's
unit inversion (prime k+1 ⟹ the class is complete: all coprime pairs) is class-completeness
at prime level. The LRC's interval core {1..k} is itself the contiguous arc of the speed
circle — contiguity extremality is why near-APs and interval cores keep appearing as the
covering extremizers.

**(b) Circulant tournaments = the class circle with an orientation choice.** Arcs class by
difference d; a circulant tournament is a choice of one orientation per class (the symbol
set S); everything spectral (THM-125 block-diagonalization, Paley Gauss sums) is Fourier on
the class circle. New this session: the arborescence count evaluates as a product over the
class circle — τ_r(Paley_p) = p^{(p−3)/2}((p+1)/4)^{(p−1)/2} (refereed p = 3, 7, 11), a
closed form living entirely at the torsion points.

**(c) The staircase's range classes are the CUT-OPEN parallel classes.** Tiles (x,y) class
by range r = x−y — the difference invariant — but the base path (Rédei's spanning tree)
cuts the circle open, so the class structure carries a boundary coordinate d = min(y−1, n−x)
(distance from the cut). Referee (rangeclass_parallel_kps_S128c33.py, n = 5, 6): the wiggly
self-loop rate is constant on (r, d) cells — but honesty first: each (r, d) cell is exactly
one grid-reflection mirror orbit (the two distances sum to n−1−r), so the constancy is
SYMMETRY-FORCED, not a new law. The genuine content is quantitative: the class coordinate
DOMINATES (r = 2 skip-1 tiles: SL ≈ 2.3% vs r ≥ 3: ≈ 7% at n = 6) while the cut coordinate
is a small correction (0.2pp) — the cut-open circle is NEAR-CIRCULANT, with the skip-1
class as the special "adjacent-transposition" class (the S20er anisotropy, now read as a
class phenomenon, not a tile phenomenon).

**(d) The mod-8 ladders.** E4 levels (step-8 progressions, THM-854 no-holes) and the
Milgram/Gauss-sum probes (signatures mod 8) live on the 8-torsion circle; the "level"
coordinate of the metagraph is a class coordinate for the Z_8 action on quadratic phases.

## What the circle "represents abstractly"

It is the MODULI OF PAIR-SLOPES: a pair {a,b} on a cyclic vertex set has a well-defined
slope/midpoint class, and the circle parameterizes those classes. Crossing (interleaving),
resonance (LRC), and flip-interaction (staircase) are all functions of slope classes only —
that is the descent. The repeated discovery is that OUR EXTREMAL OBJECTS ARE INTERVALS in
this moduli: Guy's split, the interval core {1..k}, the corridor's single class, the
contiguous covering families. "Contiguity beats scattering" is the one theorem-shape
appearing in four costumes.

## Where it points next (named, honest)

1. **The 2-page/hypercube ground state** (THM-917's note): page codes live on 2^{C(n,2)},
   crossings are a quadratic form, conjectured ground value Z(n). The quadratic form is
   circulant-structured by (b); a Fourier/torsion analysis of that energy is the natural
   attack, and it is exactly the shape of our Milgram/discriminant machinery.
2. **Class-level transfer law**: THM-923's arborescence transfer F(u,v) descends, on a
   circulant tournament, to a class function on the circle — the transfer flow between
   rotation classes is a circulant Markov chain; its spectrum is computable in closed form.
   (Small named computation.)
3. **LRC composition**: the remaining named LRC(14) items (opus's level-5 wall; THM-518
   Weyl large-spread glue) are both statements about which residue classes a large-spread
   family can occupy — parallel-class language; worth one session of translation, since
   the corridor law shows single-class corridors produce one-line formulas.
