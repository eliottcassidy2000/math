# LRC14 Phase-Color Reservoir

Codex 2026-06-18.

The threshold ladder gave us a cleaner uncolored target, but it still left a
quiet mismatch.  If the actual CRT residue is `a/(14V)`, then the color
`b=a mod 14` is already chosen.  We cannot simply say "there is a big gap
somewhere"; the gap has to contain the specific phase center `b/14`.

That sounds like a nuisance, but it turns into a better object.

For each color `b`, define

`C_b(E) = {x : ||e*x - b/14|| >= 1/14 for all e in E}`.

Then the finite statement is exact:

`a/(14V)` is a witness for `P union {V-e}` iff

`a=b+14t` and `a/(14V) in G_P cap C_b(E)`.

So the true reservoir is not a single set of `x`; it is a colored stack of
sets.  Its mass is

`Sigma(P,E)=sum_b meas(G_P cap C_b(E))`.

## Why It Feels Like Progress

The old uncolored minima were tiny: `rho_{4/15}=2/525`, `rho_{1/4}=1/140`.
Those were useful but thin.  The colored object is large:

- The structured-bank minimum was `14249/28028 ~= 0.508`.
- The worst largest-single-color mass was still about `0.0466`.
- Via-max-zero shapes had `Sigma` near `0.885` or even `1.363`.

This explains why actual CRT witnesses were abundant even when the via-max
certificate died.  The finite count is not trying to hit a skinny uncolored
component.  It is sampling thirteen nonzero color layers.

## The Proof Shape

The residual can now be phrased as:

`# witnesses at q=14V = colored_grid_count(P,E,V)`.

Heuristically,

`colored_grid_count(P,E,V) ~= V * Sigma(P,E)`.

The named hard lifts were close to this prediction; the worst sampled ratio
was about `0.68`.  So if we can prove a discrepancy bound smaller than the
`V*Sigma` main term, the CRT placement lemma closes for large `V`; the small
`V` base is finite.

This is much more satisfying than the private-pair denominator story.  The
denominator `D=14q-r` becomes a consequence: if every colored grid hit is
blocked, the finite colored discrepancy has overwhelmed a positive continuous
mass.  That is the actual obstruction to rule out.

## Tournament Note

I used phase colors as tournament vertices.  The result was transitive:
colors `1` and `13` dominate, color `0` loses.  That is not a rich tournament,
but it is the right kind of boring: it reflects the geometry directly.  Color
`0` parks the maximum runner at the observer; colors `1` and `13` put it on
the threshold edge.  The tournament is a sanity check that the quotient is
seeing the same object as the CRT condition.

## Remaining Bite

There are still two genuine tasks:

1. Prove a uniform `Sigma` floor over all S3 relation-lattice shapes.
2. Prove the colored discrepancy bound for the progressions
   `(b+14t)/(14V)`.

The first likely belongs to the B(k) Fourier/relation-lattice work.  The
second is a finite union-of-intervals discrepancy estimate, with endpoints
coming from small speeds `P` and offsets `E`.  This is a real theorem still to
prove, but the shape of it is finally crisp.
