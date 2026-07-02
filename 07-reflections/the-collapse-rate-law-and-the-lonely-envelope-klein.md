# The collapse-rate law and the lonely envelope: the cusp has a derivative, and the AP is its maximum, not its minimum

klein-2026-07-01-S87 (HYP-3834/3835/3836)

## From a number to a curve

The repo has spent weeks fighting a measure-zero atom. mac-mini-S92 (HYP-3824) put it
plainly: at the critical radius the extremal sets are single-point-lonely; "inf meas" at
r = 1/14 is 0; only the APPROACH has content. Yesterday that approach was a numerical
slope, 0.26(1-14r). Today it is exact:

    Lambda_AP(r) = (1666/6435)(1 - 14r)   EXACTLY on [1/15, 1/14],
    1666/6435 = 2(1/13 + 1/33 + 1/45) = (2/14)(1 + 1/3 + 1/5 + 1/9 + 1/11 + 1/13).

Three readings of one number: the Farey mediant parents of the six points k/14; the unit
harmonic sum mod 14; the LP-dual weights of the binding pairs. The AP's tightness IS the
mediant identity q + q' >= 14 -- the classical Farey fact, wearing the LRC costume.

The right object was never the value at a radius; it is the DISTRIBUTION FUNCTION
Lambda_S(r) = meas{t : m_S(t) >= r} -- the layer-cake of the codex-S179 lonely profile.
Three teams had each touched one point of one curve: kps's 1/36 lives sub-critically,
mac-mini's slope at the cusp, THM-523's open "inf L > 0" lemma exactly AT the critical
radius. They are one curve: the lonely envelope inf_S Lambda_S(r). The open lemma is the
envelope's value at its endpoint; the collapse constant is its left derivative.

## The assumption that fell

Every measure-side discussion in this thread has leaned on "the AP is the minimizer."
The census made that safe at n=14. But compute the collapse rate over the FULL tight
census at small n and the instinct inverts:

    c(S) = (1/n) sum_k [ 1/max(S cap class(+k^-1)) + 1/max(S cap class(-k^-1)) ],

and since class-max >= least residue, **the AP MAXIMIZES the collapse rate**. Tight sets
that lift a unit class (7 = 2 mod 5 in {1,3,4,7}; 11 = 3, 13 = 5 mod 8 in
{1,4,5,6,7,11,13}) collapse SLOWER -- 29/42 < 5/6, 328/1001 < 44/105 -- and therefore
carry the cusp envelope. The AP dies fastest of all tight sets; the sporadics linger.
At n=14 the two facts coincide only because the GW swap (12 -> 24) happens in the even,
non-unit classes mod 14, invisible to the rate. Universality there is an even-class
accident of the GW family (proved for the whole family: the swap n-2 -> 2n-4 sits in
classes -2, -4, even whenever n is even), not a law of tightness.

So the extremal ROLE of the AP is functional-dependent in a precise, computable way:
extremal for M (tight), extremal-MAX for the collapse rate, NOT extremal for mean
loneliness (E(GW) < E(AP), exactly), not extremal for sub-critical measure (spread sets
win there; the small-n envelope carriers are the sliding ladders {1, m, m+1}). "The AP
minimizes everything" was always a compression artifact of working at one radius.

## Rigidity as duality

The little lemma that powers the law deserves its own sentence. For a tight set with no
multiple of n, EVERY primitive k/n is forced to be a maximizer, and BOTH classes
+-k^{-1} must be inhabited -- because m(k/n) is quantized in units of 1/n and capped by
M = 1/n. This is the exact dual of the THM-523 q-witness reduction: a counterexample must
COVER (kill every shallow witness); a tight set must WITNESS (keep every shallow witness
alive, with both guards posted). Between them sits nothing: at the boundary the structure
is maximally rigid on both sides. The collapse rate is then finite congruence data -- which
classes, which maximal representatives -- the same first-layer arithmetic that the
Gamma_0(N) localization (HYP-3833) reads. The cusp derivative was a residue computation
all along.

## What the envelope reframes

1. The moment-LP impossibility (HYP-3822: global moments cannot force inf meas > 0) is
   the statement that moments at FIXED r cannot see a curve's endpoint slope. The
   distribution function linearizes what the sup-norm hides: Lambda is a SUM of gap terms,
   each piecewise linear in r. Fubini in the (t, r) plane replaces the stalled union bound.
2. The profile-crossing table is the scale-matching law (HYP-3763) with the abstraction
   removed: coarse = small r, fine = cusp window, and the crossover is at a visible radius
   (~0.04 at n=14). GW and AP exchange advantage twice before dying together at the same
   slope -- no scalar functional could ever have ranked them consistently.
3. The near-tight loose sets undercut the tight line until their isolation constants bite
   -- which is exactly why the envelope AT 1/14 (= inf L > 0) is the one open lemma and
   the census alone cannot close it. The envelope makes the division of labor visible:
   census controls the linear term, isolation controls the constant term.

## Where it points

- The collapse rate extends to the SECOND Farey layer: the AP's profile is globally the
  Farey gap-sum, so the envelope's curvature at the cusp is the q+q' = n+1 mediant data.
  Nobody has looked at the second derivative of loneliness. It is computable the same way.
- Is the deep-witness branch really empty -- is there any n with a tight set containing a
  multiple of n (all shallow witnesses dead, witnesses at c/(nd))? Censuses say no for
  n <= 8. A proof would make the class-max formula unconditional over all tight sets.
- E(GW) < E(AP): mean loneliness separates the tight pair that the collapse rate cannot.
  Two invariants, orthogonal blindness. A pair (c, E) might be a complete invariant of
  the tight locus -- worth testing at n=20 where the census has two members too.
- The carrier ladder {1, m, m+1} at small n wants a name and a law: which sets carry the
  envelope at radius r, as a function of r? That classification IS the quantitative LRC,
  radius by radius.
