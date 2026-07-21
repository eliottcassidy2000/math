# GMC(2) as a GIT nullcone, and the three cancellation theorems

boxeph-2026-07-20-S188 (HYP-8600, THM-1800).

The moment functional of the circular pair is the Fischer/Bargmann pairing;
one-sidedness is Hilbert-Mumford instability for the hyperbolic torus; so
GMC(2) says: analytic instability (all moments vanish) = algebraic
instability (a 1-PS kills you). Kempf-Ness shape.

The tournament shadow: three theorems with one shape —
1. REDEI: h is odd. Transitive tournaments survive mod 2; everything else
   pairs off under path reversal.
2. VANDERMONDE: prod(x_i - x_j) = signed tournament sum. Transitive
   tournaments survive exactly; everything else cancels under cycle
   reversal (verified n <= 5; the naive canonical pairing needs care at
   n >= 4 — the instability of lex-first choices under the very reversal
   they index is itself a recurring theme: cf MISTAKE-199's race lessons).
3. GMC(2) (conjecture): one-sided P survive the moment pairing trivially
   (all moments vanish); the conjecture is that everything two-sided is
   SEEN (some moment nonzero) — intransitivity cannot hide.
In all three, an invariant pairing kills exactly the intransitive/
two-sided/non-triangular part. Tournaments are in/transitivity itself:
they are the universal witnesses of what such pairings cancel.

Where this could bite mathematically: Kempf-Ness gives a VARIATIONAL
characterization of GIT-(semi)stability via moment-map norms; if the
analytic moment-nullcone question can be put in Kempf-Ness form (the
E-pairing as a moment-map square), the two-sided nonvanishing would follow
from a convexity/properness argument — a genuinely different route than
the resurgence ladder. Next session's shape: write E[P^m] as the m-th
Fischer power against the coherent state; ask whether log E[e^{tP}] is
convex along the gauge orbit (the Kempf-Ness function analog). The S180
trap (gauge-invariance of critical values) says the NAIVE version fails;
the Kempf-Ness function is built from the NORM, not the critical values —
the distinction that time was missed.
