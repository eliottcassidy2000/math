# Open Questions Resolved in S116

## OPEN-Q-031: RESOLVED (approximate, not exact)

arg(lambda_c)/pi = 0.692717, ln(2) = 0.693147. Difference 4.3e-4.
Cannot be exact: tau is algebraic, ln(2) is transcendental (Lindemann-Weierstrass).

## OPEN-Q-030: KEY LEMMA PROVED

**Lemma:** If arc a->b is in no transitive triple in tournament T on n vertices,
then score(a) = 1 and score(b) = n-2.

**Proof:** For every c != a,b, the triple {a,b,c} must be a 3-cycle (since it's
not transitive by hypothesis). The only 3-cycle with T[a][b]=1 is a->b->c->a.
This requires T[b][c]=1 and T[c][a]=1 for ALL c. Therefore:
- b beats every vertex except a: score(b) = n-2.
- a loses to every vertex except b: score(a) = 1. QED.

**Verified:** 100% at n=4 (24/24), n=5 (160/160), n=6 (1920/1920).

**Consequence for Q-030:** The non-core arc a->b connects the unique
score-1 vertex to the unique score-(n-2) vertex. Since these vertices
are unique (by the score constraint), there is at most ONE non-core arc
in any tournament. This, combined with the transitive core analysis,
should give sim_H in {0,1} for all n >= 4.
