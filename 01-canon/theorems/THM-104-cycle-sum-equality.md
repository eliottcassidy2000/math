# THM-104: Cycle Sum Equality for Edge-Sharing 3-Cycles

**Status:** PROVED (algebraic)
**Filed by:** kind-pasteur-2026-03-08-S43

## Statement

If two directed 3-cycles C_1 and C_2 in a tournament share a directed edge,
then for any TT-cocycle z (i.e., z in Z_1 with z(a,b)+z(b,c)=z(a,c) for all
TT triples (a,b,c)):

  sum_{C_1}(z) = sum_{C_2}(z)

where sum_C(z) = z(a,b) + z(b,c) + z(c,a) for cycle C = a->b->c->a.

## Proof

Let C_1 = (a->b->c->a) and C_2 = (a->b->d->a) share directed edge a->b.
Then c->a and d->a (from the cycle directions).

For the edge between c and d (exactly one direction in a tournament):

**Case 1: c->d.**
- TT(b,c,d): b->c->d, b->d => z(b,c) + z(c,d) = z(b,d)
- TT(c,d,a): c->d->a, c->a => z(c,d) + z(d,a) = z(c,a)
- Therefore: z(b,c) - z(b,d) = -z(c,d) = z(d,a) - z(c,a)
- So sum_C1 - sum_C2 = [z(b,c)+z(c,a)] - [z(b,d)+z(d,a)]
                      = [z(b,c)-z(b,d)] - [z(d,a)-z(c,a)] = 0

**Case 2: d->c.**
- TT(b,d,c): b->d->c, b->c => z(b,d) + z(d,c) = z(b,c)
- TT(d,c,a): d->c->a, d->a => z(d,c) + z(c,a) = z(d,a)
- Therefore: z(b,c) - z(b,d) = z(d,c) = z(d,a) - z(c,a)
- So sum_C1 - sum_C2 = 0

## Corollary

All 3-cycles in a connected component of the 3-cycle adjacency graph
have the same TT-cocycle sum.

## See Also
- THM-105 (dominant vertex forcing)
- THM-107 (at most 1 free component)
