# Berggren edge transport: how the three tree children move the Collatz guard (a,b,k)

**Status: PROVED (the root-coordinate form of the Berggren maps, the complete
transport table (a,b,k) -> (a',b',k') for all six child/orientation cases, the
multiplier-keeping classification, the complete list of sporadic legal children
for multiplier 3 and parameter +-1, the parent table on the three THM-3756 cones,
the identification of the inverse-fibre braid as the Berggren word B1^((4^j-4)/3) B2)
+ FINITE-EXACT (census of all coprime odd pairs s<=2000 against all guards a<=15,
|b|<=15; hypotenuse-<=10^6 triangle counts; 1849203 consecutive orbit triangles)
+ CITED (THM-3756, odd_square edges/synthesis notes, arithmetic braids)
+ REFUTED (the brief's "depth <= 1 above every k=1 edge": one witness, 3->5)
+ SCOPE (negative multipliers).  Nothing here bears on Collatz convergence or
cycle completeness; no literature-priority claim.  Adversarially audited
2026-09-22 ([audit script](../../04-computation/experiments/collatz_mod6_20260921_berggren_edge_transport_audit.py),
[audit output](collatz_mod6_20260921_berggren_edge_transport_audit.out)): every
key number was recomputed by independent code (Berggren's triple matrices
instead of Euclid algebra, direct edge enumeration instead of the numpy pair
grid, Euclid (m,n) enumeration for the hypotenuse count, an exact integer
linear solve for section 7) and agrees; four statements of the first version
were wrong or overstated as written and are corrected in place (section 10):
"the parent of a k=1 edge is never a (3,+-1) edge" (three accidental
exceptions, the sporadic list read upward), the "b<=-3 / b>=3" rows of
Corollary 2.2 (necessary only; smallest witnesses b=-7 and b=9), "the fibre is a
subset of the B1-ray" (true for the fibre minus its k=1 element), and the
section 7 "PROVED sketch" (now a proof: ratio bound plus exact solve).**  Lane `berggren_edge_transport`
of session `collatz-mod6-20260917` (machine `mac-mini`), wave of 2026-09-21/22.

## Inheritance and concept board

Inherited without re-proof: the odd-square encoding of a directed odd edge
x -> y as the primitive triple with roots (s,t)=(max(x,y),min(x,y)) and the
fact that the Berggren parent of the U_1 edge 7->11 is the pair (7,3), which
is not a U_1 edge ([odd_square edges, section 1 and 3](odd_square_20260921_edges.md),
[odd_square synthesis, section 4](odd_square_20260921_synthesis.md), with the
edge-count law C_edge sqrt(X), C_edge=0.5078...); the affine Berggren children
L,M,R on the odd roots (q,d) with determinants 1,-1,1, and the unique inverse
descent on three disjoint cones
([THM-3756, (25) and (28)](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md));
the inverse-fibre braid R_b(x)=4x+b and the inverse guard (I1)
([odd_square inverse, section 2](odd_square_20260921_inverse.md),
[arithmetic braids collatz](arithmetic_braids_20260917_collatz.md) section 2);
the Gaussian-squaring/Pell reading of triples (THM-3341) is only cited, not used.
The session synthesis is [collatz_mod6_20260917_synthesis.md](collatz_mod6_20260917_synthesis.md);
the sibling semicircle lane [pythagorean_semicircle](collatz_mod6_20260917_pythagorean_semicircle.md)
uses the same triangle chart in the angle coordinate and is not repeated.
"The brief" below is the prompt that launched this lane; it is not a
repository file, so the statements attributed to it (the b=1, k=1 table, the
"depth <= 1" claim, the Euclid example (9,2)->(20,9)->(64,53), the orbit count)
are UNCITED-RECOLLECTION as to wording and are verified here on their
mathematical content only (audit A11 recomputes the orbit count).

Concept board. Closest proved mechanism: THM-3756's affine children and cone
descent, here written on the roots as B1:(s,t)->(s+2t,t), B2:(s,t)->(2s+t,s),
B3:(s,t)->(2s-t,s). Canonical hostile: a k>=2 Collatz edge, whose three
children all leave the multiplier 3 by identity (x-type children shift a by
+-2^(k+1), the y-type child acquires the odd cofactor 2^(k-1)+3); the only
accidental exceptions are the sporadic parents 7->5 and 3->1 of section 3.
Corrected near miss: the opus statement "neither 3->7 nor 7->3 is a U_1 edge"
is correct as stated; the near miss is to read it as "the parent is illegal,
full stop", whereas the parent pair (7,3) is the halving edge 7->3 of the
multiplier-1 family E_{1,-1,1} and, in the opposite orientation, the 3x+5 edge
3->7, with the affine parameter law b''=(x+3b)/2. Least-used sidecar: the
B1-ray {(x,y): t=y fixed}, on which the inverse fibre of a target y, apart
from its k=1 element, sits at the B1-heights 2^(k0-1)(4^j-1)/3.

| Object | Kept by the bridge | Lost / decisive test |
|---|---|---|
| Directed edge x->y of E_{a,b,k} | coprimality, both parities, primitivity | orientation and the guard (a,b,k) |
| Berggren child of the root pair | one endpoint (source for x-type, target for y-type) | the other endpoint is replaced |
| Transport table | multiplier shift 2^(k+1), parameter sign, exponent | which children are E-edges at all |
| Inverse fibre of y | a B1-ray with t=y (all elements with k>=2) | the k=1 element is off the ray and needs one B2 |
| Sporadic legal children | five readings, all with x<=7 | none beyond the Berggren root cluster |

Notation. E_{a,b,k} = {(x,y) odd positive, coprime, x!=y : ax+b = 2^k y}, a odd
positive, b odd, k>=1; the exponent k is the exact 2-adic valuation of ax+b
because y is odd. "(3,b,k)" abbreviates E_{3,b,k}. Root pair (s,t) with s>t.
Script: `04-computation/experiments/collatz_mod6_20260921_berggren_edge_transport.py`;
all numbers below are in its output `collatz_mod6_20260921_berggren_edge_transport.out`
(claims S1-S11 there).

## 1. The Berggren maps in root coordinates (S1)

**PROVED.** With Euclid parameters (m,n), m>n>=1, coprime, opposite parity, and
roots s=m+n, t=m-n (so m=(s+t)/2, n=(s-t)/2), the three children
(2m-n,m), (2m+n,m), (m+2n,n) are

```text
B1 (s,t) = (3m-n, m-n) = (s+2t, t),
B2 (s,t) = (3m+n, m+n) = (2s+t, s),
B3 (s,t) = (m+3n, m+n) = (2s-t, s).                                (1)
```

*Proof.* Substitute m=(s+t)/2, n=(s-t)/2: 3m-n = (3s+3t-s+t)/2 = s+2t and
m-n = t; 3m+n = 2s+t and m+n = s; m+3n = 2s-t and m+n = s. QED. The matrices
[[1,2],[0,1]], [[2,1],[1,0]], [[2,-1],[1,0]] have determinants 1, -1, 1 and
preserve oddness, coprimality and s>t>=1 (checked on all 2211 children of the
primitive Euclid pairs with m<=60). This is THM-3756 (25) with (q,d)=(s,t) and
L=B1, M=B2, R=B3 (CITED); the point of restating it is that B1 fixes t, while
B2 and B3 move the old s into the new t.

## 2. The complete transport table (S2, S3)

Let x -> y be an edge of E_{a,b,k}, so y=(ax+b)/2^k. Each child pair contains
exactly one of x, y. Call the child *x-type* if it contains the source x, and
*y-type* if it contains the target y.

**Theorem 2.1 (PROVED; symbolic identity in x, a, b, K=2^k, verified on 20233
edges with a<=15, |b|<=15, x<400).**

```text
case A, y>x, (s,t)=(y,x):
  B1: (y+2x, x)   x-type   edge x -> y+2x   in  E_{a+2^(k+1), b, k}
  B2: (2y+x, y)   y-type   edge 2y+x -> y   with  a u + b = (2^k + 2a) y
  B3: (2y-x, y)   y-type   edge 2y-x -> y   with  a u - b = (2a - 2^k) y
case B, y<x, (s,t)=(x,y):
  B1: (x+2y, y)   y-type   edge x+2y -> y   with  a u + b = (2^k + 2a) y
  B2: (2x+y, x)   x-type   edge x -> 2x+y   in  E_{2^(k+1)+a, b, k}
  B3: (2x-y, x)   x-type   edge x -> 2x-y   in  E_{2^(k+1)-a, -b, k}          (2)
```

An x-type child is always an E-edge when its multiplier is positive (source x
kept, exponent k kept, multiplier shifted by +-2^(k+1), parameter b or -b; the
sign case 2^(k+1)-a<0 is the SCOPE item below). A y-type child keeps the target y and
the multiplier a; it is an E-edge, namely of E_{a, +-b, log2 M}, exactly when
M = +-2^k + 2a is a positive power of two; otherwise it is only a
*generalized edge* a u +- b = 2^j m y with an odd cofactor m>1.

*Proof.* x-type: u = 2x + y (case B, B2) gives 2^k u = 2^(k+1) x + ax + b, i.e.
(2^(k+1)+a)x + b = 2^k u; u = 2x - y gives (2^(k+1)-a)x - b = 2^k u;
u = y + 2x (case A, B1) gives (a + 2^(k+1))x + b = 2^k u. The multiplier is odd
and the exponent is exact because u is odd. y-type: u = x + 2y gives
a u = ax + 2ay = (2^k y - b) + 2ay, i.e. a u + b = (2^k+2a) y; u = 2y - x gives
a u - b = (2a - 2^k) y. A y-sourced identity a'y + b' = 2^j u is
impossible: multiplying a u = (alpha 2^k + 2a) y - alpha b by 2^j and comparing
coefficients of y would give a a' = 2^j (alpha 2^k + 2a), even on the right and
odd on the left. The unordered child pairs agree between the cases
({x, y+2x}, {x+2y, y}, {2y-x, y} or {2x-y, x}), as they must. QED

The brief's table for b=1, k=1, x=3,7,...,31 is printed in S2: B1 gives x->y+2x
in (7,1,1), B2 gives (4x+1)->y in (3,1,3), B3 gives (2x+1)->y in (3,-1,2); e.g.
7->11 has children (25,7), (29,11), (15,11).

**Corollary 2.2 (PROVED, S3; which children keep the multiplier).** An x-type
child never keeps a (a+2^(k+1)>a, and 2^(k+1)-a=a would need a=2^k). A y-type
child keeps a and is an E-edge iff

```text
k=1:   M = 2(a+1) [alpha=+1]  or  M = 2(a-1) [alpha=-1]  is a power of two,
       i.e. a = 2^j - 1, resp. a = 2^j + 1;
k>=2:  alpha=+1 never (odd part 2^(k-1)+a >= 3);
       alpha=-1 iff a = 2^(k-1)+1, and then M=2, k'=1.                  (3)
```

For a=3 the exhaustive list for k<=12 is (k,alpha,M) in {(1,+1,8), (1,-1,4),
(2,-1,2)}. Hence the multiplier-3 children of a multiplier-3 edge are exactly:
from k=1 with y>x, B2 -> (3,b,3) [this is the inverse braid R(x)=4x+b] and
B3 -> (3,-b,2); from k=1 with y<x (needs x<-b, so b<=-3; b=-3,-5 give no such
edge and the smallest witness is b=-7, the edge 3->1 of (3,-7,1) with B1 child
(5,1) = 5->1 in (3,-7,3); audit A4), B1 -> (3,b,3); and from k=2 with y>x
(needs x<b, so b>=3; b=3,5,7 give none and the smallest witness is b=9, the
edge 1->3 of (3,9,2) with B3 child (5,3) = 5->3 in (3,-9,1); audit A4),
B3 -> (3,-b,1). For
parameter b=+-1, S4 shows y>x iff k=1 (99998 edges with odd source <=10^5,
no violation; proof: k=1 gives y-x=(x+b)/2>0 for x>=3; k>=2 gives
2^k y = 3x+b < 4x <= 2^k x for x>=3; x=1 only produces the excluded fixed edge),
so **no child of a k>=2 edge of (3,+-1) has multiplier 3 by identity** (PROVED).

**SCOPE.** 2^(k+1)-a<0 (x-type, case B) needs b<0 with |b| > (a-2^k)x, e.g.
the edge 3->1 of E_{9,-23,2} (27-23=4=4*1) has B3 child (5,3) with -3+23=20=4*5, multiplier -1; and M=-2
(a=2^(k-1)-1, alpha=-1) would be multiplier -a. Negative multipliers are outside
the families and are not counted anywhere below.

## 3. The multiplier-3 sub-language: generic children and five sporadic readings (S5, S6)

Identity readings are not the only way a child pair can be a (3,+-1) edge: the
pair may accidentally satisfy 3p+b'=2^j q for another (b',j) or in the other
orientation. This is the question actually asked ("does the tree restricted to
legal Collatz edges have depth <= 1 above any k=1 edge?").

**Theorem 3.1 (PROVED + FINITE-EXACT).** Let x->y be an edge of (3,b,k), b=+-1,
and let (p,q), p>q, be one of its three Berggren children. Then (p,q) is a
(3,b',j) edge, b'=+-1, in some orientation exactly in the following cases:

* generic (every k=1 edge, i.e. every x>=3 with 3x+b = 2 mod 4): the B2 child
  (4x+b, y) read as (4x+b)->y in (3,b,3), and the B3 child (2x+b, y) read as
  (2x+b)->y in (3,-b,2);
* sporadic, exactly five readings:

```text
E1  3->5   (3,+1,1):  B3 child (7,5)   also reads  5->7   in (3,-1,1)
E2  7->5   (3,-1,2):  B2 child (19,7)  reads       19->7  in (3,-1,3)
E3  7->5   (3,-1,2):  B3 child (9,7)   reads       9->7   in (3,+1,2)
E4  3->1   (3,-1,3):  B1 child (5,1)   reads       5->1   in (3,+1,4)
E5  3->1   (3,-1,3):  B3 child (5,3)   reads       3->5   in (3,+1,1)         (4)
```

*Proof.* (i) Orientation and exponent. If p->q is a (3,b',j) edge then
2^j q = 3p+b' < 2^j p forces j>=2; if q->p then 3q+b' = 2^j p > 2^j q forces
j=1. (ii) Ratios. The children have p/q = rho+2 (B1), 2+1/rho (B2), 2-1/rho (B3)
with rho=s/t>1. For B2, 3p/q in (6,9) and 2^j is within 1/q of it, so j=3 is
the only p->q exponent, and q->p would give s = b'-2t < 0. For B3, 3p/q in
(3,6) with q>=3 (q=1 would force p>=3>2q), so j=2 only; q->p with j=1 gives
s = 2t+b'. For B1, q->p gives 2s = b'-t <= 0, so only p->q occurs.
(iii) Each surviving reading is a linear equation in x after y=(3x+b)/2^k.
B2, p->q, j=3: 2s = 3t+b'. Case A (k=1, s=y): b'=b, the generic identity.
Case B (k>=2, s=x): (2^(k+1)-9)x = 3b+2^k b'; k=2 gives x=-3b-4b', positive
only for b=b'=-1, x=7 (E2); k=3 gives 7x in {+-5,+-11}; k>=4 has
2^(k+1)-9 > 2^k+3 >= |RHS|. B3, p->q, j=2: 2s = 3t-b'. Case A: b'=-b, generic.
Case B: (2^(k+1)-9)x = 3b-2^k b'; k=2: x=4b'-3b, i.e. x=7 for (b,b')=(-1,1)
(E3) or x=1 for b=b'=1 (excluded fixed edge); k=3 and k>=4 as before.
B3, q->p, j=1: s = 2t+b'. Case A: x = b-2b', i.e. x=3 for (b,b')=(1,-1) (E1) or
x=1 (excluded). Case B: (2^k-6)x = 2b+2^k b'; k=2 gives x=-b-2b' with y even or
x=y; k=3 gives x=b+4b', i.e. x=3, y=1 for (b,b')=(-1,1) (E5) or y even; k>=4:
x=1 is impossible and x>=3 has 3(2^k-6) > 2^k+2. B1, p->q: (2^j-6)t = 3s+b'.
Case A (k=1): (2^(j+1)-21)x = 3b+2b', whose only positive solution is x=1
(j=3, b=b'=-1), excluded. Case B: 3x(2^j-6-2^k) = 2^k b' - (2^j-6)b. For
k=2 the solutions are x=1 (excluded) only; k=3, j=4 gives 6x = 8b'-10b, i.e.
x=3, y=1 for (b,b')=(-1,1) (E4); k=4 has no solution; for k>=5 the cases
j<=k-1, j=k, j=k+1, j>=k+2 are each excluded by comparing 3|2^j-6-2^k| with
2^k+2^j-6 (j=k+1 leaves only x=1, which fails). So k<=4 and j<=k+2 in every
case. The script solves every reading exactly for k,j<=40 (well beyond the
bound) and finds precisely the generic identities and (4); an independent brute
force over all 4165 edges of (3,+-1) with max(x,y)<=5000 finds 3332 = 2*1666
generic legal children (two per k=1 edge, 1666 k=1 edges) and the same 5
sporadic readings. QED

Note that E2 and E3 are the *generic* children of the k=1 reading 5->7 of the
pair (7,5); the sporadic content of (7,5) is E1 alone.

**Theorem 3.2 (PROVED; the legal sub-forest).** Restrict the Berggren tree to
pairs that are (3,+-1) edges in some orientation. Every k=1 edge has exactly two
legal children, both k>=2 edges and hence leaves, except that the child (7,5) of
3->5 is itself the k=1 edge 5->7 with two leaf children (19,7), (9,7). A k>=2
edge has no legal child except 7->5 (via 5->7) and 3->1 (children (5,1) and
(5,3)). Consequently:

* the brief's statement "depth <= 1 above any k=1 edge" is **REFUTED** by the
  single witness 3->5 -B3-> (7,5)=5->7 -B2-> (19,7)=19->7 (depth 2);
* corrected: 3->5 is the only k=1 edge of depth 2, every other k=1 edge has
  depth exactly 1, and the Berggren root cluster
  {(3,1),(5,1),(5,3),(13,5),(7,5),(19,7),(9,7)} (a tree of depth 3 rooted at
  the Berggren root (3,1)=3->1) is the unique legal component with more than
  three pairs. The sub-language of Berggren words that keeps (3,+-1) is
  {B2, B3} applied once to a k=1 edge, plus the finitely many words inside the
  root cluster.

Census (S6) of all 1664 legal pairs with max(x,y)<=2000, testing each child's
legality directly: 666 pairs with min k=1, all with exactly 2 legal children;
498 (k=2), 249 (k=3), 125 (k=4), 63, 31, 15, 8, 4, 2, 1, 1 (k=5..12) pairs
with 0 legal children; and the single pair (3,1) (k=3) with 2 legal children.
Audit A6 repeats this at max(x,y)<=4000 (3330 legal pairs, 1332 with k=1 and
two legal children each, (3,1) again the only k>=2 pair with children) and
computes the connected components of the legal sub-forest under parent links
inside the window: sizes {1: 832, 2: 500, 3: 497, 7: 1}, the single size-7
component being the root cluster (sizes 1 and 2 are window truncation: k>=4
pairs are isolated, and k=1 pairs near the boundary have one or both children
above 4000).

## 4. Parent transport on the three cones (S7)

**PROVED.** THM-3756's inverse maps on the roots are B1^-1(s,t)=(s-2t,t) for
s/t>3, B2^-1(s,t)=(t,s-2t) for 2<s/t<3, B3^-1(s,t)=(t,2t-s) for 1<s/t<2
(s/t=2 is impossible for odd s; s/t=3 only at the root (3,1)). For a (3,b) edge
with b=+-1 the cone depends on k alone (verified on all edges with
max(x,y)<=2000: 666 in C3 with k=1, 499 in C3 with k=2, 249 in C2 with k=3,
125, 63, 31, 15, 8, 4, 2, 1, 1 in C1 with k=4..12, plus the root):

```text
k=1  (y/x in (1,2)):  parent (x, (x-b)/2)   = halving edge x -> (x-b)/2 of E_{1,-b,1};
                      forward: (1,-b,1) -B3 (case B)-> (2^2-1, b, 1) = (3,b,1).
                      Other orientation: (x-b)/2 -> x is a (3,b'',1) edge with
                      b'' = 2x - 3(x-b)/2 = (x+3b)/2.
k=2  (x/y in (1,2)):  parent (y, (x+b)/2)   = k=1 edge (x+b)/2 -> y of E_{3,-b,1}   (B3^-1)
k=3  (x/y in (2,3)):  parent (y, (x-b)/4)   = k=1 edge (x-b)/4 -> y of E_{3,b,1}    (B2^-1, the braid)
k>=4 (x/y > 3):       parent (x-2y, y) with 3(x-2y)+b = (2^k-6) y,
                      odd cofactor 2^(k-1)-3 >= 5: never a (3,.,.) edge by identity.   (5)
```

*Proof.* k=1: 2t-s = 2x-(3x+b)/2 = (x-b)/2 and x-b = 2*(x-b)/2 is the
E_{1,-b,1} identity; the B3 row of (2) with a=1, k=1 returns (2^2-1,b,1).
k=2: 2y-x = (x+b)/2 and 3(x+b)/2 - b = (3x+b)/2 = 2y. k=3: x-2y = (x-b)/4 and
3(x-b)/4 + b = (3x+b)/4 = 2y. k>=4: 3(x-2y)+b = 2^k y - 6y. The cone
inequalities follow from y/x = (3+b/x)/2 for k=1 and x/y = 2^k x/(3x+b) for
k>=2 with x>=3. QED

This settles the opus statement in both directions: the parent of 7->11 is the
pair (7,3), which is the multiplier-1 edge 7->3 (7-1=6=2*3) and, in the other
orientation, the 3x+5 edge 3->7 (9+5=14=2*7, b''=(7+3)/2=5). The parent of a
Collatz k=1 edge is always the halving edge x -> (x-b)/2 and is never a (3,+-1)
edge *by identity*; accidentally it is one in exactly two cases, namely 3->5
(parent (3,1) = 3->1 in (3,-1,3)) and 5->7 (parent (5,3) = 3->5 in (3,1,1)),
which are the sporadic readings E5 and E1 of section 3 read upward (audit
fix: the first version said "never", which is false as written; audit A7 lists
the exceptions, and the sporadic list of Theorem 3.1 proves there are no
others). The parents of k=2 and k=3 edges are the k=1 edges of their own fibre
or of the mirrored parameter; the parents of k>=4 edges are generalized edges
with an odd cofactor by identity, with the single accidental exception 5->1,
whose parent is the root (3,1) (E4 read upward).
For general odd a the parent pair is still unique; it is an identity pre-image
under (2) only when it reads as (a-2^(k+1),b,k) [B1 or B2 of the appropriate
case], (2^(k+1)-a,-b,k) [B3, case B] or (a,+-b,k') with +-2^k'+2a=2^k; the
symbolic x-type inverse identities are verified in S7 for k<=5.

## 5. Census against all guards a<=15, |b|<=15, and the hypotenuse count (S8, S9)

**FINITE-EXACT.** Among the 405432 coprime odd pairs (s,t), t<s<=2000, exactly
54755 (13.51%) are E_{a,b,k} edges for some a<=15, |b|<=15 (odd), k>=1, in some
orientation, with 57837 directed realizations in total; 52860 pairs have a
single realization, 1415 have two, and the histogram tail reaches one pair with
40 realizations. The a=3 row has 833 realizations at b=-1 and 832 at b=1
(1664 pairs = 833+832-1: the pair (7,5) carries two b=-1 readings, 7->5 and
5->7; no pair carries both a b=+1 and a b=-1 reading, audit A8); per multiplier the totals are
12685, 10582, 8253, 6804, 5817, 5049, 4516, 4131 for a=1,3,...,15. For every
realization all 173511 children were checked: the identities (2) hold on all
of them; for 37529 children the predicted child guard lies inside the census
window and is found in the child's realization list; 73450 predictions fall
outside the window (a'>15, |b'|>15 or s'>2000); 62532 y-type children have M
not a power of two (generalized edges).

Hypotenuse count. There are 159139 primitive triangles with (s^2+t^2)/2 <= 10^6.
Of these, 507 are 3x+1 edges (507 directed, all distinct), against the opus law
C_edge*sqrt(X) = 0.5078*1000 = 507.8; 506 directed 3x-1 edges on 505 triangles
(the 5<->7 cycle is one triangle); the union over b=+-1 is 1012 triangles from
1013 directed edges (the difference is again the single triangle (7,5) with
its two b=-1 readings; the b=+1 and b=-1 triangle sets are disjoint); the union over
|b|<=5 has 406, 337, 506, 507, 338, 404 directed edges for b=-5,-3,-1,1,3,5
and 2483 distinct triangles, 1.560% of all primitive triangles at this height.
The near-equality 507 vs 507.8 is the opus O(log X) error term at X=10^6, not a
new law.

## 6. The inverse-fibre braid is a Berggren word (S10)

**PROVED.** Fix an odd target y not divisible by 3 and parameter b. The inverse
fibre is x_k = (2^k y - b)/3 over the k of one parity class (2^k y = b mod 3),
and x_{k+2} = 4x_k + b = x_k + 2^k y. Since B1^m(s,t) = (s+2mt,t):

```text
k>=2 (x_k > y):   (x_{k+2}, y) = B1^(2^(k-1)) (x_k, y);
k=1  (y > x_1):   (x_3, y)     = B2 (y, x_1).                                  (6)
```

Hence the word from the k=1 triangle to the k=2j+1 triangle of the same fibre is
B1^((4^j-4)/3) B2 (B2 first; exponent 0+4+16+...+4^(j-1)), and from k0>=2 to
k0+2j it is B1^((2^(k0+2j)-2^k0)/6). The fibre minus its k=1 element (whose pair (y,x_1) has t=x_1, off the ray)
is a subset of the B1-ray {(x,y): t=y} at B1-heights m_j = 2^(k0-1)(4^j-1)/3;
R is a tree step (B2) only once, at k=1, and afterwards a power of B1 whose
exponent quadruples (4, 16, 64, ...). Verified on 8002 fibre steps (y<=2001, k<=12, b=+-1). The
brief's example in Euclid coordinates, (9,2)->(20,9)->(64,53), is the root chain
(11,7)->(29,11)->(117,11), i.e. 7->11 (k=1), 29->11 (k=3), 117->11 (k=5);
(29,11)->(117,11) is B1^4 and not a child (the children of (29,11) are (51,11),
(69,29), (47,29)), which is exactly what the brief observed.

This is the sidecar that THM-3756 does not have: its outer-rank fibre fixes
s (the shell B+C=s^2), whereas the Collatz inverse fibre fixes t=y and is a
B1-ray with a lacunary set of heights.

## 7. Consecutive orbit triangles (S11)

**FINITE-EXACT.** Over all odd starts <=10^5 under U_1 there are 1849203
consecutive orbit-triangle pairs (recomputed independently in audit A11) and 0
of them are Berggren parent/child in either direction. **PROVED (audit A11;
the first version only sketched this).** Write r for the ratio max/min of an
edge pair. For b=+1 and x>=3, k=1 gives r=y/x=(3+1/x)/2 in (3/2,5/3] and
k>=2 gives r=x/y=2^k x/(3x+1) in [0.3*2^k, 2^k/3) (checked on x<=10^5 in
A11). A B2 child has ratio 2+1/rho in (2,3) and a B3 child 2-1/rho in (1,2),
so the child exponent is at most 3 (j>=4 forces r>=4.8), and the parent ratio
rho=1/(r2-2) resp. 1/(2-r2) lies in (3/2,5/2], (2,3] or [5/4,3/2), so the
parent exponent is at most 3 as well. A B1 child (s+2t,t) keeps t as its
smaller coordinate, so both pairs would share their minimum: with k>=2 that
forces z=x+2y, ratio >2, impossible for j=1 (ratio <=5/3) and impossible for
j>=2 (then y is the larger coordinate of (y,z)); with k=1 it forces z=x<y,
contradicting the retained smaller coordinate. Hence every parent/child
incidence between consecutive orbit triangles has k,j<=3, and the exact
integer linear solve of all 3x2x2 cases for k,j<=60 has the single solution
x=y=z=1, k=j=2 (the pair (1,1) is its own B3 child), which no orbit forms. QED

## 8. Typing the bridge

Source: directed edges of the Collatz families E_{3,+-1,k} (more generally
E_{a,b,k}). Target: the Berggren tree of coprime odd pairs (s,t), s>t. Map:
edge -> root pair (max, min). Preserved: coprimality (gcd(x,y)=gcd(x,b)=1),
oddness of both roots, primitivity of the triple, and one endpoint per child
(the source for x-type children, the target for y-type). Lost: the orientation
and the guard (a,b,k) (a pair can carry up to 40 guards in the window of S8).
Sidecar: the transport table (2)-(3), the parent table (5) and the braid word
(6). Test: the child pairs of every k>=2 Collatz edge, all of which leave the
multiplier 3 by identity; and the five sporadic readings (4), which show that
the legal sub-forest is not closed under any tree operation beyond one step,
except inside the root cluster.

What this adds to THM-3756 and the odd_square note: THM-3756 gives the tree and
its cones on all pairs; the odd_square note gives the edge encoding, the angle
bounds, the C_edge law and one illegal parent. This lane adds (i) the exact
guard transport along every tree edge, in both orientations, including the
identification of the Collatz k=1 parent as a halving edge and of the opus
parent as the affine parameter shift b''=(x+3b)/2; (ii) the proof that the
multiplier-3 language is generic {B2,B3} once above k=1 plus a finite root
cluster, with the complete sporadic list; (iii) the braid R as B2 then powers
of B1 along the t=y ray. None of this is a descent argument: the only integer
quantity that decreases along a legal tree step is the Berggren height of the
pair, and Theorem 3.2 says legal steps do not chain.

## 9. Reproduction block

```text
cd <worktree>
python3 04-computation/experiments/collatz_mod6_20260921_berggren_edge_transport.py \
    > 05-knowledge/results/collatz_mod6_20260921_berggren_edge_transport.out
python3 -O 04-computation/experiments/collatz_mod6_20260921_berggren_edge_transport.py | diff - \
    05-knowledge/results/collatz_mod6_20260921_berggren_edge_transport.out   # identical
python3 04-computation/experiments/collatz_mod6_20260921_berggren_edge_transport_audit.py \
    > 05-knowledge/results/collatz_mod6_20260921_berggren_edge_transport_audit.out   # ALL AUDIT CHECKS PASSED
```

Runtime a few seconds, RAM well under 1 GB; sympy for the symbolic identities,
numpy for the pair censuses; explicit `raise` only; the run ends with
`ALL CHECKS PASSED`. Every number quoted above appears in the `.out` or, for
the audit figures (A4 witnesses, A6 component census, A11 ratio intervals), in
the audit `.out`; the audit script runs in about one second and is identical
under `python3 -O`.

## 10. Audit record (2026-09-22)

Verdicts per claim (audit items A1-A11 in the audit `.out`). S1 root form:
CONFIRMED by the triple matrices on 4141 pairs. S2 transport table: CONFIRMED
by solving the child's guard from the child pair on 631320 children (a<=31,
|b|<=31). S3 classification (3): CONFIRMED for all odd a<=129, k<=14. S4
orientation = valuation: CONFIRMED (x<=2*10^5). S5 sporadic list: CONFIRMED by
an independent exact solve (k,j<=60) and brute force to max(x,y)<=20000 (16665
edges, 6666 with k=1, 13332 generic children, the same 5 sporadic readings).
S6 sub-forest: CONFIRMED, with the component census above. S7 parents:
WEAKENED as first written ("never a (3,+-1) edge") and corrected: three
accidental legal parents, 3->5, 5->7 and 5->1, all in the root cluster. S8
census: CONFIRMED number for number; "overlap" wording sharpened. S9 hypotenuse
counts and C_edge: CONFIRMED. S10 braid word: CONFIRMED on 10983 fibre words
(y<1000, k0<=11, j<=3) by direct B1 counting; "fibre is a subset of the ray"
sharpened to "fibre minus its k=1 element". S11: CONFIRMED and upgraded from a
sketch to a proof. Citations: THM-3756 (25) and (28) match the maps and cones
used here (checked against the canon file); the odd_square edges/synthesis/
inverse sections and arithmetic_braids section 2 say what is attributed to
them. Attributions to "the lead" (now "the brief") are not repository
citations and are marked UNCITED-RECOLLECTION as to wording.

Edits made: status line (audit sentence, "lead" -> "brief"); inheritance
paragraph (provenance of the brief); concept board (hostile qualified "by
identity", near miss reworded so that the opus statement is not called wrong,
sidecar and table row qualified "apart from the k=1 element"); Theorem 2.1
("always" qualified by the sign SCOPE); Corollary 2.2 (smallest witnesses
b=-7, b=9); Theorem 3.2 (component census); section 4 (accidental legal
parents, k>=4 exception 5->1); section 5 (overlap wording, twice); section 6
(ray statement, "quadruples"); section 7 (proof); reproduction block; this
section. The explorer's script and `.out` were not changed (all checks
correct; rerun identical under python3 and -O).

## Stopping boundary / next question

Stopped at: the complete transport and parent tables for all (a,b,k), the exact
multiplier-3 sub-language for b=+-1 (generic plus five sporadic readings), and
the braid-as-word law. Not done: the sporadic classification for |b|>=3 (the
S8 census has the data; the linear-equation method of Theorem 3.1 applies
verbatim, with the y>x cases at k=2 no longer vacuous), and a statement of which
guards (a,b,k) with a>3 have generic legal children of the same multiplier
(Corollary 2.2 gives the rule a=2^j+-1 at k=1, a=2^(k-1)+1 at k=2).

Next question: the y-type child of a k>=2 Collatz edge is the generalized edge
3u+b = 2(2^(k-1)+3) y. Is the set of generalized edges with a fixed odd
cofactor m (here m = 2^(k-1)+3) closed under some Berggren word, i.e. is there
a transport table on the (a,b,k,m) guards in which the Collatz edges are the
m=1 slice and the tree acts with bounded depth? The transport identities in (2)
already hold verbatim with 2^k replaced by 2^k m, so the question is only which
words return to m=1.
