# One-contact sources: an exact divisor sieve and a forbidden strip near R₊=2H

**Scope:** one outcontact u at the apex v, with absolute targets H
and R=R₊. This is not a restriction on witnesses allowed to change
their number of contacts, and it does not settle scalar Chan–Pak.

Set d=R−H and e=2H−R=H−d. Assume d>0. [The finite-support
identity](ORDINAL_DEFECT.md) says d is exactly the number of full
puncture paths ending at u. Hence u is in the terminal strong
component A of the puncture. Put m=H(G), where G is the ordinal
prefix before A, with H(empty)=1. This prefix beats both A and v.
It contributes a Cartesian product factor to all three counts:
(H,R,d)=m·(H',R',S), where H',R' are the outputs from A with the
one-contact apex, and S counts full paths of A ending at u.
In particular m is an odd divisor of gcd(H,d). Conversely every
realized scalar m is a physical prefix factor.

## Near the upper boundary

Let a=|A|, h=H(A), P be its exposed start response for {u}, and
t be the outdegree of u in A. Since R'=h+P and d'=S, we have

```
e/m = (h−S) + (P−S).
```

Each term is nonnegative. The support {u} contributes h(A−u)
to P, at least S by deleting u from paths that end there. If
A is nontrivial, h−S≥a−1 because every other vertex can end a
full path (rotate a Hamiltonian cycle). More is true for the
start term. Every nonempty subset of the t outneighbors can
follow u as a path, giving at least 2^t−1 extra supports. For
each of the a−1−t other vertices take one shortest path from u.
These have distinct supports: distance from u increases by
one at each step, and the endpoint is determined by its unique
maximal distance in the support. None belongs to the outneighbor
subsets, since its endpoint is not an outneighbor. A complement
always supplies at least one path. Thus

```
e/m ≥ 2(a−1) + 2^t−1−t ≥ 2(a−1)     for nontrivial A.
```

It follows that R≤2H. At equality A is a singleton, and
**(N,2N) occurs with one contact iff N occurs as a scalar**:
append an apex/contact edge to any scalar witness. That
construction is ordinal, not generally strong.

At any fixed *positive* e the allowed H form a finite set.
Indeed a≤1+floor(e/(2m)) and S≤(a−1)!, so

```
H = e + mS ≤ e + m·floor(e/(2m))! ≤ e + floor(e/2)!.
```

For the last step put k=floor(e/(2m)); a≥3 gives k≥2, and
(mk)!/k!≥m for m≥1. Therefore every fixed positive offset
from the boundary is eventually forbidden, even without fixing
the order. Equivalently, in any unbounded sequence with positive
deficit, e grows at least on the order of log H / log log H.
The boundary itself is precisely the scalar spectrum, so the
two behaviors are genuinely different.

The [standalone replay](replay_one_contact_upper.js) classifies
the strip 1≤e≤13. At m=1 the proof restricts to strong A of
order at most seven. For any m>1 the first nontrivial deficit
is at least 4m, so m=3 is the only additional option in the
strip, and it adds no distinct pair. The exact pairs (H,R) are

| e | Pairs |
|---:|---|
| 4 | (5,6) |
| 7 | (9,11) |
| 8 | (9,10) |
| 10 | (11,12) |
| 11 | (15,19) |
| 12 | (13,14), (15,18) |

Thus (63,126) occurs, but no (63,126−e) with 1≤e≤13 occurs
with one contact. In any non-boundary one-contact realization
at H=63 the scalar prefix factor must be one: every proper
nontrivial factorization 63=m·H' uses 7 or 21 on one side,
and m=63 would force H'=1 and e'=0. This makes the 7/21/63
distinction a structural constraint on the choice of source,
not just a coincidence of residues.

## A bound working from the other side

For arbitrary d>0 with nontrivial A, delete u and write D=A−u.
Let d'=d/m and R'=R/m. Then d' is the number of full paths in D
ending at an inneighbor of u. If D is strong, at most d' vertices
are inneighbors of u; every other vertex is a contact from u.
All their nonempty subsets contribute to the new P_D, so
R'=h_D+K_D+P_D≥2^max(|D|−d',0). The identity for R' follows
by stripping u from each subset path counted by P.

If D is not strong, let F be its first strong component.
Because A is strong, u beats at least one vertex of F.
Therefore P_F≥1, and the ordinal transfer gives
P_D≥P_F·2^(|D|−|F|). Also h_F divides d' by factoring a
full D path at the component boundaries. Since |F|≤h_F,
this again gives the bound

```
|A| ≤ 1 + d' + floor(log₂ R').
```

At prescribed H,R test each odd divisor m of gcd(H,d).
Check m is a scalar count and compare actual marked terminal
components up to the smaller of the two bounds. For m=1 use
the empty prefix. For m>1 the prefix check itself is finite:
remove its singleton strong components, use that each other
component has order at most its count, and use that the sum
of those counts is at most their product m. This is an exact
witness-or-rejection procedure for the one-contact branch, not
a claim that the resulting witness always exists.

The divisor factorization and two size bounds originated in a
background pass. The strengthened outdegree term, fixed-positive-
offset obstruction, and complete strip through 13 are derived
here. Run **node replay_one_contact_upper.js** to replay the strip.
