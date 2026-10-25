# Consecutive-tail repairs in every endpoint layer

**PROVED:** let k>=4 and n>=2k+3. Every closed blue marked path
`P=(0,1,...,n-1)` whose incoming endpoint ports form the consecutive tail

```
N^-_T(0)={n-k,n-k+1,...,n-1}
```

has three successive directed black K moves to the literal labelled target
`K_P(T)`. Thus its endpoint scores are `(n-1-k,k)`. There is no additional
middle, module, clone, strongness or isomorphism restriction. General
endpoint layers and global black-core connectivity remain **OPEN**.

The k=4 case is the [tail theorem](TAIL.md). The same ordered-graph surgery
extends with its two literal anti-reversal witnesses. Its only candidate
boundary obstruction is a uniquely forced circulant tournament, and a
uniform second surgery repairs that entire exceptional family.

## 1. Two complementary chord directions

Put `a=0,b=1,x=2,d=3,e=n-k,z=n-1`,
`U=(4,...,e-1)` and `R=(e,...,z)`. Blue reflection gives

```
N^+_T(z)={0,1,...,k-1},
N^+_B(a)={b} union {e,e+1,...,z},   B=K_P(T).             (1)
```

If x->4, use

```
Q1=(a,b,x,U,R,d),
Q2=(b,x,U,a,R,d),
Q3=(b,x,d,U,a,R).                                       (2)
```

All blocks are nonempty in the stated range. Q1 uses x->4 and z->d;
Q2 inserts a between e-1 and e, using the complements of a->(e-1)
and e->a; Q3 restores x->d and d->4 after two complements. Boundary
pair cancellation gives `E(Q1) xor E(Q2) xor E(Q3)=E(P)`.
Q1 reversal equates true `T(n-2,a)` with false `T(d,x)` if blue.
Q2 reversal at `S=K_Q1(T)` would equate false `S(z,b)` with true
`S(d,x)`. These prove its first two words black for every k>=4.

For Q3 at B, the product `g=r_Q3 r_P` has the following action:

```
g(i)=i+1       for 0<=i<k,
g(i)=i+2       for k<=i<=e-3,
g(e-2)=0, g(e-1)=e,
g(i)=i+1       for e<=i<=n-2,
g(z)=k+1.
```

This follows by composing the positive cyclic rotations on the first and
last e positions. For odd n it is one n-cycle. If Q3 were blue, B would
therefore be regular; (1) forces `n=2k+3`. For even n the a-orbit is

```
{0,1,...,k,k+2,k+4,...,e-2},
```

of size n/2. Only b from the outneighbour set (1) lies in it. Orbit
regularity would force `(n/2-1)/2=1`, or n=6, impossible here. Thus Q3
is black except for a possible source at the single odd order 2k+3.

If instead 4->x, put `V=(4,...,n-2)` and use

```
(a,x,d,V,z,b), (d,a,x,V,z,b), (d,a,b,x,V,z).              (3)
```

These are the inherited ell=1 odd-prefix words with literal ports retained.
Their legality follows from a->x, a->d, z->b and 4->x; their edge XOR
is E(P). The first reversal equals that of a cyclic cut. If it were blue,
shift by two would be an automorphism of T. At odd n this forces n=2k+1.
At even n, the even-label orbit contains exactly floor(k/2) incoming
neighbours of a, so its regularity forces
`n=4 floor(k/2)+2<=2k+2`. Both are below the stated range.
After the first move, a has outdegree k+1 and z indegree k; the second
reversal pairs them and cannot be anti. The final reversal product with
r_P has two disjoint four-cycles, so cannot be a tournament automorphism.
Thus (3) is a three-black detour throughout the stated range.

## 2. The unique candidate at n=2k+3

In the exceptional odd order the permutation g in section 1 is the cycle

```
C=(0,1,...,k,k+2,k+3,...,2k+2,k+1).
```

If Q3 is blue, g-invariance and the complete a-row (1) determine every
arc of B uniquely. In the cyclic coordinates given by C, it is

```
B(u,v)=1 iff (pos(v)-pos(u)) mod (2k+3) belongs to
D={1,k+2,k+3,...,2k+1}.                                (4)
```

D contains exactly one of every nonzero pair {t,-t}, so (4) is a tournament.
In these coordinates r_P acts as `t -> 2k+1-t`, making P blue. Every
consecutive P arc has difference 1 or k+2, hence is directed in B.
The source `T=K_P(B)` has the prescribed tail ports, closed P and x->4.
Thus the exception exists, is unique as a labelled source for each k,
and must be repaired rather than dismissed by a score argument.

## 3. A uniform repair for the whole exceptional family

Write

```
A=(0,...,k), u=k+1, v=k+2, w=k+3, s=k+4,
Z=(k+5,...,2k+2).
```

At the source T defined by (4), use

```
C1=(A,w,u,v,s,Z),
C2=(A,w,u,Z,v,s),
C3=(A,u,Z,v,w,s).                                      (5)
```

Here Z is nonempty for every k>=4. The new C1 joins k->w,w->u,v->s
are true in T: their B cyclic differences are 2,k,2, respectively,
which do not belong to D; they are all nonpath pairs of P. C2's new
joins u->first(Z) and z->v are supplied by complementing the true T
arcs first(Z)->u and v->z. Their B differences are k-2 and k, again
outside D. C3 restores the original path joins k->u, v->w and w->s
after two complements; every other join is retained. Direct boundary
cancellation gives the same E(P) XOR identity.

All three words have uniform literal anti-reversal failures:

* C1 at T maps `a->v` to `(k-1)->z`. Their values are one and zero,
  since v is an ordinary outgoing neighbour of a and k-1 is an outgoing
  port of z.
* C2 at `S1=K_C1(T)` maps `a->w` to `w->s`. The first is one because
  w->a was complemented; the second is zero because the original path
  w->s was omitted by C1 and complemented.
* C3 at B maps `a->3` to `z->s`. The first is zero by (1); the second
  is one, since its cyclic difference is k+5, which belongs to D for
  every k>=4 (at k=4 it equals the upper endpoint 2k+1).

Thus (5) repairs the uniquely forced candidate for every k>=4. Together
with (2) and (3), this proves the all-k theorem at every n>=2k+3.

## 4. Scope and verification

This is a family of exact legal surgeries, not a quotient connectivity
argument. Every labelled vertex, original arc, path-edge set and actual
anti-reversal remains present. The missing general step is to treat
nonconsecutive port sets; the order threshold is sufficient and is not
claimed sharp.

Run `python general_tail_check.py` and `python -O general_tail_check.py`.
The checker retains local arc-set/bit-row evaluators, tests both chord
directions for k4..20 and n2k+3..2k+22, reconstructs the exceptional
circulant family directly for k4..80, and checks every displayed witness.
It also records that merely translating the old six-vertex numerical
certificate to the last six labels fails already at k=5; the successful
repair retains the moving junction k and the full Z block.
