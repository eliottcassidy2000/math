# Closed-middle endpoint-three repairs

Status: **PROVED** lemmas with the explicit guards below. This is a partial
structural reduction, not a full endpoint-three theorem or a bridge claim.

## 1. Inheritance and retained coordinates

The closest mechanism is the [endpoint-two closed-middle repair](../endpoint_two_repairs_20261020g/theory.md).
Its common-final anti-automorphism test and exact three-K path ledgers survive,
but its small neighbourhoods gain one member. The hostile boundary is an
actually blue word, not the failure of a sufficient score test. The least-used
sidecar is the permutation product of two actual anti-reversals. Our live
board is: mirrored exceptional ports; path-edge XOR; small exact neighbourhoods;
permutation orbit capacity; and positional failure families.

Write `P=(a,b,W,c,z)`, `W=(x,...,y)`, for a closed blue Hamiltonian path.
All degrees refer to T unless a subscript says otherwise. Assume
`d(a)=n-4,d(z)=3`. There is a two-element middle set I such that

```
N-_T(a)={z} union I,   N+_T(z)={a} union J,   J=r_P(I).
```

The sets I and J may overlap. The standard path forbids b in I and c in J.
In this note the middle is closed: `c->b`. Set
`epsilon=[b->y]=[x->c]`; the equality follows from the blue reflection.
We use the standard P coordinates `a=0,b=1,x=2,...,y=n-3,c=n-2,z=n-1`.
No module hypothesis is imposed on W.

For a directed path Q, `K_Q(T)=T^op xor E(Q)` complements the nonpath
pairs. Three K moves whose path-edge XOR is E(P) end literally at K_P(T).
The inherited ledgers below have this identity, independently of all arc
orientations. Blue means the actual positional reversal is an anti-automorphism.

## 2. A common final word is black above order ten

**Lemma (PROVED).** Suppose `c notin I`. At `B=K_P(T)`, the directed word

```
F=(a,z,b,W,c)
```

is black for every n>=11.

**Proof.** The first two joins are complemented closing pairs; `b notin J`
follows from `c notin I`. In B,
`N+(a)={b,z} union I`, so d_B(a)=4. If F were blue, its endpoint c
would have B-indegree4 and therefore T-outdegree4. Since c loses to a,
beats z, and beats b, its middle outneighbours have the form `{b,t1,t2}`.
The preceding path vertex y is not one of them. Hence
`N-_B(c)={y,b,t1,t2}`. As `r_F(z)=y`, anti-reversal gives

```
r_F({b} union I)={b,t1,t2}.
```

The vertex `r_F(b)=w=P[n-4]` is not b. Thus some member of I must be w.
Since P remains blue at B, `g=r_F r_P` would be an automorphism of B.
Its explicit action is

```
g(0)=n-3, g(1)=0, g(2)=n-1,
g(i)=i-2 for 3<=i<=n-2, g(n-1)=n-2.
```

Following the parity strings gives one n-cycle for odd n, or two n/2-cycles
for even n. In the odd case B would be regular, forcing n=9. If n=2m,
the orbit containing a is `A={0,1,3,5,...,n-3}`. It contains b but neither
z nor w=n-4. Therefore a has at most two outneighbours in A. Transitivity
on A makes B[A] regular, so `(m-1)/2<=2` and n<=10. Both cases contradict
n>=11. QED.

The same argument works for endpoint scores `(n-1-k,k)` with k>=3:
I has size k-1, B(a) has degree k+1, and F-blue forces w in I. The odd
exception is n=2k+3, and the even bound is n<=4k-2. Therefore F is black
whenever n>4k-2. This bound is sufficient; no sharpness claim is made here.

## 3. The epsilon-zero branch is fully repaired

**Theorem (PROVED).** If `c notin I` and epsilon=0, every n>=11 source
has the exact three-black K detour

```
Q1=(c,W,b,z,a),
Q2=(z,c,W,b,a),
Q3=(a,z,b,W,c).
```

**Proof.** The joins c->x, y->b and b->z are available; successive
complementation gives the inherited legal ledger. Its edge XOR is E(P).
After Q1, z and a have degrees n-4 and2 respectively, so Q2 is black by
its endpoint sum n-2. Section2 makes Q3 black at its other endpoint B.

For Q1, c beats b,x,z. Unless d_T(c)=3, its endpoint sum exceeds n-1.
If Q1 were blue, therefore `N+_T(c)={b,x,z}`. Its reversal sends x to z
and b to `t=P[3]`. Matching N+(c) with N-(a) forces `I={x,t}`.
The product `g=r_Q1 r_P` would be an automorphism of T, with action

```
g(0)=2, g(1)=0, g(i)=i+2 for 2<=i<=n-5,
g(n-4)=1, g(n-3)=n-1, g(n-2)=3, g(n-1)=n-2.
```

It is one n-cycle for odd n, and two n/2-cycles for even n. Odd regularity
would force `n-4=(n-1)/2`, hence n=7. For n=2m, the a-orbit is
`A={0,1,2,4,6,...,n-4}`. Of a's three inneighbours `{z,2,3}`, exactly
2 belongs to A. Thus regularity of T[A] forces `(m-1)/2=1`, so n=6.
Consequently Q1 is black already for n>=8. QED.

## 4. Epsilon one away from the second port

Assume `epsilon=1` and `c,x notin I`. The inherited words

```
Q1=(c,b,z,a,W),
Q2=(z,c,b,a,W),
Q3=(a,z,b,W,c)
```

are directed and have edge XOR E(P). The required a->x join is exactly
the additional guard x notin I. By section2, Q3 is black for n>=11.

**First-word lemma (PROVED).** Q1 is black at n>=8. Its product
`r_Q1 r_P`, if an automorphism, has cycle lengths n for odd n, two
cycles of length2m at n=4m, and lengths m,m,m+1,m+1 at n=4m+2.
This is the explicit permutation identity proved in the inherited note,
independent of the endpoint degree. Every even n gives even permutation
order, impossible for a tournament automorphism. Odd regularity forces
`n-4=(n-1)/2`, hence n=7.

**Second-word reduction (PROVED).** At n>=11, Q2 is black unless

```
I={P[4],P[5]}  or  I={P[4],y}.
```

**Proof.** Put S=K_Q1(T). If Q2 were blue, its endpoint equality gives
`d_T(y)=n-5`, because d_S(z)=n-4 and d_S(y)=n-2-d_T(y).
Thus y has T-indegree4. Since c,x notin I implies b,y notin J,

```
N-_S(z)={b} union J.
```

Let w immediately precede y. The Q1 move retains w->y, so
`N+_S(y)=N-_T(y) minus {w}`. This set contains b, and also a whenever
y notin I. The reversal of Q2 sends b to P[n-5] and a to P[n-6].
Neither image is b for n>=11. Anti-reversal must therefore put
P[n-5] in J, hence P[4] in I. If y notin I, it must also put P[n-6]
in J, hence P[5] in I. If y in I, the pair is instead {P[4],y}.
These exhaust the alternatives. QED.

**Exceptional-pair exclusion (PROVED, n>=11).** Neither displayed pair
actually permits a blue Q2 at n>=11.

**Proof.** Let `t=P[n-10]`, `s=P[n-7]`, and write a directed-pair value
as T(v,w) in {0,1}. Under the Q2 reversal, the images needed below are

```
r(z)=y, r(t)=P[5], r(a)=P[n-6], r(s)=x.
```

Hence blue Q2 would imply

```
S(z,t)=S(P[5],y),  S(a,s)=S(x,P[n-6]).
```

All four pairs are nonpath pairs of Q1 at n>=13. Blue P also gives
`T(x,P[n-6])=T(P[5],y)` and `T(a,P[9])=T(t,z)`. Substituting the
nonpath complements into the two equations gives

```
T(a,P[9])=1-T(a,s).
```

In both exceptional port pairs, P[9] and s are outside I, distinct from
z, for n>=13. Both values in this last equation are therefore1, a
contradiction.

At n=11, use Q2's anti-equations on the pairs `(c,P[5])` and `(b,P[5])`.
Their right sides are `(a,P[7])` and `(a,P[6])` respectively. All four
are nonpath pairs of Q1, and a beats P[6],P[7]. Consequently
`T(c,P[5])=T(b,P[5])=1`, contradicting blue P, which identifies
`T(b,P[5])` with `T(P[5],c)`.

At n=12, Q2's anti-equations on `(z,b)`, `(c,a)`, `(P[2],P[4])`
have right sides `(P[7],P[9])`, `(P[6],P[8])`, `(P[3],P[5])`.
The pair b->z is retained by Q1, while c,a and all three right sides
are nonpath pairs; a->c is true. Hence these equations imply
`T(7,9)=1`, `T(6,8)=0`, and `T(2,4)=T(3,5)`. Blue P identifies the
last two values with T(7,9) and T(6,8), a contradiction. This handles
both exceptional port pairs without enumerating their interiors. QED.
Thus the three words in section4 are all black at n>=11 whenever
c,x notin I.

## 5. A second-position port needs only one excluded pair

Suppose epsilon=1, `x in I`, and `I!={P[2],P[3]}`.
Write `P=(a,b,x,U,y,c,z)`, with U nonempty. At every n>=9 the following
is an exact three-black K detour:

```
Q1=(a,b,x,U,z,y,c),
Q2=(a,x,U,z,y,b,c),
Q3=(a,x,U,y,b,c,z).
```

**Directedness and target.** The sole changed boundary join from the
endpoint-two proof is last(U)->z. It is present exactly when P[3] is
not in I, which is our excluded-pair guard. The joins x->a and z->y
are available, and all other required retained or twice-complemented
joins are as in the inherited proof. The three path-edge sets XOR to
E(P), so the target is exactly B=K_P(T).

**First two blackness witnesses.** The Q1 reversal pairs (a,z) with
(x,c), whose source values are respectively false and true. The Q2
reversal pairs (a,b) with (x,c); at that source the first was retained
and is true, while the second was complemented and is false. Both
contradict the anti-automorphism equation.

**Final-word orbit proof.** Test Q3 at B. If blue, `g=r_Q3 r_P` would
be an automorphism fixing a and z. On the middle positions it has one
cycle of length n-2 when n is odd and two cycles of length (n-2)/2 when
n is even. This follows directly from

```
g(a)=a, g(z)=z, g(1)=2,
g(i)=i+2 for 2<=i<=n-5,
g(n-4)=1, g(n-3)=n-2, g(n-2)=3.
```

Because g fixes a it must preserve the three-element middle set
`N+_B(a) intersect M={b} union I`. Such a set cannot be a union of the
displayed orbits at n>=9: the odd orbit has size at least7, and in the
even case each orbit has size at least4. This proves Q3 black. No
assumption on any arc inside U enters.

## 6. The adjacent pair of incoming ports

Suppose epsilon=1 and `I={P[2],P[3]}`. Put `e=P[n-4]` and
`U=P[3:n-4]`, so U is nonempty at n>=8. The two outgoing exceptional
ports are exactly e,y. The following ledger, proposed independently by
the bounded-template lane, inserts z immediately before both:

```
Q1=(a,b,x,U,z,e,y,c),
Q2=(a,x,U,z,e,y,b,c),
Q3=(a,x,U,e,y,b,c,z).
```

It is an exact three-black detour at every n>=9. The only joins changed
from section5 are last(U)->z and z->e, both available: last(U) is neither
e nor y. The remaining joins use x->a, b->y, c->b and retained path
arcs. The new edge pairs last(U)z and ze each occur twice, while the
surviving last(U)e and ey are original P edges, so the edge XOR is E(P).

For Q1, reversal pairs `(a,y)` with `(b,c)`: a->y is true and b->c is
false. For Q2, reversal pairs `(a,b)` with `(x,c)`; at that source the
first is true and the second false, as before. Q3 is literally the same
relocated-b word as in section5, so its orbit proof applies unchanged.
This also explains why treating the two outgoing ports as one port would
have broken the old ledger's directedness.

## 7. Closed-middle theorem and remaining boundary

**Theorem (PROVED).** At every n>=11, every closed-middle endpoint-three
source with `c notin I` has an exact three-black K detour to K_P(T).
This combines sections3--6, with no additional strongness assumption
on the partner and no restrictions on overlaps between I and J.

## 8. The paired last-position and second-position ports

**Theorem (PROVED).** If the middle is closed and I={x,c}, there is an
exact three-black detour for every n>=9. The epsilon-zero branch below
already works for n>=6; the epsilon-one branch is section5, whose proof
never needed c notin I.

For epsilon=0, put `D=P[3:n-3]`, allowing D empty at n=6. The
bounded-template lane supplied the following ledger and its last arc
witness:

```
Q1=(a,b,x,D,z,y,c),
Q2=(b,x,D,z,y,a,c),
Q3=(b,x,D,y,a,c,z).
```

The required source boundary arcs are `c->a,b,x,z`, `z->a,b,y`,
`a->y`, `x->z`, and the P path. Q1 skips y and uses z->y; Q2 gains
y->a and a->c by complementing their nonpath pairs. The last(D)y and
cz pairs (or xy when D empty) are complemented twice and hence return
for Q3. The surviving path-edge XOR is exactly E(P): the added pairs
last(D)z,zy,ya,ac occur twice, with x replacing last(D) when D is empty.

Since c beats a,b,x,z, Q1's endpoint excess is d(c)-3>=1. At the
second source, Q2's endpoint sum is
`(n-1-d(b))+(n-2-d(c))=n-2`, because P-blue gives d(b)+d(c)=n-1.
Finally Q3 reversal pairs `(b,c)` with `(x,z)`. At its source, b->c is
false, having been complemented twice; x->z is true, having been
complemented twice if D is nonempty or retained twice if D is empty.
Thus Q3 is also black. All arguments retain arbitrary D.

## 9. A movable cut repairs every nonadjacent last-position pair

**Theorem (PROVED).** Suppose the middle is closed, `I={h,c}`, and
`2<=h<=n-4`. Every n>=6 source has an exact three-black K detour.
There is no guard on the chord x,c or on any interior pair.

Put `k=n-1-h` and split the middle before k:

```
L=(b,P[2],...,P[k-1]),  R=(P[k],...,y).
```

Here L has at least two vertices and R is nonempty. The extra outgoing
ports of z are b,k. Use

```
Q1=(a,L,z,R,c),
Q2=(L,z,R,a,c),
Q3=(L,R,a,c,z).
```

**Directedness.** Q1 uses last(L)->z, valid because last(L) is neither
b nor k, followed by z->k. It otherwise follows P. Its K move supplies
y->a, because h!=y, and a->c, because c->a in T; all other Q2 joins
are retained. The two complements restore last(L)->k and c->z for
Q3, while its other joins survive Q2.

**Target.** The added pairs last(L)z, zk, ya, ac each occur twice.
The surviving boundary pairs are ab, last(L)k, yc, cz, and all internal
L and R path pairs occur three times. Thus the three path-edge sets
XOR to E(P), giving exactly B=K_P(T).

**First word.** The source pair c->b is true by the closed-middle
hypothesis. Under Q1 reversal it maps to y->a, which is false because
h!=y. This proves Q1 black, without a degree estimate.

**Second word.** Q2 has endpoints b,c. At its source their degree sum is
`(n-1-d_T(b))+(n-2-d_T(c))=n-2`, using d_T(b)+d_T(c)=n-1. It is black.

**Final word.** Its reversal pairs `(b,c)` with `(x,z)`. At the third
source b->c is false, since it was omitted and complemented twice.
The original x->z is true because k>=3, and this pair was retained
twice if k=3 or complemented twice otherwise. The anti equation fails.
This short witness, identified by the companion probe lane, removes
the threshold required by the following independent orbit test.

At B the word is independent of the cut: `Q3=(b,x,...,y,a,c,z)`.
If it were blue, `g=r_Q3 r_P` would be an
automorphism of B. Its explicit standard-position action is

```
g(0)=1, g(1)=2, g(i)=i+2 for 2<=i<=n-5,
g(n-4)=0, g(n-3)=n-2, g(n-2)=n-1, g(n-1)=3.
```

For odd n it is an n-cycle, forcing regularity and n=9 because
d_B(a)=4. At even n=2m it has two m-cycles, with a-orbit
`A={0,1,2,4,...,n-4}`. Of `N+_B(a)={b,z,h,c}`, at most b,h lie in A.
Thus regularity of B[A] requires `(m-1)/2<=2`, hence n<=10.
This separately proves Q3 black for n>=11 without using the x,z arc.

The movable cut is essential: placing z immediately before its actual
extra outport k repairs the precise failed boundary incidence, while
the final word and its permutation obstruction stay unchanged.

**Degree-free boundary corollary (PROVED).** More generally let n>=6 and
P=(a,b,x,...,y,c,z) be any closed blue path, with

```
c->b,   z->b,   a->y,   d_T(z)>=3.
```

There is an exact three-black detour of the same form, without any
specified endpoint degree and with arbitrarily many additional ports.
Choose k to be the first middle outneighbour of z after b. It exists
because z beats a,b and has degree at least3. Blue P makes a->y
equivalent to x->z, so k is at position at least3. By its minimality,
the predecessor of k beats z. Use the same L,R split before k.
Blue P also gives c->a from z->b. Thus every directedness argument
above still applies. The three blackness proofs use only c->b, a->y,
the reflected equality d(b)+d(c)=n-1, and x->z. Additional exceptional
ports do not occur in any of those tests. This proves the corollary
directly; the endpoint-three orbit bound is not used in it.

**Combined frontier (PROVED reduction).** Sections3--9 repair every
closed-middle endpoint-three source at n>=11 except the adjacent last
pair `I={y,c}`. That case is completed by root R3 and lift Theorem E2;
see the [full theorem and exhaustive partition](README.md).

## 10. Reproduction and precise finite scope

`theory_verify.py` has its own tournament construction, path predicate,
anti-reversal predicate, and K operation; it imports no inherited checker.
It exhausts all 12,96,640,7680 marked endpoint-three presentations at n6,7,8,9
and replays exactly the stated branch guards. Small legal words may be
blue below the theorem thresholds, and are counted rather than erased.
It also constructs 6072 deterministic larger sources at n11--60, with no
strong-partner or mixedness filter, and checks the permutation orbit
identities through order200. There are 7031 guarded ledger replays in
total. The degree-free corollary is independently checked on every
standard-blue presentation of orders6--9, with 4090 qualifying sources
and last-vertex degrees3--6 represented. Normal and optimized Python
both pass. These controls are separate from the
all-order proofs, and do not extrapolate their conclusions.
