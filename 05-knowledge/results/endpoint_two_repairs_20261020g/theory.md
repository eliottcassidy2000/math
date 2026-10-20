# Closed-middle repairs with two exceptional endpoint incidences

Status: **PROVED** symbolic three-black ledgers and their stated guards;
**FINITE-EXACT** bounded replay. All paths below are marked. No change of
Hamiltonian presentation or unrecorded isomorphism is free.

## 1. The retained ports and current board

The inherited [cyclic singleton extension](../kernels_games_detours_20261020f/detours.md)
handles closed blue paths with endpoint scores `(n-2,1)`. The
[two-block repair](../strong_blue_games_20261020f/detours.md) retains internal
endpoint-preserving paths rather than assuming their existence. Their near
misses are respectively discarding an exceptional boundary incidence and
using a fixed four-vertex reroute at larger block sizes. The least-used
sidecar here is the product of two actual positional reversals.

The live concepts are: marked blue reflection; the two exceptional ports;
literal three-K ledgers; degree-two neighbourhoods; and permutation orbit
obstructions. The anchor is the strong mixed-to-mixed blue case, not an
inference of black connectivity from failed local templates.

Let `P=(a,M,z)` be closed and blue with `d(a)=n-3,d(z)=2`. Write `r` for
its positional reversal. There is exactly one middle vertex u with `u->a`,
and exactly one middle vertex v with `z->v`; necessarily `v=r(u)`. Every
other middle vertex lies on `a->w->z`. Thus T is the cyclic extension of
the blue middle M with exactly these two boundary pairs reversed.
The path condition excludes u=first(M) and v=last(M). The ports may
coincide only at the fixed middle position. This is a full coordinate
description, not an assumed module.

For the rest of this note M is closed. Write

```
P=(a,b,W,c,z),   W=(x,...,y),   c->b,
epsilon=[b->y]=[x->c].
```

The equality follows from the blue middle reflection. Initially assume
`u!=c`, so `v!=b`. Every ledger keeps the entire arbitrary W path, or a
specified subpath of it; no interior arcs are erased.

**Main theorem (PROVED).** For every n>=8, every such marked blue path
with a closed middle and `u!=c` has a three-black K detour to its exact
labelled partner `K_P(T)`. Neither endpoint tournament needs an additional
strongness hypothesis. The second-position-port branch is valid already
at n>=7.

## 2. A common final word is black at n>=8

Let `B=K_P(T)`. Whenever `u!=c`, the word

```
F=(a,z,b,W,c)
```

is directed in B: a->z is the complemented closing arc, z->b is the
complement of b->z, and the middle path survives.

**Lemma (PROVED).** F is black at B for n>=8.

**Proof.** In B, `N+(a)={b,z,u}` and d(a)=3. If F were blue, its
endpoint c would have indegree3 in B and therefore degree3 in T. Since
u!=c, c beats z but not a, so it has exactly two outneighbors in M:
`{b,t}`. Here t is neither y (which precedes c) nor c. In B this gives
`N-(c)={y,b,t}`. Reversal of F sends a to c and z to y. Thus

```
{r_F(b),r_F(u)}={b,t}.
```

At n>=6, b is not fixed by r_F. Consequently `r_F(b)=t` and
`r_F(u)=b`, so `u=t=w`, where w is the vertex at position n-4 of P.

Since P is blue at B as well, `g=r_F r_P` would be an automorphism of B.
Its cycle structure is one n-cycle when n is odd, and two n/2-cycles when
n is even. To verify this directly, identify P with `0,...,n-1`; then

```
g(0)=n-3, g(1)=0, g(2)=n-1,
g(i)=i-2 (3<=i<=n-2), g(n-1)=n-2.
```

Following the two descending parity strings gives the asserted cycles.
For odd n, B would be regular; d_B(a)=3 would force n=7. For even
n=2m, the orbit containing a is

```
A={0,1,3,5,...,n-3}.
```

It contains b but not z or u=w=n-4. Therefore a has exactly one
outneighbor inside A. An automorphism transitive on A makes B[A] regular,
so `(m-1)/2=1`, forcing n=6. Both possibilities are excluded at n>=8.
QED. This proof does not assume the whole partner B is strong.

## 3. Closed middle with epsilon=0

Use the inherited cyclic ledger

```
Q1=(c,W,b,z,a),
Q2=(z,c,W,b,a),
Q3=(a,z,b,W,c).
```

Because epsilon=0 supplies c->x and y->b, and v!=b supplies b->z,
the three words are directed at their successive K sources. The same
boundary XOR as in the inherited proof gives the exact target K_P(T).

Here c beats b,x,z and not a, so d_T(c)>=3. Q1 has endpoint score sum
`n-3+d_T(c)>=n`, and is black. After Q1, z and a have scores n-3 and
1 respectively, so Q2 has endpoint score sum n-2 and is black. Q3 is
black by section 2, tested at its other endpoint B. This proves this
branch for n>=8. The sharper inherited neighbourhood argument also gives
n>=7; the present main theorem does not need that sharpening.

## 4. Closed middle with epsilon=1 away from the second port

Suppose additionally `u!=x`. Use

```
Q1=(c,b,z,a,W),
Q2=(z,c,b,a,W),
Q3=(a,z,b,W,c).
```

The only old cyclic-module join affected by the new ports is a->x in
Q1; u!=x retains it. The other required joins are c->b, b->z, z->a,
and the retained W path. Thus both first words are directed, the inherited
edge XOR is unchanged, and Q3 is directed at B. The composite ends at B.

### First word: an explicit permutation obstruction

If Q1 were blue, `g=r_Q1 r_P` would be an automorphism of T. Its cycle
lengths are

| n | Cycle lengths |
|---|---|
|odd|n|
|4m|2m,2m|
|4m+2|m,m,m+1,m+1|

This is an elementary permutation identity. For n>=8, follow

```
0 -> n-5, 1 -> n-3, 2 -> n-2, 3 -> 1,
4 -> n-1, 5 -> 0, i -> i-4 (6<=i<=n-3),
n-2 -> n-4, n-1 -> n-6.
```

The descending residue-class strings give the displayed lengths; n=6,7
are checked directly from the words. If n is even, at least one cycle
has even length, impossible for a tournament automorphism: an even-order
automorphism has a nontrivial involutory power, which swaps the endpoints
of some directed edge. If n is odd, the automorphism would be transitive
and T regular, forcing `n-3=(n-1)/2`, hence n=5. Thus Q1 is black for
every n>=6 in this endpoint layer.

### Second word: its two outneighbors retain the exceptional port

Let `S=K_Q1(T)` and suppose Q2 were blue. Its endpoint-score equality
forces `d_T(y)=n-4`, because `d_S(z)=n-3` and
`d_S(y)=n-2-d_T(y)`. Since v!=b and v!=y (the latter is u!=x),

```
N-_S(z)={b,v}.
```

Write w for the vertex immediately preceding y in W. Because epsilon=1,
both b and w beat y in T. The exact value d_T(y)=n-4 says y has exactly
three inneighbors in T.

If u!=y, the third is a, giving `N+_S(y)={a,b}`: Q1 retains w->y and
complements the other two incoming pairs. Anti-reversal of Q2 would
therefore require `{r_Q2(b),r_Q2(v)}={a,b}`. But for n>=8,
r_Q2(b) is the P vertex at position n-5, neither a nor b.

If u=y, let the third inneighbor be t; now `N+_S(y)={b,t}`. Since b
is not fixed, anti-reversal forces `r_Q2(b)=t` and `r_Q2(v)=b`, hence
t=v. Yet u=y gives v=x, while r_Q2(b) has P position n-5, different
from x at position2 for n>=8. Again a contradiction. Thus Q2 is black.
Q3 is black by section 2. This proves the branch.

## 5. The second-position exceptional port has its own local repair

Let n>=7, u=x, and epsilon=1. Put

```
P=(a,b,x,U,y,c,z),   U nonempty.
```

The relevant arcs are c->b, b->y, x->c, x->a and z->y, together with
P's path arcs, z->a, and `a->U->z`. Use the new ledger

```
Q1=(a,b,x,U,z,y,c),
Q2=(a,x,U,z,y,b,c),
Q3=(a,x,U,y,b,c,z).
```

**Directedness.** Q1 uses last(U)->z and z->y. Its K move complements
x->a, b->y and c->b, supplying a->x, y->b and b->c in Q2; all other
Q2 joins are retained. The two nonpath complements of last(U)->y and
c->z restore their original directions for Q3. Its other joins survive
Q2. Hence every word is directed, for an arbitrary interior U.

**Exact target.** Internal U edges occur three times. The added edges
last(U)z, zy, ax, yb and bc occur twice. The surviving boundary edges
are ab,bx,x first(U),last(U)y,yc,cz, exactly those of P. Thus the path
edge XOR is E(P), and three K moves end literally at K_P(T).

**Blackness from three actual arc witnesses.** No degree estimate is used:

1. Q1 reversal maps the pair (a,z) to (x,c). In T, a->z is false while
   x->c is true.
2. Q2 reversal maps (a,b) to (x,c). At its source, a->b was retained by
   Q1 and is true, while x->c was complemented and is false.
3. Write t=first(U). Q3 reversal maps (a,b) to (t,z). After Q1,Q2,
   a->b is false, while t->z is true: that pair was retained twice when
   U is a singleton, and complemented twice otherwise.

Each contradicts the necessary anti-automorphism equation. All three
moves are black. The nonempty-U boundary n>=7 is explicit; when U is
empty the last witness is absent, so this proof does not cover n=6.

## 6. Checks, failure boundaries, and remaining work

`theory_verify.py` replays every closed standard-blue middle at total
orders8..10, with every admissible u except the last middle vertex. It
checks actual directedness, actual positional blackness, endpoint scores,
and exact labelled targets. It separately checks every second-port source
at n=7, the two permutation identities through order200, and deterministic
larger samples. These tests validate the stated proofs; none supplies an
all-order existence inference by extrapolation.

`theory_probe.py` is the provenance search for section5: it keeps a
two-vertex U block, exhausts all 7! block permutations against all16
order-eight sources in that branch, and records successful ledgers. It
does not replace the arbitrary-U proof, and its block model is not a
claim that U is a module.

Reproduction:

```
python 05-knowledge/results/endpoint_two_repairs_20261020g/theory_probe.py
python 05-knowledge/results/endpoint_two_repairs_20261020g/theory_verify.py
python -O 05-knowledge/results/endpoint_two_repairs_20261020g/theory_verify.py
```

An additional independent integration receipt, `theory_unified_audit.py/json`,
replays the full case partition from the companion proofs at orders9 and10:
3,072 and 28,672 marked presentations respectively. No strong-partner or
mixedness filter is imposed. This includes all512 order-nine presentations
where the two exceptional ports coincide at the reflection-fixed vertex.
Every selected ledger is checked for all three black moves and its literal
labelled target, in both normal and optimized Python. This is a finite
cross-implementation audit of the combined theorem; the all-order claims
still come from their individual proofs.

```
python 05-knowledge/results/endpoint_two_repairs_20261020g/theory_unified_audit.py
python -O 05-knowledge/results/endpoint_two_repairs_20261020g/theory_unified_audit.py
```

The boundary `u=c` is outside this note's main theorem and must keep its
actual mirrored port `v=b`; none of sections2--4 may silently be applied
there. The companion last-port and block-lift proofs handle its two chord
branches. Open middles are handled by the companion synthesis. Failure of
any one ledger is not evidence of a bridge or of disconnected black support.
