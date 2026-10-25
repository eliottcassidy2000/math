# A cyclic-prefix repair at the actual third port

**PROVED, n>=11.** This is a sufficient three-black detour for an arbitrary
actual third port. The complete endpoint-four layer and all-order black-core
connectivity remain **OPEN**. All vertices remain literal, so no clone merger,
isomorphism quotient, or port identification is used.

## Inheritance and exact scope

The inherited [second-port/prefix packet](../../ports_owners_f1_20261025/endpoint.md)
closes `I={2,y,c}` for n>=11 and supplies a guarded growing-prefix repair.
The [odd-prefix packet](../../endpoint_four_frontier_20261021/README.md)
handles certain early chords. The new operation moves the whole prefix ending
immediately before the **actual** port h; its useful chord can be near the far
end of the word. The nearest prior ledger is Case A of the second-port packet.
The canonical hostile is the order-seven pure-blue target in
[the low-order audit](../../shape_ports_20261025/tournaments.md); it prevents an
unfiltered small-order endpoint-four theorem. A failed fixed cut is not such a
hostile. The least-used sidecar is the exact four-point endpoint neighbourhood
under the product of two positional reversals.

The board is: actual exceptional ports; cyclic prefix exchange; path-edge XOR;
four-port orbit obstruction; intermediate anti-reversals; and legal guards.
The map from the ordered graph chart to its re-ordered chart preserves every
arc through the trijection's affine bit formula, rather than preserving only
an unlabelled graph or a score statistic.

Let `P=(0,1,...,n-1)=(a,b,x,3,...,e,y,c,z)` be directed, closed and blue in T,
where `a=0,b=1,x=2,e=n-4,y=n-3,c=n-2,z=n-1`. Assume

```
n>=11,  2<=h<=e,  I={h,y,c},  c->b,
N^-_T(a)={z} union I,  N^+_T(z)={a,b,x,n-1-h},
3->c,  (n-h)->x.                                           (G)
```

The two new guards are actual arc conditions; neither can be omitted.
Put `W=(1,...,h-1)` and `H=(h,...,e)`, both nonempty. Then

```
Q1=(z,a,W,H,y,c),
Q2=(z,a,H,W,y,c),
Q3=(a,H,W,y,c,z)                                           (L)
```

are directed black words in their successive tournaments, and their three
K moves end at the literal labelled tournament `K_P(T)`.

## Proof of legality and target

Q1 is a cyclic cut of P. Write `S=K_Q1(T)` and `B=K_P(T)`.
The only new Q2 joins are `a->h`, `e->b`, and `(h-1)->y`.
The first reverses the incoming port `h->a`. P-blue gives
`T(b,e)=T(3,c)=1` and `T(y,h-1)=T(n-h,x)=1`; both pairs were
nonpath pairs of Q1, so the other joins hold at S. All within-block joins
are retained. Q3 removes `z->a` and restores `c->z`, the latter having been
complemented twice. Its other joins are retained by Q2.

Q1 differs from P only by replacing `cz` with `za`; Q3 differs from Q2 by
the reverse replacement. Therefore `E(Q1) xor E(Q2) xor E(Q3)=E(P)`.
Three complements and that path-edge identity give exactly B.

## Q1 is black already for n>=11

If P and its cyclic cut Q1 were both blue, the product of their reversals
would shift the labels by two modulo n and be an automorphism of T.
For odd n this action is transitive, so T is regular and
`n-5=(n-1)/2`, forcing n=9. For even `n=2m`, the even labels form one
transitive orbit, whose induced tournament is regular. The incoming
neighbours of a inside that orbit are c and possibly h: z and y are odd.
Thus `(m-1)/2<=2`, forcing n<=10. Both contradict n>=11.
This improves the generic cyclic-cut threshold specifically because both
y and c are retained as literal incoming ports.

## Q2 is black

Its positional reversal r sends a to y. It sends e to

```
w=2h-1       if 2h<=n-3,
w=2h-n+3     if 2h>n-3.
```

The Q1 path retains `e->y`, so `S(e,y)=1`. Except when n is even and
`h=(n-2)/2`, w is an ordinary vertex distinct from b and every member of
I. Consequently `S(a,w)=0`, contradicting the Q2-blue equality
`S(a,w)=S(e,y)`.

In the exceptional case `h=(n-2)/2`, the reversal sends 3 to `n-6` and
c to z. For even n>=14, `n-6` is outside the original outgoing ports
`{a,b,x,n/2}` of z. Both relevant pairs are complemented by Q1. The guard
`3->c` gives `S(3,c)=0`, while `S(z,n-6)=1`, another contradiction.

Only n=12,h=5 remains within the theorem. Here

```
Q2=(11,0,5,6,7,8,1,2,3,4,9,10).
```

If Q2 were blue, `S(1,9)=S(0,8)` would imply
`T(1,9)=T(0,8)=1`, since both are complemented nonpath pairs.
P-blue then gives `T(2,10)=1`, hence `S(2,10)=0`.
But Q2 reversal equates this with `S(11,7)=1`: vertex 7 is not a z-port.
Thus Q2 is black also in the sole exceptional small case.

## Q3 is black

Check at its other endpoint B. P remains blue there and
`N^+_B(a)={b,h,y,c,z}`. If Q3 were blue, `g=r_Q3 r_P` would be an
automorphism fixing both a and z, and hence preserving the exact four-point
set `N^+_B(a) minus {z}={b,h,y,c}`.

For `3<=h<e`, both `g(y)=Q3[h]` and `g(c)=Q3[h+1]` lie in the block
`H,W`, which is the cyclic order h,h+1,...,e,1,...,h-1 on {1,...,e}.
They are consecutive in that cyclic order. The only port values inside
the block are {1,h}. For `3<=h<e` those values are not consecutive in
either direction, so the two images cannot both be ports. If h=2,
`g(c)=4` is not a port (n>=11). If h=e, `g(y)=e-1` is not a port.
Every case contradicts invariance of the four-point set. This proves L.

## Genuine extension and remaining boundary

For h=3 the second new chord is simply `y->x`. Thus at every n>=15 the
previously unresolved `h=3, c->x, y->x` branch is now repaired whenever
`3->c`. This does not depend on the inherited growing-prefix condition
`c->4`. The checker stores a source satisfying `3->c,4->c` and rejecting
every applicable inherited odd-prefix chord; its detour is consequently
a real extension of those sufficient cases, rather than a new description
of a previously repaired source.

At `h=e`, the chord `(n-h)->x` becomes `4->x`; the old odd-prefix cut
already handles that chord at n>=15. Thus L does **not** close the hard
`h=e, x->4` boundary. The subsequent companion [tail theorem](TAIL.md)
does close the entire `h=e` family at n>=11, by a different ledger and a
uniquely forced order-eleven exception. For general h, the remaining sources must additionally
fail at least one of `3->c` and `(n-h)->x`. The overlapping long-prefix
continuation in the sibling packet may reduce this further; no claim is
made here that all surviving chord patterns are impossible.

The order-ten target in the historical prompt is already closed by the
[audited endpoint-three packet](../../endpoint_three_order10_20261021/README.md)
and the independent all-layer census recorded in
[THM-4732's result packet](../../black_detours_closed_blue_lines_klein_dictionary_20261020.md).
There is no new order-ten census here. The complete endpoint-four layers
at orders 11..14 remain open; this theorem and the tail theorem repair
specified subsets.

## Reproduction and controls

Run `python check.py` and `python -O check.py` from this directory.
The checker exhausts the stated guarded sources at orders 9 and 10 as
diagnostics, and replays deterministic controls at n11..60 for every h.
Arc-set and bit-row evaluators independently check intermediate directedness,
actual anti-reversals, and literal target. It checks the proof witnesses,
the n12 exceptional identity, the newly covered source, and both omitted-guard
hostiles. Small diagnostic success does not lower the proved n>=11 scope.
