# The whole three-port tail is repaired from order eleven

**PROVED:** every closed blue endpoint-four marked path with
`I={n-4,n-3,n-2}` has a literal three-black detour for every n>=11.
The middle may be open or closed. This closes the entire `h=n-4` branch
of the inherited large-order residual. It does not close every other
endpoint-four port set or prove global black-core connectivity.

## Exact statement and inheritance

Use `P=(0,1,...,n-1)=(a,b,x,d,U,e,y,c,z)`, where
`a=0,b=1,x=2,d=3,e=n-4,y=n-3,c=n-2,z=n-1` and
`U=(4,...,n-5)`. Assume P is closed, directed and blue and

```
N^-_T(a)={e,y,c,z},   N^+_T(z)={a,b,x,d},   n>=11.          (1)
```

The inherited [odd-prefix repair](../../endpoint_four_frontier_20261021/README.md)
handles the chord `4->x` from order fifteen; its ell=1 proof extends to
order eleven under (1), as shown below. The new complementary `x->4`
ledger came from the ordered-graph operation compiler: retain every unknown
interior edge inside U, and ask which changes of block boundaries have the
correct path-edge XOR. Its three blackness arguments retain the literal
ports and the actual positional involutions.

No clone contraction is used. Consequently clones, shared endpoint
neighbourhoods and all unrelated interior arc choices remain exact.

## 1. The new x->4 repair

If `x->4`, use

```
Q1=(a,b,x,U,e,y,c,z,d),
Q2=(b,x,U,a,e,y,c,z,d),
Q3=(b,x,d,U,a,e,y,c,z).                                  (2)
```

Q1 replaces the consecutive segment `x,d,4` by `x,4` and moves d after z;
the new joins are the assumed x->4 and the outgoing port z->d. At
`S=K_Q1(T)`, Q2 inserts a between `n-5` and e. Both new joins are
complements of original a->(n-5) and e->a. Q3 restores x->d and d->4
after two complements. All other joins are retained from Q2.

The changed unordered pairs cancel twice; equivalently direct cancellation
of the displayed block boundaries gives
`E(Q1) xor E(Q2) xor E(Q3)=E(P)`. Thus the target is exactly `B=K_P(T)`.

Q1 is black: its reversal would equate `T(c,a)=1` with `T(d,x)=0`.
Q2 is black: at S its reversal would equate `S(z,b)=0` with
`S(d,x)=1`. The first pair reverses the original outgoing port z->b;
the second reverses the original path edge x->d, which Q1 omitted.

For Q3, check blackness at B, the other end of its own K move. If it
were blue, `g=r_Q3 r_P` would be an automorphism of B. Its explicit action,
obtained by reading the two displayed reversals, is

```
g(0)=1, g(1)=2, g(2)=3, g(3)=4,
g(i)=i+2                 for 4<=i<=n-7,
g(n-6)=0, g(n-5)=n-4, g(n-4)=n-3,
g(n-3)=n-2, g(n-2)=n-1, g(n-1)=5.
```

For odd n this is one n-cycle. Hence B would be regular, but
`N^+_B(a)={b,e,y,c,z}` has size five, forcing n=11. For even n there
are two n/2-cycles, and the one containing a is

```
{0,1,2,3,4,6,8,...,n-6}.
```

Only b of the five outneighbours of a lies in this orbit. The induced
orbit tournament would therefore be regular of outdegree one, forcing
n/2=3, or n=6. This is impossible in the stated range. Thus all three
words in (2) are black for n>=12.

## 2. Order eleven reduces to one labelled tournament

At n=11 the candidate automorphism is the cycle

```
g=(0 1 2 3 4 6 7 8 9 10 5).
```

An invariant tournament under a transitive cycle is uniquely determined by
the complete outgoing row of one vertex: transport every ordered pair by
a power of g until its first vertex is a. Here the row is already fixed
as `N^+_B(a)={1,7,8,9,10}`. Its transport is consistent and gives exactly
one B, and hence exactly one source T. The source bit rows are

```
(126,220,472,976,1952,1102,1936,1825,1571,1063,15).        (3)
```

This source really has Q3 blue, so the exception is retained rather than
hidden in a score argument. It has the alternative three-black certificate

```
(0,1,2,3,4,7,5,6,8,9,10),
(0,1,2,3,4,7,5,9,10,6,8),
(0,1,2,3,4,5,9,10,6,7,8).                               (4)
```

Direct intermediate arc-set and bit-row replay proves (4). Its discovery
used the inherited labelled first-success search in
`endpoint_three_repairs_20261020h/probe_search.py`; that search stopped
after three tested first/last path pairs. The proof here uses only the
displayed certificate and the transitive-cycle uniqueness argument, not
the search's enumeration of Hamiltonian paths.

Every other source at n=11 has Q3 black, so (2), with (4) replacing it for
the single source (3), proves the x->4 branch at every n>=11.

## 3. The complementary 4->x branch

For `4->x`, put `V=(4,...,c)` and use the inherited ell=1 words

```
R1=(a,x,d,V,z,b),
R2=(d,a,x,V,z,b),
R3=(d,a,b,x,V,z).                                       (5)
```

Their legality uses a->x, a->d and z->b from (1), followed by
complementing 4->x. Their boundary pairs have the same XOR identity as P.
The first reversal is that of a cyclic cut. If it were blue along with P,
the shift-by-two automorphism would be transitive for odd n, forcing n=9.
For even n its even-label orbit has at most two incoming neighbours of a,
forcing n<=10. Thus R1 is black for n>=11.

After R1, R2 reversal pairs a with z, but a has outdegree five and z has
indegree four. Hence R2 is black. At B, the product of the R3 reversal
with r_P consists of two disjoint four-cycles on the first and last four
positions, fixing the other positions. This cannot be a tournament
automorphism: a suitable power is an involution exchanging a pair and
reversing its arc. Thus R3 is black. This proves (5) from order eleven.

Together the two chord directions exhaust the tail family (1).

## 4. Remaining frontier and verification

Combining this theorem with the inherited n>=15 partition and the
[arbitrary-port cyclic-prefix repair](THEOREM.md), a surviving endpoint-four
source must have closed middle, `I={h,y,c}`, `h` different from both x and e,
and `c->x,y->x`. It must reject every applicable old/new prefix guard,
including the new simultaneous guard `3->c` and `(n-h)->x`.
These are exact necessary conditions for being outside the proved bank;
they do not exhibit a disconnected black component.

`tail_check.py` replays both generic ledgers, independently reconstructs the
unique order-eleven exception from its cyclic symmetry, replays its repaired
certificate, and exhausts all reflected tail interiors at n=9,10,11 as
bounded controls. Its larger controls keep every h=e literal and test both
directions of x versus4. The declared theorem rests on the proofs above.

Run `python tail_check.py` and `python -O tail_check.py`; receipts must agree.
