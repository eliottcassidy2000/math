# Moving both reflected ports: a uniform open end-pair repair

October 20, 2026. **PROVED**, with the quantifiers below. This is an actual
three-black detour to a literal labelled partner, not a path relocation or
an isomorphism inferred from endpoint degrees.

## The open end-pair theorem

Let n>=7 and let P be closed and blue with endpoint scores (n-4,3).
Write

```
P=(a,b,x,W,y,c,z),       I={y,c},       J={b,x},
```

where W is nonempty. Suppose the middle is open, b->c. No internal
chord or module condition on W is required. The following three K words
are directed and black and end exactly at K_P(T):

```
Q1=(a,b,c,z,x,W,y),
Q2=(b,c,z,x,W,a,y),
Q3=(b,x,W,a,y,c,z).
```

Here K_Q complements all nonpath pairs of Q. A word is blue precisely
when its positional reversal is an anti-automorphism.

**Directedness.** Q1 uses the open chord b->c, the path arc c->z,
the exceptional arc z->x, and the remaining displayed path. At its
target, last(W)->a and a->y hold because a->last(W) and y->a were
complemented. These are the only new joins of Q2. In the third source,
b->x and y->c have each been complemented twice; c->z has been
retained twice. The other Q3 joins were retained by Q2. Thus all three
words are directed for every internal orientation of W.

**Exact target.** Added pairs bc,zx,last(W)a,ay occur twice. The remaining
boundary pairs ab,bx,x first(W),last(W)y,yc,cz occur once; the internal
W edges occur three times. The path-edge XOR is exactly E(P). Three
K moves therefore end at the literal labelled tournament K_P(T).

**First move.** Its reversal sends a to y and c to the P vertex at
position n-5. The source arc a->P[n-5] is true because that vertex is
neither exceptional port nor z. The image arc c->y is false because
y->c is a P edge. Hence Q1 is black. This positional formula uses n>=7;
at n=6 the image of c is z, so this witness must not be exported there.

**Second move.** Its reversal sends a to c and b to y. At the source
K_Q1(T), a->b was retained and is true, while y->c was complemented and
is false. These are the two arcs that an anti-automorphism would equate.

**Third move.** Its reversal sends b to z and c to x. After the first
two moves, b->c is still true because it was retained twice, while
x->z remains false because the actual exceptional arc z->x was
retained twice. Again the anti equation fails. Neither of these last
two witnesses uses a score estimate.

The three witnesses prove the theorem. In particular this repairs the
open end-pair branch even when the earlier sufficient guards
b->last(W) and x->last(W) fail.

## Inheritance and the decisive change

The inherited endpoint-two B ledger moved a immediately before z in
its second word. That join is wrong when c is the additional incoming
port. The new word instead moves a between W and y, using the actual
complemented port y->a. The final XOR is recalculated and the actual
reversal checked at each step. This avoids the documented error of
reusing a singleton's positional reversal after changing a block.

The smallest discovery controls were the two order-eight inherited-bank
residuals

```
(30,68,56,178,162,195,157,7),
(30,76,56,176,226,195,141,7).
```

Both admit the displayed template, with W=(3,4). Their endpoint-excess
triples are (1,0,-1) and (1,-1,0): each contains a balanced step, so
degrees alone would miss precisely the reason this repair works.

The independent reflected-port construction and exact small controls
are recorded in [probe_refined.py](probe_refined.py). The branch proof
is uniform in W and supplies the all-order quantifier. Finite agreement
is a control, not the justification for extrapolation.

## Scope and retained information

The source is a marked closed blue path, the map is the displayed triple
of K operations, and the preserved predicate is exact labelled target
equality together with blackness of every move. The sidecar is I={y,c}
and J={b,x}, plus the open-middle arc b->c. Forgetting those incidences
would destroy the legality argument. This theorem does not classify the
closed end-pair branch, other endpoint deficits, or global black components.
