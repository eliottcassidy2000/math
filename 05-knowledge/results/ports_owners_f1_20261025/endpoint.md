# Exact endpoint-four repairs at the second port and along a growing prefix

**PROVED:** the closed endpoint-four family with incoming middle ports
`I={x,y,c}`, for every n>=11. A further guarded prefix theorem handles
arbitrary h with `I={h,y,c}`, `c->x`, `c->(h+1)`, and n>=2h+8.
**FINITE-EXACT:** independent arc-set and bit-row controls below.
Full endpoint-four and all-order black-core connectivity remain **OPEN**.

## Inheritance and current target

The requested endpoint-three layer `(n-4,3)` is already closed for every
n>=7: the [structural theorem](../endpoint_three_repairs_20261020h/README.md)
handles n>=11, and the [audited order-ten packet](../endpoint_three_order10_20261021/README.md)
joins the earlier finite n7..9 head. The unfiltered n>=6 extension is false:
source rows `(6,4,8,51,39,7)` have pure-blue partner `C3 -> C3`.
MISTAKE-646 corrects that quantifier and the duplicate-discovery claim.

The closest current mechanism is the
[endpoint-four odd-prefix frontier](../endpoint_four_frontier_20261021/README.md).
Its exact remainder includes `h=x,c->x`, `h=n-4`, and `c->x,y->x`;
its sufficient cases can overlap. Our new complementary ledgers close the
first entire family. The canonical hostile is the pure-blue order-six
target, not a failed fixed template. The corrected near miss is the
reflection-fixed-pair error: c->b and c->x do not determine x versus y.
Neither new ledger assumes the direction of x versus y. The least-used
sidecar is the composed positional reflection on the unchanged interior:
it turns a shifted prefix into a test of two literal endpoint incidences.

The board is: actual exceptional ports; positional anti-involutions;
path-edge XOR; growing prefix blocks; exact intermediate arc sets;
and the distinction between a failed template and a disconnected component.
All claims concern labelled tournaments with one marked path. No clone
contraction or isomorphism quotient is performed.

## 1. Setting

Let `P=(0,1,...,n-1)=(a,b,x,3,...,y,c,z)` be a closed blue Hamiltonian
word in T, with `d_T(a)=n-5`, `d_T(z)=4`. Write r_P(v)=n-1-v.
The exceptional set I is defined by

```
N^-_T(a)={z} union I,     N^+_T(z)={a} union r_P(I).
```

We treat a closed middle, c->b. A K move complements all unordered pairs
outside the directed path word. A word is blue precisely when its actual
positional reversal r is anti: T(u,v)=T(r(v),r(u)) for every ordered pair.
Put B=K_P(T). Blue/black status of a word agrees across its own K edge.
Thus the final word may be tested at B. For any three directed words,
`E(Q1) xor E(Q2) xor E(Q3)=E(P)` proves that their literal target is B.

## 2. The second exceptional port: I={x,y,c}

The incoming frontier already handles x->c for n>=7 (section 3 there).
Suppose c->x. Put d=3 and split by the actual arc d versus c.

### Case A: d->c

Write D=(3,...,n-4). Use

```
Q1=(z,a,b,x,D,y,c),
Q2=(z,a,x,D,b,y,c),
Q3=(a,x,D,b,y,c,z).                                  (A)
```

Q1 is a cyclic cut of P, hence directed. At S=K_Q1(T), a->x is the
complement of x->a. If e=n-4, blue symmetry gives b->e iff d->c, so
e->b holds at S. Also y->b iff c->x, hence b->y holds at S. These are
exactly the new joins of Q2. Q3 restores c->z after two complements;
all its other joins are retained. Pair by pair the added za,ax,eb,by
occur twice and the remaining odd-multiplicity pairs are exactly E(P).

Each word has a literal anti-reversal contradiction:

* In T, Q1 reversal maps c->a to y->z. The former is true and the latter
  false, since x in I forces z->y.
* At S, Q2 reversal maps z->x to b->c. The former is false: Q1 does not
  retain the original true z->x. The latter is true: c->b is a nonpath
  arc of Q1 and is complemented.
* At B, Q3 reversal maps a->y to d->z. The former is true after the
  original y->a is complemented; the latter is false after d->z is
  complemented. Here d is not one of z's original ports b,x,y.

Thus (A) is a three-black detour. These arguments actually hold from
n>=7 in this branch; n>=11 is the uniform theorem's range.

### Case B: c->d

Write U=(4,...,n-4). Use

```
Q1=(b,x,a,d,U,y,c,z),
Q2=(b,x,a,d,c,z,U,y),
Q3=(a,b,x,d,c,z,U,y).                                (B)
```

Q1 uses x->a and a->d. At S its new Q2 joins d->c and z->4 are
complements of c->d and 4->z. Q3 restores a->b, x->d, and the
other omitted original path joins after two complements. The added
ax,ad,dc,z4 cancel in pairs; all remaining odd pairs are E(P).

For Q1, suppose both r_P and r_Q1 were anti. Their product would be an
automorphism rotating the first three positions and the last three,
and fixing e=n-4. In particular it sends a to b and fixes e. But
T(a,e)=1 while T(b,e)=T(d,c)=0. This is impossible.

For Q2 and n>=12, its reversal sends a to n-5 and n-8 to z.
At S, the arc a->(n-8) is false and z->(n-5) is true: both underlying
pairs were complemented by Q1, and the relevant vertices are ordinary
interior vertices, distinct from d and the endpoint ports. These two
bits would have to agree under Q2 reversal. At n=11, use instead
S(b,z)=1 versus S(z,y)=0; Q2 reversal interchanges b,y and fixes z.

At B, Q3 reversal maps a->e to b->y. The first is false and the second
true, because T(a,e)=1 and T(b,y)=T(x,c)=0. Thus Q3 is black.
This proves (B) for every n>=11, regardless of the x,y chord.

Combining the two cases with the inherited x->c ledger proves the full
closed `I={x,y,c}` endpoint-four family at every n>=11.

## 3. A growing-prefix repair

The second case has a useful variable-length extension. Let

```
I={h,y,c},    h>=2,    n>=2h+8,    c->b, c->x, c->(h+1).
L=(2,...,h), d=h+1, U=(h+2,...,n-4).
```

Then use

```
Q1=(b,L,a,d,U,y,c,z),
Q2=(b,L,a,d,c,z,U,y),
Q3=(a,b,L,d,c,z,U,y).                                (C)
```

The last vertex of L is the actual incoming port h, so Q1 is directed.
The first vertex of U is not one of z's outgoing ports b,x,n-1-h;
the bound n>=2h+8 guarantees this. The direction and XOR checks of (B)
therefore apply verbatim to the four boundary pairs ah,ad,dc,z(h+2),
with the full interior L retained. In particular no arbitrary vertex is
silently substituted for h.

The first word's position permutation rotates the prefix 0,...,h by
one place and fixes all other positions. Its reflection product rotates
that prefix and its disjoint reflected suffix. It fixes v=n-h-2 while
sending a to b. However T(a,v)=1 and T(b,v)=T(d,c)=0. Hence Q1 is black.

Q2 reversal sends a to u=n-h-3 and v=n-h-6 to z. The order bound places
u,v in U; neither is an exceptional port or a Q1 path neighbour of a
or z. Thus S(a,v)=0 whereas S(z,u)=1, proving Q2 black.

The last word has the same endpoint witness as (B): at B it sends
a->(n-4) to b->y, with values zero and one. This proves (C).

The proof uses the *actual* port h to choose its prefix. Its boundary
n>=2h+8 is a sufficient positional guard, not a proved optimal threshold.
The result advances some `c->x,y->x` cases beyond the second-port family;
it does not close h=n-4 or force the new chord c->(h+1).

## 4. Exactness, hostile probes, and remaining boundary

The checker constructs every reflected middle completion at n=7..11 for
the stated second-port, closed-middle, c->x branch, and independently
replays the two ledgers using arc sets and bit rows. At n=11 it checks
all16,384 sources (8,192 per chord direction). The n8..10 head has16,
128 and1,024 sources, all repaired. At n7, two Case A sources pass;
the two Case B words have an illegal join because U is empty. This
is recorded without extending the uniform n>=11 theorem.
The larger controls use1,120 deterministic second-port sources at
n11..80 and1,156 growing-prefix sources, including34 at the exact
boundary n=2h+8 for (C).
The source, actual intermediate directedness, actual anti predicate,
literal labelled target, and every proof witness are checked separately.
Normal and optimized runs must agree; the receipt gives exact counts.
No SCC, mixedness or isomorphism filter selects the sources.

Two hostile controls retain the guards' meaning: exchanging the two
complementary d,c ledgers produces a failed direction in each branch;
the order-six pure-blue target remains a genuine scope obstruction.
The proof keeps arbitrary clone configurations because every original
vertex and arc remains present. No weighted-core merger is assumed to
commute with K; that would be a different claim.

After combining this note with the incoming endpoint-four frontier,
the n>=15 residual can no longer have h=x. The remaining sources lie in
`h=n-4` or `c->x,y->x` with h!=x, and must fail both the incoming
odd-prefix guards and the new growing-prefix guard. This is a reduction
by sufficient detours, not a witness of disconnected black components.

Reproduction from the repository root:

```
python 05-knowledge/results/ports_owners_f1_20261025/endpoint_check.py
python -O 05-knowledge/results/ports_owners_f1_20261025/endpoint_check.py
```
