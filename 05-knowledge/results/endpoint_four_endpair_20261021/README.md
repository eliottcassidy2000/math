# One more endpoint-four subbranch and genuine final-word hostiles

21 October 2026. **PROVED** for the displayed branch; **FINITE-EXACT**
for four blue-template hostiles. This is an addendum to
endpoint_four_repairs_20261021, not a full endpoint-four theorem.

Let P=(a,b,x,W,y,c,z) be closed and blue, with n>=15,
d_T(a)=n-5, d_T(z)=4, and I={h,y,c}. Assume

    c->b, c->x, x->y,    h not in {x,n-4}.

Then there are three black K moves to the exact labelled K_P(T):

    Q1=(a,x,W,y,c,z,b),
    Q2=(z,W,y,x,b,a,c),
    Q3=(b,z,W,y,x,a,c).                                (A)

This is the inherited endpoint-three E2 x->y ledger with its additional
port retained. The new work is checking both extra blackness exceptions,
not claiming the inherited word list as new.

## Directedness and exact target

The first word needs a->x and z->b, provided by h!=x and c in I.
For Q2, z->first(W) at K_Q1(T) is the complement of first(W)->z;
this is supplied exactly by h!=n-4. Its new joins y->x, x->b, b->a,
a->c are the complements of x->y, b->x, a->b, c->a. All other joins
are retained. The two newly required joins b->z and x->a for Q3 are
complements of retained Q1 arcs; all other Q3 joins survive Q2.
The added pairs ax,zb,z-first(W),yx,ac occur twice; every surviving
pair is exactly a P pair. Thus the three path-edge sets XOR to E(P).
The endpoint is literal, not an isomorphism class.

## Three actual blackness arguments

Q1 reversal compares first(W)->z, which is true, with x->c, which is
false. Thus Q1 is black. Write S=K_Q1(T) and X_v=T(a,v).

If h!=5, Q2 reversal compares S(a,5), false, with S(x,3), the retained
true P arc. If h=5, assume Q2-blue and use these four ordered-pair
bit equations, alternating Q2-blue at S and P-blue at T:

    (a,4) -> (b,3) -> (n-4,c) -> (7,z) -> (a,n-8).

The Q2 equations have respective right-hand sides zero and one after
accounting for which pairs Q1 retained; the P equations have zero.
Hence X_4 xor X_(n-8)=1. At n>=15 with I={5,y,c}, both vertices 4 and
n-8 are outside I union {z}; both bits are one. This contradiction
proves Q2 black in its only additional exceptional port case.

Check Q3 at the opposite endpoint B=K_P(T), where blue/black is the
same as at its own source. If h!=3, its reversal compares B(a,3),
false, with B(x,z), true. If h=3, it instead compares B(n-4,z), true,
with B(a,5), false: the reversal sends z to a and n-4 to 5. The first
bit is true because n-4 is the reflected extra outgoing port of z in
T; the second is false because a->5 was a nonpath arc in T. Thus Q3
is black in all cases. This proves (A).

The additional repaired region leaves the following frontier after the
base packet: closed I={h,y,c}, with h in {x,n-4}, or c->x and y->x.
No claim about black-core connectivity follows.

## Why the opposite chord case cannot simply reuse its old final word

When c->x and y->x, the inherited E2 list is

    Q1=(a,x,W,y,c,z,b),
    Q2=(W,a,x,y,c,z,b),
    Q3=(W,a,b,x,y,c,z).                                (B)

The receipt gives four exact marked sources: n=15 with h=3,6 and
n=18 with h=3,9. In each source, all three words of (B) are directed,
the XOR target is exact, but Q3 is actually blue in its source and
at B. Their full adjacency rows are saved in witnesses.json and
independently replayed by verify.py. These examples concern this
particular word, not existence of another three-black detour or a
metagraph bridge. No assertion that these are the only hostile orders
is made.

## Reproduction and provenance

    python3 verify.py > receipt-normal.json
    python3 -O verify.py > receipt-optimized.json
    cmp receipt-normal.json receipt-optimized.json

The checker uses directed arc sets, no inherited imports, verifies the
complete blue predicate at every step, checks the literal endpoint,
and checks the four-link affine identity across n15..200. It samples
all allowed h at each n15..40 with deterministic reflected interiors,
plus the stored hostile sources. Finite controls support, not replace,
the displayed uniform proof. Sources are pinned in SOURCE-MANIFEST.json;
only MANIFEST.sha256-listed files plus the manifest are publication
payload. No repository or remote was changed by this lane.
