# Odd-prefix cuts repair the four endpoint-four final-word hostiles

21 October 2026. **PROVED guarded all-order repairs; FINITE-EXACT hostile
replays.** This bounded continuation leaves the frozen endpoint-four base
and end-pair addendum unchanged. It does not prove full endpoint-four or
black-core connectivity.

## Result first

All four previously stored blue-final-word hostiles are repaired on the
same labelled sources by one uniform cut at vertex 6. There is no new
source selection, isomorphism transport, or change of the target.

The structural mechanism is a family of odd-prefix cuts: a local chord
supplies legality, while the last word's reflection product has even
order on two disjoint endpoint blocks. A second short argument closes
the entire remaining h=x, x->c branch. The genuinely open remainder is
stated below rather than counted as a failed connectivity theorem.

## Coordinates and inheritance

Let P=(a,b,x,3,...,y,c,z), numbered 0,...,n-1, be a closed blue
Hamiltonian word in T, where n>=15, d_T(a)=n-5 and d_T(z)=4. Put

    a=0, b=1, x=2, y=n-3, c=n-2, z=n-1,
    N^-_T(a)={z,y,c,h},
    N^+_T(z)={a,b,x,k},  k=n-1-h.

Thus I={y,c,h}; h is a middle vertex other than b,y,c. Blue means that
the actual word reversal is anti. K_Q complements every nonpath pair.
Three directed words whose unordered path-edge XOR is E(P) have literal
labelled target K_P(T).

The inherited source is endpoint_three_repairs_20261020h/lift.md,
section 5, which proves odd-prefix blackness in its endpoint-three
port setting. This note retains the additional port y, checks the
changed exact neighbourhood, and applies the mechanism to the new
hostile sources. The mechanism itself is inherited, not claimed new.
The base endpoint_four_repairs_20261021 section 2 supplies the
four-port relocated-b obstruction used in the last theorem below.
All dependencies and frozen prior packets are hash-pinned.

## 1. The guarded odd-prefix theorem

Let ell be a positive odd integer, put t=ell+3, and assume

    n>=max(15,2(ell+3)),   t->x,
    h not in {2,ell+2}.

When ell>=3 require additionally h not in {4,ell-1}. At ell=1 there
are no additional exclusions. Write

    D=(3,...,ell+2),    U=(ell+3,...,c).

Then the following three K words are directed, black, and end exactly
at K_P(T):

    Q1=(a,x,D,U,z,b),
    Q2=(D,a,x,U,z,b),
    Q3=(D,a,b,x,U,z).                                      (1)

There is no condition on c versus b, c versus x, or x versus y.
In particular this theorem applies inside the exact residual frontier,
not only inside a previously solved chord case.

### Directedness and path-edge XOR

Q1 needs a->x, supplied by h!=2, and z->b, supplied by c in I.
Its other joins are P arcs. In S=K_Q1(T), the new join last(D)->a
is the complement of a->last(D), supplied by h!=ell+2; and x->t is
the complement of t->x. The a->x join is retained. These are exactly
the new joins of Q2. For Q3, the original arcs a->b and b->x are
complemented twice and restored; all other joins are retained by Q2.

The internal D and U edges appear three times. The added pairs ax,
zb, last(D)a, and xt appear twice; the surviving boundary pairs are
ab,bx,x-first(D),last(D)t,c z. Thus the three path-edge sets XOR to
E(P), proving the exact endpoint.

### First word

Q1 has the same positional reversal as the one-step cyclic cut
(b,x,3,...,c,z,a): only its first and last positions are exchanged.
The inherited sharp cyclic-cut theorem, endpoint_three root.md R1,
proves this reversal black at n>4d_T(z)-2=14. Hence Q1 is black.
The sharp order-14 limitation belongs to that test, not a bridge claim.

### Second word

At S the exact neighbourhood is

    N^-_S(z)={a,c,x,k}.

If ell=1, Q2 reversal pairs a,z. But d_S(a)=5 and d^-_S(z)=4;
these cannot be paired by an anti-automorphism. Thus Q2 is black.

If ell>=3, Q2 reversal sends 4 to z and a to n-ell. The first
vertex 4 is inside D, and S(4,a)=1 because h!=4 and Q1 did not
retain a4. The vertex n-ell belongs to U and is distinct from a,c,x.
It equals k only if h=ell-1, which was excluded. Consequently
S(n-ell,z)=0. The actual reversal therefore equates a true arc with
a false one, proving Q2 black. The size bound on n ensures all listed
positions lie in the stated blocks, including at equality n=2(ell+3).

### Final word

Work at B=K_P(T), where P is still blue and Q3 has the same blackness
as at its own source. The position permutation q for Q3 rotates the
first L=ell+3 labels left by three, fixing the rest. Its reflection
product is

    r_Q3 r_P = q r_P q^(-1) r_P.

The supports of q and r_P q^(-1) r_P are the first and last L
positions, disjoint because n>=2L. On either block its nontrivial
cycle lengths are L/gcd(L,3). Since ell is odd, L is even and this
quotient is even. The product thus has even order. An even-order
permutation cannot be a tournament automorphism: its involution power
swaps a pair and reverses that pair's arc. Both reversals cannot be
anti, so Q3 is black. This proves (1) uniformly in the interior.

## 2. One cut repairs all four exact hostile sources

Take ell=3, so t=6 and D=(3,4,5). The theorem reduces to

    n>=15, 6->x, h not in {2,4,5}.

The same three words therefore work for every such source:

    (a,x,3,4,5,6,...,c,z,b),
    (3,4,5,a,x,6,...,c,z,b),
    (3,4,5,a,b,x,6,...,c,z).                              (2)

The frozen hostile rows have (n,h)=(15,3),(15,6),(18,3),(18,9).
Every one satisfies 6->x and the port exclusions. Their old final
word was genuinely blue, not merely score-balanced. The checker
verifies that fact again and then replays (2) in the same labelled
source: three directed black words, identical literal K_P(T) target.
This gives an actual repair of the obstruction, not just a different
coordinate description.

The case ell=1 is also useful: t=4 requires h not in {2,3} and
4->x. Larger odd prefixes cover further chord patterns. The theorem
does not assert that a qualifying early chord always exists.

## 3. The entire h=x, x->c branch

Suppose the middle is closed, I={x,y,c}, and x->c. Put
D=(3,...,n-4). Then use

    Q1=(a,b,x,D,z,y,c),
    Q2=(a,x,D,z,y,b,c),
    Q3=(a,x,D,y,b,c,z).                                  (3)

The predecessor n-4 beats z because its reflection 3 is not in I;
z->y holds because x is in I. These supply Q1. At its target,
a->x, y->b, and b->c are complements respectively of x->a,
b->y=x->c, and c->b. These supply the new joins of Q2. The omitted
pairs (n-4),y and c,z are complemented twice and return for Q3;
all other joins are retained. The path-edge XOR is E(P).

Q1 reversal would equate false a->z with true x->c. Q2 reversal
would equate retained true a->b with complemented false x->c.
The last word is the relocated-b word F from the base packet. If it
were blue at B, r_F r_P would fix a,z and act by one n-2 cycle or
two (n-2)/2 cycles on the middle. The middle outneighbour set of a
in B is {b,x,y,c}, of size four. It cannot be a union of the odd
middle orbits at n>=7; even cycles cannot occur in an automorphism.
Thus Q3 is black. This proof shows that the earlier c-not-in-I
bookkeeping restriction was unnecessary for this specific ledger.
It does not authorize dropping its x->c guard.

## 4. Exact remaining scope after this pass

Combine the frozen base and end-pair addendum with this note. Any source
not yet covered must be a closed middle with I={h,y,c} and satisfy one
of the following:

- h=x and c->x;
- h=n-4;
- c->x and y->x.

In addition, it must fail every applicable odd-prefix guard in section 1
(otherwise that theorem already repairs it). The rows can overlap.
This is an exact reduction by the proved sufficient cases, not a claim
that any surviving source is a bridge or lacks longer black detours.
No census of these remaining tournaments was attempted.

This pass stops at a new all-order repair mechanism, a complete repair
of the four known hostile sources, and a sharper stated frontier. The
unresolved task is finding another guaranteed chord or a new surgery
when all short odd-prefix cuts are obstructed, particularly h=x,c->x.

## 5. Verification and reproduction

    python3 verify.py > receipt-normal.json
    python3 -O verify.py > receipt-optimized.json
    cmp receipt-normal.json receipt-optimized.json

The checker constructs sources directly, retains arbitrary reflected
interiors, tests every word in its actual intermediate tournament,
checks the full anti predicate and exact final arc set, and replays all
four stored hostile rows. Permutation controls cover every admissible
odd prefix at n15..200; these controls are separate from the uniform
permutation proof. Normal and optimized receipts must agree.
The previous packets remain byte-for-byte unchanged. Publish only files
listed by MANIFEST.sha256 plus that manifest. No repository or remote
mutation was performed by this lane.
