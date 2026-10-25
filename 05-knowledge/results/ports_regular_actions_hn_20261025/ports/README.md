# Nonconsecutive ports: the first missing tail vertex and a protected-port family

2026-10-25. **PROVED** structural detours, with **FINITE-EXACT** controls.
All vertices and intermediate arcs remain literal. No clone contraction,
isomorphism quotient, or middle-orientation filter is used. Full black-core
connectivity and all nonconsecutive port sets remain **OPEN**.

## Inheritance and exact setting

The closest mechanism is the
[consecutive-tail repair](../../trijection_hourglass_ports_20261025/endpoint/GENERAL_TAIL.md)
and its [affine compiler](../../trijection_hourglass_ports_20261025/trijection/README.md).
The canonical hostile is an illegal join hidden by a correct-looking block
picture. The corrected near miss is moving the insertion predecessor into
the very block being removed; this really occurs at n11 below. The less-used
sidecar is the full outneighbour set inside an orbit of a composed reversal.

Let P=(0,1,...,n-1) be a directed, closed, blue word in T. Thus
r_P(i)=n-1-i is an anti-automorphism and z=n-1 points to a=0. Write
b=1,x=2,y=n-3,c=n-2. A K move keeps the arcs of its directed word and
reverses every other arc. A word is black iff its actual positional reversal
is not anti in its current source. Put B=K_P(T); P is also blue in B.

If N^-_T(a)={z} union I and k=|I|+1, then

    N^+_T(z)={a} union r_P(I),
    N^+_B(a)={b,z} union I.

Three words have literal target B exactly when their unordered path-edge
sets XOR to E(P). Blackness can be checked at either endpoint of a word's
own K edge. Every construction below checks these distinct obligations.

## 1. The entire first nonconsecutive endpoint-four tail

**Theorem A.** For every n>=12, the port set

    I={n-5,n-3,n-2}

has a three-black detour to B. The middle may be open or closed.

Put h=n-5 and R=(h,...,z). We use a complete four-case partition:

1. If 4->2, use the inherited odd-prefix ledger with ell=1.
2. If 2->4 and 6->2, use that ledger with ell=3.
3. Otherwise 2->4 and 2->6. Split by 3 versus5 as below.

The inherited bound is n>=max(11,2ell+3). Its excluded port values for
ell1 are h=2,3, and for ell3 are h=2,4,5. Since h=n-5>=7, both applications
are valid, with no further orientation assumptions.

### Case 3a: 3->5

Set V=(5,...,h-1) and use

    Q1=(a,b,x,3,V,R,4),
    Q2=(b,x,3,V,a,R,4),
    Q3=(b,x,3,4,V,a,R).                            (A)

Q1 removes4 and places it after z. Its new joins are 3->5 and z->4;
the latter is an endpoint port. Q2 inserts a between h-1 and h, using
the complements of a->(h-1) and h->a. Q3 restores3->4 and4->5 after
two complements. Boundary-pair cancellation gives the E(P) XOR identity.

Its three blackness witnesses are:

* Q1 at T: c->a has value1 but its reflected reverse4->2 has value0.
* Q2 at S=K_Q1(T): z->b has value0 but its reflected reverse4->2 has
  value1, because the original2->4 was not a Q1 edge and was complemented.
* Q3 at B: a->c has value1 but its reflected reverse2->6 has value0.
  Here r_Q3(a)=6 and r_Q3(c)=2; 2->6 was complemented by K_P.

### Case 3b: 5->3

Set U=(6,...,h-1), which is nonempty from n12, and use

    Q1=(a,b,x,U,R,4,5,3),
    Q2=(b,x,U,a,R,4,5,3),
    Q3=(b,x,3,4,5,U,a,R).                          (B)

Q1 replaces the block3,4,5 by the join2->6 and places4,5,3 after z.
The new joins z->4,4->5,5->3 are known. Q2 again inserts a between
h-1 and h. Q3 restores2->3,3->4 and5->6 after two complements;
4->5 was retained throughout. All other joins are retained. Added
nonpath boundary pairs occur twice, proving the E(P) XOR identity.

Its blackness witnesses are:

* Q1 at T: a->2 has value1 and its reflected reverse4->3 has value0.
* Q2 at S: 2->4 has value0 and its reflected reverse6->5 has value1.
  Neither pair is a Q1 edge, and their original directions are2->4,5->6.
* Q3 at B: exactly the a->c versus2->6 contradiction used in (A).

These cases exhaust all arc choices and prove Theorem A. Both branches
have the same final word (1,...,h-1,0,h,...,z); its actual anti-reflection
test has become a single boundary-chord test. This is why the extra
sixth vertex in the repair matters.

### Honest small-order boundary

At n11, h=6, the U block in (B) is empty and h-1=5 has been removed.
The proposed Q2 insertion is then after2 rather than after5. Q3 contains
5->a, which remains false after two complements, and the path-edge XOR
also fails. This is a hostile to this ledger, not a black-connectivity
counterexample. The checker supplies an actual source satisfying all
case3b guards and exhibits this failed join. No n11 claim is made here.

## 2. A general rotation lemma for insertion at the first port

For 4<=h<=n-4 put Q=(1,...,h-1,0,h,...,n-1). The product
g=r_Q r_P is the composition of positive one-step rotations on the first
and last h positions, first applying the last rotation. Let q=n-h.

If 2h<=n, these supports are disjoint (or touch without sharing a vertex),
and a's orbit is exactly0,...,h-1. If 2h>n, write j=2h-n:

* For odd j, g is one n-cycle. At j=1 this is the ordinary cyclic shift.
  At j>=3 one can read its single cycle as
  0,1,...,q,q+2,...,h-1,h,...,z,q+1,q+3,...,h-2.
* For even j, g has two n/2-cycles. The orbit containing a is
  {0,1,...,q,q+2,...,h-2}.

These formulas follow directly by composing the two rotations; they also
prove the asserted coverage and cycle lengths. Since j and n have the
same parity, the overlap alternatives can be stated in terms of n.

Suppose h=min I. Every member of N^+_B(a) other than b lies at least h.
Thus in the disjoint case a has exactly one outneighbour in its h-orbit;
if g were an automorphism its induced regular tournament would force h=3.
In the even overlap case the same argument forces n/2=3. Both are excluded.
In the odd overlap case transitivity forces n=2k+3, because B(a) has
outdegree k+1. Consequently Q is black whenever h>=4,n>=8 and
n is different from2k+3. This conclusion uses h=min I, not merely |I|.

## 3. Arbitrary nonconsecutive sets with two retained end ports

**Theorem B.** Let k>=4,n>=4k-1, |I|=k-1, h=min I>=5 and

    {n-4,n-2} subset I.

Then the marked closed blue P has a three-black detour. All other members
of I are arbitrary, including shared/reflected positions.

If 2->4, put U=(4,...,h-1), R=(h,...,z), and use

    (a,b,x,U,R,3),
    (b,x,U,a,R,3),
    (b,x,3,U,a,R).                                  (C)

The same exact surgery as the inherited consecutive-tail theorem is legal:
z->3 follows from n-4 in I; a->(h-1) and h->a follow from h=min I.
The first two black witnesses are c->a versus3->2 and, after Q1,
z->b versus3->2. They hold because c in I and2->4 is the removed join's
replacement. The final word is covered by section2, since n>=4k-1>2k+3
for k>=4. Its edge XOR is E(P).

If 4->2, use the inherited ell1 ledger. It is legal since a->2,a->3 and
z->b hold. Its first reversal agrees with a cyclic-cut reversal. If blue,
the shift-by-two permutation would be an automorphism of T. For odd n it
is transitive and forces n=2k+1. For even n the even-label orbit contains
at most k-1 incoming neighbours of a, since z is odd. Regularity forces
n<=4(k-1)+2=4k-2. Both contradict n>=4k-1. The second word pairs
outdegree k+1 with indegree k, and the final reversal product has two
disjoint four-cycles. Thus both are black, exactly as in the inherited
proof. This proves Theorem B without making the ports consecutive.

The h>=5 condition is load-bearing for this surgery. At h4 its insertion
predecessor h-1=3 has been moved to the end; the proposed last word then
needs the false join3->a. The checker retains a source witnessing this.
This failed implication is the same block-overlap error as the n11 case
of Theorem A, not a mathematical obstruction to another detour.

## 4. A lower threshold for one missing tail position in higher layers

**Theorem C.** For k>=5 and n>=2k+5, the one-hole port set

    I={n-k-1} union {n-k+1,...,n-2}

has a three-black detour. This has k-1 interior ports; its missing vertex
inside the terminal interval is n-k. The k4 version lacks the required
n-4 port and is handled separately by Theorem A.

Use the two ledgers of Theorem B with h=n-k-1. For2->4, all legality and
blackness arguments still apply, since h>=k+4 and n>2k+3. For4->2,
only the cyclic-cut capacity needs improvement. At even n, its even-label
orbit contains k/2-1 incoming neighbours when k is even, and (k+1)/2
when k is odd. Hence blue would force respectively n=2k-2 or n=2k+4.
At odd n it would force n=2k+1. All are below n>=2k+5. The other two
blackness arguments are unchanged. This proves Theorem C.

## 5. Frontier and reproduction

At endpoint four and n>=15, the former h=n-5 residual is now completely
removed, regardless of its middle state. Together with the inherited
tail and second-port theorems, the large-order bank remainder has

    3<=h<=n-6, I={h,y,c}, c->x, y->x,

closed middle and failure of every applicable old/new prefix guard.
These are necessary conditions for lying outside the proved bank, not
an existence or disconnection claim. Theorem B covers other higher-layer
port sets but does not silently enlarge this endpoint-four statement.

Run `python check.py` and `python -O check.py`. Its complete affine
universes, literal checks, and hostiles are recorded in the receipts.
The structural proofs give the unbounded conclusions; exact controls do
not infer them from a finite interval. The exploratory `probe.py` imports
the preceding compiler only for discovery and is not the independent
verification path used by `check.py`.
