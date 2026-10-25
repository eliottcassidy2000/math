# The local cyclic contradiction removes the protected-port size threshold

2026-10-25. **PROVED:** the protected-port theorem below, using the inherited
consecutive-tail repair for one precisely identified exception. **FINITE-EXACT:**
`boundary_check.py`. The earlier [Theorems B and C](README.md) remain valid;
their sufficient size bounds are superseded by this stronger scope. Their
proofs and receipts are retained unchanged. Full black-core connectivity is
still **OPEN**.

## Theorem

Use the literal marked-path setting of [README.md](README.md): P=(0,...,n-1)
is directed, closed and blue in T, a=0,b=1,x=2,z=n-1,c=n-2,
`N^-_T(a)={z} union I`, `k=|I|+1`, and `B=K_P(T)`. Suppose

    n>=9, h=min I>=5, and {n-4,n-2} subset I.

Then T has a three-black detour to the exact labelled target B. No restriction
on the other ports, their gaps, the middle tournament, shared reflected
endpoints, or clone identities is imposed. These assumptions imply k>=3.

In particular, the one-hole family

    I={n-k-1} union {n-k+1,...,n-2}

is covered for **every k>=5 and n>=k+6**. This is the exact h>=5 range for
this family. The k=4 one-hole family lacks the required port n-4 and still
uses the separate Theorem A; it is not silently included here.

## 1. The first cyclic move is blocked by its own legality guard

If 4->2, use the ell=1 words

    L1=(0,2,3,4,...,z,1),
    L2=(3,0,2,4,...,z,1),
    L3=(3,0,1,2,4,...,z).                         (1)

The earlier edge ledgers already prove their legality and exact target from
`a->2`, `a->3`, `z->1`, and `4->2`. The only unnecessarily weak step was the
degree bound for L1. Its positional reversal sends 2 to z and 4 to n-3.
Thus an anti-reversal would equate

    T(2,4)=0,
    T(n-3,z)=T(a,2)=1,                           (2)

where the second equality is the original blue reflection r_P. This is a
literal pair contradiction at every n in the theorem. In terms of composed
reversals, shift by two would send `a->2` to `2->4`, contradicting the guard.
There is no cyclic degree exception to repair in this branch.

The other two inherited blackness arguments do not need a size bound in k:
after L1, a has outdegree k+1 and z has indegree k, while r_L2 exchanges them.
For L3, the product with r_P has two disjoint four-cycles (n>=9 keeps them
disjoint), impossible for a tournament automorphism. This completes (1).

## 2. The only remaining exception is the already repaired consecutive tail

If 2->4, put U=(4,...,h-1), R=(h,...,z) and use

    Q1=(a,b,x,U,R,3),
    Q2=(b,x,U,a,R,3),
    Q3=(b,x,3,U,a,R).                            (3)

Legality and XOR are exactly those of Theorem B. The retained ports n-4 and c
give `z->3` and `c->a`; h=min I gives `a->(h-1)` and `h->a`. The first two
black witnesses are unchanged:

* In T, `T(c,a)=1` versus `T(3,2)=0` under r_Q1.
* In S=K_Q1(T), `S(z,b)=0` versus `S(3,2)=1` under r_Q2.

For Q3 it suffices to work at B. Its word is `(1,...,h-1,0,h,...,z)`.
The [rotation lemma](README.md) computes `g=r_Q3 r_P` exactly. Disjoint
supports force a regular h-vertex orbit in which a has just one outneighbour,
contradicting h>=5. Even overlapping supports similarly force an orbit of
size n/2 to have outdegree one, which would give n=6. Odd overlapping
supports are transitive, so Q3 can be blue only if

    n=2k+3.                                      (4)

At (4), h<=k+3 follows simply because I has k-1 members in [h,n-2].
If h<=k+1 the supports are disjoint, already excluded. This leaves two cases:

* **h=k+2.** Here n=2h-1. The reversal r_Q3 fixes a and exchanges b and z.
  It would equate `B(a,b)=1` with `B(z,a)=0`. The first is a retained path
  edge; the second is the complement of the closing edge z->a. Thus Q3 is
  black. This is also a direct proof for the one-hole n=2k+3 boundary.
* **h=k+3.** The interval [h,n-2] has exactly k-1 vertices. Therefore
  `I={k+3,...,2k+1}` is the entire consecutive tail. Its inclusion of n-4
  forces k>=4. This is precisely the unique cyclic family covered by
  [GENERAL_TAIL.md, sections 2--3](../../trijection_hourglass_ports_20261025/endpoint/GENERAL_TAIL.md).

For completeness, the last citation supplies more than a numerical
connectivity assertion. If Q3 is blue, B is the unique cyclic tournament in
the coordinates

    C=(0,1,...,k,k+2,k+3,...,2k+2,k+1),
    B(u,v)=1 iff pos(v)-pos(u) mod n is in {1,k+2,...,2k+1}.

With A=(0,...,k), u=k+1,v=k+2,w=k+3,s=k+4 and Z=(k+5,...,2k+2), its
three literal repair words are

    (A,w,u,v,s,Z),
    (A,w,u,Z,v,s),
    (A,u,Z,v,w,s).                               (5)

The cited proof establishes their joins, their E(P) XOR and three explicit
anti-reversal failures for every k>=4. Thus (3) works unless its final word
is blue, and the exact remaining source has (5). This finishes the theorem.

## 3. Failure boundary and controls

The h>=5 condition is still load-bearing for these ledgers. At h=4, vertex
3 is moved to the end by Q1, so the putative Q3 requires the false join
3->a. The source n=15, I={4,11,13}, with 2->4 is retained as a feasible
hostile; the failed join is forced, not merely unproved. This invalidates
that word substitution, not the existence of another black detour.

`boundary_check.py` verifies all sources in each affine cell, rather than
sampling their middle orientations. It uses the shared signed equation
solver in `check.py` but supplies the new complete guard split and the
conditional exceptional repair. It independently replays literal arcs for
model sources, so the universal certificates and direct executions have
different representations. No source is discarded because a test fails.

Its declared universe is every protected-port set for n=9..17, both chord
guards; every one-hole set for k=5..24 and n=k+6..2k+7, both guards; and
the unique exceptional consecutive-tail source for every k=4..40. It also
checks the explicit local witness (2), the h=k+2 endpoint witness, and the
h=4 failed-join hostile. Run `python boundary_check.py` and
`python -O boundary_check.py` from this directory.

The research move is precise: test a cyclic symmetry against the guard
that made the path legal before relying on orbit-degree capacity. Here
the guard already rules the symmetry out. This is a scope strengthening,
not a correction of a false previous theorem or a proof of unrestricted
black-core connectivity.
