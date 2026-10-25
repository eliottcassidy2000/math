# The h=3 forward branch: move the cut to the later port

2026-10-25. **PROVED:** every marked closed blue endpoint-four source with
I={3,n-3,n-2}, 1->3, and n>=12 has a literal three-black detour. The proof
retains all middle arcs and labels. **FINITE-EXACT:** the affine, literal,
degree and small-hostile controls in `check.py`. Global black-core
connectivity is not asserted.

## Inheritance and exact scope

The nearest mechanism is the
[unprotected-port complementary split](../../unprotected_ports_cancellation_20261025/ports/README.md).
It removed the n-4 port but left h=3 because moving vertex 2 also removed
the intended insertion predecessor h-1. The
[two-C3 endcap repair](../../unprotected_ports_cancellation_20261025/ports/ENDCAP_NEXT.md)
subsequently closed h=3 with 3->1. It showed that the later port y can
serve as an insertion point while the middle remains literal.

The canonical hostile here is the old h=3 forward template: its second
move requires the false join 1->0, its third requires the false join 2->0,
and its edge XOR is wrong. The corrected near miss is using the first
incoming port as though its predecessor survives every extraction. The
least-used sidecar is the induced degree sum on a permutation orbit, when
global degree sums do not constrain both orbits separately.

The concept board is retained middle / later insertion port / two endpoint
degree changes / reflection-product cycles / induced orbit degrees / literal
edge XOR. The Anchor is the remaining h=3 branch. The Niche is a sparse
degree-change identity. The Wildcard is an even orbit whose single reversed
arc gives an exact degree total.

Let P=(0,...,n-1) be directed, closed and blue in T. Write

    a=0, b=1, x=2, y=n-3, c=n-2, z=n-1,
    N^-_T(a)={3,y,c,z}, T(b,3)=1, B=K_P(T).

K_Q retains the directed edges of its word Q and reverses all other pairs.
Blue means that Q's actual positional reversal is anti in its source.
Since P is blue, r_P(v)=z-v is anti in both T and B. In particular

    N^+_T(z)={a,b,x,n-4}, N^+_B(a)={b,3,y,c,z}.       (1)

**Theorem.** For every n>=12, all such sources have the following three
successive directed black K moves to the exact labelled target B.

## 1. One unconditional three-word ledger

Set U=(3,...,n-4), R=(y,c,z), and use

    Q1=(a,b,U,R,x),
    Q2=(b,U,a,R,x),
    Q3=(b,x,U,a,R).                                 (2)

The extracted vertex is still x=2, but a is inserted immediately before
the later incoming port y. Its predecessor n-4 has not been extracted.

Q1 uses b->3, z->x and original path joins. Q2's new joins (n-4)->a and
a->y are the complements of a->(n-4) and y->a. Q3 restores b->x and
x->3 after two complements and retains the other joins. The added pairs
b3,zx,(n-4)a,ay each occur twice, and the odd-multiplicity pairs are
exactly E(P). Hence

    E(Q1) xor E(Q2) xor E(Q3)=E(P).                   (3)

This proves directedness and the literal target. No Hamiltonian word has
been changed for free; all three intermediate words are checked below.

Q1 is black by a direct closing-edge witness. Its reversal would equate
T(a,z)=0 with T(b,x)=1.

## 2. Q3: its actual orbit has two outgoing ports

At the target B, Q3 has word (1,...,n-4,0,y,c,z). Let g3=r_Q3 r_P.
Directly composing the two reversals gives

    0->1->2->3;
    g3(v)=v+2 for 3<=v<=n-6;
    g3(n-5)=0;
    n-4->y->c->z->4.                                (4)

These are an explicit extension of the earlier insertion rotation
calculation to insertion at n-3; no out-of-range use of that lemma is
needed. If n is odd, (4) is one n-cycle. If Q3 were blue, g3 would be a
transitive automorphism of B. Its outdegree at a is five, forcing n=11.

If n is even, the orbit of a has size n/2 and is

    {0,1,2,3,5,7,...,n-5}.

Exactly b and 3 from (1) lie in it. A cyclic automorphism of that induced
tournament would make every outdegree two, forcing n/2=5, hence n=10.
Neither alternative occurs for n>=12. Thus Q3 is black at B, equivalently
at its own preceding source.

## 3. Q2: two endpoint degree changes control the odd orders

Let S=K_Q1(T), and compare it to B. Their pair difference is

    F={12,23,13,2z}.                                (5)

Write beta_v=d^+_B(v). Only two degrees change from B to S:

    d^+_S(v)=beta_v+delta_v,
    delta_2=-1, delta_z=+1.                         (6)

This also follows because P and Q1 have the same first vertex and their
last vertices are z and 2. Suppose, for contradiction, that Q2 is blue,
and put g=r_Q2 r_P. The two anti equations imply

    beta_(gv)-beta_v=-delta_(r_Pv)-delta_(gv).       (7)

The nonzero increments are -1 at v=0,n-4, and +1 at v=n-3,n-2.
Also beta_0=5. The permutation itself is explicitly

    g(0)=3, g(1)=4, g(2)=5;
    g(v)=v+4 for 3<=v<=n-8;
    g(n-7)=0, g(n-6)=y, g(n-5)=c, g(n-4)=z;
    g(y)=1, g(c)=2, g(z)=6.                          (8)

For odd n, it is a single cycle. All progressions in the next table step
by four; isolated endpoints are written separately.

| Order | Cycle starting at 0 | Sum of offsets from beta_0 |
|---|---|---:|
| n=4t+1 | (0,3,7,...,n-6,y,1,4,8,...,n-5,c,2,5,9,...,n-4,z,6,10,...,n-7) | 0 |
| n=4t+3 | (0,3,7,...,n-4,z,6,10,...,n-5,c,2,5,9,...,n-6,y,1,4,8,...,n-7) | -n |

In the first row, the t positions from 3 through y have offset -1,
the t positions from 2 through n-4 have offset +1, and all others have
offset zero. In the second row the increments in (7) give total offset
-n. Equivalently the first segment 3,...,n-4 contributes -t, the segment
z,6,...,n-5,c contributes -2(t+1), and 2,5,...,n-6,y contributes -(t+1).

The mean outdegree of an n-vertex tournament is (n-1)/2. The first row
would therefore force n=11, and the second n=9. Both contradict n>=12.

## 4. Q2: the two even residue classes

If n=4t+2, (8) has a cycle

    (0,3,7,...,n-7).

Only v=0 on this cycle has a nonzero increment in (7), and it is -1.
A degree cannot decrease by one upon returning to the same vertex.
Thus Q2 is black in this residue class.

Now let n=4t. The g-orbit of a is the 2t-element set with cyclic order

    C=(0,3,7,...,n-5,c,2,5,9,...,n-7).               (9)

Its complement is r_P(C), and is also a g-orbit. Among the four pairs
in (5), only 23 is internal to C; none is internal to r_P(C). In B the
arc 2->3 is an original retained path edge. Thus the induced tournament
S[C] differs from B[C] only by reversing 2->3, while S[r_P(C)] equals
B[r_P(C)]. Let alpha_v be the outdegree inside B[C]. Combining the two
anti maps on these two sets gives

    alpha_(gv)-alpha_v = -1 if gv=3,
                         +1 if gv=2,
                          0 otherwise.             (10)

Here r_P maps C to its complement and r_Q2 maps that complement back to
C; this is an induced-degree equation, not an assumed automorphism of
B[C]. Exactly t successive vertices, from 3 through c in (9), have
degree alpha_0-1; the other t have degree alpha_0. By (1), the outgoing
neighbours of a within C are precisely 3 and c, so alpha_0=2.

The degree sum from (10) is therefore 4t-t=3t. An actual tournament on
2t vertices has degree sum t(2t-1). Equality forces t=2, or n=8, contrary
to n>=12. This settles the last case and proves Q2 black. Together with
sections 1--2, it proves the theorem.

## 5. Exact controls and retained failure boundaries

The checker quantifies every source completion in the stated affine
system for n=12..100, plus 120 and 160; it does not enumerate Hamiltonian
paths or isomorphism classes. It separately replays two literal models
per order, checks (4),(8), the odd degree sums, the nonclosing degree
cycle, and the induced orbit equation through order 500.

Two genuine lower-order failures are retained: at n=8 the second word
can be blue, and at n=10 the third can be blue. These are template
failures, not connectivity counterexamples. The checker also preserves
the original h=3 insertion error, including its two forced false joins
and wrong XOR. No failed source is silently filtered from the theorem.

Run `python check.py` and `python -O check.py` from this directory. All
checks use explicit exceptions. The proof concerns literal matrices and
vertices throughout, so arbitrary middle arcs and existing clone identities
remain present. Together with the preceding h=3 backward theorem, this
closes the entire h=3 endpoint-four family for n>=12; it makes no claim
about the separately assigned h=4 residual or the global conjecture.
