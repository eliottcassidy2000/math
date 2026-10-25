# The h=4, 3-to-1 residual has a three-black detour

2026-10-25. **PROVED:** the entire endpoint-four source family below has
a literal three-black detour. **FINITE-EXACT:** complete affine source
cells and independent literal replay corroborate the proof. Full
black-core connectivity and the h=3 residual are not asserted here.

## Inheritance and the changed cut

The inherited [unprotected-port theorem](../../unprotected_ports_cancellation_20261025/ports/README.md)
handles every h>=5 and the h=4,1->3 branch. Its backward ledger fails at
h=4 because its first new join is the false arc 0->4. This is the
canonical hostile and corrected near miss: retaining the same block
boundary cannot fix that ledger merely by a different blackness proof.

The new cut moves the whole block 1,2,3,4 to the end, making the first
new join 0->5. The least-used sidecar is the pair of actual nonports
6 and n-4: they give a uniform middle-word anti-reversal contradiction.
No degree estimates, source symmetry quotient, clone contraction, or
Hamiltonian-path census is needed. The concept board is the movable cut /
literal path-edge XOR / one-row port data / actual anti-reflection.

## Exact source and theorem

Let n>=12 and let P=(0,1,...,n-1) be a directed, closed blue path in a
tournament T. Thus r_P(v)=n-1-v is an anti-automorphism and n-1 points
to 0. Write a=0, b=1, x=2, y=n-3, c=n-2, z=n-1, and suppose

    N^-_T(a)={4,y,c,z},        T(3,1)=1.

Equivalently, the interior incoming ports are I={4,y,c}. Let V=(5,...,z).
Then the following three words are directed and black at their respective
sources, and their successive K moves end at B=K_P(T):

    Q1=(a,V,b,x,3,4),
    Q2=(x,a,V,b,3,4),
    Q3=(x,a,b,3,4,V).                                  (1)

Here K_Q retains exactly the directed path edges of Q and reverses all
other arcs. Every vertex and arc in (1) remains labelled throughout.

## Legality and literal target

The new joins of Q1 are a->5 and z->b. The first follows because 5 is
not an incoming port. The second follows by blue reflection from c->a.
All other Q1 joins are original P edges.

In S=K_Q1(T), Q2 introduces x->a and b->3. Both corresponding pairs
were omitted by Q1, and their T directions were a->x and 3->b. All
remaining Q2 joins are retained from Q1. Q3 introduces a->b and 4->5;
both original P edges were omitted by Q1 and Q2 and so have been reversed
twice. Its other joins are retained from Q2.

For a word Q let E(Q) be its unordered path-edge set. Put

    X={01,45,05,1z},        Y={12,23,02,13}.

These are actual unordered pairs. Direct comparison gives

    E(Q1)=E(P) symmetric_difference X,
    E(Q2)=E(P) symmetric_difference X symmetric_difference Y,
    E(Q3)=E(P) symmetric_difference Y.

Consequently E(Q1) XOR E(Q2) XOR E(Q3)=E(P). Three K moves therefore
flip exactly the same pairs as K_P. Together with the separately proved
legality, this proves the literal target B, not merely an isomorphism.

## Three uniform blackness witnesses

Write r_i for the actual positional reversal of Q_i. A word is blue at
its source exactly when it preserves every arc after reversing direction.

**First word, at T.** Its reversal sends the ordered pair (a,3) to the
reversed pair (5,4): r_1(3)=5 and r_1(a)=4. But

    T(a,3)=1,        T(5,4)=0,

since 3 is not a port and 4->5 is an original path edge. Thus Q1 is black.

**Middle word, at S.** Its reversal sends (3,z) to (6,a): r_2(z)=6 and
r_2(3)=a. Blue reflection in T gives

    T(3,z)=T(a,n-4)=1.

The equality uses the actual P reflection; the final value uses the fact
that n-4 is not in {4,y,c,z}. Also T(6,a)=0 because 6 is not a port.
Neither pair {3,z} nor {6,a} occurs in Q1, so

    S(3,z)=0,        S(6,a)=1.                         (2)

This contradicts the proposed anti-reversal and proves Q2 black for all
n>=12. There are no exceptional residue classes or critical degree cases.

**Last word, checked at B.** The actual reversal sends (a,3) to
(n-4,c). Both tested pairs are nonpath pairs of P. Therefore

    B(a,3)=0,
    B(n-4,c)=1-T(n-4,c)=1-T(1,3)=1.

The last equality uses 3->1, and the middle equality uses the original
blue reflection. Thus Q3 is black at B. Blue/black status is the same
at both ends of its K edge: positional reversal preserves E(Q3), so
complementing its remaining pairs preserves the anti test. Q3 is
therefore black at its actual preceding source as required. QED.

## Failure boundary and scope

The old first join 0->4 remains false throughout this source family.
The new construction repairs it by starting at 5 and preserving the
extra vertex 4 until the last surgery. The exact port set matters: if
5 is added as an incoming port, the first new join 0->5 is false. The
checker includes that feasible outside-scope hostile; it is a failure
of this word, not a connectivity counterexample.

Combined with the inherited forward h=4 branch, this closes all h=4
endpoint-four sources at n>=12. In the inherited frontier the remaining
first-port value is h=3, subject to its previously stated guards. This
does not claim full black-core connectivity.

## Reproduction

Run `python check.py > normal.json` and `python -O check.py > optimized.json`.
The affine engine checks all completions in every source cell n=12,...,120,
with no middle-orientation filter. A separate arc-set implementation
replays three seeded literal models per order, checking legality, every
positional anti equation, and the exact final labelled target. Reversal
formulas and the edge-set rectangle are also checked through n=1000.
No Hamiltonian words are enumerated beyond the three specified words.
