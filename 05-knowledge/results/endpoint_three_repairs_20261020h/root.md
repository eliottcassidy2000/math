# Endpoint-three repairs: orbit capacity and open middles

Research continuation, October20,2026. **PROVED** for the explicitly
stated branches. Together with the companion branches these give the
[complete endpoint-three theorem for n>=11](README.md).

## Inheritance and working coordinates

The closest proved mechanism is the complete
[endpoint-two theorem](../endpoint_two_repairs_20261020g/README.md).
Its hostile is the full-tiling suffix whose only exit from its prescribed
first vertex is the prescribed last vertex. Whole-block transport repairs
that obstruction. The corrected near miss is exporting a singleton's
positional-reversal identity to a longer block. The least-used sidecar
is the actual exceptional incoming set, together with its reflected
outgoing set and their possible overlap.

Our board is: exceptional port sets; path-edge XOR; actual anti-involutions;
orbit score capacity; block boundary joins; and operation-specific debt.
The anchor is the endpoint-three layer. The niche is the sharp cyclic-cut
capacity bound. The wildcard is whether the recovered 3-4-3 and 4-2-4
decompositions carry any of these operations, rather than just equal counts.

Write K_Q(T)=T^op xor E(Q), and R_Q(T)=T xor E(Q). A Hamiltonian word
Q is blue iff its positional reversal r_Q is anti. Let

    P=(a,M,z),  d_T(a)=n-1-k,  d_T(z)=k,  z->a,

be blue and closed. There is an exceptional incoming set I in M of size
k-1 such that

    N^-_T(a)={z} union I,    N^+_T(z)={a} union J,
    J=r_P(I).

Conversely a blue middle and these mirrored boundary incidences determine
the whole source. The path condition requires first(M) not in I. No
disjointness of I and J is assumed. For endpoint three, k=3 and |I|=2.

## R1. A sharp all-defect cyclic-cut bound

**PROVED.** For k>=2, every cyclic cut of P obtained by moving its first
vertex to the end is black if n>4k-2. The bound is sharp for every k>=2.
The same conclusion holds for any directed word with that positional
reversal.

If the cut were blue, composing its reversal with r_P would give the
automorphism shifting positions by two modulo n. For odd n it is
transitive, forcing n-1-k=(n-1)/2 and n=2k+1.
For n=2m its two parity orbits have size m. The induced orbit tournament
containing a is regular, so a has (m-1)/2 outneighbors inside it. Its
opposite orbit contains z, which beats a. Hence

    2m-1-k <= (m-1)/2 + m-1,

which forces m<=2k-1, or n<=4k-2. This includes the obstruction from
even m, since no even-order tournament is regular.

For sharpness, put m=2k-1 and use vertices e_i,o_i for i modulo m.
Within each parity class orient i->j when j-i lies in {1,...,k-1}.
Orient every cross pair e_i->o_j except o_(i-1)->e_i. Then

    P=(e_0,o_0,e_1,o_1,...,e_(m-1),o_(m-1))

is closed and blue. Reflection sends e_i to o_(m-1-i) and o_i to
e_(m-1-i), preserving the required anti equations both within and
between parity classes. Simultaneous cyclic shift of both indices is
an automorphism, so the one-step cyclic cut of P is blue too. Every
even vertex has score (k-1)+(m-1)=3k-3=n-1-k; every odd vertex has
score (k-1)+1=k. The endpoint-three boundary is therefore a genuine
ten-vertex obstruction to this particular blackness argument. It is
not a metagraph bridge.

## R2. Open-middle detours away from the penultimate port

**PROVED, arbitrary k>=2.** Write P=(a,b,W,z), x=first(W), c=last(W),
and let y precede c. Suppose b->c, y not in I, and n>4k-2.
Then the following three K moves are directed and black, and their
literal target is K_P(T):

    Q1=(b,W,z,a),
    Q2=(z,W,b,a),
    Q3=(a,z,W,b).

The reflected exclusion gives x not in J, so x->z.
Together with b->c these are the exact guards of the inherited cyclic
ledger; its path-edge XOR is E(P). Q1 is black by R1. Q2 has endpoint
score sum n-2 and is black.

To prove Q3 black without a degree shortcut, work at B=K_P(T), where P
is still blue. The product r_Q3 r_P is exactly (a c)(b z), fixing all
other vertices. If Q3 were blue this nontrivial involution would be a
tournament automorphism, impossible: any transposed pair has its arc
reversed. Thus Q3 is always black whenever it is a directed word,
independently of endpoint scores or the location of c in I. In the
requested k=3 layer, the threshold is n>=11.

## R3. Last-port two-block repair with a second exceptional port

**PROVED.** In the endpoint-three layer, write

    P=(a,b,x,D,c,z),  I={c,h},  x->c.

Assume n>=11 and h is neither x nor the P vertex at position n-4.
Then the three R words

    Q=(a,x,D,c,z,b),
    S=(a,b,x,c,reverse(D),z),
    V=(b,z,D,c,x,a)

are directed black moves to R_P(T). Equivalently, Q,reverse(S),V are
three black K moves to K_P(T).

Use pivot b, left block (a), and right block (x,D,c,z) in the inherited
two-block theorem. After reversing the displayed right-block path, its
endpoint-preserving repair is (x,c,reverse(D),z). The first chord is
assumed. Its final join first(D)->z holds because J={b,r_P(h)} and
h!=P[n-4], the reflection of first(D)=P[3]. The remaining new source
join a->x is retained by h!=x.

Q has the cyclic-cut reversal and is black by R1. S has endpoint sum
n-2. At the other endpoint of the last R move, R_P(T), z->a remains
true and b->x is false because it was a P arc. Reversal of V would
equate these two arcs. Thus the last move is black irrespective of the
second exceptional port. The actual port exclusions are not dropped.

## R4. Complete penultimate-port block repair

**PROVED.** Let P=(a,b,x,D,c,z) have endpoint scores (n-4,3), with
b->c, I={y,h}, y=last(D), c not in I. If n>=11, the following are three
black K moves to K_P(T), for every admissible second port h:

    Q1=(a,b,c,z,x,D),
    Q2=(b,c,a,z,x,D),
    Q3=(b,x,D,c,a,z).

Legality and the exact XOR are the inherited whole-block penultimate
ledger. They only require b->c, z->x, and a->c; these follow from the
stated ports. The internal tournament of D is arbitrary.

For Q1, if it were blue, g=r_Q1 r_P would be an automorphism. Its
cycles have lengths: n for odd n; two cycles of n/2 for n=2 modulo4;
four cycles of n/4 for n=0 modulo4. This follows by tracing

    g(0)=n-6, g(1)=n-5, g(2)=0, g(3)=1,
    g(4)=n-2, g(5)=n-1, g(i)=i-4 (6<=i<=n-3),
    g(n-2)=n-4, g(n-1)=n-3.

In the odd case regularity forces n=7. In the n=2 modulo4 case the
orbit containing a is the even positions; neither z nor y belongs to
it. Thus a has at most one inneighbor within that orbit, forcing its
size at most3, or n<=6. For n=0 modulo4 its orbit is
{0,2,6,10,...,n-6}, again excluding z and y. The same argument forces
its size n/4 at most3, hence n<=12. So Q1 is black for n>=13. At n=12,
its only possible exceptional h values are2,6, the other members of a's
orbit. But Q1 reversal sends a to y and P[n-5] to c. Unless h=P[n-5],
it would equate the true arc a->P[n-5] with the false arc c->y. At
n=12 this forces h=7, contrary to h in {2,6}. The odd case already
excludes n=11. Thus Q1 is black throughout the stated range.

At S=K_Q1(T), N^-_S(z)={a,c,r_P(h)}. If r reverses Q2, then
r(z)=P[n-6] and r(a)=P[n-5]. If h is neither4 nor n-6, the
second vertex is not r_P(h), and a->r(z) in T.
That pair is absent from Q1, so r(z)->a in S, while r(a)->z is false.
This contradicts antisymmetry. It remains to handle h=4 or h=n-6.

For n=12,13 or n>=15 put t=P[n-10], w=P[n-9], s=P[n-5]. Anti-reversal
of Q2 would give

    S(a,t)=S(P5,s),    S(z,w)=S(P4,P[n-6]).

All four pairs are nonpath in Q1. Blue reversal of P gives

    T(P5,s)=T(P4,P[n-6]),    T(w,z)=T(a,P8),

where T(u,v) denotes the zero/one arc indicator. Combining these exact
equalities forces T(a,t)=1-T(a,P8). Both a-arcs are true for either of
the two exceptional h choices in these orders, a contradiction.

At n=14 the anti equations for Q2 position pairs (0,8) and (1,9) give
T(1,6)=T(3,11) and T(12,7)=T(2,10), since the pairs are nonpath in
Q1. P's reversal gives T(1,6)=1-T(12,7) and T(2,10)=T(3,11), an
immediate contradiction independent of the ports.

At n=11 and h=4 or5, S(z,8)=S(z,7)=1: those middle vertices beat z
in T and neither pair is retained by Q1. Q2 anti equations at position
pairs (0,7),(1,7) would give S(1,5)=S(9,5)=1. P reversal fixes5 and
exchanges1,9, while both pairs are complemented by Q1. Hence these
last two indicators sum to1, a contradiction. This proves Q2 black
in every case without enumerating tournaments.

Finally, at B=K_P(T), Q3 reversal would equate the true path arc a->b
with the false arc z->x. Thus Q3 is black, without a further guard.
