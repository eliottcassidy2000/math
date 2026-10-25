# Removing the n-4 port: the two terminal ports suffice away from two degrees

2026-10-25. **PROVED** structural three-black detours; **FINITE-EXACT** controls
in `check.py`. Full black-core connectivity remains **OPEN**. The two critical
degree cases below are failures of this ledger, not connectivity obstructions.

## Inheritance and scope

The closest proved mechanism is
[the protected-port theorem](../../ports_regular_actions_hn_20261025/ports/NEXT_BOUNDARY.md).
It requires n-4 and n-2 to be incoming interior ports. The present construction
replaces n-4 by n-3 and splits on the first chord 1 versus 3. The canonical
hostile is a feasible source on which a legal intermediate word is genuinely
blue. The corrected near miss is mistaking a finite absence of such sources
for an unrestricted claim: the general family has two explicit degree loci.
The least-used sidecar is the degree change between two path-complement moves
with the same first vertex and different last vertices. It has support two,
even though their edge difference has support four.

Anchor: endpoint-four nonconsecutive ports. Niche: sparse degree changes and
signed reflection products. Wildcard: antipodal pairs inside even permutation
cycles. The concept board is literal path words / first port / endpoint degree
defect / affine anti equations / exact edge XOR / unchanged clone labels.

Let P=(0,...,n-1) be directed, closed and blue in T. Set
a=0,b=1,x=2,y=n-3,c=n-2,z=n-1. Let

    N^-_T(a)={z} union I, k=|I|+1, h=min I, B=K_P(T).

Thus N^+_T(z)={a} union r_P(I), where r_P(v)=z-v. Every K move retains
its word's directed path edges and complements all other unordered pairs.
Blue means its actual positional reversal is anti. All vertices remain
literal; no clone contraction, relabelling quotient, or middle filter occurs.

**Main theorem.** If n>=12 and

    I={h,n-3,n-2}, 5<=h<=n-4,

then T has three successive directed black K moves to the literal target B.
No closed-middle or additional chord assumption is required. In the inherited
large-order endpoint-four remainder, only h=3 or h=4 survives this theorem.
The h=4,1->3 branch is additionally closed by section 1a below.

**General theorem.** Suppose n>=12, k>=4, h>=5 and {y,c} subset I. The same
conclusion holds whenever

    k is odd, or n is neither 2k+1 nor 2k+3.                 (G)

It also holds for the exact consecutive tail at n=2k+3 by the inherited
tail theorem. In particular n>=2k+4 is always covered. The main theorem
follows with k=4, whose only excluded degree values are 9 and 11.

## 1. The chord 1->3: remove 2, move the cut, restore 2

Write U=(3,...,h-1), R=(h,...,z). Use

    A1=(a,b,U,R,x),
    A2=(b,U,a,R,x),
    A3=(b,x,U,a,R).                                      (A)

The first new joins are b->3 and z->x, supplied by the chord and y in I.
The second move inserts a between h-1 and h. The requisite joins are the
complements of a->(h-1) and h->a. The third restores b->x and x->3 after
two complements, retaining the other joins. The nonpath boundary pairs
occur twice; direct edge cancellation gives

    E(A1) xor E(A2) xor E(A3) = E(P).                     (1)

A1 is black without a degree calculation: its reversal sends (a,z) to
(b,x) in reverse order. T(a,z)=0 whereas T(b,x)=1.

For A2, work in S=K_A1(T). Its reversal sends z to 3. If 2h<=n+2 it
sends a to w=n-h+2, whose P-reflection is h-3. If 2h>=n+4 it sends a to
w=n-h+3, whose P-reflection is h-4. Both reflected values lie below h,
so w is not a z-outport. The pairs a3 and zw were omitted by A1; w is
neither c nor x. Hence

    S(a,3)=0, S(z,w)=1,                                 (2)

which contradicts an anti-reversal. The intervening case 2h=n+3 fixes a.
Since d^+_S(a)=k+1, blue would force n=2k+3 and h=k+3. Then cardinality
forces I to be the entire consecutive tail {k+3,...,2k+1}.

A3 has word (1,...,h-1,0,h,...,z). The
[inherited rotation lemma](../../ports_regular_actions_hn_20261025/ports/README.md)
and the subsequent boundary refinement show that it is black unless this
same exact consecutive-tail case occurs. In detail, disjoint supports would
force h=3, even overlapping supports would force n=6, and odd overlapping
supports require n=2k+3. At that order h<=k+3; h<=k+1 is disjoint, while
h=k+2 has B(a,b)=1 versus B(z,a)=0 under the reversal fixing a. Only h=k+3
remains. The
[all-k consecutive-tail theorem](../../trijection_hourglass_ports_20261025/endpoint/GENERAL_TAIL.md)
already gives a three-black detour for that entire family. Thus (A) needs
no new cyclic exception repair.

## 1a. The first omitted boundary: h=4 with 1->3

The forward ledger (A) also works when h=4, n>=12 and {y,c} subset I,
with arbitrary additional ports. Its legality still holds: the insertion
predecessor is 3, which was retained when vertex 2 moved. A1 has the same
closing-edge witness. The preceding A2 witness does not extend verbatim:
here w=c and the edge c->z was retained, so S(z,c)=0. Instead its reversal
maps (3,4) to (y,z), and

    S(3,4)=1, S(y,z)=0.

The first is a retained original path edge; the second complements the
original T(y,z)=T(a,2)=1. A3 has two disjoint four-cycles in its reflection
product with r_P, impossible for a tournament automorphism. Thus all three
moves are black. No claim is made here for h=4 with the opposite chord.

## 2. The chord 3->1: a different three-block ledger

Put V=(4,...,z). Use

    D1=(a,V,b,x,3),
    D2=(x,a,V,b,3),
    D3=(x,a,b,3,V).                                      (D)

D1 uses a->4 (h>=5), z->b (c in I), and retained path joins. D2 introduces
x->a and b->3, the complements of a->x and 3->b. D3 restores a->b and
3->4 after two complements and retains its other joins. Each nonpath pair
again occurs twice, and the edge XOR is E(P).

D1 is black: its anti-reversal would equate T(a,x)=1 and T(4,3)=0.
For D3 work at B. Its anti-reversal would equate

    B(a,3)=0, B(n-4,c)=1.                               (3)

The second bit uses T(n-4,c)=T(b,3)=0 from the original blue reflection
and the current chord. Both tested pairs are nonpath pairs of P.

Only D2 requires a longer argument. The next section proves it can be blue
only when k is even and n is 2k+1 or 2k+3. This proves (G).

## 3. The degree defect and the exact two critical orders

Let S=K_D1(T), and compare it to B. Their flipped-pair difference is

    F={01,34,04,1z}.                                    (4)

Write b_v=d^+_B(v) and d^+_S(v)=b_v+delta_v. All changes cancel except

    delta_3=-1, delta_z=+1.                             (5)

For example the two pairs at a remove its edge to 1 and add its edge to 4.
The description also follows directly from the endpoint degrees of K moves:
D1 and P have the same first vertex, and their last vertices are 3 and z.

If D2 is blue, let r be its reversal and g=r r_P. The two anti equations give

    b_(g v)-b_v = -delta_(r_P v)-delta_(g v).            (6)

The increments are -1 at v=0,n-5, +1 at v=n-4,n-3, and zero elsewhere.
Also b_0=k+1.

For odd n, g is a single cycle. Its exact cycles and the sum of its degree
offsets from b_0 are as follows; all displayed progressions step by four:

| Order | Cycle of g | Sum of offsets |
|---|---|---:|
| n=4t+1 | (0,4,...,4t,1,5,...,4t-3,2,6,...,4t-2,3,7,...,4t-1) | -n |
| n=4t+3 | (0,4,...,4t,3,7,...,4t-1,2,6,...,4t+2,1,5,...,4t+1) | 0 |

For the first row, 0 has offset zero; the nonzero multiples of four other
than z have offset -1; z and the residue-1 progression have offset -2;
the residue-2 progression has offset -1; residue 3 has offset zero.
For the second row, the nonzero multiples of four have offset -1, the
residue-2 progression except z has offset +1, and all others have offset
zero. This proves the table directly from (6).

The average tournament outdegree is (n-1)/2. The two rows therefore force
respectively

    n=2k+1 when n=4t+1, or n=2k+3 when n=4t+3.          (7)

In both cases k=2t is even. These are necessary conditions, not merely
upper bounds on the size of a finite search.

## 4. Even orders: an unchanged antipodal arc cannot reverse

Still under the assumption that D2 is blue, the exact arc form of the two
anti equations is

    B(gu,gv) = B(u,v) xor 1_F({r_Pu,r_Pv})
                            xor 1_F({gu,gv}).          (8)

For n=4t, g has two cycles of length 2t, on the even and odd labels. On
each cycle the step is +4 until it switches to the other residue of the
same parity. The only F pairs internal to these cycles are 04 and 1z;
both are adjacent cycle pairs. For n=4t+2, the four g-cycles are the
four residue classes modulo four, of lengths t+1,t+1,t,t. Again 04 and
1z are the only internal F pairs, and they are adjacent. The other two
F pairs join different cycles.

At n>=12, choose an even cycle of g with length L>=4. For n=4t take
either cycle. For n=4t+2 take a residue-0 cycle if t is odd and a
residue-2 cycle if t is even. Its antipodal pair {u,g^(L/2)u} and its
entire pair orbit avoid F. Their r_P images also avoid F: the identity
r_P g r_P=g^(-1) maps a cycle to a cycle and antipodal pairs to antipodal
pairs; such pairs cannot be adjacent when L>=4.

All right-hand flip terms in (8) therefore vanish along this orbit.
After L/2 steps g swaps the two endpoints of an unchanged tournament
arc. That is impossible. Thus D2 is black at every even order in scope.
Together with (7), this proves the general and main theorems.

## 5. Hostiles, exact universe, and frontier

The two odd critical degrees are genuine blue-word possibilities. The
checker retains sources with D2 blue at n=13,I={5,7,9,10,11} (k=6),
and n=15,I={6,8,10,12,13} (k=6). Both satisfy all new source assumptions
and the chord 3->1. Their ledgers are legal and have the correct target,
but D2 cannot be certified black. This is not a disconnected-component
claim, and no source was excluded based on a numerical outcome.

The h>=5 guard is also real for (D): if h=4 then a->4 is false at its
first new join. For example n=12,I={4,9,10},3->1 is a feasible source
with that forced failure. This explains why the h=3,4 endpoint-four
remainder is not closed by (D). Combining section 1a, the remaining
endpoint-four source cells are h=3, and h=4 with 3->1, subject to the
previously inherited frontier guards. The backward h=3 ledger is legal
and its first two moves are black, but its final word can genuinely be
blue; the checker retains this further failed-template control.

Run `python check.py` and `python -O check.py` in this directory. The
checker imports the previous exact signed-equation engine, supplies these
new ledgers and source cells, and separately replays literal arc sets.
It checks the complete endpoint-four family for n=12..40 and every h;
all port subsets containing y,c with h>=5 for n=12..17 when (G) applies;
and the degree/permutation/antipodal formulas for n=12..300. Exact
consecutive tails are routed to the inherited protected-port certificate,
not removed from the source universe. All quantifiers are supported by
the proofs above; finite checks are corroboration.
The added h=4 forward check covers every extra port subset at n=12..17.
