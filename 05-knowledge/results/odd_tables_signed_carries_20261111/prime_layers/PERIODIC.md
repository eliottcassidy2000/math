# A sevenfold paired actuator with an arbitrary-prefix preparation

**PROVED IN NOTE:** an actual H13/C3 paired rewrite under multiplier7,
its additive defect loop, and a preparation for every finite prefix.
**PROVED BY COMPLETE FINITE CERTIFICATE:** the multiplier5 graph with
these weights has only two cyclic states, both zero-drift self-loops.
**FINITE-EXACT:** named literal tournament and arithmetic controls.
**OPEN:** an equality-preserving correction at an arbitrary original
target. This extends the [prime-layer compiler](PROOF.md) and complements
the [general phase-bank theorem](../actuators/PROOF.md).

## 1. The actual target equation

Retain the inherited H13 pair and second C3 pair:

    L(B)=13g(B)+32B+6(B mod2)+18,      L(B)+2,
    N(C)=3g(C)+12C+10,                N(C)+2.         (1)

Here `g(0)=g(1)=1`, `g(2q+1)=g(q)`, and
`g(2q)=g(q)+q+(q mod2)`. The two H13 charts use core
`(18,8,3,5,14)` with cuts6 and3 and right address B. The C3 charts use
core `(2,4,1)` with cut3/right address `2C+1`, and cut1/right address
`4C+3`. These are literal constructors from the inherited
[residue-phase note](../../h13_paired_descent_20261111/square_graceful/RESIDUE_PHASE.md).

**Sevenfold identity.** Put

    B=64n+36,          C=224n+135=32(7n+4)+7.

Then the difference for **both** pairs is

    E(n)=L(64n+36)−N(224n+135)
        =13g(n)−3g(7n+4)−n−3(n mod2)−80.           (2)

To verify it symbolically, binary recursion gives

    g(64n+36)=g(n)+54n+32,
    g(32x+7)=g(x)+3x+(x mod2).

Substitute x=7n+4. The source constant is1586 and target constant1666;
the coefficient difference of n is−1. This proves (2) with the actual
charts retained. In particular the constant is−80, not a freely chosen
offset that could be adjusted independently of the tournaments.

## 2. An exact loop, including its integer multiplication carry

The state `(j,u,v)=(4,−1,−3)` represents the nonconstant part
`13g(n)−3g(7n+j)+u n+v(n mod2)`. Pulling it through the four bits
of suffix1001 returns that state and adds4 to the constant. Equivalently,

    E(16n+9)=E(n)+4.                               (3)

For a direct derivation, use

    g(16n+9)=g(n)+6n+4,
    7(16n+9)+4=16(7n+4)+3,
    g(16x+3)=g(x)+3x+(x mod2).

The last identity at x=7n+4 gives the target increment `21n+12+(n mod2)`.
Substitution in (2) cancels the n and parity changes and leaves4.

Define

    n_m(q)=16^m q+9(16^m−1)/15.

Then, for every q,m>=0,

    E(n_m(q))=E(q)+4m.                             (4)

The ordinary-binary source suffix is1001. Multiplication by7 sends its
periodic phase `9/15=3/5` to `1/5`, with integer carry4:

    7·9 = 4·15+3.

The corresponding target suffix is0011. This is the reverse phase
direction of the multiplier3 example. The common denominator5 is the
finite multiplicative clock, while j=4 is its retained translation.
Discarding j would replace g(7n+4) by the wrong function g(7n).

This is a concrete use of the prime-layer method: multiplication by7
is the same finite phase as multiplication by2 modulo5, but their
integer carry and source scaling differ. The period relation selects
the phase; retaining the integer affine action supplies (2)-(3).

## 3. Every prefix admits a positive number of loop repetitions

**Theorem.** Given any integer prefix p>=0, form

    q1=16p+7,       q2=32p+15.

Exactly one of E(q1), E(q2) is divisible by4. Choose its q and put

    m=−E(q)/4,       n=n_m(q),
    B=64n+36,        C=224n+135.                    (5)

Then m is a positive integer and

    L(B)=N(C),       L(B)+2=N(C)+2.                 (6)

Thus both resulting H13 targets have the specified literal C3 witnesses.
The original prefix p is retained, followed by either0111 or01111,
then m copies of1001, and finally100100 in B. This family occupies the
even H13 address class `B=36 mod64`, unlike the earlier odd terminal0101
family. Its lower/upper target residues are1/0 modulo3.

**Proof.** For x>=0, direct recursion yields

    E(4x+3)=13g(x)−3g(7x+6)−25x−3(x mod2)−104,
    E(8x+7)=13g(x)−3g(7x+6)−50x−132.              (7)

Every g value is positive and odd, with `g(x)<=x+1−(x mod2)`.
Consequently both expressions in (7) are strictly negative: the first
is at most `−12x−16(x mod2)−94`, and the second at most
`−37x−13(x mod2)−122`.

Now take x=4p+1. The two arguments in (7) are q1,q2, and subtraction
gives

    E(q2)−E(q1)=−100p−50=2 mod4.                   (8)

Both defects are even, so exactly one is divisible by4. Negativity
makes m in (5) positive, and (4) gives E(n)=0. Equation (2) supplies
both actual witnesses. QED.

A shorter named seed is q=13, where g(13)=5 and g(95)=3, giving E(q)=−40.
Ten copies of1001 yield

    n=14953358137753,
    B=957014920816228,
    C=3349552222856807,
    H pair=(41199492341138701, 41199492341138703).

The actual C3 right addresses are6699104445713615 and13398208891427231.
These large targets are consequences of a short exact construction,
not an expanded tournament census. No smallest-seed claim is made.

The [new C3 rerooting theorem](../tournament/PROOF.md) applies to every
exact equality L(B)=N(C), so it applies here too. It supplies a beneficial
paired reset with the specified root change and subsequent all-Z language.
One must retain its changed full-state difference; it is not automatically
a pure-H gap under arbitrary root-changing words.

## 4. How the finite search found an actual channel

For any odd multiplier t, consider

    F(n)=13g(n)−3g(tn+j)+u n+v(n mod2)+w,
    0<=j<t.

For a pulled-back bit a∈{0,1}, write `ta+j=2J+r` and δ=J mod2.
The exact next state and constant increment are

    j'=J,
    u'=2u+13(1−a)−3t(1−r),
    v'=13(1−a)−3(1−r)(1−2δ),
    w'−w=(u+v)a−3(1−r)(J+δ).                     (9)

This is an affine carry search, not a search over H-value ranges. Every
closed nonconstant-state loop must stay in

    −13<=u<=3t,    v∈{−3,0,3,10,13,16}.            (10)

Indeed outside the u interval the recurrence moves strictly farther out,
while every transition places v in the displayed list. For t=3, the +3
case is absent; it is necessary at general odd multipliers. Thus a finite
graph captures every possible additive binary loop with these weights.

The first direct alternatives were actual channels, not arbitrary states:

    L(2n+1)−N(5n+j)
       =13g(n)−3g(5n+j)+4n+46−12j;
    L(2n)−N(7n+j)
       =13g(n)−3g(7n+j)−7n+13(n mod2)+8−12j.

None of their five or seven initial states reaches a cycle in its finite
return graph. This stopped that specific phase choice, but did not stop
the search for a different physical phase.

To realize a selected cyclic state, compare dyadically translated
actual addresses

    L(2^s n+a),       N(2^k(tn+j)+r).

Let x=rev_s(a), y=rev_k(r), with s,k>=1. The inherited block formula
`g(2^s n+a)=g(n)+(2^s−1−x)n+(1−x mod2)(n mod2)+g(a)−1`
gives

    u=45·2^s−15t·2^k+3t−13−13x+3ty,
    v=13(1−x mod2)−3(1−y mod2)(1−2(j mod2)).       (11)

For the t7 state `(j,u,v)=(4,−1,−3)`, choose

    s=6, k=5, x=9, y=28,
    a=rev_6(9)=36, r=rev_5(28)=7.

Equations (11) give the exact selected state. Evaluating the inherited
chart constants then gives w=−80, producing (2). This is the missing
realization step between a formal phase loop and an actual tournament
rewrite. The prepared prefix in section3 pays the resulting constant.

## 5. A complete multiplier-five stopping certificate

For t=5, the graph (9)-(10) has870 vertices and864 retained edges. Its
only cyclic states are

    (j,u,v)=(0,2,10), with the zero-bit self-loop;
    (j,u,v)=(4,0,0), with the one-bit self-loop.

Both have constant increment0. All other retained edges strictly decrease
the integer rank saved in `five_rank_certificate.json`. Consequently
**no multiplier5 binary actuator with weights13 and−3 has nonzero
additive drift**, at any word length, in this exact fixed-carry language.

This is a finite proof certificate for an all-length statement. Its
coverage is justified by (10), rather than by imposing a maximum word
length. It does not prohibit other weight pairs, other source operations,
or return laws that change the nonconstant state rather than returning it.
The obstruction also has a clear phase interpretation: multiplication
by5 destroys the nonzero denominator5 phases used by the successful
four-phase bank. The complete graph certifies that another denominator
does not silently restore an additive loop in this particular language.

## 6. Reproduction and the remaining target gate

```bash
python3 -B 05-knowledge/results/odd_tables_signed_carries_20261111/prime_layers/periodic.py
```

The script reconstructs the complete multiplier5 graph and its decreasing
rank certificate, the two direct-channel reachable graphs, and the
multiplier7 loop including its constant drift. It checks (2)-(3) on
n=0..255, the full preparation on p=0..63, and all four literal source
charts at n=0,1 by independent subset path counting. The canonical
constructor is read from its Git blob and its hash is recorded. Large
prepared targets use the all-address constructor theorem and identities,
not a claim of literal subset enumeration at their full orders.

The connection contract is: synchronized finite phase -> retained integer
carry -> coefficient loop -> actual chart-address realization -> prepared
constant cancellation. Each arrow preserves the paired target equality
it actually supplies. Padding changes the address B and therefore the
original target. The remaining Chan–Pak obligation is a way to enter or
use such a loop while retaining an arbitrary supplied target, not merely
to construct more exact families with prefixes of every type.
