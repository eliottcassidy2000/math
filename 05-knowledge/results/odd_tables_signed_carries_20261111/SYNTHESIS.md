# Modular carry clocks and paired tournament source resets

**PROVED IN NOTES and FINITE-EXACT as stated below. Chan–Pak remains OPEN.**
The new simultaneous reset reaches every binary prefix in two different
address channels. A root change then preserves both realized counts while
improving each source for all later permitted continuations. The modular
table connection supplies the exact clock and the integer carry behind
these operations.

## 1. Four phases give an infinite family of exact actuators

Retain the binary function

    g(0)=g(1)=1,
    g(2n+1)=g(n),
    g(2n)=g(n)+n+(n mod2).

For i=1,2,3,4, let a_i=i^(-1) mod5 in{1,2,3,4} and d_i=[i≤2].
The four integer potentials

    P_i(n)=5g(n)−a_i n+5d_i(n mod2)

obey

    P_i(16n+3i)=P_i(n)+2.

Their four blocks are0011,0110,1001,1100. For any odd multiplier m
coprime to5, write mi=5j+i' with1≤i'≤4. The synchronized address map

    m(16n+3i)+j=16(mn+j)+3i'

transports one phase to the other. If integer weights a,b satisfy a≡b
mod5, an explicit integral linear correction makes
`a g(n)−b g(mn+j)` increase by exactly `2(a−b)/5` on every repetition.
For the H13/C3 weights13 and3, the drift is4 at **every** odd multiplier
coprime to5. The translation j is essential.

This explains an appearance of five by a literal operation. The finite
phase bank stores the recursive pattern; an unbounded integer stores
the remaining defect. Neither coordinate is discarded.
[General phase theorem and complete ternary return graph](actuators/PROOF.md).

The [general phase-clock theorem](prime_layers/PHASE_CLOCKS.md) extends
this to every odd modulus. Bit reversal gives an exact integrality test;
doubling and negation preserve drift. If k is the order of2 modulo p,
the sum of the p−1 phase drifts is

    k(p−χ4(p))²/8.

Within the entire Fermat-shaped family2^m+1, only3 and5 have uniform
phase drift; within2^k−1, only3 and7 do. Both families have integral
coefficients at every size, so integrality and uniformity are distinct.
Clock11 supplies another integral uniform clock, with drift18, and an
actual infinite H25/C3 paired equality with drift36. These are proved
family boundaries, not a classification of all primes by a finite scan.

## 2. Two physical paired rebases and a uniformly beneficial root change

The actual H13 pair and C3 pair have counts

    L(B)=13g(B)+32B+6(B mod2)+18,       L(B)+2;
    N(C)=3g(C)+12C+10,                 N(C)+2.

The core, apex cuts, and literal flank addresses are retained in the
[tournament proof](tournament/PROOF.md). Two exact channels now have
arbitrary-prefix zero-defect preparations.

**Multiplier three.** For

    D(B)=L(B)−N(3B),
    F(n)=13g(n)−3g(3n)+n+13(n mod2)+20,

we have `D(16n+5)=F(n)` and `F(16n+3)=F(n)+4`. Two short parity
preparations choose a negative multiple of4 after any supplied prefix;
repeating0011 cancels it. Both H13 targets then have actual C3 witnesses.

**Multiplier seven.** A different literal pair of charts gives

    L(64n+36)−N(224n+135)
      =13g(n)−3g(7n+4)−n−3(n mod2)−80=E(n),
    E(16n+9)=E(n)+4.

For any prefix p, exactly one of16p+7 and32p+15 gives a negative
multiple of4. Repeated1001 cancels it. The source phase is3/5 and the
target phase is1/5, with integer carry4.
[Sevenfold construction and exact preparation](prime_layers/PERIODIC.md).

**Root reset.** Whenever L(B)=N(C), B≥1,C≥2, retain the last right-flank
vertex as the new root of each C3 witness. Write C0 for C with its leading
one removed, and ε for its second leading bit. The new lower state is

    (H,E,U,d)=(N(C),N(C0)−3ε,N(C)+3ε,N(C0));

the upper state adds2 to all four coordinates. Both preserve the current
target and strictly reduce the last three coordinates of

    q=(H,U,U−d,H+U−E)

relative to the H13 witnesses. Thus each remains cheaper under every common
O/Z/J continuation and fixed right-root amalgam. The proof works for both
address channels, without a leading-carry exception.

The smallest displayed reset is

    (133,13,141,13) → (133,25,133,25),
    (135,13,141,13) → (135,27,135,27).

Their exact q savings are(0,8,20,20) and(0,6,20,20). The new pair preserves
its unit gap along every all-Z continuation. General continuations change
the gap to2g(address), so shared ports require a further operation.
A literal asymmetric birth restores common ports(81,453,81) at counts
453 and399, giving a pure-H width54. This is a specific compatible pair,
not permission to subtract54 from an arbitrary target.

The [arithmetic content of the response](tournament/CONTENT.md) is an
additional retained resource. The133 branch keeps content1, while the135
branch's reset has content27, forcing all its permitted future counts to
be divisible by27. Smaller future counts alone do not imply better
coverage. Indeed appending and marking a global sink always attains the
absolute q minimum(H,H,0,H), but has content H. The useful source-selection
problem therefore retains both q and arithmetic content, as well as the
actual core and root.

**Composite clock341.** The H25 pair at odd address2n+1 has counts
`P(n)=25g(n)+244n+186` and P(n)+2. Against N(19n), its defect D19 has
a shifted phase F19 satisfying

    D19(1024n+3)=F19(n),
    F19(1024n+3)=F19(n)+4.

The repeated block0000000011 has periodic fraction1/341, and its target
block has fraction19/341. The common period10 synchronizes the factors
11 and31. A choice between the seven-bit preparations128p+27 and128p+28
makes F19 a negative multiple of4 after **every** prefix p. Repeating the
block cancels the defect. The concrete meeting at n890 gives targets
231421/231423, independently counted on all four literal graphs.
[All-prefix H25 reset and its response bounds](tournament/COMPOSITE_CLOCK.md).

The same last-flank root choice also improves the inherited repeated-9
H13/C3 family in residues0/2. Thus its earlier unfavorable apex-response
comparison is resolved by changing the retained root; its original-target
padding limitation remains.

## 3. Prime and composite table layers retain different information

The central block of an odd multiplication table depends only on its
modulus, not on primality. With m=4n+1 its entries are n and3n+1, each
twice; with m=4n−1 they are n and3n−1, each twice. Each row sums to m.
The inherited mod-eight decoder distinguishes the3n+1 and3n−1 sheets
and their edge clocks. Every nonzero row is reversible for a prime
modulus; for a composite modulus, the row's gcd with the modulus records
its information loss. The central rows remain reversible in both cases.

The new [prime-layer compiler](prime_layers/PROOF.md) replaces the generated
multiplier group `<2,3> mod m`, for gcd(m,6)=1, by the rank-two relation
lattice

    L_m={(u,v):2^u3^v=1 mod m}.

Three integers specify its triangular basis. Each prime-power digit adds
one linear constraint over the prime field. Coprime factors combine by
lattice intersection over an explicit shared quotient. At1729 the last
gluing quotient is C6, so the actual multiplier group has216 elements
instead of the local product1296. This is a concrete meaning of a shared
six-phase constraint.

The full affine action additionally retains translation. The3-power part
of an arbitrary odd modulus has an exact contraction/fibre law, completing
the description without treating3 as invertible there. The signed legal
ancestor construction realizes these finite actions from a supplied hub;
it does not assert descent from every prescribed source.

## 4. Infinite signed trajectories can be stored by one exact recurrence

Take an actual finite odd Collatz trajectory on either sign sheet and
choose an odd observer modulus separately at every marked time. The new
[observer-preserving lift](history/PROOF.md) computes the **least** period
that preserves all those residues, every valuation in the finite head,
and the entire tail after a selected splice. Its infinite source family
satisfies one explicit affine recurrence.

Choosing each observer to be that state's original table modulus makes
every marked table in every positive lift composite. For the plus route
3→5→1, the exact family begins

    n_0=3,       n_(k+1)=4096n_k+2275,
    14563 → 21845 → 1.

The marked new table moduli are58253=13·4481 and87381=21·4161.
The first valuation remains1; the final valuation is4+12k. The same
construction on a minus-sheet seed retains that seed's tail rather than
changing its basin. Thus prime/composite table factors and the actual
signed dynamics can be varied with precisely stated retained information.

## 5. Local repair now has an exact entry condition

The older square-sum port constructions and the strategy-cube heredity
argument suggested preserving a whole boundary response while moving a
block. This gives a complete H25 fixed-port grammar:

    (A,2^t C+x) → (2^t A+c,C),
    304c+25rev_t(x)=147(2^t−1).

It preserves E,U,d,K and graph order. The first solution is t7,c61,x80.
Every nontrivial move strictly lowers both branch-count sets, even with
independent branch toggles. A separate all-context entry-gap theorem
shows that its small lost count cannot be repaid by the permitted positive
continuations at that fixed state. The lesson is an exact boundary for
this operation: a usable repair must also change the source interface.
[Complete port classification and entry gap](history/PORT_REPAIR.md).

Incoming work [on target-preserving ready returns](../berggren_target_return_20261111/reset/TARGET_RETURN.md)
supplies another actual baseline. Two native arc changes on its overshoot
source give the all-length count2673·2^k+645, exactly two below the desired
2673·2^k+647. The marked indegree-one root survives. The local search and
full response are recorded in [NATIVE_ZERO.md](history/NATIVE_ZERO.md).

This yields a universal **repeatable strong zero-digit operation**. If a
root and its sole predecessor each have indegree one, two births through
the predecessor preserve that interface and strongness while mapping

    (h,F,s) → (3h,3h+2s,h).

It retains all old arcs. In debt coordinates b=F−h, this map is
`(h,b,s)→(3h,2s,h)`, while the ordinary two births give
`(3h+2b,h−s,h+b)`. The latter matrix has determinant1 and inverse
`h=h'−2s', b=3s'−h', s=h'−2s'−b'`; the zero-digit map instead forgets
the old debt b. Thus there is an exact two-operation arithmetic machine
on actual strong sources, with its retained and lost information visible.
The universal zero digit is proved. The desired nonzero +2 repair remains
open, including after the narrowly specified local tests.

The final incoming checkpoint supplies a sharper stopping condition:
a one-vertex increase by2 occurs only at1→3 or3→5. Thus replacing the
5991 source cannot make a single newborn reach5993. It also supplies a
complete finite decision procedure for each fixed unit-offset ray
multiplier, with the unbounded odd-slope existence question still open.
[Incoming proofs and their precise interfaces](INCOMING.md).

## 6. The remaining proof obligation

The arbitrary-prefix constructions reach a zero defect by appending
digits. They preserve the supplied prefix but change the original target.
Once a zero is reached, the reroot preserves that target and improves its
future response. Joining those two stages at an **arbitrary prescribed
original count** remains unresolved.

A completion must select an actual smaller source and an operation
returning exactly the requested odd N, with compatible ports and a
well-founded decrease. The new finite clocks, integer debts, source
resets and native repairs make several pieces of that requirement
constructive. They do not yet supply the selection law for every N.
