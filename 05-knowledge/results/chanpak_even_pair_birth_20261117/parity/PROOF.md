# Even-level birth, odd recovery, and a parity-complete source profile

**PROVED all-order parity reductions and smooth-chain classification.** No assertion that all odd counts occur. Counts refer to ordinary tournament Hamiltonian paths, not automorphism orbits.

## 1. Recovering true first order from an even-level atlas

For an attained h let E(h) and O(h) be its minimum positive even and odd tournament orders. Adding a singleton SCC preserves H, so

    |E(h)-O(h)|=1.

If the first even level is E=2k, the true first order is exactly one of E-1 and E. It is E-1 **iff some H=h tournament at order E has a singleton SCC**. Deleting that singleton proves one direction; padding an odd-order E-1 witness proves the other. Otherwise every order-E witness has singleton-free core, and true birth is E.

The quantifier 'some' is indispensable. A chosen even-level witness with no singleton does not certify even first birth. For h=15, the strong interval chain IC(2,3) has order6 and has no singleton. But a regular five-vertex tournament already has H=15, so E(15)=6 while the true birth is5. There are also order-six witnesses obtained by padding those order-five sources. A marked even-level atlas that stores only one representative can erase the decisive information.

If h is genuinely even-first at E, its first odd realization is E+1, but every odd recovery containing an even-first witness as an induced subtournament is obtained by singleton insertion. The same-count insertion theorem in `birth-atlas-now/PROOF.md` establishes this precisely. In particular a **strong** odd recovery cannot be an induced one-vertex extension of that even-first source; it must reorganize the core.

Strong odd recovery may not exist: H=5 first appears at order4, while Busch's floor at5 is9 and increases thereafter. In contrast H=33 first appears at order6 and has a strong odd recovery at7 (the extreme-back-edge tournament with2^5+1 paths). This contrast is structural, not a padding ambiguity.

## 2. A two-state first-birth recursion

The correct source state before padding is a parity pair, not a single minimum. For d>1 let s_e(d) be the minimum strong order of parity e realizing d, or infinity. Let b_e(h) be the analogous singleton-free minimum allowing an ordered sequence of nontrivial SCCs. Set b_0(1)=0,b_1(1)=infinity for the empty core. Then

    b_e(h)=min_(d|h,d>1, eta in{0,1})
              [s_eta(d)+b_(e-eta mod2)(h/d)].

This is min-plus divisor convolution with a Z/2 parity state. Only strong divisors are searched, through their exact Busch caps. For h>1 the full minima are

    E(h)=min(b_0(h),b_1(h)+1),
    O(h)=min(b_1(h),b_0(h)+1).

For h=1 separately E=2,O=1. This recovers the true first order and states exactly what an even-level atlas lacks. The parity-specific strong minima cannot be replaced by the global strong minimum: a later opposite-parity strong class can change b_e even though it does not change the ungraded birth.

Multiplicity is obtained by the same first-SCC recursion, summing products over all minimizing divisor/parity choices. Arbitrary binary parenthesizations would overcount. Converse reverses and converses the ordered SCC word, as in `birth_task/factor/PROOF.md`.

## 3. A coefficientwise obstruction to even-order symmetry

Every self-converse even-order tournament on n vertices satisfies

    P_z(T)=z^n mod4, coefficient by coefficient in Z[z].     (1)

In particular H(T)=1 mod4. Hence if h=3 mod4, every even-order realization is non-self-converse; its isomorphism-class fiber is a disjoint union of distinct converse pairs. This applies in particular to genuine even-first fibers and to padded even recovery of an odd-first value.

Proof of(1). A tournament automorphism group has odd order: a nonidentity involution would swap a pair of vertices and reverse their arc, which is impossible; Cauchy's theorem excludes an even group order. An anti-automorphism has order2m with m odd, because its square is an automorphism of odd order. Its mth power is an involutive anti-automorphism. Such an involution has at most one fixed vertex, and at even order none.

Use the inherited odd-cycle expansion of P_z: a disjoint collection of nontrivial odd directed cycles has weight2^k z^c, where k counts those cycles and c counts all cycles and remaining singleton vertices. The anti-involution acts on these collections. A non-fixed collection with k>=1 is paired with another equal monomial, giving a multiple of4. A fixed collection cannot contain an individually invariant odd cycle: its odd vertex set would contain a fixed vertex of the involution. Thus all nontrivial cycles occur in pairs, k is even and at least2, and their contribution is again divisible by4. Only the empty collection survives, contributing z^n.

At odd order the involution has one fixed vertex v. Modulo4, the only extra fixed collections consist of a single invariant odd cycle through v. Thus a more precise odd profile is

    P_z(T)=z^n+2 sum_(C invariant odd directed cycle through v)
                       z^(n-|C|+1) mod4.                  (2)

This explains exactly which central-cycle coordinate becomes available when moving from even to odd symmetry. It is not a construction of an odd recovery at the same H.

## 4. The entire 3/5-smooth interval-chain grammar has one parity

Moon's interval chain with parts p_i>=2 has

    H=product_i (2^(p_i-1)+1), order=1+sum_i p_i.

Suppose H=3^a5^b, a,b>=0, a+b>0. Then **every** part belongs to{2,3,4}. Here is an elementary proof, avoiding any assumption about Catalan-type conjectures. If2^r+1 is3/5-smooth and r is odd, it is divisible by3 and not5, so equals3^u. If u is odd, v₂(3^u-1)=1, forcing(r,u)=(1,1). If u is even, the two factors3^(u/2)-1 and3^(u/2)+1 must be powers of2 differing by2; they are2,4, giving(r,u)=(3,2). If r is even, the value is not divisible by3 and therefore equals5^v. For v odd, v₂(5^v-1)=2, giving(r,v)=(2,1); for v even the same adjacent-power argument would require5^(v/2)=3, impossible.

Let c be the number of parts4. The counts of parts2,3,4 are respectively

    a-2c, b, c,  0<=c<=floor(a/2).

All chains therefore have exactly

    n=2a+3b+1,             n mod2=b+1 mod2.                 (3)

Thus no amount of rearranging or splitting the smooth chain atoms can recover the opposite order parity at the same count. Opposite-parity recovery requires leaving the entire interval-chain grammar, not simply adding another chain word. This is a scoped impossibility, not an obstruction to other strong tournaments (h=15 already supplies a counterexample to that stronger claim).

Moon-chain recognition gives distinct ordered compositions distinct classes, and converse reverses the composition (the triangle boundary is individually harmless). Hence the exact total chain-class count is

    sum_c (a+b-c)! / ((a-2c)! b! c!).                       (4)

For each c put(A,B,C)=(a-2c,b,c). Its self-converse contribution is zero unless at most one of A,B,C is odd; otherwise it is

    floor((A+B+C)/2)! /
       (floor(A/2)! floor(B/2)! floor(C/2)!).               (5)

Summing(5) gives the self-converse count; the remainder forms converse pairs. In particular a and b both odd forbid every self-converse chain. By(3), these are exactly the smooth chains with even order and H=3 mod4, matching the general obstruction(1).

For H=45(a=2,b=1), the five chains all have order8, with one self-converse class and two converse pairs. For H=15(a=b=1), the two chains have order6 and form one converse pair, while the genuine order-five birth has two self-converse classes. The even-level symmetry can therefore differ completely from the true first-birth symmetry.

## Provenance and replay

Moon's chain product, recognition and converse law are inherited from THM-4744, which explicitly credits Moon(1972) and related extremal results. The old smooth-number/Catalan discussions suggest the factorization argument, but the elementary classification above is supplied in full. SCC/parity minimum formulas refine the earlier first-birth and plateau-core packages. The modulo4 consequence H=1 for even self-converse tournaments is already recorded in THM-791; the proof here retains the entire slot polynomial and its odd-center counterpart.

The replay checks literal small smooth chains and all distinct word orders, their H and self-converse counts, coefficientwise congruence for self-converse tournaments through order5, the h15 even-level trap, and the h33 strong odd recovery. These controls are finite supplements to the all-order proofs.
