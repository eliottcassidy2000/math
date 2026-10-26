# Source-fixed shape switching: the exact refill defect

2026-10-26. **PROVED elementary statements:** exact center-switch cocycle,
commutator/primitive-word classification, source pullback, and an unbounded
ordinary-integer refill family. **FINITE-EXACT:** a seven-bit refill on the
actual marked orbit, complete bounded block controls and literal replay.
**OPEN:** an arbitrary-shape global rank on the root-27 orbit and all-prefix
nonhalting. No claim of native-machine formalization, Collatz convergence,
or novelty for affine conjugation or word commutation is made.

## Result first

There is an exact missing term in a sum of fixed-shape budgets. At macro
index 86 of the marked source `z=7,q=1` (the same guarded process whose
first odd-`x` state is `(a,x)=(6,27)`), the legal shape block `P=100`
consumes three bits of its periodic fuel, taking it from 5 to 2. Changing
the reference shape to `Q=101000` then refills that fuel to **9**. Every
one of these macros is legal on the actual source; this is not a random
cylinder or an unrelated root. The extra seven bits are an exact equal-
valuation cancellation against the shape-switch defect.

The defect has a closed formula. Its zero locus consists precisely of
blocks that are powers of one common primitive word; this is the one case
where fixed-shape fuel glues without a reset. A fixed pair of different
centers also admits arbitrarily large refills at positive integer sources.
Thus `periodic fuel + a shape-dependent constant` cannot be a universal
nonincreasing rank for legal switching across all sources. That obstruction
is deliberately weaker than nonexistence of a richer root-27-specific rank.

## 1. Inheritance, exact source and guards

The pinned inherited anchor [source/04-README.md](source/04-README.md),
sections 1–5, gives

    T(z)=(3z+c)/2, c=z mod 2,       q'=q+3c-1.

A native macro requires `q>=1`, and is blocked precisely when `q=1,c=0`.
An ordered full ternary tree has grammar `0 | 1 Tree Tree Tree`. Deleting
the final leaf of any nontrivial tree gives a nonempty **balanced block**
P: every prefix has credit change at least zero and the total change is
zero. Its length L is a positive multiple of three, and its first bit is
1. Any concatenation of such blocks is balanced. Appending a final 0
forms a tree, but that last 0 is the blocked symbol, not an executed macro.

For a word w of length n define

    C_w = sum_(i<n) c_i 2^i 3^(n-1-i),
    F_w(z) = (3^n z+C_w)/2^n.

The inherited cylinder theorem is exact:

    w is the first n bits of z  iff  3^n z+C_w = 0 mod 2^n.
    C_(uv)=3^|v| C_u+2^|u| C_v.

Hence all intermediate parities and guards must be replayed; a rational
formula without its cylinder is not an executed macro. A balanced block
whose cylinder holds is legal at every initial `q>=1` and leaves q unchanged.
This does not imply the entire root-7 itinerary decomposes into a fixed
bank of balanced blocks.

The physical initial segment is

    (z,q)=(7,1)->(11,3)->(17,5)->(26,7)->(39,6).

At the third cut `x=z+1=27,a=q-1=6`. The odd-`z` and odd-`x` sections are
different, so their indices and root functions are not identified. The
actual first 110 macros used below are all legal by literal replay.

Closest proved mechanism: inherited first-mismatching-leaf theorem and
fixed-shape fuel. Canonical hostile: its nineteenth root-27 reset, where
`a+v_3(x-1)` falls from 25 to 22. Corrected near miss: the 344-shape bank
proves infinitely many repetitions of each listed shape, not arbitrary
shape changes. Least-used sidecar: the **signed, ordered affine defect**,
not just its valuation. Board: original source / chronology / tree shape /
periodic center / cancellation defect / clock-preserving merging.

## 2. Exact switch cocycle

For a balanced block P write

    L_P=|P|, A_P=3^L_P, B_P=2^L_P, D_P=A_P-B_P,
    N_P(z)=D_P z+C_P,
    r_P=-C_P/D_P,          h_P(z)=v_2(N_P(z)).

The integers D_P and C_P are odd and positive. At positive integer z,
N_P(z)>0, so h_P is finite. Since D_P is odd, h_P(z) is exactly the
2-adic distance exponent `v_2(z-r_P)`. Moreover P is realized at z iff
`h_P(z)>=L_P`. Keeping the same block consumes precisely its length:

    h_P(F_P(z))=h_P(z)-L_P,                         (1)

whenever P is legal.

For another balanced block Q define the **ordered center defect**

    Delta_(P,Q)=C_Q D_P-C_P D_Q.

**Theorem 1 (switch cocycle).** If P is legal at positive integer z and
`y=F_P(z)`, then

    B_P D_P N_Q(y)
       = A_P D_Q N_P(z) + B_P Delta_(P,Q).         (2)

If Delta is nonzero, put `t=v_2(Delta)`, `u=h_P(z)-L_P`. Then

    u != t  => h_Q(y)=min(u,t),
    u =  t  => h_Q(y)>=t+1.                        (3)

In the second line the surplus is the valuation of the sum of the two
normalized odd terms; it is not bounded merely by t. If Delta=0, then
`h_Q(y)=u` exactly.

**Proof.** Substitute `y=(A_P z+C_P)/B_P`, multiply through by `B_P D_P`,
and collect the coefficient of z. Equation (2) follows directly.
The two summands on its right have valuations `h_P(z)` and `L_P+t`.
Unequal valuations give their minimum; equal valuations cancel at least
once. The left side has valuation `L_P+h_Q(y)` because D_P is odd. For
Delta=0 there is one nonzero term. This proves all branches. The integer
guard was used to assert y is an actual state, not merely a rational value.

A useful exact ledger is therefore

    h_Q(y)=h_P(z)-L_P + switch_gain(z;P,Q),
    switch_gain = h_Q(y)-[h_P(z)-L_P].

The gain can be negative or positive. Replacing it with a nonnegative
shape count, dropping its sign, or treating it as independent of z loses
the cancellation branch that matters.

## 3. Zero defect, order and exact clone limits

The same integer is a literal failure of affine commutation:

    F_Q(F_P(z))-F_P(F_Q(z))
        = -Delta_(P,Q)/(B_P B_Q).                  (4)

It is independent of z, but its effect on a valuation depends on the
source-dependent normalized odd part in (2).

**Theorem 2.** For nonempty balanced P,Q, the following are equivalent:

1. Delta_(P,Q)=0;
2. r_P=r_Q;
3. their affine maps commute;
4. PQ=QP as literal bit words;
5. P and Q are positive powers of one common primitive bit word.

**Proof.** The first three equivalences are the definition and (4).
Commutation says `C_(PQ)=C_(QP)` at the same total length. At a fixed
length the carry map is injective: its residue `-C_w 3^(-|w|) mod2^|w|`
uniquely determines the itinerary word by the cylinder theorem. Hence
commutation gives PQ=QP; the converse is immediate. The common-word
claim follows by cancellation: if `|P|<=|Q|`, PQ=QP forces Q to begin
with P; write Q=PR, cancel P and obtain PR=RP. Induct on total length
until a common primitive remains. Conversely powers commute. This is a
finite word argument, not an assumed infinite periodic orbit theorem.

This is a sharply scoped merger rule. It identifies a shared periodic
center, **not** the same executed transition. For example P=100 and PP
share r=-9/19, but their positive-source maps differ and have durations
3 and 6. A clocked affine interface can substitute P for Q on every
source only if both lengths and carries agree, hence P=Q literally.
Different lengths have different slopes `(3/2)^L`; at the same length
carry injection recovers the word. Ordered tree shape is therefore not
a disposable clone coordinate. A broader contextual merger needs an
explicit weaker observable and its lost-duration/source sidecar.

## 4. Pull every comparison back to the original source

Let s be the fixed positive integer source and let w be an actually
executed prefix of length n, with state `z_n=F_w(s)`. Pull the periodic
center for a prospective block P back through the prefix:

    rho_(w,P)=(2^n r_P-C_w)/3^n.

Its denominator is odd. Exact affine substitution gives

    h_P(z_n)+n = v_2(s-rho_(w,P)).                 (5)

This is a genuinely source-preserving representation: s never changes.
If P is executed and retained as the reference shape, then

    rho_(wP,P)=rho_(w,P).

The total precision on the right of (5) stays fixed while n grows by
L_P and local fuel drops by L_P. When the shape is changed at the same
prefix w,

    v_2(rho_(w,Q)-rho_(w,P))
          = n+v_2(Delta_(P,Q))                    (6)

for nonzero Delta. Thus the refill is exactly a three-point 2-adic
cancellation at the original source. It cannot be removed by pretending
the source is independently restarted at every tree.

Equations (5)–(6) extend to every finite legal prefix and every balanced
reference block without a fixed-shape restriction. They give the needed
source-preserving **carrier**, not a well-founded global rank: the pulled
centers change under switching, and arbitrarily much new source precision
can appear at an equal-depth cancellation.

## 5. A refill on the actual root-27 process

Take P=100 and Q=101000. Their exact constants are

    (L_P,C_P,D_P)=(3,9,19),
    (L_Q,C_Q,D_Q)=(6,351,665),
    Delta_(P,Q)=684,        v_2(Delta)=2.

Literal iteration from `(z,q)=(7,1)` reaches at macro time 86

    z=10952205219944621, q=47.

The next bits begin `1001010001011101000`; hence P followed by Q is
legal with q returning to 47 after each block. After P,

    y=36963692617313097,
    h_P(z)=5,       h_P(y)=2,       h_Q(y)=9.

The two terms in (2) have equal valuation 5; their exact sum has
valuation 12, giving h_Q(y)=12-3=9. Switching contributes **+7** beyond
the remaining P-fuel. This defeats naive monotonicity even at the actual
source, without claiming that every possible shape-only weighting is
refuted by this one witness.

The checker separately reconstructs the inherited old-section hostile
`x=57971215,a=21,k=2,v_3(x-1)=4` at odd-return index 19, where charge
falls from 25 to 22. It verifies the repaired Antihydra identity

    a_i = 6 + sum_(l=1..i)(2v_3(o_l)-1),

never including o_0=27. The new cocycle does not repair that old charge;
it identifies another exact refinement loss that any proposed repair
must also pay.

## 6. Unbounded refill at a fixed pair of shapes

**Theorem 3.** With the fixed P,Q above, for every integer m>=1 there
is a positive integer source z at initial q=1 such that P is legal,
`y=F_P(z)`, and

    h_P(z)=5,              h_Q(y)=6m.              (7)

The process then realizes `P Q^m 0` as a tree certificate, whose last
0 is blocked. It executes exactly `3+6m` native macros before that halt.

**Proof.** The word `P Q^m 0` is a full tree: its blocks have nonnegative
prefix credit and zero total credit, followed by one final leaf.
Its binary cylinder is a nonempty residue class modulo `2^(4+6m)`.
Choose its least residue plus that modulus, ensuring a positive integer
source in the physical arithmetic sector. At y the word `Q^m0` is
realized. The inherited pumped criterion, or a direct geometric-sum
calculation, gives `h_Q(y)=mL_Q=6m`. Since this exceeds t=2, (2) forces
`h_P(z)-L_P=t`, proving h_P(z)=5. All proper tree-prefix q values stay
positive. The last leaf is the first blocked symbol, not another native
macro. The macro count is the word length minus one.

For clarity, the geometric-sum criterion used here is exact:
with A=3^L, B=2^L, C=C_Q and D=A-B,

    3^|Q^m0| y+C_(Q^m0)
       = 3/D * (A^m N_Q(y)-C B^m).

This is divisible by `2^(Lm+1)` iff `v_2(N_Q(y))=Lm`; unequal valuations
cannot reach the required modulus, while equal valuations yield one
further cancellation of odd parts.

**Corollary (limited rank obstruction).** No finite constants w(P),w(Q)
make `R(z,P)=h_P(z)+w(P)` nonincreasing on *every* legal transition
`(z,P)->(F_P(z),Q)` over positive integer sources. Indeed (7) would
require `6m+w(Q)<=5+w(P)` for all m. This also obstructs any uniformly
bounded shape-only correction to h.

This is not a root-27 counterexample to survival. The family varies the
ordinary integer source and initializes q=1; it is not asserted to lie
on the marked source with its inherited credit. Nor does it rule out a
rank with unbounded source/history-dependent coordinates, a different
nonlinear statistic, or a rank tailored solely to the marked ray.

## 7. What remains, and where defect/clone machinery can really attach

The inherited exact all-prefix target remains

    sum_(j<J) ell_j <= 3J,  for every J on the odd-z ray from 7,
    ell_j=v_2(3z_j+1).

Equivalently every finite tree word w fails the source-7 cylinder:

    v_2(7*3^|w|+C_w) < |w|.

The present construction handles **arbitrary ordered shape switches**
without changing that source, but supplies no bound on their cancellation
surplus. For a bank-based strategy the exact next obligation is to prove
an amortized upper bound on the accumulated switch gains along the one
actual legal trajectory, with coverage of that trajectory and its native
guards established separately. An infinite bank chosen after observing
the whole unknown halting tree simply restates the original target.

The affine defect (4) is not the Frobenius orientation defect from the
supersingular-j/tournament result. One is a signed rational obstruction
to commuting ordered affine words; the other counts unoriented-pair
orbits whose reversal forces orientation failures. No arithmetic map
identifying them is supplied here, so no moonshine/genus theorem is
imported into root-27 survival by shared terminology. The valid reusable
lesson is to preserve the native operation and optimize only the declared
observable. Likewise the exact merger test above is a concrete source-
and-clock interface, not a claim that abstract tournament clone mergers
act on the native macro process.

The sibling modular-carry lane supplies a signed-port permutation and
orientation-defect specialization, not an identified arithmetic Frobenius
action, and should be cited independently rather than conflated with (4).

## 8. Closed executable universe and evidence

Run from this directory:

    python3 -B check.py > normal.json
    python3 -O -B check.py > optimized.json
    cmp normal.json optimized.json

No load-bearing assertion disappears under optimization. The checker uses
only exact integers/Fractions and imports no earlier research code.

- All 1,023 binary words of lengths 0..9: independent direct and recursive
  carries, cylinder residues, and same-length injection.
- All 71 nonempty balanced tree-prefix blocks from trees with 1..4 internal
  vertices; all 5,041 ordered pairs, including 87 common-center pairs.
- 10,082 exact source-lift cocycle tests with literal guarded block replay.
- The source-7 first 110 macros, 111 source prefixes and all 71 blocks:
  7,881 original-source pullback identities, including same-center transport.
- The actual macro-86 seven-bit refill, plus the inherited odd-return-19
  charge failure and all 20 corrected returned-state v3 counter checks.
- Nine exact controls of the all-m refill construction, at
  m=1,2,3,4,5,8,16,32,64, with 837 legally executed macros in total;
  each of the nine final blocked leaves is read without executing its update,
  and attempted execution is separately required to raise the native guard.
- Positive/hostile merger control: P and P^2 have the same periodic center
  but different exact duration and affine action.

The paired receipts are byte-identical PASS. The unbounded claims are
proved above, not inferred from these caps. Source pins retain the original
63f52936 snapshot and the repaired THM-4771 at incoming2807ee; the incoming
commit does not change the inherited root-27 anchor or ancestry files.
