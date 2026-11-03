# Arithmetic clocks and exactly observable boundary data

2026-11-03. **PROVED** elementary statements CB1–CB6 below, with the explicit
S168 seed and gate lift inherited as finite-certified inputs. **FINITE-EXACT**
small tables have two counting/linear-algebra implementations. No external
priority claim, autonomous unmarked tournament dynamics, or spectrum
completeness theorem is asserted.

## Inherited object and operation scope

The current source is [S168 cellular, Theorems 1–6](../session168_computational_structures_20261103/cellular/NOTE.md),
merged at `34815c52f7`. It defines six labelled six-vertex tournaments

    B_h, h=0,...,5: masks 153,8345,8337,8369,177,185.

Bits are in lexicographic pair order; bit one means smaller label points to
larger label. These form a literal arc-flip hexagon. Their classes depend
exactly on h mod3; the involution `(0 4)(1 2)(3 5)` carries B_h to B_(h+3).
Every seed is strong, with H=29 and labelled-path-slot polynomial

    f(r)=r^6+12r^4+16r^2.

For a ring of m sites there are two tracks a,b over F3. The actual rule is

    R(a,b)=(b,a+b+Sb),       (Sb)_i=b_(i+1),
    R^(-1)(c,d)=(d-c-Sc,c).

Use the four inherited layers: a+=b+Sb; a+=b; b=a-b; a=a-b.
Each controlled single-trit change lifts from a physical h in Z/6Z by
adding the unique representative of target-source in {-1,0,1} modulo6.
This changes at most one literal arc. The physical encoding is the ordinal
sum of 2m blocks B_h, so has order12m, H=29^(2m), and full profile f^(2m).
All intermediate gate states retain that profile. The block labels, spatial
neighbor relation, and layer schedule remain specified data.

S168 proves logical and physical **round** order8m at m=3^a. We extend its
state classification and boundary analysis; we do not rebrand that theorem
as a new result. The earlier S167 single-track I+S clock has order2m and is
a different update rule.

## CB1. Every source's exact period, including its physical lifts

Let m=3^a and identify the shift ring with

    A=F3[z]/(z^m-1)=F3[u]/u^m,   u=z-1,
    M(z)=[[0,1],[1,1+z]].

For a nonzero v=(a,b) in A^2 put

    nu=min(ord_u(a),ord_u(b)),    d=m-nu.

Then every physical lift of v has the same least round period as v, namely

    p(v)=8*3^ceil(log_3 d).                              (1)

The zero logical state has period1; all 4^m physical states over it are fixed.
For each d=1,...,m, exactly 8*9^(d-1) logical states have depth d. Hence the
number of logical cycles at this exact depth is

    3^(2d-2-ceil(log_3 d)).                              (2)

There are 4^m times as many physical cycles of each depth and period.

**Proof.** The inherited calculation is

    M(1)^4=-I,    M(1) has order8,
    M(z)^8=I+uB(u),    B(0)=[[1,2],[2,2]], det B(0)=1.

M(1) has irreducible characteristic polynomial T^2+T-1 over F3. Its
identification with multiplication by a primitive element of F9 shows every
nonzero vector has period8. The first nonzero u-coefficient of v therefore
forces any period to be a multiple of8.

B is invertible over A. If t=3^j c with 3 not dividing c, then
`(I+uB)^t-I` is u^(3^j) times an invertible matrix: the remaining constant
matrix is c B(0)^(3^j). Thus it kills v exactly when 3^j>=d. This proves (1)
for logical states, including minimality, not just a return bound.

Write q=3^ceil(log_3 d). From M^4=-I+uD and characteristic three,

    M^(4q)v=-v,

since u^q annihilates v. Each of the four intermediate gate maps is linear.
The second half of this state's p(v)-round gate itinerary is therefore the
negative of the first half. Negation reverses every signed unit increment
in its physical lift. All increments cancel cell by cell, so every lift
returns after p(v) rounds. Projection excludes an earlier physical return.
For v=0 every gate has zero increment, including on nonzero deck bits.

The number of pairs divisible by u^(m-d) is9^d; subtract9^(d-1). Divide by
the period for (2). Over each logical cycle, 4^m physical starts above a
fixed logical point belong to different same-period cycles. QED.

For m=3 there is one logical fixed point, one 8-cycle and thirty 24-cycles.
The physical census is 64 fixed points, 64 8-cycles and1920 24-cycles.
The all-state physical census covers all6^6=46,656 labelled configurations.
A spatially constant nonzero source has depth1 and period8, even on a large
ring. The impulse has depth m and period8m. Maximum period is not the
period of every nonzero source.

## CB2. One fixed vertex attachment reads all six literal states

To B_h adjoin a vertex w with the fixed boundary rule

    1 -> w,        w -> {0,2,3,4,5}.

For h=0,1,2,3,4,5 the resulting Hamilton-path counts are respectively

    (109,115,101,77,73,81).                              (3)

They are all different. All six extensions are strong: the old block is
strong, w is reached from1, and w reaches another old vertex.

**Finite proof and formula.** For an arbitrary tournament T with marked q,
let s_q count Hamilton paths starting at q. Let K_q count ordered covers
by two paths, the first nonempty and ending at q, the second possibly empty.
The same one-vertex attachment has exactly

    H(T)-s_q+K_q                                         (4)

Hamilton paths. If w is first, the old path must not start at q, giving
H-s_q. Otherwise the prefix must end at q; cutting at w gives exactly the
K_q covers. For the six seeds, K_1=(83,89,75,52,49,56) and
s_1=(3,3,3,4,5,4). With H=29, (4) proves (3). `check.py` counts subset paths
and covers; `audit.py` separately examines all7! orders of each extension.

Among the64 one-vertex orientation cuts, exactly16 give six distinct
counts, a FINITE-EXACT bank in `results.json`. Cut2 above is one such choice.
Within the model of unchanged old arcs and scalar H readout, zero added
vertices cannot separate the seeds; one added vertex is therefore minimal.
This is not minimality among all allowed probes or all tournament gadgets.

For a chosen block j in an ordinal sum of L blocks, make every earlier
block point to w and w point to every later block. The result is the same
ordinal sum with B_h replaced by its seven-vertex extension. Its count is

    29^(L-1) times the corresponding entry of (3).         (5)

This gives a literal boundary query at every designated cell. It uses a
fresh copy of the state, or a query followed by explicit deletion of w.
The queried tournament itself no longer has the original profile.

**Why the mark is necessary.** The deck involution swaps labels1 and2.
Moving the probe along with that involution would produce isomorphic
extensions and would not read the deck bit. Here label1 is held fixed.
This is boundary information supplied by a retained physical port.

## CB3. One boundary history is a complete logical observer, at every width

For every m>=1, set ell(a,b)=a_0. The first2m observations

    ell(v), ell(Rv), ..., ell(R^(2m-1)v)                  (6)

determine v uniquely. They are an invertible linear map from F3^(2m) to
F3^(2m); fewer than2m scalar trit observations cannot be injective on all
logical states. The same impulse v=(delta_0,0) is a cyclic vector for R.

**Proof.** As matrices on both tracks,

    M-M^(-1)=I+S,    hence S=M-M^(-1)-I.

Also ell M reads b_0. Thus ell S^j and ell M S^j recover both tracks at all
sites. These Laurent-polynomial expressions in M span the full dual.
M is invertible; by its characteristic polynomial, negative powers are
polynomials in M, and powers reduce to the first2m. This proves (6).
For the impulse, Mv=(0,delta_0), and S^j produces both spatial bases by
the same argument, proving cyclicity. The fewer-observations bound follows
from dimension or the cardinality3^(2m). QED.

Apply the boundary attachment CB2 at the a_0 block at each round. Formula
(5) decodes its physical h, and reducing h modulo3 gives (6). Consequently
these actual Hamilton-count queries recover the entire logical state.
The dimension argument alone bounds only trit-valued observations. The
following carry calculation proves the sharper physical-reader boundary.

### The deck bit is a partial look-ahead carry, with an exact census

Let y_t=a_0(t) in F3, h_t the literal Z/6Z state of that cell, and let
sigma:F3->{-1,0,1} be the signed representative. Summing the three gates
that update a_0 in one round gives

    h_(t+1)-h_t = sigma(y_(t+2)-y_t)+sigma(y_(t+1))
                   -sigma(y_(t+2))                 (mod6).       (9)

Indeed, with A=a_0,B=b_0,C=b_1, those increments are
sigma(B+C), sigma(B), and -sigma(A+B+C). The next boundary trits are
y_(t+1)=B and y_(t+2)=A+B+C, proving (9).

A physical history of length L records y_0,...,y_(L-1) and the initial
deck bit. All its increments except the last are then predetermined by
(9). The last increment tests y_L. When y_(L-2)=0 it supplies no information;
otherwise it partitions the three possible y_L into a singleton and a
pair. The two increment values differ by3 modulo6, so this is exactly the
additional deck information. It cannot recover all three values.

For 2<=L<2m, CB3 says the first L+1 trits are independent coordinates.
Summing the one, two, two possible carry outcomes for y_(L-2)=0,1,2 gives
exactly

    2 * 3^(L-1) * (1+2+2) = 10*3^(L-1)                (10)

distinct physical histories. There are6 histories at L=1 and2*9^m at
L=2m. Consequently **exactly2m consecutive round-boundary queries are
necessary and sufficient in the worst case**, even for this six-valued
physical reader, to recover the
logical state and observed deck bit. At L=2m-1 the two unresolved values
in a carry-pair give distinct logical sources with identical histories.
This is specific to the inherited four-layer lift and the fixed a_0 port,
not a lower bound for all possible controllers, probes, or adaptively
spaced query times.

### One-lap physical obstruction is a boundary winding parity

For any logical orbit of period p, at any width, put

    W_j(i)=sum_(t=0)^(p-1) sigma(y_i(t+j)-y_i(t)),  j=1,2.

These are integer multiples of3. Summing (9) telescopes the last two
terms, so the a_i deck change is W_2(i)/3 mod2. The single b_i update in
each round has increment sigma(y_i(t+2)-y_i(t+1)); its deck change is
W_1(i)/3 mod2. Thus the full one-lap deck vector is known from these two
boundary winding parities at each site, without reconstructing all the
intermediate tournament flips. The physical period is p iff every entry
is zero, and2p otherwise.

There is a still smaller parity test: W_j(i)/3 mod2 is the parity of the
number of t for which the unordered pair {y_i(t),y_i(t+j)} is {0,1}.
The closed word's three undirected edge counts have the same parity,
since each vertex has even total incidence. Their sum has that parity;
it is also W_j(i) modulo2, equal to W_j(i)/3 modulo2. This is an exact
necklace-crossing invariant for this physical lift, not a general formula
for arbitrary gate schedules.

### Exact remaining ambiguity

Under the operation alphabet consisting only of this fixed clock and this
fixed a_0 query, two physical initial configurations have identical infinite
query histories iff they have the same entire logical state and the same
physical a_0 state. The other2m-1 deck bits are invisible forever.

Indeed (6) recovers the logical state and the first query recovers physical
a_0. Conversely, changing any other deck bit by3 commutes with every lifted
gate and leaves all logical controls and a_0 unchanged. This proves both
directions. Each future-equivalence class therefore contains exactly
2^(2m-1) physical states; the number of distinguishable classes is2*9^m.

One initial CB2 query at each remaining cell completes recovery of the full
physical source. A decoder supplied with the logical state and the a_0
deck bit must still retain2m-1 bits to reconstruct arbitrary physical states.
If queries at all ports are allowed, the equivalence relation becomes
equality. If only the invariant profile is queried, every state is equivalent.
The operation alphabet, not the number alone, determines the memory cost.

## CB4. An exact logical clock prescription for every ring length

The characteristic and minimal polynomials of R on m sites coincide:

    P_m(T)=(T^2-T-1)^m-T^m.                              (7)

Write m=3^a b with 3 not dividing b. Factor P_b over F3 into distinct
monic irreducibles f_i with multiplicities e_i, and let o_i be the
multiplicative order of a root of f_i in its finite field. Then

    ord(R)=lcm_i(o_i) * 3^(a+[8 divides b]).              (8)

Every e_i is1, except that T^2+1 occurs with multiplicity2 exactly when
8 divides b. Equation (8) is a finite algorithm for each width, not a
width-independent numerical bound or a primality test.

**Proof.** Block elimination gives
`det(TI-M)=det((T^2-T-1)I-TS)=P_m(T)` for a cyclic shift S.
CB3's cyclic vector makes this polynomial minimal as well.
Characteristic three gives P_m=P_b^(3^a).

For nonzero T put r(T)=T-1-T^(-1); P_b=T^b(r(T)^b-1).
Every spatial b-th root is simple because3 does not divide b. Multiplicity
can increase only where r'(T)=1+T^(-2)=0, i.e. T^2=-1. At such a point
r=-1-T, r^2=-T, r^4=-1, so r has exact order8. The second derivative of r
is nonzero, giving multiplicity exactly2 when8 divides b. Thus the only
repeated irreducible is T^2+1 in precisely this case.

For any irreducible f, its root order o is prime to3 and f divides T^o-1
with multiplicity1. Raising to3^j multiplies the multiplicity by3^j.
The least order modulo f^e is therefore o*3^ceil(log_3 e). Taking the lcm
over all primary factors, whose maximal exponent is3^a or2*3^a, proves (8).
QED.

The same reasoning applies to a particular state's annihilator Q_v dividing
P_m: use only the primary factors and multiplicities present in Q_v. A
physical lift can then require one or two logical laps, by the inherited
deck-parity test; (8) alone does not decide that extra lap.

| m | Logical operator order | Impulse physical period, where checked |
|---:|---:|---:|
|1|8|8|
|2|8|16|
|3|24|24|
|4|40|80|
|5|328|328|
|7|2920|2920|
|8|120|240|
|9|72|72|
|11|968|1936|
|16|4920|not needed for this table|

These finite entries have polynomial and independent dense-matrix checks.
The 5-site clock includes a degree8 factor of root order164, giving
lcm(8,164)=328. At11 sites the remaining degree20 polynomial splits into
four degree5 factors, whose combined root order is242; thus the clock is
lcm(8,242)=968=8*11^2. Its impulse has a nonzero winding obstruction,
so the physical period is1936, unlike the undoubled odd width5. These are
explicit field-order and boundary-carry calculations, not
inferences from the numerals5 and11. The 8-site extra factor3 is ramification
at the critical points of r, not the physical double-cover bit.

## CB5. Compatible finite clocks need not give a finite return at the limit

The quotients F3[[u]]/(u^(3^a))^2 carry the preceding clocks compatibly.
On the inverse limit F3[[u]]^2, M(1+u) has no nonzero periodic vector.

**Proof.** If8 does not divide t, M(1)^t-I is invertible, hence so is
M(1+u)^t-I over the power-series ring. If t=8*3^j c with3 not dividing c,
the proof of CB1 writes this difference as u^(3^j) times an invertible
matrix. Multiplication by u^(3^j) is injective on F3[[u]]^2, so again no
nonzero vector is fixed. Every positive integer t is covered. QED.

The single compatible impulse has period8*3^a in every finite quotient but
no finite return in the limit. The bonding map is reduction modulo u^m;
in spatial coordinates it folds a larger ring by summing coefficients in
each congruence class, not by deleting the extra sites. Such logical
projection has not been realized as a physical port-preserving tournament
move here. This provides a precise obstruction to transporting finite
modular recurrence into an integer termination assertion.

## CB6. The integer schedule clock remains separate

For any fixed sequence of L invertible gates, adjoin a program counter
j mod L and apply the corresponding gate at each tick. If its macrostate
has least period p under the whole sequence, the state with counter has
least tick period exactly Lp. A return of the counter requires a multiple
of L ticks; those returns apply precisely the whole-round map.

Our sequential schedule has L=4m gate attempts, including idle attempts.
Thus a source of CB1 period p has program-counter period4mp. In particular
the impulse has8m rounds, but32m^2 scheduled ticks. The number of actual
nontrivial arc flips is a different, source-dependent clock bounded by that
tick count. Parallel layers also give a different cost convention.
This does not correct S168's period theorem: that theorem explicitly counts
rounds. It prevents a later consumer from silently equating the clocks.

## Consequence for the research targets

The complete path profile is constant on arbitrarily long literal periodic
computations, while an explicitly marked one-vertex context separates their
local states. Consequently a descent certificate depending only on that
profile cannot be strictly decreasing on these moves. A useful completion
must retain its source, permitted probe/controller, and boundary response.
This is a concrete information obstruction and a working observer, not a
proof of Collatz convergence or of all odd H values except7 and21.
