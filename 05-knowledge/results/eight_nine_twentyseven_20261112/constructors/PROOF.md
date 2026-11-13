# Actual escapes from the asymmetric eight-vertex divisor-27 source

**PROVED:** the literal module-response identities, a universal two-arc
modulo8 selector with sharp relative bounds, their marked extensions,
the sharp count-preserving reduction27→3, and the
all-length primitive rooted return below. Source statistics are finite exact
certificates; unbounded conclusions use the stated substitution and birth
bijections. **OPEN:** target-directed coverage and Chan–Pak. No new claim
that the displayed individual Hamilton counts were previously unrealized.

## 1. The inherited source and two physical response directions

Use the source in
[Session293's exact context-ideal theorem](../../session293_cycle_returns_20261113/control/MODULE_IDEAL.md),
with zero-based outgoing rows

    C=(56,117,129,22,100,76,141,59).

It is asymmetric, H(C)=513, its ordered-cover gcd is9, and its exact
uniform-context gcd is27. We accept the inherited asymmetry assertion;
no automorphism audit is repeated. Let C_k(T) count ordered spanning
nonempty k-path covers and define the slot polynomial

    F_T(x)=sum_(k≥1) C_k(T) binom(x,k).

For integer x≥1 this counts spanning directed-path covers in x labelled
slots, allowing empty slots. In particular F_T(1)=H(T), and
F_(TT_d)(x)=x^d. The exact source polynomial is

    F_C(x)=x^8+34x^6+232x^4+246x².                    (1)

Now make either of two **single-arc changes**:

    C2=(56,53,129,22,100,76,143,59),  reverse1→6;
    C4=(56,117,137,18,100,76,141,59), reverse3→2.

Both remain strong. Their full polynomial responses are

    F_C2−F_C=2x^4=2F_(TT4),
    F_C4−F_C=4x²=4F_(TT2).                           (2)

Thus H(C2)=515 and H(C4)=517. Their exact uniform-context divisors are5
and1, respectively. The +2 move removes the3 obstruction while exposing
a5 obstruction; the +4 move is primitive already for the whole uniform
context language. Neither claim is inferred merely from its current H.

The complete finite ordered-cover certificates are:

| source | (C_1,...,C_8) | raw gcd | uniform gcd |
|---|---|---|---|
| C | (513,6102,32508,99432,187200,216000,141120,40320) | 9 | 27 |
| C2 | (515,6130,32580,99480,187200,216000,141120,40320) | 5 | 5 |
| C4 | (517,6110,32508,99432,187200,216000,141120,40320) | 1 | 1 |

The uniform gcd is `gcd_k (k−1)! C_k`, by the inherited exact theorem.
Converting these complete lists in the binomial basis proves (1)–(2),
not just equality at finitely many arbitrarily chosen colours.

## 2. Universal positive responses, including exterior marked states

Let E be **any fixed uniform tournament exterior**: every exterior vertex
sees every vertex of the substituted module in the same direction. Keep
all exterior arcs and these polarities fixed. Then (2) gives the actual
all-order identities

    H(E[C2])−H(E[C])=2H(E[TT4]),
    H(E[C4])−H(E[C])=4H(E[TT2]).                     (3)

Proof: compress every maximal module run in a spanning path to a slot.
The resulting exterior skeleton coefficients a_k(E) do not depend on the
module. The count is `sum_k a_k(E) C_k(module)`. Polynomial equality (2)
is equivalent to the same coefficientwise equality of the complete C_k
vectors, proving (3). Every graph on either side is a literal tournament.
No signed or fractional graph is required to realize the left-hand move.

The same proof applies to ordered path-cover observables with **all marked
vertices exterior to the module**. In particular, for any retained exterior
root, each coordinate of `(h,s,q,F)` obeys (3): h is H, s counts Hamilton
paths starting at the root, q is H after its deletion, and F counts ordered
spanning covers (P,Q), where P is nonempty and ends at the root and Q may
be empty. For q simply
delete the exterior root first. For the other coordinates compress runs
inside the marked covers. Likewise the H,U,W interface of the offset-tail
constructors is transported when its marks are exterior.

This gives exact positive response vectors through every continuation that
keeps the old source module uniform. It is an actual addition-like
operation, but the added quantity is the corresponding **same exterior**
applied to TT4 or TT2. It is not an independently chosen integer.

There is also a sharp unit-gap boundary. The first difference in (3) is2
exactly when E[TT4] is transitive. In a uniform context that means a
transitive prefix, followed by the module, followed by a transitive suffix.
If the exterior is nonempty, the resulting E[C] is not strong. Therefore
a strictly larger strong uniform wrapper of this pair has gap at least6.
The second pair always has gap at least4, and at least12 in a larger strong
uniform wrapper. These statements use the inherited fact that a tournament
has exactly one Hamilton path iff it is transitive. They identify the
precise renewal boundary of these two local repairs.

### A two-arc selector for every odd residue modulo8

Apply **both** disjoint arc changes, giving the strong tournament

    C24=(56,53,137,18,100,76,143,59),     H(C24)=495.

Its full polynomial is `x^8+34x^6+226x^4+234x²`, its ordered covers are
`(495,5994,32292,99288,187200,216000,141120,40320)`, and its raw and
uniform contents both equal9. The mixed response is exact:

    F_C24−F_C2−F_C4+F_C=−8(x^4+2x²)
                        =−8F_(C3→TT1).             (3a)

Thus the switches are coupled, but their mixed carry is divisible by8 in
every uniform context. Set

    h=H(E[C]),  A=H(E[TT4]),  B=H(E[TT2]).

All three numbers are positive odd integers. For switch bits a,b∈{0,1},
the four literal outputs have counts

    h(a,b)=h+2aA+4bB−8ab(A+2B).                     (3b)

Their values are therefore

    h,   h+2A,   h+4B,   h−6A−12B.                 (3c)

Since A and B are odd, these are **all four odd residues modulo8, exactly
once**, for every fixed exterior. The graph order, every exterior arc, and
every exterior/module incidence remain unchanged. Only the two named
internal arcs are chosen. The three-coordinate response `(h,A,B)` retains
the entire carry, and also gives the exact weighted balance

    H(E[C24])+3H(E[C2])+3H(E[C4])=7H(E[C]).          (3d)

No inference of an unconstrained arithmetic operation from the numeral7
is needed: this identity is a consequence of the literal two-arc response.
The same linear relations hold coordinatewise for exterior marked
path-cover observations; the modulo8 conclusion specifically uses the
oddness of the scalar Hamilton counts.

**Sharp relative cost.** Every one of the four chosen counts lies in

    (495/513)h ≤ h(a,b) ≤ (517/513)h.                (3e)

Both constants are attained by the empty exterior, so they cannot be
improved uniformly. To prove the bound put L_d=H(E[TT_d]). Increasing d
admits a path injection: insert the new last transitive clone immediately
after the previous last clone. Its exterior neighbourhood is identical,
and its only possible subsequent module vertices would have larger labels,
so this preserves every Hamilton path and is injective. Consequently

    L8≥L6≥A≥B>0,
    h=L8+34L6+232A+246B≥267A+246B.

In particular h≥513B, and
`h≥171A+342B` because A≥B. Hence
`6A+12B≤18h/513`, proving the lower bound. For the upper bound,
`4B≤4h/513` and `2A≤4h/513`, the latter following from h≥267A.
The unchanged branch lies inside the same interval. This is a universal
low-digit selector within about−3.51% to+0.78% of its source count.

This theorem gives residue selection near any source containing the
declared uniform C module. It does not assert that such source counts
meet every target interval, or that controlling a residue determines the
remaining quotient. Those are additional exact-target obligations.

## 3. A count-preserving native reset reaches the best possible divisor

The directed Hamilton path

    P=(1,0,4,5,6,2,7,3)

in C gives a native operation: complement every arc except the seven
spine arcs of P. The resulting strong tournament is

    N=(212,137,186,97,170,195,22,72).

Its ordered-cover list is

    (513,6294,33804,102696,190800,217440,141120,40320),

so H(N)=513, and both its raw and uniform-context gcds equal3. Thus this
literal source change removes **two powers of3 while preserving H**.
Its polynomial is

    F_N=x^8+36x^6+238x^4+238x²,
    F_N−F_C=2x²(x²−1)(x²+4).                       (4)

The reduction27→3 is optimal among **all** count-preserving replacements,
even those of a different vertex order. The inherited slot-polynomial
theorem gives integer monomial coefficients, F_T(0)=0 and

    F_T(−x)=(−1)^|T| F_T(x).

See [Session139, all-moduli, section0](../../session139_incidence_moduli_and_stability_20261102/all_moduli/PROOF.md)
and the direct permutation/cut proof in
[Session137, width, section4](../../session137_overnight_proof_frontiers_20261102/width/PROOF.md).
Modulo3 every integer is0,1 or−1. Therefore if3 divides H(T)=F_T(1),
then3 divides every integer value of F_T. Binomial inversion implies3
divides every C_k(T), hence also the exact uniform-context ideal. Every
H513 replacement has this unavoidable factor. N attains the lower bound.

There is a substantive connection to the older five-vertex native reset
in [RESET_JUNCTION J2](../../coverage_machine_induction_20261103/RESET_JUNCTION.md).
That actual H15 pair has response `E(x)=2x(x²−1)`, and (4) is

    F_N−F_C=(x³+4x)E(x).                            (5)

Thus the eight-vertex pair packages a response in the positive polynomial
cone generated by that smaller native operation: one TT3-padded response
plus four TT1-padded responses. The multiplier x³+4x is **not** a tournament
polynomial: it would require a three-vertex tournament with H5. Equation(5)
is a relation between responses; it does not silently turn the formal
sum of five experiments into one padding module. The actual pair C,N is
the physical realization of its left side.

The lost3-adic thickness can be observed at every length. Put k≥0 and wrap
`C→TT_k` cyclically with two new vertices a,b, using old→a→b→old. The
inherited cyclic-wrapper identity H(G(T))=H(T)+F_T(2) gives

    H(G(C→TT_k))=513+7128·2^k,
    H(G(N→TT_k))=513+7320·2^k.                       (6)

The first count has exact3-adic valuation3, and the second exact valuation1,
for every k. Their difference is192·2^k. These are strong actual graphs.
Thus the new smaller ideal is witnessed by an explicit unbounded family,
rather than only by a signed gcd combination of hypothetical contexts.

## 4. Breaking the last factor: an all-length primitive rooted return

Uniform incidence cannot remove the final3 at fixed H513. The rooted
language can deliberately breach it. Take the converse of C, retaining
physical root5. Its literal rows and rooted state are

    C^op=(198,136,122,225,139,147,50,68),
    z=(h,s,q,F)=(513,77,129,711).                    (7)

Use the inherited actual root-child commands from
[ROOTED_RESETS R1](../../two_port_return_induction_20261103/ROOTED_RESETS.md):
add x with `root→x→all other old vertices`; R retains the old root, and
N marks x instead. The first child already sees C nonuniformly. Every
subsequent graph and root is specified by the program.

For every k≥1 the word `R^k N R` returns a **balanced primitive state**

    (b_k,t_k,t_k,b_k),
    b_k=1746·2^k+241,
    t_k=582·2^k−17,
    b_k−3t_k=292,       gcd(b_k,t_k)=1.              (8)

These graphs have k+10 vertices and are strong. They preserve the original
eight vertices as an induced subtournament. Their first counts are3733,
7225,14209. No primality assertion is made.

Proof: the inherited birth maps are

    R(h,s,q,F)=(h+F−s,q,q,2F−q),
    N(h,s,q,F)=(h+F−s,h−s,h,h+F).

For (7), put a=h−s=436, d=F−q=582 and delta=d−a=146.
The direct iteration, equivalently
[RESET_JUNCTION J3](../../coverage_machine_induction_20261103/RESET_JUNCTION.md), gives

    t=a+q+(2^k−1)d,
    b=a+3q+(3·2^k−1)d,
    gcd(b,t)=gcd(delta,2^k d+q).

This is (8). Because b and t are odd, the only possible prime in the gcd
is73. Modulo73 the vanishing condition is `−2·2^k+56=0`, or `2^k=28`.
But the complete order-nine subgroup generated by2 is

    {1,2,4,8,16,32,64,55,37},

which excludes28. This proves primitivity for **every** k, not a finite
clock sample.

The continuation clock has a precise8/9 character description:
`73=Phi_9(2)=8²+8+1=8·9+1`. On the nonzero residues modulo73, the
character `chi(x)=x^9` has kernel exactly the nine-element subgroup
generated by2 and has eight image phases. But `chi(28)=10`, and
`10^4=−1 mod73`, so the forbidden residue28 is outside that kernel.
In fact28 has order72: its ninth power has order8, while its eighth power
is16 and has order9. This identifies the subgroup and quotient that
control the actual rooted construction's guard.

Strongness is preserved by every root-child birth: the old
root reaches the newborn, the newborn reaches another old vertex, and the
old strong tournament supplies all remaining connections.

The state is reusable, not merely a count coprime to3. By
[ROOTED_RESETS R3](../../two_port_return_induction_20261103/ROOTED_RESETS.md),
forming a fresh root-child seal from a rooted graph has continuation divisor
exactly gcd(h,s,q,F). Thus every fresh seal of the returned rooted graph
has divisor1. The **existing** final root-child pair is primitive as well:
its parent is the state immediately before the last R. With
`F_k=582·2^k+129`, its coordinate gcd is
`gcd(t_k,129,F_k)=gcd(t_k,129,146)=1`. These are distinct applications of
the same interface theorem; the old child is not silently treated as a
fresh seal. The returned root also has b≡t≡1 mod3. The inherited
physical C3 balanced-return map sends it to

    (b',t')=(61b−2t,9b),

which remains primitive and satisfies3|t', 3∤b'. The standard guarded entry
bridge in [RESET_JUNCTION J4](../../coverage_machine_induction_20261103/RESET_JUNCTION.md)
then enters the existing repeatable C3/C5 controller: apply C5 unless
t'/b'=7 in P1(F11); in that exceptional case apply one more C3, which
changes the ratio to3, then C5. A zero denominator is the infinity state,
already allowed by that bridge. These are inherited graph operations with
retained ports, not formal matrix realizations claimed from scratch.

This proves an actual unbounded exit from the supplied27-obstructed source
into a reusable primitive family. It does not prove that those continuations
select every prescribed odd target; arithmetic admissibility and exact
target supply remain distinct.

## 5. Mechanism ledger and reproduction

The inheritance pass used the exact excursion-factorial ideal as the nearest
mechanism; a fixed uniform module with divisor27 as the hostile; the false
inference “raw gcd9 is the whole continuation obstruction” as the corrected
near miss; and the older named-root defect clock as the underused sidecar.
The live concepts are uniform contexts / native source changes / monomial
response directions / named root defects /3-adic thickness / closed primitive
returns. The three successful maps change different coordinates:

| physical operation | preserved | changed | required sidecar |
|---|---|---|---|
| native C→N | H513 and all vertex names | uniform ideal27→3 | actual spine and full cover response |
| one arc C→C2 or C4 | other27 arcs | H and every uniform-context response by(3) | common exterior and marks outside module |
| root-child return | induced C^op and executable port history | current H and uniform incidence | full(h,s,q,F), plus the order-nine guard |

`escape.py` computes all five complete profile lists and polynomials, verifies
the literal native path, and counts four external-root response controls
(two exteriors, two single replacements), plus the coupled selector and its
relative bounds in both exteriors. It constructs and independently counts
the complete rooted states for k=1,2,3, and records the finite order-nine
certificate supporting the all-k proof. It does not audit the repository's
global spectrum or exhaust all larger source replacements.

Run

    python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/constructors/escape.py

To emit any literal returned graph, add `--k` and a positive integer.
The exact receipt is `escape_receipt.json`. The remaining Chan–Pak gate is
still a target-directed construction, rather than merely a primitive source
or a positive response cone. Equations(3) and(8) supply explicit operations
and their guards for that next synthesis.
