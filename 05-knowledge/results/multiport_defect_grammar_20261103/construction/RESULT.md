# Exact two-birth transport and a closed strong tripling channel

Status: **PROVED algebraic operator and scoped all-order construction**;
**FINITE-EXACT** small extension cubes. Closure under both nonzero ternary
digits remains **OPEN**. No theorem ID is allocated here.

## Inheritance and scope

The literal two-birth prototype is inherited from
`../../session180_170807_returns_20261104/phases/RESULT.md`:
`B=[8,1,3,6]`, old mark `w=0`, auxiliary `a=1`, and `(h,F)=(5,7)`.
Its two outputs have `(17,21)` and `(19,21)`.

The closest general mechanism is the standard two-birth ready reset in
`../../two_entry_funnels_20261103/READY.md`, especially section 4. Our C2
operation is **exactly its converse**. The scalar reset and its failure to
remain 1/2-ready are inherited, not new claims. The refinement here retains
a second physical mark w, computes its exact joint endpoint response, and
obtains a closed simultaneous `(h,F)` tripling law.

The canonical hostile is the actual `(19,21)` child: even arbitrary two births
and arbitrary output re-marking cannot realize its required next response pair.
The near miss is the apparent `5 -> 17 -> 53` recurrence; its six fixed-mark
second-step realizations are all one isolated clone-block enlargement. The
least-used sidecar is the joint count for a path starting at a and ending at w.
The concept board is: exact old relation; two named ports; empty-cut conventions;
positive entry defects; clone-block clocks; inherited readiness obstruction.

Nothing here identifies these small cores with the saved H=170807 tournament,
its equal-count alternatives, or any other graph having matching scalars.

## 1. The subset ring and its unit

Let B have finite vertex set V. Work in the commutative square-zero ring

    R_V = Z[z_v : v in V] / (z_v^2 : v in V).

Write `[f]` for the coefficient of `z_V=product_(v in V) z_v`. The path
enumerator below is a polynomial in **vertex variables**, not the univariate
ordered-path-cover or unordered-path-cover polynomial used elsewhere:

    P = 1 + sum_(nonempty S subset V) H(B[S]) z_S.

In particular **P has constant coefficient 1**. Empty path slots are allowed
by that unit; omitting it invalidates every formula below. Define

    A = sum_S s_a(B[S]) z_S,     W = sum_S e_w(B[S]) z_S,
    E = sum_S e_a(B[S]) z_S,     Q = sum_S q_aw(B[S]) z_S.

Here `s_a` counts paths starting at a, `e_w` paths ending at w, and `q_aw`
paths starting at a and ending at w. These four polynomials have zero constant
coefficient. They also vanish on subsets missing a required port.
The case `a=w` is allowed: `Q=z_a`, since the singleton is the only simple
path starting and ending at the same vertex.

Products enforce disjoint vertex sets. Repeated P factors represent ordered
slots: their coefficients count all assignments to the indicated factor slots.
Set

    h=[P],     F=[WP],     S=[PA],
    L=[WPA],   K=[P^2 Q],  G=L+K,   J=[WE],   e_a=[E].          (1)

Thus `F=E_w(B,2)` and S is the dual two-slot start response at a. All
coordinates in (1) are nonnegative integers.

Both `sigma=S-h` and `tau=G-F` are nonnegative. For sigma, cut a Hamilton
path immediately before a. This injects Hamilton paths into the two-slot
covers counted by S. For tau, take a two-slot cover counted by F and cut the
path containing a immediately before a. If a lies in the second path, this
gives a term of L; if a lies in the first path ending at w, it gives a term of
K. The slot types distinguish the cases, and restoring the removed join
recovers the input. This is an injection, including the empty-cut cases.
Consequently the deficits count extra disconnected joins, not signed formal
cancellations.

## 2. The actual births and general operator

Append a physical vertex x with `x -> a` and `v -> x` for every old `v != a`.
Call the resulting tournament X. Two possible second births are:

* **C2:** append y with every original vertex `v in V` beating y, and `y -> x`.
* **C4:** append y with `w -> y`, and y beating every other vertex of X.

These names identify recipes, not guaranteed offsets on arbitrary inputs.
All old arcs and old vertex names are retained. In out-row bit encoding, with
old labels `0,...,n-1`, x=n, y=n+1 and birth cuts listing vertices beating the
newborn, the first cut is `((1<<n)-1) xor (1<<a)`. The C2 second cut is
`(1<<n)-1`; the C4 second cut is `1<<w`.

For a nontrivial strong B (necessarily order at least 3), X and both final
cores are strong: X has an incoming old incidence and the outgoing arc to a;
C2 has `B -> y -> x -> a` and the old strongly connected B; C4 has an incoming
arc from w and outgoing arcs to old vertices. The singleton boundary is kept
separate: X need not be strong, although the polynomial identities still hold.

Put `U=P-E+PA`. The path and old-w endpoint polynomials of X are

    P_X = P + z_x U,        W_X = W + z_x P Q.                (2)

Indeed a path containing x either ends at x, in which case its old prefix
cannot end at a, or has x immediately followed by a. Hence

    H(X)=h-e_a+S,       F_w(X,2)=F-J+G.                      (3)

Write `s_X=s_w(X,1)` and `D=H(X-w)`. The two complete operators are

    H(C2)=h+2S,             F_w(C2,2)=F+2G;                 (4)
    H(C4)=H(X)-s_X+F_w(X,2),
    F_w(C4,2)=2F_w(X,2)-D.                                 (5)

For (4), the complete enumerators are

    P_C2=P+z_x U+z_y P+z_x z_y(P+2PA),
    W_C2=W+z_x P Q+z_x z_y P Q.                            (6)

Extracting the full coefficient of `P_C2` and `W_C2 P_C2` gives (4).
For (5), a path through y starts at y, or enters y immediately after w.
Thus `P_C4=P_X+z_y(P_X-S_w,X+W_X P_X)`. A path ending at w either omits y
or starts at y, excluding the singleton-w continuation. Therefore
`W_C4=W_X+z_y(W_X-z_w)`. In the product, terms containing w twice vanish;
the remaining coefficient is `2F_w(X,2)-H(X-w)`.

These are exact operator laws. In particular C2 is an all-clock digit d iff

    S=h+d/2,        G=F,                                   (7)

and C4 is such a digit iff the two expressions in (5) equal `3h+d,3F`.
The full response, rather than h alone, is necessary for an all-clock lift.

## 3. Transporting a to x closes a strong zero-digit channel

After C2, keep the old physical w and designate x as the new auxiliary.
Its start and joint polynomials are

    A_x'=z_x(1+A)+z_x z_y A,       Q_xw'=z_x Q.

Using (6), coefficient extraction gives

    S_x'=h+2S=H(C2),
    L_x'=F+2L,       K_x'=2K,
    G_x'=F+2G=F_w(C2,2).                                  (8)

The unit in `1+A` is essential: it counts the singleton path x.
Equations (4) and (8) give the exact transported map

    (h,F,S,G) -> (h+2S, F+2G, h+2S, F+2G),                (9)

or, in nonnegative entry-defect coordinates,

    (h,F,sigma,tau) -> (3h+2sigma, 3F+2tau, 0, 0).         (10)

Thus the class of actual marked strong cores with `S=h,G=F` is closed under
C2 followed by this physical auxiliary transport. Every C2 image lies in it.
On this class the exact action is `(h,F) -> (3h,3F)`.

This iterates without an order bound: after k further steps the core has
`n+2k` vertices and response `(3^k h,3^k F)`. At each step x is a fresh named
vertex, while every old vertex and old arc remains. The graph relation itself
grows and is retained; this is not replacement by a scalar-equivalent graph.

For the standard zero-ready attachment `Z_s(B,w)`, attach each newborn below
the retained prefix and retained outer root. Since the core w is unchanged,
the old outer root keeps its old star and sole inport. The full counts satisfy

    H(Z_s)=h+2^s F,
    H(Z_s after entry)=3H(Z_s)+2sigma+2^(s+1)tau,
    H(Z_s after each subsequent channel step)=3H(Z_s).     (11)

This is a source-preserving, strong, arbitrary-depth **zero-letter channel**.
It does not produce either nonzero digit after entry.

## 4. Why the natural auxiliary choices do not rearm the digit

Let `R=E_w(B-a,2)` when `a!=w`, and R=0 when `a=w`. Directly expanding the
old-a and new-y start/joint polynomials gives

    S_a'=4S-H(B-a),      G_a'=4G-R;
    S_y'=H(C2)+h-e_a,    G_y'=F_w(C2,2)+F-J.               (12)

The old-a start polynomial used here is
`A+z_x(A-z_a)+z_y A+z_x z_y A`; its joint-to-w polynomial remains Q.
The new-y start polynomial is `z_y+z_x z_y(1+A)` and its joint-to-w polynomial
is `z_x z_y Q`. These give (12) with no asymptotic assumptions.

There are strict debts on a nontrivial strong input:

    F-J >= e_w(B,1) > 0,       F-R >= e_w(B,1) > 0.        (13)

For the first inequality, `[WE]` counts only covers whose second slot ends
at a, whereas the empty-second-slot covers already contribute `e_w(B,1)`.
For the second, insert a into the second path of a cover on B-a. A vertex can
always be inserted into a directed path while preserving its relative order:
prepend if possible, append if possible, otherwise insert at the first change
from an edge pointing into a to an edge pointing out of a. The image has a in
a nonempty second slot, and deletion of a recovers the input. Empty-second-slot
Hamilton paths on B are additional. The case a=w has R=0 directly.

To see `e_w(B,1)>0` without a counting assumption, take a longest directed path
ending at w. An omitted vertex either inserts into it or is dominated by every
path vertex. If every omitted vertex is dominated, the path's vertex set
dominates its nonempty complement, contradicting strongness.

If the entry step had `G=F`, then `F'=3F`. Keeping old a makes
`G_a'-F'=F-R>0`. Choosing new y makes `G_y'-F'=F-J>0`. Therefore either choice
incurs positive endpoint debt when used for the next C2 step. Even when its H
offset happens to be 2 or 4, the matching all-clock F identity fails.

There is a stronger inherited obstruction covering **every** auxiliary, not
just these natural choices. Converse transforms C2 exactly into the standard
two-birth reset. By `two_entry_funnels_20261103/READY.md`, section 4, the reset
image of any strong input of order at least 4 has no 1-ready or 2-ready vertex,
under arbitrary re-marking. Dualizing, **no auxiliary in any such C2 image has
`S_a(2)-H` equal to 1 or 2**. Consequently C2 cannot produce either next
nonzero digit at any auxiliary, even before imposing the F guard.

The order boundary matters. A C3 input has a C2 image of H=9, whose y can give
H=31 at the next C2 application; equation (13) still prevents the corresponding
same-clock F tripling. A singleton source gives C3 and more than one zero-ready
start vertex. For nontrivial strong input, x is the unique zero-ready start
vertex in C2: every old vertex retains an old out-neighbor and additionally
beats y, while y deletion leaves strong X; x deletion leaves `B -> y`.

## 5. The apparent 17-to-53 chain is a single clone-clock coincidence

Keep the first almost-sink x and replace y by a nonempty transitive block Y
of length m, with `B -> Y -> x`. For every m>=1,

    h_m=h+2^m S,       F_m=F+2^m G.                       (14)

If x ends a Hamilton path, all Y vertices precede it in one terminal block,
giving h. Otherwise x is followed by a, and the Y vertices can be assigned
independently to the block before x or the terminal block, giving `2^m S`.
For a two-slot cover ending at w, x in the second terminal block gives F.
Otherwise the two possible old-core arrangements counted by L and K each
offer two Y slots, giving `2^m G`.

For the prototype `(h,S,F,G)=(5,6,7,7)`, this is

    (h_m,F_m)=(5+6*2^m, 7+7*2^m).

The six fixed-mark H=17 to H=53 extensions are precisely the six transitive
orders on the old y and its two newborn clones: m=1 becomes m=3. An increase
by two clones gives `h_(m+2)=3h_m+2` iff `2^m S=2h+2`; here that forces m=1.
Thus these actual witnesses do not define an unbounded +2 channel.

## 6. Exact finite cubes and audit

`probe.py` checks all 1105 labeled marked cases through order 4, including weak
tournaments, arbitrary w/a, coincident ports, empty cuts and the singleton.
It checks (2)--(5), (8), (12), and the positive entry defects against direct
subset enumeration. It also enumerates all ordered two-birth extensions of
the actual prototype and its two actual children. The strong final core is
required; weak intermediate cores are permitted, so this is broader than a
stepwise strong birth menu. Every old arc is retained.

| Actual source | Cube size | Fixed w: +2 / +4 | Arbitrary output mark: +2 / +4 |
|---|---:|---:|---:|
| B, `(5,7)` | 512 | 2 / 4 | 6 / 12 |
| displayed C2, `(17,21)` | 8192 | 6 / 0 | 12 / 0 |
| displayed C4, `(19,21)` | 8192 | 0 / 0 | 0 / 0 |

Each hit must have `(h',F')=(3h+d,3F)`. Moving-mark hits are response-pair
certificates only: attaching them under an existing outer root would require
explicit old/new inport-star changes. No phase preservation is inferred.
Both a deletion recurrence and an independently directed forward recurrence
agree on every full menu; Python separately checks every positive literal row.
There is no all-order birth-only impossibility claim here.

`amplifier.py` records the first four exact zero-channel states
`(17,21),(51,63),(153,189),(459,567)`. Every old relation is checked. Its empty
nonzero-digit guards are finite controls for the inherited all-order C2
obstruction, not its proof. C4-recipe absence is claimed only on those listed
states; other old-arc changes or larger construction schedules remain open.

Reproduction from the repository root:

```text
python 05-knowledge/results/multiport_defect_grammar_20261103/construction/probe.py
python -O 05-knowledge/results/multiport_defect_grammar_20261103/construction/probe.py
python 05-knowledge/results/multiport_defect_grammar_20261103/construction/amplifier.py
python -O 05-knowledge/results/multiport_defect_grammar_20261103/construction/amplifier.py
```

The strongest next target is a literal repair of the nonnegative defect pair
`(sigma,tau)` that enables a nonzero digit and preserves the enlarged named
boundary for the following step. This is more specific than finding another
small core with the same `(h,F)`. The closed zero channel, the two-birth hostile
cube, and the all-order re-marking obstruction supply independent positive and
negative controls for that repair.

## 7. A nonempty +4 preparation locus and a typed next target

The locus `sigma=2,tau=0` is not merely a formal condition. An actual strong
five-vertex core provides it:

    B=[20,1,2,7,14],       w=1, a=2,
    (h,F,S,L,K)=(9,21,11,9,12).

The C2 child is `[116,97,66,103,110,4,32]`, retains every old arc and has
`(h',F')=(31,63)`. It therefore realizes the all-clock +4 law by C2 itself.
In contrast, the inherited four-vertex `(5,7)` prototype supplies the
`sigma=1,tau=0` preparation and hence the C2 +2 law.

`entry_guards.py` checks that no `sigma=2,tau=0` guard occurs among all 403
labeled strong marked `(B,w,a)` cases through order 4. Across the 180 saved
marked representatives in Session180's small bank, it checks 1044 auxiliary
choices and finds six +4 guards. This second universe is the saved
representative bank, **not** a complete tournament census through order 6.
No representative in that bank has both guard types at alternative auxiliaries.
The actual rows and responses, including positive outputs, are saved in
`entry_guards.json`; the smallest found +4 witness has order 5.

A precise sufficient preparation problem for this **C2-based** grammar is:
given a reachable actual marked core `(B,w)` in the zero channel, construct,
for each d in {2,4}, a legal named source repair to `(~B,w)` with

    H(~B)=H(B),       E_w(~B,2)=E_w(B,2),
    some actual auxiliary a satisfying sigma=d/2, tau=0.

Applying C2 then produces the required digit and returns to `sigma=tau=0`
at its new x. A reachable class closed under these preparations and the
zero channel would give an iterated three-letter grammar. This is a **typed
sufficient criterion**, not a claimed reduction of every possible grammar or
a proof that such repairs exist. C4 and larger schedules may solve the original
problem without satisfying this particular preparation condition.

Reproduction of this additional exact probe:

```text
python 05-knowledge/results/multiport_defect_grammar_20261103/construction/entry_guards.py
python -O 05-knowledge/results/multiport_defect_grammar_20261103/construction/entry_guards.py
```
