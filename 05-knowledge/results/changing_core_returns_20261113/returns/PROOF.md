# Exact target inversion with a retained start channel

**PROVED IN NOTE:** the arbitrary-core target interface, start/defect
propagation, and the all-source resource formulas below. **FINITE-EXACT:**
the first thirty ray targets, the bounded-core exclusion, and the stated
shadow search. **OPEN:** a uniform selector, a repeatable ready return, and
Chan--Pak. These are endpoint constructors, not a proved path of
count-preserving native edits between arbitrary sources.

## 1. Inheritance, scope, and concept board

The nearest proved mechanism is the six-coordinate compiler in
[funnel_six_state.py](../../berggren_target_return_20261111/reset/funnel_six_state.py),
using the complete unique-start grammar of
[two_entry_funnels](../../two_entry_funnels_20261103/READY.md).
The incoming [multicore obstruction](../../multicore_metagraph_20261113/FAMILIES.md)
excludes homogeneous inflation from every ready881 source to1691, and
excludes the pure405*2^k+71 ray from every one-clock outer-tail family.
Those exclusions are inherited, not rerun here.

The scalar binary recurrence and section inverse are also inherited. See
[central compiler, section2](../../central_modular_chanpak_20261110/compiler/PROOF.md)
and, more directly,
[recursive seed switching, S1](../../recursive_seed_switching_20261103/selector/NOTE.md).
The latter already proves an exact binary-section inverse and a
square-root-scale bound. Our interval-pruned implementation is an
application to a different physical target interface, not a claim of a
new inversion principle. The exact periodic-word drift criterion is in
[paired descent, section4](../../h13_paired_descent_20261111/construction/PROOF.md);
the irregular words below do not evade that criterion by a claimed
periodic selector.

The canonical hostile is k17 below. The corrected near miss is extrapolating
twenty-nine successful unit-ready constructions to a uniform ray. The
least-used sidecar is the start multiplicity S(1), which must accompany
the six-coordinate compiler once unique-start cores are abandoned.

The board is: exact target inversion (anchor); the retained start channel
(niche); arbitrary-source resource profiles (wildcard); physical
indegree-one roots; bounded-core failure; and target-neutral marked
preparation. A scalar target hit is not a controller-entry certificate.

## 2. An arbitrary strong marked core and its exact target

Let C be a nonempty strong tournament with marked vertex u. Set

    (h,e,s,p,f,t)=(H(C),E_C(1),S_C(2),P_C(2),E_C(2),S_C(3)),
    b=S_C(1),    G=s+p.

Here P_C(j) counts ordered spanning covers by j directed paths with empty
bins allowed; E and S require the first bin to be nonempty and end or
start at u. In particular b is an actual number of Hamilton paths,
not an indegree.

For an integer n>=0, let w be its binary expansion read low bit first;
n0 gives the empty word. Add a transitive tail dominated by C, and a
new vertex x which beats u, loses to every other vertex of C, and loses
to tail position i exactly when w_i=1. Every nonempty word ends1. This
is the inherited literal funnel constructor. Call the result V.

The full six-coordinate state of V is given by the inherited compiler.
The extra coordinate obeys the particularly simple identity

    S_V(1)=b.                                             (1)

Indeed a Hamilton path beginning at x must enter C at u, traverse one
complete path of C starting at u, and then traverse the tail in its
unique order. Entering the tail before covering C would make return to
C impossible. This is a bijection, including the empty-tail case.
The funnel is strong when C is strong of order>=3; the singleton core
is also allowed with a nonempty word, which creates a directed cycle.
The singleton with an empty word is excluded from strong constructors.

Now add y with x->y and y beating every other vertex of V. Then y has
exactly one in-neighbor x. The inherited birth identity gives

    H(T)=H(V)+E_V(2)-b,
    (H(T),S_T(y),H(T-y),E_T(y;2))
      =(N,H(V)-b,H(V),N+b).                              (2)

Consequently the actual indegree-one readiness defect is

    H(T-y)-S_T(y)=b.                                     (3)

Indegree one by itself does not imply unit readiness. Equation(1) lets
the whole construction carry arbitrarily large defects faithfully.

Write g(0)=g(1)=1 and use the inherited recurrence

    g(2a+1)=g(a),      g(2a)=g(a)+a+(a mod2).

The tail coefficients are L_2(w)=n+1, B_1(w)=g(n), and
D_1(w)=1-(n mod2). Substituting the full funnel compiler into(2) gives

    N=G n+h g(n)+e(n mod2)+(G-e-f-b).                    (4)

This is the complete Hamilton count of the literal graph, not a selected
path subclass. The observer tuple is

    (h,e,G,f+b),

with b retained separately to recover its readiness defect. Equal
six-coordinate states without S(1) are not sufficient for this use.

### Exact inverse and its finite search bound

For A,B>0 and C>=0, solve

    A n+B g(n)+C(n mod2)+D=N,     n>=0.                  (5)

The two parity sections have coefficients

    even: (A,C,D) -> (2A+B,B,D),
    odd:  (A,C,D) -> (2A,0,A+C+D),                       (6)

while B stays fixed. Both branches are taken. No global monotonicity is
assumed: g has large downward jumps at suitable adjacent integers.

The elementary bounds 1<=g(n)<=n+1 give the necessary interval

    max(0,ceil((N-B-C-D)/(A+B))) <= n <= floor((N-B-D)/A). (7)

The upper bound for g follows by induction, using the stronger
g(n)<=n for odd n>=1 in the even step. The implementation prunes empty
intervals and directly checks intervals of width at most4. At every
section A at least doubles and D never decreases. Hence the upper
endpoint eventually becomes at most0 and every branch terminates.
Equations(6) partition all nonnegative integers, so the algorithm
returns every solution of(5), with no external word-length cutoff.

For comparison with the older inverse, put
N0=max(0,N-B-D) and C0=max(C,B). At depth j, the width in(7) is at most
B*N0/A_j^2+C0/A_j. Thus A_j>=max(sqrt(B*N0),C0) makes the width at most2.
This recovers a square-root-scale section count up to coefficient factors;
it is not a polynomial-in-log-N coverage algorithm.

## 3. Thirty exact ray values and a real defect-three escape

For N_k=405*2^k+71, k1..30, the initial bank comprises every unique-start
rooted strong core through order10, plus the singleton, deduplicated only
by the four scalar observer coordinates. It has1596 keys. Exact inversion
finds unit-ready sources at27 targets; the misses are k8,17,21.

Deeper cores recover two more unit-ready sources:

|k|target|core order/rank|chronological core words|outer address|total order|
|--:|--:|--:|---|--:|--:|
|8|103751|11/1561|1,111,empty,empty,1|3|15|
|21|849346631|12/1001|101,101011|7957|27|

All earlier27 actual graphs and marked roots are in
[inverse_receipt.json](inverse_receipt.json); these two are in
[deepening_receipt.json](deepening_receipt.json). The first30 therefore
have29 unit-ready witnesses in this packet. This count is not a statement
that k17 has no unit-ready source.

For k17, exact inversion over every unique-start rooted core through
order17, plus the singleton, finds no solution at any outer word length.
The explicit universe has1,346,269 rooted classes. This is a bounded-core
exclusion, not an all-order exclusion of ready sources.

The start channel gives a small concrete repair of the target. Begin with

    C0=(16,9,3,5,14),   marked vertex1,
    six-state(C0)=(15,3,19,108,19,75),   S_C0(1)=3.

An empty-tail funnel produces

    C=(48,9,35,37,46,2),   marked vertex5,
    six-state(C)=(31,12,34,328,89,183),   S_C(1)=3.

Use outer address n145623, the low-bit-first word

    111010110001110001.

Its parent state is

    (5319651,31,9576017,90429050940,47764583,92907885679).

The final birth produces an actual26-vertex strong tournament, root25,
unique predecessor24, and complete marked state

    (H,S,H-delete,E2)=(53084231,5319648,5319651,53084234).   (8)

The defect is exactly3. Complete outgoing masks, physical indices, and
both core profiles are in the deepening receipt. The two small core
profiles and their start multiplicities are independently counted from
their literal graphs. The large graph's full count follows from the
proved arbitrary-word compiler and(1)-(4); no 26-vertex subset-DP run is
claimed.

This fills the thirtieth **indegree-one target**, while explicitly
changing the readiness class. It does not turn the ray into a closed
unit-ready recursion. In particular, swapping the final root and its
predecessor as physical names preserves H but gives marked state

    (53084231,31,31,95529135).

For this state Delta=F-H=42444904 and d=F-q=95529104. The R^kNR return
has t_k=2^k*d+q-Delta, which is never0 modulo3. Thus that particular
endpoint-neutral preparation fails the native controller guard for
every k; the failure is arithmetic, not merely a failed short search.

A uniform replacement of C by a strong unit-start core with the same target
response(4) would require observer(31,12,362,92), hence E2=91 rather
than89. None exists among unique-start cores through order7. The
inherited strong-order lower bound in
[ready-return obstruction](../../ready_return_pentagon_20261103/OBSTRUCTION.md)
is H>=45 from order8 onward, so this is an all-order exclusion of that
specific H31 observer replacement. It does not forbid a different
observer hitting one selected target.

## 4. A closed interface with genuinely unbounded defect

The sibling [preparation theorem](../preparation/PROOF.md) studies
C3[TT_m,{x},B], oriented TT_m->x->B->TT_m. Let B be any nonempty
tournament with profile(h,p,t)=(H(B),P_B(2),P_B(3)). Mark a resource
vertex with l predecessors and r successors, l+r=m-1, m>=2.
Write A=2^l3^r and D=2^r. In addition to the sibling H and E2 formulas,
the exact inverse interface is

    H=(2^m-1)(p-h)+2h,
    b=2^r h,
    E1=2^l(p-h)-[r=0](p-2h),
    S2=(2A-D)p+(4D-2A-[l=0])h,
    P2=2[(3^m-2^m)t+(3*2^m-3^m-1)p+(1-2^m)h],           (9)
    E2=2(2^r3^l-2^l)(t-p)+3*2^l p-2*2^l h-[r=0](p-h).

Thus just three source evaluations determine all coefficients of(4)
for this entire family. Applying(1)-(3) makes the final indegree-one
defect b=2^r h unbounded as r increases. Nothing here asserts that a
chosen exact target is accepted by this family.

**Proof of P2.** More generally, place x in one of j labelled path bins.
Its bin has a prefix and suffix; each other bin is a B-run followed by
a resource run. There are j+1 B-slots and j+1 resource slots. Before
the two local exclusions, the count is (j+1)^m P_B(j+1).
The invalid prefix has nonempty B-prefix and empty resource-prefix,
giving j^m(P_B(j+1)-P_B(j)). The invalid suffix has empty B-suffix and
nonempty resource-suffix, giving ((j+1)^m-j^m)P_B(j). Their intersection
is (j^m-(j-1)^m)(P_B(j)-P_B(j-1)). Inclusion-exclusion, multiplied by
the j choices of the x-bin, gives

    P_C(j)=j[((j+1)^m-j^m)P_B(j+1)
        +(3j^m-(j+1)^m-(j-1)^m)P_B(j)
        +((j-1)^m-j^m)P_B(j-1)].                         (10)

Empty slots are allowed throughout. Setting j2 proves the P2 formula;
j1 independently recovers H.

**Proof of S2.** In the distinguished root-start bin, resources preceding
the root cannot occur. Split according as x is in this bin or the other
bin. For j general bins set
A_j=j^l(j+1)^r, D_j=(j-1)^l j^r, C_j=(j-2)^l(j-1)^r.
When x is in the distinguished bin, the count is
A_j P_B(j)-(A_j-D_j)P_B(j-1). When x is in one of the other j-1 bins,
the same two prefix/suffix exclusions as above give

    (j-1)[(A_j-D_j)P_B(j)
      +(3D_j-A_j-C_j)P_B(j-1)+(C_j-D_j)P_B(j-2)].

For j2, use P_B(0)=0 and 0^0=1; the resulting coefficients are exactly
S2 in(9). Finally, a full path ending at the root has two possible
B-slots with the second nonempty and2^l resource assignments. If r0,
the assignment with empty resource-prefix forbids both B-slots being
nonempty. This gives E1 in(9). The remaining formulas are the sibling
theorem and its direct root-start bijection. QED.

The compiler extension is useful even though the cheap hostile probe
B=(6,12,8,1), every m2..20 and root phase, has no inverse hit at k17.
Changing the defect interface does not make every parameter family
cover every target.

## 5. Verification and relation to the preparation lane

Run the following from the worktree root, then repeat each with `-O -B`
in place of `-B`:

```text
python -B 05-knowledge/results/changing_core_returns_20261113/returns/inverse.py --core-order 10 --k 30
python -B 05-knowledge/results/changing_core_returns_20261113/returns/deepening.py --max-core 17
python -B 05-knowledge/results/changing_core_returns_20261113/returns/growing.py
python -B 05-knowledge/results/changing_core_returns_20261113/returns/wrapper.py
```

The inverse controls compare3553 exact cases with brute integer inversion
or an independent tail automaton. The arbitrary-core controls use all
strong labelled four-vertex cores, every mark, and four tail words, and
check full H and actual defect by literal subset paths. The wrapper
controls compare all675 contexts formed from labelled B of orders1..4,
widths2..4, and every root phase, against the inherited ordered-bin
counter. All checks use explicit exceptions and survive optimized Python.
The graph receipts retain root indices and complete adjacency masks.

An independent audit of the sibling all-B preparation proof confirms
that full paths use at most two B-runs and its marked two-covers at most
three. Hence the K1,K2,K3 skeleton argument and coefficient identification
from all transitive B are valid. Expansion gives its claimed coefficients
(2A+2D-I,4A-D-I,2(A-D)). Its resultant and prime-phase proof also checks.
Our literal subset counter, without importing that lane's counter,
independently verifies H199 before and after preparation, root-start
counts5 and10, and the RNR output H4267,start1041,deletion1041.
That output has balanced defect0, not unit defect1.

The source-to-target connection is exact: an actual marked core and
binary tail map to an actual strong indegree-one tournament with full
count(4) and defect(3). Scalar inversion forgets all other root/port
geometry; the source masks, word, marked root, and b restore it.
The finite witnesses do not provide a uniform next-source selector or a
target-preserving native return. Those remain the specific open steps.
