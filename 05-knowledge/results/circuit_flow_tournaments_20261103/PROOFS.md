# Composable tournament gates, hidden returns, and an exact clock

**PROVED in the specified construction language.** The compiler specializes
[THM-3121, path-cover walk-content substitution](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
The two actual gates are inherited from [Session176](../session176_tournament_information_20261104/spectrum/PROOF.md).
The new calculations below concern their complete operators, order of
composition, a sharp modular quotient, and a concrete neutral reset. No
literature-priority claim or complete Hamilton-spectrum claim is made.

## 1. The interface and its compiler

Let `F_A(r)` count covers of a tournament A by r labelled directed path
slots, allowing empty slots. If `pc_k(A)` counts unordered nonempty k-path
covers, then

    F_A(r) = sum_k pc_k(A) (r)_k.

Thus F is an integer polynomial, `F_A(1)=H(A)`, and `F_A(0)=0` for nonempty
A. For the empty tournament F is identically one. A transitive m-vertex
input has `F(r)=r^m`. We retain literal input vertices, internal orientations
and their names; a numerical polynomial alone does not reconstruct them.

Fix a tournament Q with a module vertex v and t other vertices. Substitute
A at v. The following finite compiler gives the **whole** output profile.
Partition the t exterior vertices into an unordered family of nonempty
words. Order these words canonically by their minimum label. Write ell for
the number of words; assigning them to distinct labelled slots has
`(r)_ell` possibilities. Their internal order need not be a directed path:
an intervening A-run can repair an incorrectly directed exterior adjacency.

A gap between a,b is available for an A-run exactly when `a->v->b`.
The initial gap is available when `v->b`, the final gap when `a->v`.
If `a` does not beat `b`, the gap must be available and occupied; otherwise
discard the pattern. Let g count available gaps and f count required
occupied gaps. Slots with no exterior vertex supply another `r-ell`
available gaps. Inclusion-exclusion gives the contribution

    (r)_ell sum_(i=0)^f (-1)^i binom(f,i) F_A(r+g-ell-i).       (CF1)

Sum over the exterior-word patterns. This is a finite shift operator with
integer-polynomial coefficients, determined solely by the actual context.

**Proof.** Each global slot path has a unique sequence of exterior vertices
and maximal intervening A-runs. Empty/occupied gap legality is exactly the
listed condition. Conversely choose such exterior words, assign slots and
place the covering A-paths into the available gaps with all required gaps
occupied. Concatenation reconstructs the unique original cover. Excluding
empty required gaps yields (CF1). For integers r below ell the falling
factorial vanishes. Prove first for r>=t, then extend as a polynomial
identity; negative intermediate evaluation arguments cause no boundary
assumption. This includes empty A. The operator is linear in the complete
profile, and has bounded shift range for fixed exterior size. ∎

This is an executable version of the older run kernel, not a replacement
for it. In particular its size depends on the exterior context, and it is
not a polynomial-time algorithm for arbitrary tournaments.

## 2. Two strong gates, now with their complete return information

Use outgoing masks

    L=(16,17,3,7,12), U=(16,25,3,5,12)

and substitute at vertex 0. Their sole differing edge is between exterior
vertices 1 and 3. They share the cycle `0->4->3->2->1->0`; nonempty
substitution therefore gives strong tournaments. Writing f for F_A, (CF1)
gives

    (Lf)(r) = -r(r-1)(2r-1) f(r-1)
              +r^2(r^2-2) f(r) +r(r+1)(2r+1) f(r+1),

    (Uf)(r) = -2r(r^2-r+1) f(r-1)
              +r^4 f(r) +2r(r^2+r+1) f(r+1).                 (CF2)

These are actual contexts, not arbitrary algebraic instructions. As checks,
at r=1 they give `H(L[A])=6F_A(2)-H(A)` and
`H(U[A])=6F_A(2)+H(A)` for nonempty A, as in S176. The full difference is

    ((U-L)f)(r) = -r(r+1) f(r-1)+2r^2 f(r)-r(r-1) f(r+1).   (CF3)

Scalar unit control retains only (CF3) at one. For a transitive m-module,
the ordinary difference is 2, whereas the same triangle wrapper from S159
changes the difference to

    8*2^m - 2*3^m - 4.

It is 6,10,6,-38 for m=1,2,3,4; it is strictly negative for every m>=4.
The latter follows since `2*3^m>8*2^m` for m>=4. This is a new explicit
specialization of the sign-reversal mechanism already proved in
[apex-wedge PROFILE](../apex_wedge_constructor_20261103/PROFILE.md).

## 3. Circuit order measures three disjoint paths

For every nonempty tournament A,

    H(L[U[A]]) - H(U[L[A]])
       = -12 (F_A(3)-3F_A(2)+3F_A(1))
       = -72 pc_3(A).                                      (CF4)

Indeed `L at 1 =6 eval_2-eval_1` and `U at 1=6 eval_2+eval_1`.
Their commutator is `6 eval_2(U-L)-eval_1(U+L)`, which by CF3
is `-12 eval_3+36 eval_2-36 eval_1`. The third forward difference
at zero counts three ordered nonempty path slots and equals `6 pc_3`.

Consequently the two nesting orders have the same H **exactly when
|A|<=2**. For |A|>=3 a Hamilton path can be cut into three nonempty
paths, so pc_3 is positive. Existence of that path follows from the usual
vertex insertion argument; no spectrum-completeness assumption is used.

This is a concrete distinction from the older
[S36 ordered-join cocenter](../session36_information_curl_contexts_20261023/root/CURL_COCENTER_AND_LIE_BOUNDARY.md):
H annihilates commutators of *ordered joins*, because H multiplies there.
Here the bracket is between *module-context operators*, and its H image
is a nonzero occupied-run count. The two algebras have different products.

## 4. A physical neutral reset with a later readable effect

Start with the same named singleton and apply two gates. In the following
table words are written **inner to outer**:

| Word | Order | H |
|---|---:|---:|
| LL | 9 | 589 |
| UL | 9 | 611 |
| LU | 9 | 611 |
| UU | 9 | 637 |

The UL and LU graphs retain the same nine labels and differ on exactly
two exterior edges. Their full profiles differ by

    F_UL(r)-F_LU(r)=6r(r-1)(r+1)(2r^2+3).                 (CF5)

Their two-slot counts are 14884 and14488. One common L attachment reads
the distinction:

    H(ULL)=88693, H(LUL)=86317, difference2376.

Thus the observer completion depth for this specified distinction is
exactly one: H alone cannot see it; one legal common context does.

The inherited [BC2 exterior-flip compiler](../boundary_carry_returns_20261103/PROOFS.md)
turns the two flips into **four legal native complements**, with six helper
births and six helper deletions, three helpers live at once. All nine old
names, their induced source vertex, and its incident star survive. Each
completed flip is strong. The supplied native Hamilton spines have the same
predecessor/successor frame at the old source root. `native_reset.json`
retains all cuts, spines, source/target masks and event clocks; the peak
order is12. Helper-expanded intermediate tournaments are not asserted
strong, and root-response phases are not asserted unchanged. This is a
new H-neutral reset, not an assertion of one-ready ternary return closure.

The physical source root does change its response: its `(H,start,H-delete,F)`
profiles are `(611,14,187,2928)` and `(611,24,197,2884)`. Thus this reset
lowers its F coordinate by44 while preserving H. The same standard two
births at that root produce counts6467 and6379, differing by88. Both are
literal source-preserving extensions, independently counted. Their large
defects2317 and2273 show exactly why this is useful response control but
not yet a ready-source reset.

## 5. Exactly how much linear information d future gates can read

Fix d>=0. Consider all functions `A -> H(w[A])` where w has at most d
L/U gates, including the empty word. Their rational linear span is exactly

    span{F_A(1), F_A(2), ..., F_A(d+1)},

and has dimension d+1 on actual nonempty tournament inputs.              (CF6)

The upper bound follows from the three shifts in CF2 and F_A(0)=0.
For the lower bound use the d+1 contexts `1,L,L^2,...,L^d`. The coefficient
of F_A(d+1) in H(L^d[A]) is

    alpha_d=product_(j=1)^d j(j+1)(2j+1)
           =(d+1)!(2d+1)!/2^d,                            (CF7)

which is nonzero; the observation matrix is triangular. Evaluations at
1,...,d+1 are independent on the actual transitive inputs of orders
1,...,d+1, since their matrix is `j^m`, an invertible Vandermonde matrix
times a nonzero diagonal. This proves sharpness on realizable sources.

Repeated L probes recover successively every F_A(j) by exact triangular
division. At known input order n, n probes including the unmodified H
recover the full degree-n polynomial (also F_A(0)=0). They still need not
recover the tournament: the hostile control below proves that loss.

For a *specified* depth-d word, an array of d+2 initial evaluations and
the successive three-shift updates evaluates H in O(d^2) integer
arithmetic operations, O(d) registers. This cost assumes those initial
evaluations are supplied; arbitrary input profile computation is not free.
For a singleton they are simply0,1,...,d+1. Integer bit complexity is extra.
CF6 is a linear-observer theorem, not a lower bound for unrestricted
nonlinear encodings, Boolean circuits, or named-port source reconstruction.

## 6. A sharp order-blind clock, despite noncommuting exact circuits

Now restrict to words starting at the singleton, with no other operations.
Let a and b count L and U, and let `i=a+14b mod18`. For every nonempty word,

    H(w) mod8  = 3 if a is odd, 5 if a is even,
    H(w) mod27 = 2^(i+12i^2) mod27.                         (CF8)

These determine H mod216 by CRT. The empty word has H=1 and is separate.
This is a **complete finite induction certificate**, not extrapolation
from bounded word lengths. `clocks.json` supplies five reachable evaluation
vectors modulo8 and162 modulo27, every transition and a reaching word.
Both sets contain `F_singleton(r)=r`, and every L/U successor remains in
the displayed set under CF2. The exterior-word compiler independently
checks each transition.

Modulo27 the output values in clock order i=0,...,17 are

    1,11,22,8,7,14,10,2,4,26,16,5,19,20,13,17,25,23.

They are distinct, equal `2^(i+12i^2)`, and every vector with output at i
has successors at i+1 and i+14. Modulo8 the three output classes have
transitions `1 -> (3,5), 3 -> (5,3), 5 -> (3,5)`. Induction proves CF8.
No use is made of a finite-height search claiming to cover infinite words.

The combined deterministic H-mod216 observer has **exactly19 states**:
one empty state and18 nonempty clock states. All are reachable, and all
have different current outputs, proving minimality. The L/U transitions
commute on the nonempty clock, and agree in either order from the empty
state. Therefore every rearrangement of a gate word preserves H modulo216.

Moreover216 is the **largest positive integer** with this universal
order-blindness property. The two differences

    H(LLU)-H(LUL)=2376,
    H(LLLU)-H(LLUL)=1364472

have gcd216. Any modulus hiding every gate reordering divides both;
CF8 proves every divisor of216 works. Larger moduli cannot all hide order.
This statement concerns this exact grammar and initial singleton only.

## 7. Information that the profile circuit still cannot carry

The labelled 81/27 pair has the same complete path-cover polynomial.
The old construction expresses it as three cyclically substituted blocks,
with one block replaced by its converse; reversing all paths preserves
the block profile. THM-3121 therefore proves profile equality. Consequently
**every L/U word remains blind to the pair**, at every depth. The native
relations and their automorphism orders are not interchangeable.

`audit.json` independently retains the literal pair, all13 polynomial
evaluations, actual full automorphism counts, and a common nonuniform
one-vertex cut that distinguishes them. This isolates the scope boundary:
CF6 recovers exactly the uniform-context profile, not arbitrary joint-port
attachment data. Source rows and label maps remain the fallback carrier.

## 8. Consequence for the open spectrum route

The new reset and clock are tools for source selection. They neither
close the ready source family at35207 nor give a uniform inverse reaching
every required prime and7p. Indeed CF8 implies that this L/U-only grammar
has H congruent to1 or11 modulo12, a strong restriction. Its useful new
feature is controlled alteration of hidden return information at fixed H.
The next test is whether such neutral, source-preserving operations alter
the specific root response needed by a ternary continuation.
