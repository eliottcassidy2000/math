# A canonical zero-ready boundary and the new physical mark

2026-11-03. **PROVED decomposition using inherited cut and birth identities;
FINITE-EXACT identification of the two actual endpoints.** No closed ternary
return is claimed merely from this common interface type.

Let T be strong, r a vertex and b_r=E_r(2)-H(T). Then b_r=0 iff:

1. r has exactly one in-neighbor w;
2. P=T-r has at least two strongly connected components;
3. w belongs to its last component in the unique transitive SCC order.

Write P=A->B, where B is that last strong component and A is the nonempty
union of preceding components. This is a canonical marked decomposition.
Conversely any nonempty tournament A, strong B, and w in B give such a
strong zero-ready T by adding r beaten only by w.

Proof. The inherited rooted defect gives b_r>=d_r(d_r-1), so zero forces
d_r=1. For the sole in-port w the exact single-in-port identity is
b_r=start_w(P). A vertex begins a Hamilton path of a tournament exactly
when it is in the first strong component: a Hamilton cycle in a nontrivial
strong component supplies a path starting at each of its vertices, and
component paths concatenate in their forced order. Thus b_r=0 iff w is
not in the first component. Strongness of T forces w to be in the last:
a later component could otherwise not reach r or an earlier component.
Conversely, with w in the last and A nonempty, every vertex reaches w
then r, and r reaches A then all of B, proving strongness. These arguments
also prove uniqueness of the decomposition and its equality boundary.

Let P_A(2) count ordered two-slot covers of A, allowing either slot empty.
The count of the resulting tournament is exactly

    H(T)=H(A)H(B)+P_A(2)E_w(B,2).                       (Z1)

Indeed H(P)=H(A)H(B), start_w(P)=0, and the inherited one-birth formula is
H(T)=H(P)+E_w(P,2)-start_w(P). Every two-slot cover of A concatenates
uniquely with each w-ending two-slot cover of B, so
E_w(P,2)=P_A(2)E_w(B,2). No connectivity or joint-incidence condition is
dropped: all cross edges are A->B and the first B slot is nonempty.

When A=TT_s, P_A(2)=2^s and H(A)=1. Thus every fixed strong marked B
gives the all-size zero-ready family

    H(T_s)=h_B+2^s F_w, s>=1, F_w=E_w(B,2).             (Z2)

This is a sufficient two-number evaluator for this fixed attachment
operation. Full B and its named port are still required for arbitrary
complements or source changes.

## The literal new return moves the attachment port backward

For the saved H56935 source, take r=50, w=39, A={51}, and
B={1,2,4,6,7,8,14,21,22,31,32,33,39,40}. The exact fields are

    h_B=18977, F_w=18979, H(T)=18977+2*18979=56935.

For the new H170807 endpoint, switch the selected attachment port to the
existing vertex r=39. Its sole in-port is w=31. Here

    A=(50->51),
    B={1,2,4,6,7,8,14,21,22,31,32,33,40,54},
    h_B=2195, F_w=42153,
    H(T)=2195+4*42153=170807.

The old mark50 and spare6 remain physical vertices. This is a selection
of an actual named vertex, with its joint responses recomputed; it is not
a class-valued state reset. The complemented program preserved frames at
50 and6 before this final selection. Port39 did participate in a changed
edge, so its new phase is never inferred from the old marked phase.

Even arbitrary re-synchronization of these two transitive-prefix clocks
does not repeat the+2 return. For s,t>=1, the equation

    2195+42153*2^t=3(18977+18979*2^s)+2

has exactly(s,t)=(1,2). If s=t its left-minus-right exponential part is
negative. Otherwise the2-adic valuation of
42153*2^t-56937*2^s=54738 forces min(s,t)=1; direct substitution gives
only the displayed solution. Changing B or its actual incident relation
is therefore necessary for reuse, even within this common zero-ready
grammar. The general existence of reusable such B-transforms is OPEN.

The+4 equation has no solution on these two clocks: its right side becomes
54740, of2-adic valuation2. Unequal exponents force min(s,t)=2. If s=2<t,
the left side is at least109476; if t=2<s it is negative. Equal exponents
are again negative. This excludes arbitrary clock resynchronization, not
only equal increments of s and t.
