# Global fixed-receiver witness compression

**Session282. PROVED, independently reviewed.** This extends the
fixed-core receiver inverse to an explicit finite all-tournament search
boundary. It is not an efficient class enumerator or a Chan–Pak proof.

## 1. Definitions and inherited facts

Let T be a pointed tournament with root v, D=T-v, h=H(D), and receiver V,e as
in the rooted-fusion law. H(empty)=1. For a deletion vertex i define

    E_i(D)=sum_(S contains i) end_i(S) H(D\S).

Use the positive endpoint-weight and ordinal-factor formulas proved in
[Session281 receiver§2](../../session281_chanpak_20261112/receiver/PROOF.md):

    V=h+sum_(i->v) E_i(D),
    E_i(D)=(product_(a<j)Q(C_a)) E_i(C_j)
                         (product_(a>j)H(C_a)),       (1)

when C_1,...,C_m are the ordered strong components of D and i lies in C_j.
Every E_i(C_j)>=H(C_j): splitting each Hamilton path just after i injects its
paths into the defining two-slot objects. A singleton has H=1,Q=2.
Strong-component path counts multiply; a strong component of order n>=3 has
at least n Hamilton paths, since one Hamilton cycle (Camion's classical
strong-tournament theorem) has n distinct cuts. We use these known facts,
not a new proof or new priority claim for them.

## 2. Explicit bound, preserving the actual receiver

Set b(1)=0 and b(h)=h for h>1. If V=h, T has a source root and admits a
same-(H,h,V,e) witness with at most b(h)+1 vertices. If V>h then necessarily
V>=2h, and T admits a same-(H,h,V,e) witness with at most

    b(h)+2+floor(log2((V-h)/h)) vertices.             (2)

The logarithm means the largest integer s with h*2^s<=V-h; implement it by
integer inequalities, not floating rounding.

**Proof.** Write the nontrivial strong components of D as C with path counts
h_C>=3. Their product is h; their total number of vertices is at most
sum_C h_C. If there is more than one such component, a+b<=ab for a,b>=3,
so induction bounds this sum by h. If h=1 there are none. Thus all nontrivial
components together use at most b(h) vertices.

When V>h, let C_j be the LAST strong component having any vertex i beating
v. Every later component is beaten entirely by v. Therefore the full
pointed tournament is an ordinal sum T_prefix -> C_(j+1) -> ... -> C_m,
where T_prefix contains v. Delete every singleton component in this suffix.
This preserves H and h because its ordinal path factor is1. It also preserves
V and e: a path ending at v cannot visit any suffix vertex, so the suffix
lies wholly in the complementary path and contributes only its Hamilton
factor. Equivalently use the unsigned/signed version of(1).

Let s be the number of singleton components strictly before C_j. Choose one
incoming vertex i in C_j. Formula(1) and E_i(C_j)>=H(C_j) give

    V-h >= E_i(D) >= h*2^s.

Each earlier singleton contributes2 rather than its Hamilton factor1;
earlier nontrivial components have Q>=H and can only increase the bound.
Thus s is bounded by the logarithm in(2). After the deletion just performed,
there are s singleton components before C_j and at most one more at C_j.
Add b(h) nontrivial vertices and the root to obtain(2).

When V=h there are no incoming vertices, by positivity of E_i. The root is
therefore a source. Delete every singleton component of D; H=h,V=e=h remain,
and the resulting order is at most b(h)+1. The h=1 result may be the singleton
root itself. ∎

This is a physical deletion of H-neutral suffix components, not identification
of formal receiver tuples. It generally changes Q,U,d and the isomorphism
class, so does NOT preserve arbitrary two-sided fusion contexts or root moves.
It preserves exactly the receiver bank data stated above.

## 3. Consequences and the remaining gap

For any specified positive odd h and integer V, the complete possible
(H,h,V,e) fiber over ALL pointed tournaments can be decided by finite
labelled enumeration up to(2), or the source bound when V=h. If V<h or
h<V<2h, the fiber is empty. Other arithmetic values may also be impossible;
the theorem does not declare the bound sufficient for realizability.

The fixed-core inverse in [renewal](../renewal/PROOF.md) is much cheaper on
its structured universe. The global theorem explains exactly what is omitted
when only a short core menu is used: all finite strong-core patterns within
the h-vertex budget, including their ordering, incoming cuts and intervening
singleton strings. A core menu is complete only after those possibilities
are covered, not because it produced many witnesses.

The bound permits simultaneous constraints on H and e by filtering the finite
universe, but it does not give a useful V bound for a prescribed H. A huge
receiver may encode a small scalar count, and choosing which physical receiver
serves each target remains the Chan–Pak realization problem. This theorem is
an exact interface-compression principle, not an eventual-threshold theorem.
