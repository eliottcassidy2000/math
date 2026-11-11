# Bilateral length filtering preserves every arithmetic source address

**PROVED relative to the explicit F8 certificate in `../fourier/PROOF.md`:**
the bilateral all-size law, preservation of exact injectivity and balancing,
and the repeated-filter formula. **OPEN:** target coverage. This is a new
composition of two physical constructions, not a new verification of F8.

## Construction and exact state

Let T_a be any capped central-root allocation source from `PROOF.md`, with
six nonnegative multiplicities n_i and allocations 0≤a_i≤n_i. Its single
selected vertex is z. Write P=product h_i^n_i and rho_i=q_i/h_i, so

    D=H=P, W=S=E=0,
    U_a=2P product rho_i^(n_i−a_i),
    V_a=2P product rho_i^a_i.

Use the literal tournament F with outgoing masks

    (234,236,161,228,111,192,132,16).

The Fourier lane supplies its exact certificate

    h=H(F)=231, q=Q(F)=3696=16h,
    K_F(x)=231(1+x)(1+x+...+x^7).

For every integer k≥0 form the actual ordinal word

    B_k(T_a)=(F^op)^{⇒k} ⇒ T_a ⇒ F^{⇒k},

with no new marks. Empty exterior words at k=0 mean the original source.
This operation adds exactly 16k vertices. The ordinal interface gives

    H'=D'=h^(2k)P, W'=S'=E'=0,
    U'=(hq)^k U_a, V'=(hq)^k V_a,
    Q'=q^(2k)Q(T_a).                                  (1)

This proof requires no full-graph recount: adding a wholly unmarked prefix
multiplies U by its H and V by its Q; a suffix does the reverse. Their W
multipliers are both H. Converse leaves H,Q unchanged. The caps already
present in T_a keep both endpoint statistics zero, also at k=0.

## Both oriented length channels become exactly uniform

Let P^+_T(x) sum paths starting at z, weighted by the Hamiltonian count of
the complementary induced tournament, with exponent their number of vertices.
Define P^-_T similarly for paths ending at z. Thus P^+(1)=U and P^-(1)=V.

A path starting at z cannot visit an unmarked ordinal prefix. That prefix
is entirely in the complement and contributes its Hamiltonian count. An
ordinal suffix contributes its two-slot length-split polynomial K. Reversing
all paths proves the dual statement. Also K_(F^op)=K_F, since reversing paths
preserves each induced Hamiltonian count and each subset size. Therefore

    P^+_(B_k(T_a))=h^k P^+_(T_a) K_F^k,
    P^-_(B_k(T_a))=h^k P^-_(T_a) K_F^k.               (2)

For every k≥1 and every residue r modulo8, the outgoing and incoming counts
are respectively

    U'_r=(hq)^k U_a/8,
    V'_r=(hq)^k V_a/8.                               (3)

Indeed K_F modulo x^8−1 equals (q/8)S8=462S8; convolution by it sends
any residue vector to its total times462 in each coordinate. Repeating
multiplies the total by q each time, while the opposite exterior block
multiplies it by h. The displayed integers are genuine nonnegative path
counts. The k=0 case deliberately asserts no mod-eight uniformity.

## Arithmetic addresses survive phase erasure exactly

Fix k and the multiplicity vector. Multiplication of U−V by the SAME strictly
positive integer (hq)^k preserves distinctness and sign. Consequently the
filtered bank still has exactly product(n_i+1) distinct signed responses.
The balance criterion remains an iff:

    U'=V' iff every n_i is even and a_i=n_i/2 for every i.

In fact each single residue pair already retains the address: U'_r−V'_r
is the same positive scalar multiple (hq)^k/8 of U_a−V_a. An observer that
normalizes each length-residue distribution to total1 would lose this address;
the actual unnormalized path-count channels preserve it. Erasing phase is
therefore different from erasing amplitude.

All filtered allocations also retain the SAME COMPLETE unmarked cover
profile, because for every j≥1 their common original Q_j is multiplied by
Q_j(F)^(2k). No larger cover-order state is needed to repeat the filter:
the seven-coordinate carrier and the two eight-residue vectors close under
every bilateral step. The growing integer sizes and the factor h^(2k) in
the offset are real costs, not free projection or a unit-offset return.

For example the all-zero multiplicity seed is TT3 with central mark,
D=H=1 and U=V=2. One bilateral step has

    D'=H'=53361, U'=V'=1707552,
    U'_r=V'_r=213444 for each r modulo8,

on 19 vertices. Applying the same step to the six-species one-copy bank
produces 64 distinguishable sources on157 vertices with both length
channels uniform. None is balanced because all six multiplicities are odd.
For arbitrary multiplicities and k the vertex cost is exactly

    3+4n1+6n2+10n3+18n4+34n5+66n6+16k.

The existing endpoint-zero closure now gives actual outputs

    h^(2k)P + 2^t (hq)^k U_a, t≥0,

where t=0 is the bare-apex case and positive t is the selected-tail case.
This produces an explicit filtered family, not arbitrary odd-target coverage.

## Component certificate and reproduction

`construct.py` exposes `bilateral_filter(T,marks,profile,k)` to emit the
literal larger tournament and exact compiled profile. The receipt checks
k=0,1,2,3 on the actual TT3 central-root seed. It computes that seed's two
length polynomials directly and performs exact integer convolutions with
the accepted F8 polynomial, checking all sixteen residue totals against (1).
It also checks the all-zero initial counts above. These are COMPONENT
composition controls, not claimed independent subset enumeration on the
19-vertex and larger filtered graphs. The all-size proof is the unique
ordinal path/complement decomposition in (2).

    python3 -B 05-knowledge/results/session263_spin_flip_20261111/converse/construct.py

Transfer ledger: the Fourier lane's physical annihilator supplies exact
residue uniformity; the converse lane's valuation code supplies independent
amplitudes. The map is bilateral ordinal attachment, the preserved predicate
is arithmetic address and balance, the destroyed data are nontrivial mod-eight
Fourier components, and the retained sidecar is the TOTAL response and paid
source offset. The cheapest hostile is normalization to total1, which makes
all the uniform distributions identical despite their different counts.
