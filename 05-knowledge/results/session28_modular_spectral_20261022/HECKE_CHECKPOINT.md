# Rank-88 theta: modular but not a scalar Hecke eigenform

Status: CITED modular-form framework plus independently FINITE-EXACT Hecke arithmetic. This does not independently certify the minimum of an88-dimensional lattice.

For an extremal positive definite even unimodular lattice of rank88, the theta series has weight44 and starts1+0q+0q²+0q³. The four-dimensional space M44(SL2(Z)) has basis Delta^j E4^(11-3j), j=0,1,2,3. These initial conditions therefore force

    F=E4^11-2640 E4^8 Delta+1813680 E4^5 Delta²-244992000 E4² Delta³
     =1+168498000 q^4+2480127344640 q^5+6298997202432000 q^6+...

Here theta uses q^(norm²/2). The q^4 coefficient counts minimal norm8 vectors.

For the usual weight44 Hecke operator,

    [q^n]T_p F = a_(pn)+p^43 a_(n/p),

with the second term zero if p does not divide n. Since a_2=0 but [q²]T_2F=a_4=168498000, F is not a T2 eigenform. An independently computed Krylov matrix(F,T2F,T2²F,T2³F) has rank4. Consequently the smallest T2-invariant space containing F is all M44, not a one-dimensional scalar eigenline.

The characteristic polynomial factors as

    (x-(1+2^43))
    *(x³+2209944x²-15663522502656x-19976984434430705664).

The cubic is irreducible over Q and has three distinct real roots. This is compatible with the classical self-adjointness of T2 on the cusp space under the Petersson inner product. It is not evidence that all operators called Hecke, including U_p at primes dividing level, must have real spectrum.

The independent script uses a DIFFERENT basis, E4^(11-3j)E6^(2j), constructs Delta directly from q product, confirms E4³-E6²=1728Delta, and verifies T2 expansions throughq18. Its cyclic determinant in that basis is

    -218792213161134562973631210691580078125/139968 !=0.

The separate producer computation uses Delta/E4 basis and checks T2 and T3, including commutation. Agreement between distinct bases is an audit rather than two invocations of one script.

Primary source for the current existence claim: Ryotaro Sakamoto, *Two Extremal Even Unimodular Lattices of Rank88*, arXiv:2609.15633v2,27September2026, https://arxiv.org/html/2609.15633v2 . The paper reports two constructions over Q(sqrt(-23)), using Hermitian ranks4 and11 then restriction of scalars:88=2*4*11. Its existence/minimum claim is CITED here, not independently re-certified by our small modular computation. The detailed session report will distinguish this field from the separate D=-163 gluing construction.

Classical modular-form ingredients (graded ring, q expansions, Hecke operators) are documented in J.S.Milne, *Modular Functions and Modular Forms*, https://www.jmilne.org/math/CourseNotes/MF.pdf . The report's elementary noneigenform witness needs only the displayed Hecke coefficient rule; no numerical eigenvalue fitting is used.
