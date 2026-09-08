# Independent referee: fixed supports with unbounded hidden multiplicities

**Status: PASS — independent analytic audit and exact producer replay.**
The strongest fixed-source version is accepted without mathematical repair.
One surface W2, one value of h, the shear parameter lambda=1, three target
supports, and the complete pointed unit-generated Weyl module all remain
fixed. The ambient torsion multiplicities, intrinsic rational-pair degree,
and number of punctures of an ordinary original-source fibre are unbounded.
This referee is an independent analytic audit, not a second finite engine.

## 1. Exact target and inherited scope

The target is `continuing13_20260908_fixed_support_iteration`, produced in
`C:/w/continuing13_20260908_iteration/`. Its final byte pins and finite replay
are recorded in Section 7. The main theorem uses the other fixed point eta
to keep h fixed for every k. The choice of a root in the next preimage level,
which would vary h, is only a secondary variant.

The direct supplier is
[the global mutation theorem](continuing13_20260908_quintic_mutation.md),
report SHA256
`a91d2aca021335eeb182d637fc61050486ce557b021a64fcdfeb7b55cf348c4a`.
Its independently paid full charts, fibre decomposition, rational constants,
labelled principal parts and pair degree are audited in
[the combined supplier/mutation referee](continuing13_20260908_quintic_mutation_audit.md),
SHA256
`cf89c65dd95a9a3ac8d0b10591b724af27aa68777953cf0ff6e18612f9c0ba28`.
Those hypotheses are checked here for every member of the iteration family;
no conclusion is inferred merely from a rational Jacobian identity.

All response modules in this audit use the original ring C[x,t]. Global
submersion refers to both charts of W2. Fibre component and puncture counts
refer to the original affine source. The pointed unit-generated module is
distinguished from its larger ambient torsion module. The inherited
component-principal-part and canonical-connection suppliers, THM-3412
[Hamiltonian principal-part differential and Pruefer torsion arms](../../01-canon/theorems/THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md)
and THM-3770
[vertical principal-part equalizer and log-canonical dressing gate](../../01-canon/theorems/THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate.md),
retain the regular component before quotienting by the common diagonal.

The closest earlier mechanism is the fixed pointed module with hidden arms
in `continuing12_20260908_fixed_unit_hidden_arms.md`. Here the extra retained
coordinate is the entire target support set. The operation is compositional
iteration followed by a square root of the derivative and the mutation
supplier. The missing coordinate is the number of labelled components
over each retained value. The decisive hostile is a returning critical
point, which destroys squarefreeness despite finite critical-value support.

## 2. The critical orbit and the fixed source point

Fix zeta and c with zeta^2+zeta+1=0 and c^2=zeta-1, and put P(z)=z^3+c.
Both c and zeta are nonzero and zeta is not one. Direct substitution gives

    P(0)=c, P(c)=zeta*c, P(zeta*c)=zeta*c.

Thus a positive iterate of zero is never zero. This is the exact property
used in every squarefreeness and disjointness argument below.

Independently expanding the factorization gives

    P(z)-z=(z-zeta*c)(z^2+zeta*c*z-zeta^2).

Choose either root eta of the quadratic. Its constant coefficient is
nonzero; its discriminant is 1+3*zeta^2, which is nonzero. At z=zeta*c
the quadratic equals 2-3*zeta^2, also nonzero. Therefore eta is nonzero
and distinct from zeta*c. Since P(c)=zeta*c differs from c, eta is also
distinct from c. Both choices satisfy P(eta)=eta. In particular

    h=-eta/2, lambda=1

are legitimate fixed source parameters, independent of the integer k.

The uniform lambda condition requires a separate arithmetic argument.
Eliminating zeta gives c^4+3*c^2+3=0, so c is an algebraic integer.
The equation eta^3-eta+c=0 makes eta integral over the algebraic integers,
hence an algebraic integer. If sigma^2=3 and (sigma*eta)^k=-1, then
eta^(2*k)=1/3^k. A rational algebraic integer is an integer, whereas
1/3^k is not an integer. This contradiction is valid for either choice
of sigma and either admissible eta. It avoids a numerical estimate or
an unproved choice of lambda for infinitely many parameters.

## 3. Squarefreeness and all component counts, for every k

For k>=2 define

    Q_k=P^k-zeta*c,
    R_k=product_(j=0)^(k-1) P^j,
    f_k=sigma^k R_k,
    q_k=(3^k-1)/2.

Here P^0(z)=z. The chain rule gives Q'_k=3^k R_k^2=f_k^2. The critical
orbit gives Q_k(0)=0, and f_k(0)=0. The degree of f_k is q_k>=4.

Set E_j={rho:P^j(rho)=0}, including E_0={0}. If i<j and a point lies
in both levels, then P^(j-i)(0)=0, contradicting the critical orbit.
If P^j had a multiple zero, its derivative formula would force an
earlier iterate of that zero to equal zero, yielding the same
contradiction. Consequently every E_j has 3^j distinct points, the
levels are pairwise disjoint, and f_k is squarefree. In particular

    f'_k(0)=sigma^k*c*(zeta*c)^(k-2) != 0.

At the fixed point eta every iterate is eta, so

    a0_k=f_k(eta)=(sigma*eta)^k != 0,-1.

All the supplier assumptions f squarefree, f(0)=0, f(-2h)!=0 and
lambda*(lambda+f(-2h))!=0 therefore hold for every k. No further
genericity qualification is required.

For rho in E_j with j<k, its target value is

    Q_k(rho)=P^(k-j)(0)-zeta*c.

It equals a=(1-zeta)*c when j=k-1, and zero when j<=k-2. Both groups
are nonempty for k>=2. The complete finite critical-value set of Q_k
is exactly {0,a}, with numbers of critical points

    e_0=(3^(k-1)-1)/2, e_a=3^(k-1).

Every critical point has local mapping degree three because f_k is
squarefree and Q'_k=f_k^2. This also proves, without assuming a finite
experiment, the literal identities

    Q_k-a=(P^(k-1))^3,
    R_(k-1)^3 divides Q_k.

The second quotient has degree
3^k-3*(3^(k-1)-1)/2=(3^k+3)/2=q_k+2. It is coprime to R_(k-1),
and its roots are simple: a further multiple root of Q_k would be
another critical point with value zero, and all such points were
already removed with their exact multiplicity three.

The extra source component E_u has z=eta and value
b=Q_k(eta)=eta-zeta*c. The three numbers 0,a,b are pairwise distinct:
b!=0, a!=0, and a-b=c-eta!=0. No f_k-root maps to b, so the component
counts in A_k=0, grouped by target support, are exactly

    (e_0,e_a,e_b)=((3^(k-1)-1)/2, 3^(k-1), 1).

These are also the ambient arm counts. They are not the counts of all
components of the special fibres: each such fibre includes one further
irreducible regular component and has e_c+1 components. The third value
b is regular for Q_k. None of these three values is a critical value of
the global submersion T_k itself.

## 4. Applying the supplier without discarding a chart or component

Use the fixed h and lambda above and set

    u=x-h, z=u-2h+u^3*t,
    A_k=u*f_k(z), T_k=A_k+Q_k(z), G_k=1/(2*A_k^2).

The already audited identity J(A,Q)=A^3 gives J(T,G)=1. On A!=0,
G is regular and this identity implies a nonzero differential of T.
At a root of f, Q'=0 and dT=dA is nonzero. On the actual source
component u=0, the jet z=eta+u+u^3*t gives

    T_x=a0_k*(1+a0_k) != 0, T_t=0.

On the second chart x=1/r, t=-r^2-r^4*B one has z=rM, where
M=-h^2*(3-h*r)-B*(1-h*r)^3. Thus A is regular since f(0)=0, and
Q(rM) vanishes to order at least three. The tangential derivative of
T on the added divisor is -f'_k(0), which is nonzero. Every point of
both charts has now been covered.

For an arbitrary target c0, the A!=0 part of the source fibre has ring

    C[z, f_k(z)^(-1), (c0-Q_k(z))^(-1)],

by u=(c0-Q_k(z))/f_k(z) and t=(z-u+2h)/u^3. It is nonempty and
irreducible. The original fibre consists of its one irreducible closure
and exactly those A=0 components whose displayed value is c0. Since
the source function is submersive, each fibre is reduced, and distinct
components cannot intersect. This pays the regular component that makes
the number of ambient torsion arms equal to e_c rather than e_c-1.
At other targets the whole fibre is irreducible, so there is no hidden
torsion support.

The same field inverses give C(x,t)=C(T_k)(z). The induced derivation
on z is nonzero, so its rational constant field is exactly C(T_k).
This verifies the constant-field condition of the principal-part
supplier, rather than presuming it from connectedness or genus alone.

## 5. Complete principal parts and the fixed pointed module

At a root component E_rho, the simple zero of f implies
Q(z)-Q(rho)=O(A^3). Thus T-Q(rho)=A+O(A^3), giving the full
negative principal part 1/[2*(T-Q(rho))^2] and no simple pole.

At E_u, impose the actual source jet z=eta+u+u^3*t. With a0=f(eta),

    A=a0*u+f'(eta)*u^2+O(u^3),
    Q(z)-Q(eta)=a0^2*u+a0*f'(eta)*u^2+O(u^3)
               =a0*A+O(A^3).

Hence T-b=(1+a0)*A+O(A^3), and the complete negative principal part
is (1+a0)^2/[2*(T-b)^2], again with no simple pole. The regular
component has principal part zero. Both coefficients are nonzero.

At each support c in {0,a,b}, the ambient labelled quotient therefore
contains one nonzero coefficient vector at pole order two and no other
coefficient direction for the unit. It generates exactly one complete
principal-part tower. The regular component with zero principal part
ensures that this vector does not become the common diagonal class.

Write tau for the same abstract target coordinate in every member.
Chinese remainders at the three pairwise coprime polynomials
(tau-c)^2 project the unit onto its three primary parts. Multiplication
by tau-c and the canonical derivative then generate each full tower.
Consequently its exact scalar annihilator is

    ([tau*(tau-a)*(tau-b)]^2).

The pointed module is the direct sum of three fixed-support towers with
one pure double-pole distinguished vector at each support. Rescaling the
three nonzero coefficient vectors identifies the distinguished units
and commutes with both Weyl actions. This proves equality of the full
pointed unit-generated module type, not just of the scalar annihilator.

The ambient module has e_0+e_a+1=(3^k+1)/2 arms; its difference from
the three visible arms is (3^k-5)/2. The support set and unit therefore
lose an unbounded amount of component multiplicity. The proof does not
identify these differing ambient modules or their embeddings.

## 6. Intrinsic pair degree, punctures, and hostile boundaries

Q_k is monic of degree 3^k, so T_k has original t-degree exactly 3^k,
with leading coefficient u^(3*3^k). The lower-degree A_k term cannot
cancel it. Over C(T_k), the rational function
A_k=(T_k-Q_k(z))/f_k(z) has degree 3^k: numerator and denominator
are coprime over this generic field and have degrees 3^k and q_k.
The tower through A_k^2 gives

    [C(x,t):C(T_k,G_k)]=3^k*2.

Any rational mate of nonzero constant Jacobian is alpha*G_k+R(T_k)
with alpha!=0, because the rational constants have already been paid.
Its embedded pair field is the same, so this degree is intrinsic to
T_k among rational mates.

Outside {0,a,b}, the original fibre has q_k fixed missing f-roots,
3^k moving missing roots of Q_k-c0, and infinity. The moving roots
are simple and avoid the fixed ones. The puncture count is exactly

    q_k+3^k+1=3*q_k+2=(3^(k+1)+1)/2.

For each fixed k, the family of these affine curves is non-isotrivial.
An isomorphism extends to a Mobius transformation of the completions.
There are at least three fixed punctures, and the images of three
chosen fixed punctures into any one reference puncture set give only
finitely many possible transformations. For a fixed transformation the
puncture set, hence the moving root set and target value, is determined.
Thus only finitely many target values can give a prescribed isomorphism
class. The smooth completions remain P1; no varying completed genus is
asserted. The puncture count concerns the original source, not the
fibre after adding the second-chart divisor.

The c=0 hostile has P=z^3, and R_2=z^4 is not squarefree. It refutes
the implication from finite critical-value set alone to admissibility.
For k=1 the chosen Q_1 has nonzero value at zero and does not satisfy
the stated normalization; the sequence correctly starts at k=2.

The secondary E_k-source-point variant is also valid: its preimage
disjointness gives f_k(z0)!=0, its support is {0,a,-zeta*c}, and the
same integrality argument proves lambda=1 admissible. It does not
improve the main fixed-h result and is not used for the invariant
source-parameter assertion.

## 7. Producer inspection, replay, and frozen byte pins

Final producer pins:

    report 1732a019e326d16b03dcde64f424195eed4afebe6fc282bcd3edc17b73243413
    source e65e4ae328c4720510662d43e8dad2764521489fccbd858489a0bbbb2c0e0f7d
    output 345637678c2eec05006cd3bbda64bee6d9b0848e4c999dfc1f84731f5f682891
    cert   dd3b7b9141f7b9dc971ea09fe1496b43be74485e721072d286773dcee90f9d63

The final source was read in full. It uses always-active raising gates,
including under Python -O. The source is a finite exact control of k=2,3,4
over the irreducible field
Q[c]/(c^4+3*c^2+3), not a proof by finite extrapolation. The independent
all-k argument is Sections 2--6 above. No producer is imported into a
new mathematical engine, and no numerical root is used in this audit.

The finite controls cover all preimage levels through degree 81, their
squarefreeness and pairwise coprimality, the squared derivative, exact
critical groups and cubic multiplicities, the additional simple roots,
both fixed-point choices through quadratic remainders and gcds, lambda=1
admissibility, the secondary moving-h choice, and the stated hostiles.
There are 115 gates per execution. Separate normal and optimized referee
replays in private outside-repository folders both reproduced the frozen
producer stdout and regenerated certificate byte for byte, with raw LF
output and certificate bytes. The replay copied the frozen source and
executed it; it is not counted as a second independently derived engine.

Reproduction after repository integration:

    python 04-computation/continuing13_20260908_fixed_support_iteration.py
    python -O 04-computation/continuing13_20260908_fixed_support_iteration.py

Only this report is a new referee artifact. Its private replay helper and
copied producer inputs are operational files, not additional proof engines.

The accepted conclusion is a family of global polynomial first functions
with rational mates that have poles and nontrivial unit response. There
is no polynomial mate or general Jacobian-conjecture conclusion. The
general LRC(14) problem is also unaffected by this orthogonal result.
