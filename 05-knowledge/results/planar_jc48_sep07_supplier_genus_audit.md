# Independent audit: the actual finite supplier has generic genus27

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS.** This accepts the
[primary theorem](planar_jc48_sep07_supplier_genus.md), including every
polynomial A of degree m>=4, every nonzero constant b, geometric integrality,
and nonrationality of all nonzero scalar times of a nonconstant outer
Hamiltonian. It does not classify every carrier lift or different-invariant
composition. JC(2) remains OPEN.

I independently proposed the degree-four Newton polygon before receiving
the producer's Riemann--Hurwitz proof. Its vertices are (0,0),(9,3),(6,5),
(2,7), its area is59/2 and its boundary lattice count is7, giving27 interior
lattice points. This is a separate geometric consistency control; the
accepted proof below uses actual ramification and does not infer the genus
from a polygon without checking its hypotheses.

Read literally, I=p2 y3(p3-y2)(A+b y2) has y-degree seven, derivative
p2 y2[3p3 A+5(bp3-A)y2-7b y4], and two nonzero squared critical values
V1,V2. Their quadratic algebra gives an independent short norm calculation.
For its root z, with X=bp3,

    Norm(z)=-3p3 A/(7b),
    Norm(p3-z)=2p3(A+X)/(7b),
    Norm(A+bz)=2A(A+X)/7.

Multiplying p8 times their respective powers3,2,2 yields precisely
-432 p23 A5(A+X)4/(7^7 b5). Since m>3, this has degree23+9m. Independently,
one critical root has z asymptotic to -5A/(7b), giving squared critical
value of degree4+7m and leading coefficient -12500 a_m^7/(7^7 b5).
The other has z asymptotic to3p3/5, giving degree19+2m. These degrees are
different when m>3. Thus the critical values separate, and their trace
has exactly the degree used in the proof. The source's full coefficient
trace and separation identities provide additional literal checks.

The actual branch polynomial is (c2-V1)(c2-V2), not a quotient obtained
by imposing the involution y->-y on I=c. It is monic in transcendental c,
so its zeros avoid every fixed algebraic p-value where a critical root
coalesces or a critical value vanishes. At a zero, exactly one actual
critical branch has I=c, and its squared-value derivative is2c I_p.
A repeated branch root would therefore be a critical point of I at the
transcendental value c. The finite-critical-values lemma rules that out
even with positive-dimensional critical loci. Consequently all23+9m
finite nonzero branch values are simple with index two. The literal
seventh-degree discriminant in the verifier confirms the full branch
polynomial and the excluded p=0 and c=0 factors.

At p=0, substitution p=r7,y=r^-2 Y gives initial equation -bY7-c.
Its seven roots are simple and the seventh-root deck action is transitive.
This is one index-seven place and contributes six. The same fact proves
geometric integrality after extending constants: a degree-seven component
already exhausts the total y-degree. Thus applying geometric
Riemann--Hurwitz is justified, with no assumed generic irreducibility.

The infinity groups are disjoint for m>3. The D-roots have y of order
p^(3/2), the A+b y2 roots have order p^(m/2), and the remaining three roots
have order p^(-(m+5)/3). Their simple initial equations and deck actions
give respective defects1, m mod2, and3-gcd(3,m+5). The degrees2+2+3
exhaust all seven roots. This independently reconstructs

    2g-2=-14+(23+9m)+6+1+(m mod2)+3-gcd(3,m+5),

and the displayed genus formula. At the actual quartic supplier this
is -14+59+6+1=52, hence genus27. Its rational coefficients in the source
agree exactly with the completed-response supplier; both its quartic
leading coefficient and b are nonzero.

For the degree-three hostile A=-bp3, I=-bp2 y3 D2. The birational change
p=h2/w,y=h3/w gives h25=(-c/b)w11/(1-w)2. The valuations at0,1,infinity
are11,-2,-9, each coprime to25. The resulting three totally ramified
places give genus12, not the extrapolated24. This checks a genuine
failure boundary rather than merely excluding zero coefficients.

The scalar-time conclusion retains I itself, whose logarithmic first
coefficient is W=s17(A(s2)+b s6), nonzero since deg A>=4. The inherited
actual completion and fixed-input comparison give nonzero first
p-displacement at every positive multiple of a nonzero scalar time.
A rational specialization of both coordinate images would induce a
nonconstant selfmap of the smooth proper generic curve. Its genus is
at least27, so that map would be an automorphism in a finite geometric
automorphism group, contradicting those iterates. The external boundary
parameter sigma in j(sigma) is explicitly distinct from the source's
logarithmic variable s: these are scalar times in the coefficient field,
not source-dependent times. Zero time remains exempt. The finite row15
identity is unaffected by this nonrationality obstruction.

I read the complete standalone source and replayed it normally and under
optimization. Both130-gate outputs agree byte for byte with the frozen
401-byte output. Ten named controls retain the exact supplier, repeated
and vanishing A factors, degrees4..10 and nonunit b. Each finite-field
squarefree control preserves degree and avoids denominator primes; it
supports its named specialization only. The universal proof does not
extrapolate from these finite controls. The separate local chart universe
covers degrees4..12, with always-active checks under both runtime modes.

    python3 04-computation/planar_jc48_sep07_supplier_genus.py
    python3 -O 04-computation/planar_jc48_sep07_supplier_genus.py

Frozen SHA256:

    source a23ac00b69e027efca1e720b9af4ffaf285203551ebd287f05cc28a6c534c83b
    output 1d70b64d691d5aedd69246ed3d492fd0fe66060fa4dbaeeba355f439e0680b62

No mathematical correction was needed; the independent comparison added
only the scalar/source parameter distinction before status promotion.
