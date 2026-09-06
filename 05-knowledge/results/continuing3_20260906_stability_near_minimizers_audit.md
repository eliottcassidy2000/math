# Independent audit of the complete near-minimizer classification

**Status: INDEPENDENT ANALYTIC AUDIT PASS; exact independent controls.**
The new five-way equivalence, absence of nondegeneracy assumptions, local
two-atom profile, quantitative moment bridge, and mixed-dust obstruction
are accepted. No mathematical correction is required. The constant K3 is
inherited from THM-4454, not reproved or claimed new by this continuation.

Audited theorem: `continuing3_20260906_stability_near_minimizers.md`.
The mathematical input is **THM-4454**,
`01-canon/theorems/THM-4454-sharp-global-signed-root-duplication-stability.md`,
including its actual regional inequalities and strict curvature/boundary
bounds. I read the entire candidate proof before constructing the separate
verifier. The producer's bounded controls are not used as proof inputs.

## 1. Scope and compactness are correctly typed

Every finite real list has p1=p2=1 and E>0. Thus all entries have absolute
value at most one. The distances are to permutations of the indicated
zero-padded vectors; selecting the largest positive entries is correct.
Arbitrary negative entries and varying list lengths remain allowed. The
two-atom distance is strictly positive on this class, since zero distance
would contradict p1=1. No uniform first absolute moment bound is assumed.

If R tends to K3, then F_actual=(3E/4)d2^2(R-K3) tends to zero because
both multiplying factors are bounded. This is only a necessary step.
The proof then uses compactness of the two envelope coordinates, not an
unjustified compactness claim about unweighted root lists or dust mass.

The b<=z envelope has only the two inherited zeros (z,z) and (1,0).
In the b>=z region, the precise zero-set extension is valid: when the
remaining square mass c0^2 lies strictly between zero and1/3, both interval
endpoints are strictly positive by the inherited bounds, as is their
interior by concavity. At c0=0 only the equal two-atom endpoint vanishes;
at c0=z only the triple endpoint is present. This exhausts all possible
macro-coordinate limits and does not mistake c0 for an actual third root.

## 2. Both quotient degeneracies are excluded uniformly

At the one-atom alternative, q is the entire tail square mass. Its signed
cubic moment is bounded in absolute value by q^(3/2), and its fourth moment
by q^2. Direct division by q gives J->1, independently of dimension or
separate tail first masses. The resulting R limit is sqrt2-2/3>K3.

At the two-atom alternative, put v0=(a-b)^2 and q=1-a^2-b^2. The exact
top moments in the report are correct. The referee independently substitutes
q=eta*A, v0=eta*B and differentiates at eta=0, retaining both direction
coefficients. It recovers the same objective and distance linear terms.
The smooth remainder is uniform over the compact nonnegative directions
A+B=1. Signed tail moments contribute at most q^(3/2) and q^2; after division
by q+v0/2 these are o(1), with no scale relation imposed between q and v0.
The limit energy is1/4. Hence the entire profile is the stated convex
mixture of K_two and K_dust, up to a uniform o(1). Both endpoint constants
exceed K3, so this alternative is impossible for a minimizing sequence.

These checks are what permit a classification without separately assuming
the energy and distance denominators are bounded away from zero. Vanishing
unnormalized defect alone would not establish the theorem.

## 3. The discarded tail information is actually recovered

After both exclusions, every subsequential macro limit has a=b=z. The
local secant loss is a sum of nonnegative terms near that point, including
on the b>z side; its applicability here follows from the positive derivative
of f_C on the whole local interval, not a new global extension of the
b<=z envelope. Both the actual objective and the continuous secant envelope
tend to zero, so the loss does too.

The referee checks the exact individual loss identities. A negative root
pays at least f_C(b) times its square mass. A smaller positive root pays
a fixed positive multiple of its square mass times its distance below b.
The limiting taxes are strictly positive. Negative square mass therefore
vanishes. The remaining positive square mass after a,b tends to1/3; its
largest root c must approach b. This yields exactly three macroscopic
positive roots and zero remaining square mass. It is stronger than knowing
only the first two coordinates of the envelope.

The converse follows from moment continuity with the dimension-free tail
bounds. The limit energy is1/3 and the two-atom distance has a strictly
positive limit. Direct substitution gives the inherited optimum K3.

## 4. The other equivalences and the effective sidecar

M=p4-2zp3+1/3 is exactly sum r_i^2(r_i-z)^2. The band proof retains
integer multiplicity: outside [7z/8,9z/8], square mass is at most192M;
four inside roots would exceed total square mass one, while two are
insufficient when M<5/6144. Matching the resulting three roots proves
d3^2<=9600M/49. The global reverse estimate M<=(1+z)^2d3^2 is valid because
all entries have absolute value at most one. These prove the M equivalence;
the separate two-moment equivalence then follows as stated.

For complex disks, each dust entry eventually has small enough magnitude
to use the analytic logarithm near one. The remainder after its signed
linear term is bounded by the disk radius squared times total dust square
mass. No separate dust l1 bound is needed. The exact signed sum tends to
1-sqrt3, giving the stated entire-product limit. Conversely locally uniform
convergence determines coefficients near zero, and Newton identities recover
p3,p4. The referee independently reconstructs eight limiting moments from
the proposed product's ordinary coefficients. The qualitative uniform
modulus follows by the standard sequential contradiction within the already
proved varying-dimension classification.

## 5. The hostile repairs the first-moment statement

The general mixed-dust family's radical normalization is correct for every
L>=2,c>=0. Its d is positive and its scale is positive. The three nonzero
main roots ensure positive energy after p2 normalization. If c=o(sqrt L),
then d<3+c, so the displayed normalization forces the scale to sqrt3 and
all dust square mass to zero. This proves the family is near-minimizing.
The c=1 instance retains positive dust first mass; c=L^(1/4) makes both
separate masses unbounded. Only their net difference is rigid.

The final Cauchy bound means liminf N*q3 >=(sqrt3-1)^2, not convergence
of N*q3 to that number. It proves diverging dust multiplicity and no unique
schedule. The report does not assert that lost stronger conclusion.

## 6. Exact verification

The standalone referee imports no producer or repository mathematical
implementation. It checks independent directional derivatives, the
one-atom expansion, all constants and local secant identities, the band
counting bounds, eight Newton moments of the limiting entire product,
the all-parameter radical normalization and seven literal exact controls.
The universal limiting and compactness steps are the audited arguments
above; finite arithmetic is not a substitute for them.

    python -B 04-computation/continuing3_20260906_stability_near_minimizers_audit.py
    python -B -O 04-computation/continuing3_20260906_stability_near_minimizers_audit.py

The referee passes 72 always-active exact gates, with byte-identical normal
and optimized LF output. After finishing the independent path I read the
complete producer source. Its literal product controls, outward rational
interval arithmetic and all three two-atom scale regimes are consistent
with the analytic scope above. The producer's 477 checks are supporting
controls, not an inference of the theorem from a census.

Frozen review identity:

    producer source b5a2ce85df71837e5b73b3609feee0d9fcfaa14f425e22b81b563a7a3574c598
    producer output 3a0d5f4c830b62c55c48a9458f315f3ea601d68aeec40e9c634980ad2e6f8337
    producer report before promotion 5378d5699f85175074eee1037c6515f79ee92c8dfde55c3494ee722bbf51e086
    referee source ac730d76df70cad414b1fbc6b0882b9ca8c2996f32e8847bbbd475c1c938c67f
    referee output 7172844909f1602e68095bc283723c3a2c0b2fc93c8c30945a46ddab948713c5

The complete proof is accepted for promotion in this stated scope. General
LRC14 and actual Laurent noncancellation are unaffected by this
classification of the specified signed-root quotient.
