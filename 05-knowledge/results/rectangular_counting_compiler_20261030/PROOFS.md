# Proofs and exact recurrences

**R1--R4: PROVED, kernel checked by ReducedEngine.lean. R5--R7:
PROVED elementary deductions from those identities and the credited
scaled-observer law.** The rectangle identity itself belongs to the
pinned external engine, not to this packet. No global novelty claim.

## R1. Erase unused coordinates at the actual child interface

Let `Q_G=(Fin G x Fin G)^2`. Let `Rect_G` consist of quadruples
`(a,b,c,d)` with `a<b` and `c<d`. Write

    R_G=|Rect_G|=binom(G,2)^2,
    A_G=G^4-R_G, t=ceil(lambda/3), q=t+1.

The original finite choice space is `Phi=(Q_G -> Fin q)`. Its restriction
map is `rho:Phi -> (Rect_G -> Fin q)`. Extend a rectangle assignment by
zero on the complement; this gives a section `j` with `rho j=id`.
Every fibre of rho has exactly `q^A_G` members, including q=1 and the
empty-rectangle cases G=0,1. The dense representative map is injective.

In `FactoringAux.EngineGroup.nextT`, the chosen row/column network,
the comparison of their total sizes, first-good grouping, coarsening,
and final clique construction all use phi through `restrictC phi`.
Consequently

    rho(phi)=rho(psi) => nextT(T,phi)=nextT(T,psi).

This is literal equality of child tables; it is stronger than equality
of their scalar values. The Lean proof unfolds `nextT` and rewrites
the equality of restrictions. Its canonical representative uses the
existing `phiEquiv` encoding, so it connects directly to the source's
actual constructor rather than to a replacement toy model.

For any later function C on these child tables, the outputs under C
are equal. This licenses compositions and named-port observations
which receive the table. It does not license an observer also given
the raw discarded phi or its code. For example on G=2, two codes with
the same rectangle bit can have different values at `(0,0,0,0)`.
That observer explicitly defeats the proposed erasure.

Lean names: `restrict_extend`, `choiceCode_injective`,
`nextT_same_on_fibres`, in namespace `Factoring.Aux.RectangleOnly`.

## R2. One-step normalization without ring division

For every modulus m, K>=2, G=3K^8, and correctly typed clique T,
the external `value_nextT` identifies each child's value with `Z_iota`.
The inherited `FactoringAux.Level.sum_Zι` says

    sum_iota Z_iota(T)=2^(t R_G) value(T).

Combining these facts with R1 gives the reduced step

    sum_iota value(nextT(T,j(iota)))=2^(t R_G) value(T).       (1)

All identities hold in `ZMod m`, including even/composite m and m=0;
there is no division or assumption that q is a unit. Internal signed
terms retain their signs. Grouping is the source's grouping, with its
original guards and width bound.

The original summation over Phi instead multiplies (1) by `q^A_G`.
That remains its correct semantics. Choosing one representative changes
the producer; it does not recover the normalized result from a modular
value that has already lost information.

Lean names: `next_value`, `sum_next`, `choice_card`.

## R3. All-depth recurrences and exact savings in leaves

Use the source's parameters, all independent of the actual table:

    K_0=K0>=2, K_(d+1)=9 K_d^16, G_d=3 K_d^8,
    U(G,L)=floor((12G^3 L+11G^4 L+24G^4)/12),
    L'=floor((K_d U(G_d,L)+G_d^4 L)/K_d^3),
    t(L)=floor((L+2)/3), q(L)=t(L)+1.

At depth d+1 set R=R_(G_d), A=A_(G_d). Define

    E_0(L)=0,       E_(d+1)(L)=t(L) R+E_d(L'),
    P_0(L)=1,       P_(d+1)(L)=q(L)^A P_d(L'),
    Hnew_0(L)=1,    Hnew_(d+1)(L)=q(L)^R Hnew_d(L'),
    Hold_0(L)=1,    Hold_(d+1)(L)=q(L)^(G_d^4) Hold_d(L').

The reduced value Vnew is the nested sum over rectangle assignments,
with the same bottom clique value and the same child tables as R1.
Induction using (1) and distributivity proves, for every allowed T,

    Vnew_d(T)=2^E_d value(T),
    deepMult_d=P_d 2^E_d,
    Vold_d(T)=P_d Vnew_d(T),
    Hold_d=P_d Hnew_d.                                      (2)

For the last identity partition Q_G into Rect_G and its complement,
use `G^4=A+R`, and factor q^(A+R). The Lean file checks the partition
identity directly from finite-set cardinalities, not by importing an
unproved closed cardinality formula. The elementary strict-pair count
then gives R=binom(G,2)^2.

For G>=2 the ratio of the two exponents in the **per-level branch count** is

    R/G^4 = (G-1)^2/(4G^2),                                 (3)

approaching 1/4. This is neither a factor-four running-time improvement
nor a new asymptotic complexity exponent. Branches are exponentials,
and encoding, child construction, grouping search, base evaluation,
and the schedule's uniform compilation still cost time. q=1 has no
branch savings despite (3). Any practical claim must measure the
implemented algorithm in its actual regime.

Lean names: `value_eq`, `mult_factorization`, `old_value_eq_padding`,
`quad_partition`, `leaves_eq_padding`.

## R4. Named boundary assignments survive

For an arbitrary type B of named assignments and a family of same-size
cliques `T:B -> Clique`, R3 holds pointwise, hence as an equality of
whole functions:

    (b -> Vnew(T_b)) = (b -> 2^E value(T_b)).                  (4)

The same E applies because the schedule depends on size/width, not b.
Thus every fixed signed linear combination of the boundary responses
has the same scaling, and tensor contraction along declared shared
ports may be performed after this replacement. If two independent
pieces of scales 2^E and 2^F are glued by a bilinear contraction, the
resulting scale is 2^(E+F). This follows by bilinearity, not by assuming
the two pieces' port assignments are independent.

`boundary_values` formalizes (4). The tensor-contraction sentence is
an elementary consequence, not a separate formalized network compiler.
Conditional branches selecting a different schedule for each b require
retaining E(b); they cannot use the constant E version unmodified.

## R5. Exact information gained by changing the producer

Let c be an integer actual count and let a finite bank of moduli have
lcm M. The inherited exact law is that observations of Dc retain
exactly `c mod M/gcd(M,D)`. Apply it to Dold=P2^E and Dnew=2^E:

    Qold=M/gcd(M,P2^E),  Qnew=M/gcd(M,2^E),
    Qnew/Qold=gcd(M,P2^E)/gcd(M,2^E).                        (5)

If all bank moduli are odd, this gain is gcd(M,P), with Qnew=M.
For known `0<=c<=B`, zero is distinguished from all positive values
iff Q>B. This height guard remains necessary after normalization.

The source's generic G=2,lambda=4 example has q=3, R=1, A=15:
Dold=3^15*4, Dnew=4. Modulo3 the former maps both c=0 and c=1 to0,
and the latter separates them. No post-processing of the old single
reading can do so. This small case explains the mechanism, but is
explicitly not a valid point of the uniform K0=48 recursion.

At the actual uniform depth thresholds n=(Kb(48,d))^40, d=1,2,
the schedule audit gives the gains459 and513 for the specified bank.
The respective first q values have v2(q)=2 and odd-part bit lengths
3520 and56382. Thus the old odd-content obstruction occurs on the
real schedule, not just on a hypothetical G=2 grid.

## R6. Prefix-sensitive observation after normalization

The already-proved Session112 bound applies now to the normalized
semantic engine: for a genuine count `0<=c<=2^b`, all odd moduli
at most2b+8 suffice for zero detection, independently of E. For a
length-n binary witness with h fixed primary bits, b=n-h. At a full
prefix c is0 or1 and modulus3 alone suffices.

The word "genuine" matters. Unique auxiliary completion, or an
explicit multiplicity correction, is needed before substituting a
bound on primary witnesses for a bound on network assignments.
The paper's multiplication network supplies a specific unique-completion
argument; applying the mechanism to other CSPs requires its own proof.
This corollary currently concerns the reduced semantic engine. The
new loop/decoder and zero-bank implementation have not been compiled
into the external TM2 algorithm.

## R7. General weighted-fibre rule and its failure boundary

For finite X,Y, a surjection rho:X->Y and a table-valued map F with
F=Fbar composed with rho, every linear observer ell satisfies

    sum_(x in X) ell(F(x))
      = sum_(y in Y) |rho^-1(y)| ell(Fbar(y)).                (6)

This follows by partitioning X into fibres. If all fibres have the
same size p, the common factor p may be removed at the producer when
the intended target is the normalized sum. If sizes vary, they must
remain as weights. A minimal witness is fibre sizes2,1 with values3,5:
the original total11 differs from the unweighted quotient total8.
If the fibres represent actual counted objects rather than artificial
padding, even a constant p cannot be discarded from the intended
count without changing the problem.

The strongest survivor is therefore a weighted, contextual quotient.
Raw cardinality, symmetry, norm, or agreement of one scalar observer
does not establish R1's stronger hypothesis. The 81/27 tournament pair
is retained as the canonical example of this distinction.
