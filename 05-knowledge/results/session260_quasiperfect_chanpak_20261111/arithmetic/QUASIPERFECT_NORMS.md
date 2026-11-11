# Quasiperfect norm coordinates and actual cyclic tournament constructors

**Status: PROVED elementary interface theorems; FINITE-EXACT certificates explicitly scoped. Quasiperfect existence and Chan–Pak remain OPEN.** No external priority claim.

## 0. Inheritance, target, and source map

The recovered Chowla note `05-knowledge/results/ogg_triangular_chowla_20260927.md`, §5, distinguishes quasiperfect numbers sigma(n)=2n+1 from the excess-two sources in this session. It inherits the classical necessary condition that a quasiperfect n is an odd square m². We use that condition as inherited, without re-auditing its literature. Its Pell/nodal-curve material concerns a different predicate; no transfer from those integer points is assumed here.

THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`, already supplies the general positive C3 path-cover kernel. THM-3134, `01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`, retains the endpoint sidecar; THM-3202, `01-canon/theorems/THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction.md`, treats larger moving jets. The one-singleton identity below is a transparent specialization of that machinery, not a newly discovered general substitution law. The new task is to identify its exact quadratic-norm interface, prove the representation-switch obstruction, and give an actual finite target compiler.

Board: divisor-shell fixed point / ordered two-path covers / norm representations / cut-derived size bound / small-factor switching / actual graph output. Closest mechanism: C3 substitution. Hostile: a norm representation need not be a physical path-cover profile. Least-used sidecar: Q2 alongside H. The goal is constructive spectrum information, not a numerological argument for or against quasiperfect existence.

## 1. The physical quadratic form

For a nonempty tournament A let h=H(A), and define

    Q2(A)=sum_(S⊆V(A)) H(A[S]) H(A[V(A)\S]),  H(empty)=1.

Thus Q2 counts ordered spanning two-path covers allowing empty slots. It is even: exchange the slots, an involution without fixed points for a nonempty vertex set. Write

    x(A)=Q2(A)−H(A),  y(A)=H(A).

Both x,y are positive odd integers, using the inherited parity of tournament Hamilton-path counts.

**Graph identity.** For arbitrary nonempty tournaments A,B, orient the cyclic substitution A→B→z→A. Then

    H(C3[A,B,1])=x(A)x(B)+2y(A)y(B).                 (1)

Proof. If z is first or last, there are H(A)H(B) paths in either case. If z is interior, a path has the unique form A_left B_left z A_right B_right, where the displayed blocks are directed paths, A_right and B_left are nonempty, and the other blocks may be empty. The choices on A number Q2(A)−H(A), and the choices on B number Q2(B)−H(B). They are independent. This proves (1) by a literal path bijection.

In particular the equal-module construction has

    N(A):=H(C3[A,A,1])=x(A)²+2y(A)².                (2)

Every such output is3 modulo8. These are actual strong tournaments of order2|A|+1, even if A itself is not strong: the external directed triangle connects the modules.

For A=TT_t, (x,y)=(2^t−1,1), giving the infinite explicit family

    N_t=(2^t−1)²+2.

This is an infinite norm-count family. It is not claimed to be an infinite family of numbers1+2m² with integral m.

## 2. The missing graph coordinate has a strong inequality

**Cut bound.** If A has n≥1 vertices, then

    Q2(A)≥(n+1)H(A),    x(A)≥n y(A).               (3)

Proof. Cut each Hamilton path at each of its n+1 gaps, allowing either empty end. The result is an ordered two-path cover. It determines both the original path (concatenate the two lists) and its cut position, so this map is injective. Arbitrary covers need not concatenate to a directed path, which explains the possible strict inequality.

Consequences:

- x=1 forces n=y=1 and A to be the singleton.
- Identifying the distinguished representation1+2m² with the physical coordinates(1,H(A)) is impossible for m>1.
- The inequality is stronger than x≥y: a proposed norm profile has an explicit graph-size budget n≤floor(x/y).
- A cover-free scalar count H cannot recover this budget; Q2 is indispensable.

An additional content guard is

    gcd(x,y)=gcd(Q2,H).

Its square divides N(A). More structurally, every automorphism subgroup G of A acts freely on Hamilton paths and on ordered two-path covers: fixing either ordered object fixes every listed vertex. Therefore |G| divides both H and Q2 and |G|² divides N(A). For a vertex-transitive A, n divides |Aut(A)|, so n² divides N(A). This is a constructor restriction, not a theorem that every tournament realization of N has such divisibility.

## 3. What the quasiperfect equation contributes

A hypothetical quasiperfect n=m² with m odd satisfies

    sigma(m²)=1+2m²=Norm(1+m sqrt(-2)).              (4)

Equation(2) is the same quadratic form, but with a different and graph-constrained representation. The map between the two subjects preserves the integer norm, not its coordinates, prime-power factors, or existence predicate.

Because m is odd, the value in(4) is3 modulo16. A physical alternative representation must therefore satisfy x²≡1 modulo16, equivalently x≡±1 modulo8. The graph value(2) by itself may instead be11 modulo16.

For an odd prime p≡5 or7 modulo8, −2 is a quadratic nonresidue modulo p. If p divides x²+2y², then p divides x and y; iterating proves its valuation in the norm is even. Thus both(2) and(4) have even valuation at every such inert prime. This is the familiar elementary norm obstruction, applied to the actual graph coordinates. It does not distinguish the two forms by itself.

On the divisor side, write m=product p_i^e_i. Then

    product sigma(p_i^(2e_i))=1+2m².

For each inert prime l≡5,7(mod8), the sum of v_l(sigma(p_i^(2e_i))) is even. Individual divisor-sum factors need not have even valuation: cancellations of parity can occur across distinct prime-power factors. The retained sidecar is therefore the factor-indexed valuation-parity vector, not merely whether the total is a norm. Nothing here authorizes factoring each divisor-sum block separately as a norm.

### Prime-target barrier, with its exact limitation

**Theorem.** If P=1+2m² is prime and m>1, then P is not N(A) for any nonempty tournament A.

Proof. The ring Z[sqrt(-2)] is norm-Euclidean: approximate the two rational coordinates of a quotient by nearest integers; the remainder has relative norm at most1/4+2/4=3/4<1. Hence unique factorization holds. For a rational prime P represented by this norm, every element of norm P is associated to1+m sqrt(-2) or its conjugate. Units are±1. Consequently a positive representation P=x²+2y² must have(x,y)=(1,m). The cut bound(3) rules this out when m>1. QED.

**This does not obstruct a genuine quasiperfect source.** An odd prime power is deficient, since sigma(p^k)/p^k<p/(p−1)≤3/2. A quasiperfect integer therefore cannot be a prime power, and sigma(n) is already composite by multiplicativity. In fact an odd number with at most two distinct prime factors has abundancy strictly below(3/2)(5/4)<2. The prime barrier distinguishes norm representations in this restricted graph grammar, not quasiperfect existence.

## 4. Factor-three switching: positive constructions and a complete boundary

A different norm representation can repair the forbidden coordinate x=1. The identities

    (4y−3)²+2y²=1+2(3y−2)²,
    (4y+3)²+2y²=1+2(3y+2)²                       (5)

come from changing which conjugate factor of norm3 is used in Z[sqrt(-2)]. They are exact representation switches, not graph operations.

To make(5) physical requires respectively

    Q2(A)=5H(A)−3  or  Q2(A)=5H(A)+3.             (6)

**Theorem: complete physical factor-three switch classification.** Among all nonempty tournaments, (6) occurs only at the profiles

    singleton: H1,Q2=2;
    TT3:       H1,Q2=8;
    C3:        H3,Q2=12.

Thus this particular switching mechanism gives exactly m=1,5,7, and graph outputs3,51,99.

Proof. For H≥3, (3) and (6) imply n≤4+3/H≤5. The finite universe is therefore all labelled tournaments through order5, not an arbitrary cutoff. Exact exhaustive endpoint DP in norm_check.py gives only the displayed C3 case among H≥3; its certificate reports all qualifying profiles and their labelled multiplicities. For H=1 the tournament is transitive (the inherited odd-cycle expansion gives H≥3 in the presence of any triangle). Hence Q2=2^n and(6) reduces to2^n=2 or8, yielding n1 or3. This completes the all-size classification with a finite exact base certificate.

The positive controls are particularly concrete:

    H(C3[TT3,TT3,1])=7²+2=51=1+2·5²,
    H(C3[C3,C3,1])=9²+2·3²=99=1+2·7².

Both seven-vertex counts are independently checked by permutation enumeration. The second has content3 and a square factor9; it is not a primitive representation even though the distinguished representation(1,7) is primitive. Neither25 nor49 is quasiperfect: the graph construction asserts the output integer only, not the divisor-sum equation.

The cut bound explains why the attractive fixed-factor switch cannot produce infinitely many physical targets: it traps the ratio Q2/H near5 and therefore bounds graph size. Enlarging the arithmetic parameter y does not create a new graph source.

## 5. A complete finite decision procedure for this norm grammar

**Theorem: diagonal cyclic target decision.** Given a positive integer N, existence of a nonempty tournament A with N=N(A) is decidable by this finite complete procedure:

1. If N is not3 modulo8, reject.
2. Enumerate positive odd integer pairs(x,y) with x²+2y²=N.
3. For each pair enumerate all labelled tournaments on1≤n≤floor(x/y) vertices.
4. Accept exactly when H(A)=y and Q2(A)=x+y, emitting C3[A,A,1].

Sufficiency is the literal construction(1). Necessity is the cut bound and parity; therefore no omitted larger graph can rescue a failed pair. This is a mathematical completeness theorem, not an efficiency claim or a complete Chan–Pak algorithm. It gives a useful target-specific existence compiler without any unproved uniform witness-order bound.

For51 the only positive odd representations are(7,1),(1,5). The first yields TT3; the second has no size budget. For99 they are(9,3),(7,5),(1,7). The first yields C3; the others cannot have H5 orH7 within their respective size budgets. Thus the two examples genuinely use different norm coordinates and do not conceal a same-m realization.

## 6. Pell gives candidates, not graph certificates

For a proposed physical H-value y, the distinguished quasiperfect-shaped target condition is

    x²−2m²=1−2y²,   x=Q2(A)−y,   x≥|A|y.         (7)

This is a generalized Pell equation with a concrete graph-interface guard. The usual unit action generates further integral solutions from(x,m):

    (x,m) ↦ (3x+4m,2x+3m).

It preserves(7), but does not preserve the existence of A with profile(H,Q2)=(y,x+y). Example: from y1,(x,m)=(1,1), the next solution is(7,5), physically TT3; the following solution is(41,29). H1 forces a transitive graph and Q2=2^n, while x+1=42 is not a power of two, so the next Pell step is not physical. This is an exact first failure, not a merely hypothetical sidecar loss.

A promising new operation would have to implement the Pell change on Q2 while preserving or deliberately transporting H and its endpoint constraints. No such graph operation is supplied here. The actual results are the norm constructor, content/size/prime guards, complete small-switch classification, and a complete finite diagonal-grammar compiler. They give quasiperfect arithmetic a substantive but precisely bounded role without pretending to settle either open problem.

## Reproduction

Run `python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/arithmetic/norm_check.py`. The receipt gives121 literal cyclic-substitution controls, a complete factor-three profile search through order5 justified by the all-size bound, two independent seven-vertex permutation counts,500 scalar switch controls, and transitive diagonal examples through order12. The transitive examples are finite data, not a classification of all exponential/Pell intersections.
