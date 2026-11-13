# Module-response refinement of the proof forest

**PROVED with a CITED strong-core bound; FINITE-EXACT small-head certificate.
Chan--Pak remains OPEN.** The constructive comparison below is an application
of inherited frame identities. The new restrictions concern complete count
fibres and neutral padding of a fixed nonunit module.

## 1. Which smaller module could supply a prescribed prime target?

Write P_j(M) for ordered j-tuples of possibly empty directed paths covering M.
Then P_1(M)=H(M). The incoming
[core-replacement frames](../../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md)
have, for every nonempty module M with (a,p)=(H(M),P_2(M)),

    source S: 15p-a,    target B2:45p-a,    target B4:45p+a.          (1)

These identities are inherited. To supply a prescribed target N in this
language, it is necessary and sufficient to supply an actual module with

    (H,P2)=(a,(N+a)/45)   or   (a,(N-a)/45).                       (2)

Because P2 is even, the first choice requires a=-N mod90 and the second
requires a=N mod90. Cutting every Hamiltonian path at its |M|+1 boundaries
injects into ordered two-bin covers, so p>=(|M|+1)a>=2a. Thus a<=N/89 in the
minus case and a<=N/91 in the plus case. The module-count target is much
smaller than the canonical ternary source near N/3, but its P2 response is
an additional exact obligation.

For prime N, gcd(a,p)=1 is necessary: the gcd divides N and is at most a<N.
For N=7q with q prime, the gcd may instead be7; it is not automatically an
attained tournament factor. For critical N>35 one also needs gcd(a,15)=1.
A factor3 or5 of a would divide N. These are necessary tests, not sufficient
conditions for the physical response pair.

The forest can build a candidate module as an ordinal product, because

    (H,P2)(M1=>M2)=(H(M1)H(M2),P2(M1)P2(M2)).                     (3)

Choosing its factorization therefore chooses both coordinates. However,
constructing the scalar a does not certify the required P2. This is the
point at which the scalar proof forest needs a response sidecar.

## 2. Complete response fibres at H11 and H13

**Theorem C.** Every tournament T with H(T)=11 or13 has order n>=5 and

    H=11: P2(T)=100*2^(n-5),
    H=13: P2(T)=104*2^(n-5).                           (4)

Every listed order occurs. The statement includes all root choices, all
arc patterns and arbitrarily many neutral vertices, not just the displayed
representatives.

**Proof.** Strong components multiply H. Since11 and13 are prime, exactly
one nontrivial component has that count, and all other components are
singletons. For a strong component C, the published bound

    H(C)^3 >= 5^(|C|-1)

is the weaker consequence of Busch's Theorem1. Since13^3<5^5, such a component
has at most five vertices. Complete labelled enumeration through order four
has maximum H=5, so the component has exactly five vertices. The complete
order-five enumeration gives P2=100 for every H11 tournament and P2=104 for
every H13 tournament. All240 matching labelled tournaments are independently
checked by permutation enumeration. Finally (3) and P2(singleton)=2 give
(4). Adding any number of singleton components proves existence at every
n>=5. QED.

The order bound is **CITED**, not re-proved here; its formula and SCC product
statement were checked directly in [Busch's primary paper, Theorem1,
pp.3--4](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/).
The finite core-profile certificate and its full-fibre application are the
work of this note. No uniqueness of the core's isomorphism class is needed.

Convenient literal representatives are

    H11: (24,1,3,6,14),     P2=100,
    H13: (18,8,3,5,14),     P2=104.

**Consequence.** No module of count13, at any order, can make either target
in (1) a prime or a seven-prime. Its outputs are

    13*(360*2^t-1), 13*(360*2^t+1),   t>=0.

They are composite, exceed91, and have a factor13; a seven-prime divisible
by13 would have to be91. A same-H13 repair cannot remove this obstruction,
even if it changes the entire graph and uses arbitrarily many vertices.

The H11 fibre has different capability. One singleton of padding gives
(H,P2)=(11,200), with actual frame counts

    2989 -> 8989,9011.

The last value9011 is prime and its literal strong tournament is recounted.
This is a direct critical-target construction, not a fixed-digit step:
the increments relative to3*2989 are22 and44. Its canonical positive-digit
parent would be3003, not2989. The source counts11 and13 are siblings of3
under 3h+2 and3h+4, yet have different all-order module response constraints.

## 3. Fixed nonunit padding has a stronger constant-digit obstruction

The integer bin expansion used here is inherited from the path-bin calculus
and the [clone-clock P-form](../../clone_clocks_negative_collatz_20261113/CLOCKS.md).
The earlier [birth-grammar G1 argument](../../birth_grammar_repairs_20261112/PROOFS.md)
also uses independence of exponential coefficients. The following is a
specific corollary for a nonunit module, not a claim to a new spectral method.

Fix a nonempty tournament M0 of count a, and put

    M_t=M0=>TT_t,             t=0,1,2,...,

where TT_0 is empty. For a fixed rooted frame F, the inherited run expansion
and binomial inversion give

    H(F[M_t])=sum_(j=1)^J c_j P_j(M0) j^t,    c_j integers.         (5)

Indeed P_j(M_t)=P_j(M0)j^t by ordinal multiplicativity, and nonempty-bin
counts are integer binomial combinations of the P_j. One may take J to be
the number of exterior vertices plus one. All P_j(M0) are positive for
j>=1, and P_1(M0)=a.

**Theorem D (constant-term content).** Let F,G be any two fixed one-slot
frames and r any integer. If

    H(G[M_t])-r H(F[M_t])=d                            (6)

for infinitely many nonnegative integers t, then a divides d. In particular,
for odd a>1 no such pair can realize fixed digits d=2 or4 infinitely often.
The statement allows arbitrary fixed frame sizes and any number of visits
to the module.

**Proof.** Subtract d from (5) for the two frames. The coefficient of the
constant base1 is a times an integer minus d. Distinct positive exponential
bases have independent coefficients; more quantitatively, a nonzero sum of
m distinct real exponential modes has at most m-1 distinct real zeros.
To prove that fact, divide by the least exponential mode. Differentiation
leaves m-1 modes; induction and Rolle's theorem bound the original zeros.
Infinitely many integer zeros thus force every coefficient to vanish,
including the constant coefficient. Therefore a divides d. QED.

**Finite-bank version.** If a does not divide d and J is the largest run
cap of F,G, equation (6) holds at at most J-1 padding lengths. A finite bank
of frame pairs therefore covers only finitely many t, even if the selector
switches among them. Its coverage is bounded above by the sum of J-1 over
the pairs. This is a bound on possible matching parameters, not a claim
that each pair has that many physical successes.

The signed integer response lattice allows isolated unit-digit coincidences.
For the H11 representative,

    -18*P1(M_t)+P2(M_t)=-198+100*2^t

equals2 precisely at t=1. The companion also gives a four-mode integer
combination attaining2 at exactly t=0,1,2, showing that the three-zero bound
is sharp in that signed response model. Those coefficient combinations are
not asserted to be actual differences of tournament frames.

The unit boundary is a positive control. For M0 a singleton, a=1, and the
incoming fixed frames do realize digits2 and4 on M_t=TT_(t+1) for every t.
There is no conflict with their theorem. Likewise a one-time Bezout source
at a nonunit module is not excluded; the obstruction concerns continuing
through an entire neutral-padding family with finitely many fixed frames.
Changing the base module, introducing multiple independently changing slots,
or using unboundedly many frames is outside this particular conclusion.

## 4. What this changes constructively

The source selector has two distinct tasks. A scalar factorization forest
supplies smaller counts; an actual module construction must also hit the
response required by the chosen frame. At H13 that response repair is
impossible for the two incoming target frames, so the selector should change
module count or construction family. At H11 a nonunit module does give a
prime target, but padding cannot automatically convert that success into
fixed unit-digit returns with a finite frame bank.

The proposed cyclic two-run source (a1,a2)=(95,11) for N=9011 is now
**REFUTED** by the new [two-run cone, C4](../absorption/TWO_RUN_CONE.md):
95>3*11. Finite catalog absence is not the reason for this conclusion.

There is, however, a constructive acyclic alternative. The source coefficient
pair (273,0) is physically realized by adjoining an exterior Q of count273.
Take Q=C3=>B4[TT1], so H(Q)=3*91. For M=M11=>S of count11 and P2=200,

    source M=>Q has H=3003,
    target B4[M] has H=9011=3*3003+2.                  (7)

The module M is retained literally. Both graphs can use the same eighteen
vertex names: the source has six module vertices and twelve exterior vertices;
pad the fourteen-vertex target by four singleton sinks. Reassign only the
exterior roles and rewrite the required arcs. No internal module arc changes,
and it remains a uniform module at both endpoints. The receipt records this
actual batch rewire. It is not a claim that the intervening individual arc
flips have controlled counts, or that the exterior rewrite has uniformly
bounded support.

This is an instance of a general one-shot product-source rule. If a target
frame gives N from a retained module of count a, and the desired digit d has

    q=(N-d)/(3a) a positive allowed integer,

an actual Hq witness supplies the canonical source by an ordinal join. The
target is already certified by its frame. Pad the smaller exterior to use
one common vertex set, retain the module's internal arcs, and replace the
exterior context. The new source has count(N-d)/3. This does not transport
an arbitrary supplied source into that module architecture.

For M_t=M11=>TT_t the rule has two exact infinite parameter phases:

    N_plus=4500*2^t+11, d=2: t=1 modulo10,
    N_minus=4500*2^t-11, d=4: t=2 modulo10.             (8)

Indeed the respective source counts are1500*2^t+3 and1500*2^t-5. Divisibility
by11 requires2^t=2 or4 modulo11;2 has order10. In either phase q>=273 is
smaller than N and is not7 or21. Under smaller-count induction, its witness
therefore exists and the product-source rule gives a genuine fixed-digit
terminal construction. Whenever N is prime or seven-prime it supplies the
corresponding critical obligation. No infinitude of prime values is asserted.

The target family in (8) was already attained by the inherited frames; the
new connection is its canonical source with the retained nonunit module.
The exterior of count q changes with t, so this does not contradict the
finite-frame obstruction of Theorem D. At t=1 the smaller factor273 is
already explicitly attained, making (7) unconditional.

The uniform question now becomes a marked response selector with growing
physical structure. A finite set of scalar parent formulas and neutral
padding alone does not settle it. No all-prime or all-seven-prime selector
has been obtained.

## 5. Exact verification

Run

    python -B 05-knowledge/results/prime_terminal_20261113_1610/arithmetic/module_check.py
    python -O -B 05-knowledge/results/prime_terminal_20261113_1610/arithmetic/module_check.py

The universe is all1099 labelled tournaments of orders1--5. Exactly120
labelled H11 and120 labelled H13 tournaments occur, with the single respective
P2 values100 and104. All240 receive independent literal permutation counts.
The all-order theorem uses the cited core bound plus this complete head.

Further controls evaluate all four bin counts P1--P4 and literal incoming
frames on both modules with zero, one and two padding vertices; check the
unit-module boundary; and check1250 formal coefficient/digit combinations
through padding20. Explicit one-hit and three-hit positive controls prevent
the zero bound from being tested only on empty cases. These last combinations
are an arithmetic superset, not a claim of frame realizability.
