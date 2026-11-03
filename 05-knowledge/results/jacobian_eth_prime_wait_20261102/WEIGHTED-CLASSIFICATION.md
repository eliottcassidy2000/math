# Weighted Keller maps: all fibres and arbitrary inertia partitions

**PROVED here from inherited constructions; exact symbolic and independent
finite-field checks below. Self-audited, not independently agent-audited.
No claim of classifying arbitrary Keller maps or resolving JC(2).**

The inherited starting points are [THM-3438, weighted lift and degree spectrum](../../../01-canon/theorems/THM-3438-weighted-lift-keller-degree-spectrum.md),
[THM-3441, the quartic boundary](../../../01-canon/theorems/THM-3441-weighted-quartic-jelonek-components-and-inertia-parity.md),
and [THM-3448, cyclic inertia](../../../01-canon/theorems/THM-3448-weighted-keller-cyclic-jelonek-inertia-family.md).
The extension is to repeated **nonzero** roots of the auxiliary polynomial,
excluded in THM-3448 §§1–5. The new conclusions are a pointwise fibre rule,
arbitrary products of inertia cycles, and a partition lower bound for
inequivalent maps. This is not a priority claim against the literature.

## 1. Setup and the exact finite-fibre rule

Work over C. Choose a polynomial R of degree n, with p=R', and constants
b,c nonzero, satisfying

    R(0)=R(1)=p(0)=0, p(1)=-c,
    kappa=p'(1)/c != -2, a=-(1+kappa)/(2+kappa).

Write q=(wp-R)/c. For v=xy, t=x²z, u=1+v, gamma=1+av+bt, put

    alpha=u+q(u gamma)/gamma²,
    beta=c+p(u gamma)/gamma,
    F=(alpha/x², beta/x, x gamma).

THM-3438 proves polynomiality, det JF=bc, generic degree n and geometric
monodromy S_n. In particular these are composition atoms. The polynomial
degree of a coordinate is a different invariant from n.

**W1 (all target fibres within this entire family).** For target (A,B,C),
set P=BC, Q=AC² and T(w)=R(w)-Pw+cQ.

* If C != 0, the finite preimages correspond exactly to the simple roots of T.
  For D=gcd(T,T'),

      #F^(-1)(A,B,C) = n - deg D - deg rad D.             (1)

* If C=0 and p'(0) != 0, then

      #F^(-1)(A,B,0) = 3 if B²-2c p'(0) A != 0,
                     = 1 otherwise.                     (2)

* If C=0 and p'(0)=0, then

      #F^(-1)(A,B,0) = 2 if B != 0, and 1 if B=0.         (3)

These are counts over C, not counts of rational or real points.

**Proof off C=0.** The reconstruction identities are

    gamma=(P-p(w))/c=-T'(w)/c,
    x=C/gamma, y=(w-gamma)/C,
    z=gamma[gamma(gamma-1+a)-aw]/(bC²).

Every simple root reconstructs a unique finite point. Conversely C != 0
forces x,gamma != 0, so every finite point gives a simple root. For a root
of multiplicity e>1, D has multiplicity e-1 and rad D has multiplicity one;
their combined degree deletes all e copies. This proves (1). In particular
n-1 cannot occur in this chart. Multiple roots are lost to infinity,
not finite critical points.

**Proof on C=0.** The source splits disjointly into x=0 and gamma=0,x!=0.
At x=0,

    B=-c y/(kappa+2), A=h y²+b(kappa+2)z,

where h is the coefficient of v² in alpha(v,0). This triangular map has
one preimage for every (A,B). At gamma=0,x!=0, put d=p'(0). Polynomial
evaluation gives beta=c+du and alpha=u+d u²/(2c).
If d!=0, elimination gives

    (B²-2cdA)x²=c²,  u=(Bx-c)/d.

There are exactly two solutions when its coefficient is nonzero and none
otherwise; y,z then reconstruct uniquely. If d=0, instead Bx=c, u=Ax²:
there is one additional point exactly when B!=0. This proves (2)–(3).
No simplicity hypothesis on the other roots of R was used. ∎

**Exact fibre-count algorithm.** For a seed and target over Q (or another
specified exact coefficient field), the C!=0 count can equivalently be
computed as

    n-2 deg gcd(T,T') + deg gcd(D,D'), D=gcd(T,T').

This requires two univariate gcds and no algebraic root enumeration or
expansion of the three coordinate polynomials. Ordinary dense Euclidean
division uses O(n²) coefficient-field operations. This is an arithmetic
operation bound, not a claim about unit-cost arbitrary-bit integers or
about counting rational points. The C=0 branches need only the coefficient
tests (2)–(3). Thus this structured-family algorithm is unconditional and
does not need the external SAT result.

**Hostile boundary.** For THM-3438's quartic G, x=0 gives
G(0,y,z)=(87y²+6z,y,0). Its other boundary points satisfy
(B²-A)x²=1. Thus C=0 has three points when A!=B² and one otherwise.
This recovers THM-3441 and Gao's Appendix A; the entire quartic spectrum
{0,1,2,3,4} is inherited, not a new discovery. Applying (1) on C=0 would
wrongly erase both the chart information and the three-point fibres.

## 2. Repeated roots give arbitrary disjoint cycles

Factor R over C as

    R(w)=a_n w^mu (w-1) product_i (w-r_i)^s_i,
    mu>=2, r_i distinct and outside {0,1}.

The distinguished root 1 is simple since p(1)=-c. Let n=mu+1+sum s_i.

**W2 (component-generic boundary anatomy).** Around a generic point of
C=0, geometric inertia is the permutation consisting of a cycle of length
mu-1, cycles of lengths s_i, and fixed points; omit cycles of length one.
It has three finite sheets when mu=2 and two when mu>2. All other sheets
escape. The raw discriminant valuation and monogenic index are

    v_C Disc_w T = mu + sum_i(s_i-1),
    inertia defect = mu-2 + sum_i(s_i-1),
    index length = 1.                                  (4)

This is one cyclic local inertia group with the stated permutation cycle
type, not a collection of independently selectable local generators.

**Proof.** Work over an algebraic closure of C(A,B) and use C as uniformizer.
At zero, the Newton vertices (0,2),(1,1),(mu,0) give two roots of valuation
one for mu=2; for mu>2, one root has valuation one and mu-1 roots have
valuation 1/(mu-1). The latter form one cycle. The residual quadratic
when mu=2 is (p'(0)/2)W²-BW+cA, separable at the generic point.

At r_i !=0, write R(r_i+delta)=h_i delta^s_i+O(delta^(s_i+1)), h_i!=0.
The leading equation is

    h_i delta^s_i = B r_i C.

It gives one tame s_i-cycle. For s_i>1, gamma has valuation (s_i-1)/s_i,
while y=(w-gamma)/C has leading term r_i/C, so every sheet escapes.
For s_i=1, let g_i=-p(r_i)/c. Simultaneous vanishing of the leading
y and z poles would force r_i=g_i and g_i²(g_i-1)=0, hence r_i=1,
which is excluded. The designed root 1 yields the finite x=0 section.
The remaining finite sheets are exactly W1's gamma=0 solutions. This also
proves the counts without assuming any other root simple.

The zero cluster contributes mu to the discriminant valuation; a nonzero
s_i-cluster contributes 2 binom(s_i,2)/s_i=s_i-1. Differences between
clusters are units. Tame inertia contributes cycle length minus one.
The order/normalization discriminant identity makes their difference twice
the index; (4) follows. ∎

The ordinary discriminant in (P,Q) is irreducible and generically reduced:
its parametrization is (P,Q)=(p(w),(wp(w)-R(w))/c), and its tangent slope
recovers w/c. Thus exactly one pair meets generically. Over C!=0 its
pullback is one irreducible surface L=0, with generic fibre n-2 and
transposition inertia. Removing the exact C power in (4) gives its closure.

The inverse formulas show properness off V(CL); every point of L on its
generic open loses two sheets, so L is a nonproper component. W1 shows C
is a nonproper component whenever n>=4. Since the nonproper locus is
closed, these observations prove it is exactly V(C) union V(L) for n>=4.
For n=3 the inherited exceptional case has only V(L).

## 3. A partition lower bound for inequivalent primitive maps

Let p_part(j) denote the integer partition function, p_part(0)=1.
Equivalence means pre- and post-composition by polynomial automorphisms.

**W3.** For every n>=4 and every ambient dimension m>=3, there are at least

    p_part(n-3) + p_part(n-2) - 1                       (5)

pairwise inequivalent determinant-one Keller maps of generic degree n,
all with geometric monodromy S_n and hence all composition-irreducible.

**Construction.** For every partition lambda of n-3, take mu=2 and use
lambda's parts as the multiplicities s_i. For every partition lambda of
n-2 other than all ones, take one part ell>=2, set mu=ell+1, and use the
remaining parts as the s_i. Choose r_i=-i and normalize

    R(w)=-w^mu (w-1) product_i ((w+i)/(1+i))^s_i.        (6)

Then c=1 and

    kappa=-2[mu+sum_i s_i/(1+i)] < -2,

so THM-3438 applies with b=1. The first construction has three finite
boundary sheets; the second has two. Its escaping permutation has cycle
partition lambda, with ones interpreted as fixed escaping sheets.

**Why the maps differ.** Attach to each irreducible nonproper component
the pair (generic finite fibre count, geometric inertia cycle type).
Its unordered two-component multiset is preserved by polynomial source
and target automorphisms. Every construction has the same ordinary entry
(n-2, transposition). Its C entry is (3,lambda) or (2,lambda), with finite
fixed sheets added when specifying the full permutation on n letters.
Within each row, n and the nontrivial cycles recover the number of ones.
Between rows, the finite count separates them. Adding the same ordinary
entry to two multisets cannot make distinct C entries equal. This proves
pairwise inequivalence even if an equivalence tries to swap components.
Stabilization by identity coordinates preserves the two-component ledger,
generic degree, and monodromy, proving the assertion for m>3. ∎

The first lower bounds, for n=4,5,6,7,8,9, are 2,4,7,11,17,25.
They are lower bounds for maps, not a complete classification and not
counts of ambient dimensions. No claim of a finite number of classes in
a fixed degree is made.

**Concrete mixed-cycle witness.** R=w³(w-1)(w-2)³ gives c=1, kappa=0,
a=-1/2, degree seven, det JF=1, global monodromy S_7. At C=0 there are
two finite sheets and five escaping sheets, with inertia (2)(3), raw
discriminant valuation five and index one. This proves that the local
cycle pattern can be composite while the global covering is primitive.

## 4. The same positive coordinate can conceal inequivalent maps

The incoming [Session144 equivariant boundary theorem](../session144_jacobian_eth_mrb_20261103/jacobian_followup/universal_equivariant_boundary_corollary.md)
classifies the fibres of the positive-weight coordinate and proves a
composition law in the equivariant category. The following gives a precise
limit to using that coordinate as a classifier of the whole map.

**W4.** For every n>=4, at least p_part(n-3) of W3's inequivalent primitive
degree-n, determinant-one maps can be chosen with the identical component

    F_3(x,y,z)=x(1-xy/2+x²z).                            (7)

In particular even the entire positive-coordinate polynomial, not just
its degree or zero-fibre type, fails to classify the map.

**Proof.** For a partition s_1>=...>=s_r of n-3, set

    r_i=1+s_i(2^r-1)/2^i,
    R=-w²(w-1) product_i ((w-r_i)/(1-r_i))^s_i.         (8)

The positive rationals r_i-1 are strictly decreasing, hence the roots
are distinct and outside {0,1}. Their logarithmic derivative satisfies

    sum_i s_i/(1-r_i) = -sum_(i=1)^r 2^i/(2^r-1)=-2.

Thus c=1, kappa=-2[2-2]=0 and a=-1/2 for every partition.
Taking b=1 gives (7). Their C-component inertia still has the different
cycle partitions s, and three finite sheets, so W3's unrestricted
inequivalence argument applies unchanged. ∎

For n=5 the two seeds can be written

    R_1=-w²(w-1)(w-2)²,
    R_2=-(8/9)w²(w-1)(w-5/2)(w-7/4).

Both have exactly the same F_3 in (7), generic degree five and S_5
monodromy. The C-component inertia is a transposition for the first and
the identity for the second. This is not a contradiction to the incoming
coordinate-fibre theorem: that theorem classifies one coordinate, not all
coordinates jointly.

There is also an elementary degree-four separation. For R=w³-w⁴, c=1,
W1 excludes three-point fibres everywhere: C!=0 excludes n-1, while C=0
has only one or two points. At C=1, the pairs (P,Q)

    (1/8,-1/64), (0,0), (-1,-1), (0,1)

have respectively 0,1,2,4 finite preimages. In the first case
T=-(w²-w/2-1/8)²; in the third T=-(w-1)²(w²+w+1).
Thus this map has the exact spectrum {0,1,2,4}, whereas the inherited
quartic G has {0,1,2,3,4}. They cannot be equivalent under arbitrary
source/target automorphisms. The explicit separation is stronger than
merely observing different displayed formulas.

## 5. Reproduction and independent checks

From the worktree root:

```sh
python3 -B 05-knowledge/results/jacobian_eth_prime_wait_20261102/weighted_fibres.py
python3 -O -B 05-knowledge/results/jacobian_eth_prime_wait_20261102/weighted_fibres.py
```

The checker expands four maps (degrees 3,4,5,7), checks their nonzero
constant determinants, and computes full two-parameter discriminants.
An independent enumeration evaluates 7,196 actual source triples over
F_5,F_7,F_11 and checks every C=0 target against W1's finite-field version.
The complex point counts are proved above, not inferred from those finite
fields. Sixty-four root multiplicity patterns audit (1). Partitions are
enumerated independently of a coin-change recurrence through degree 14;
every seed's normalization, admissibility and distinct invariant ledger
is checked. The general conclusion rests on W3's proof, not this cutoff.
The fixed-coordinate construction (8) is checked on all partitions through
degree 14, and the four degree-four witnesses are checked by exact gcds.

Next decisive questions: can the component ledger distinguish positive-
dimensional moduli at one degree; which non-S_n global monodromies are
realizable; and what extra conditions characterize polynomial gluing of
an abstract étale cover? Those are not answered by (5).
