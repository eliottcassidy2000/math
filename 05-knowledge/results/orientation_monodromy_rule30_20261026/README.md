# Orientation, monodromy and four-state transport

2026-10-26. **PROVED elementary reformulations / CITED classical inputs /
FINITE-EXACT controls / OPEN applications**, distinguished below. Internal
review only; no independent agent audit or novelty claim.

## Inheritance and concept board

Anchor: turn the proposed geometric analogies into maps preserving a named
predicate. Niche: arithmetic exchange signs and four-vertex converse.
Wildcard: Rule 30's reset/flip/identity return maps.

The closest existing work is the [analogy atlas](../analogy_reduction_atlas_20261025/ATLAS.md),
the [known prime-tournament extension theorem](../root27_machine_bridges_20261026/anchor/LEGENDRE_SIDECAR.md),
and **THM-3456**, [left-permutive trace coordinates](../../../01-canon/theorems/THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary.md).
The canonical hostile is Rule 60: the same free trace language and Haar law,
but a constant single-seed center. The corrected near miss is that a unique
fixed point need not attract all iterates. The least-used relevant coordinate
is the *actual return map*, including resets rather than only invertible
monodromy. **THM-4047**, [left-front affine monodromy](../../../01-canon/theorems/THM-4047-rule30-left-front-affine-monodromy-clock.md),
already proves the stronger physical Rule 30 return theorem; it is inherited,
not reproved as a new discovery here.

| Object | Preserved datum | Lost datum / cheapest hostile |
|---|---|---|
| Gram matrix and norm | Length inequality | Orientation; opposite vectors saturate CS but not triangle equality |
| Probability laws in a compact product | Subsequence existence | Limit identification; one fixed point plus a two-cycle |
| Legendre symbols | Graded exchange sign | Actual edge signs after basis rephasing |
| Signed area parametrization | Oriented integral | Nonorientable gluing; a zero collapses the fibre instead |
| Logarithm and root covers | Winding and its residue | Integer winding after reduction modulo m |
| Binary return maps | Flip/reset/identity and phase | Initial boundary and the moving diagonal |

## 1. Cauchy--Schwarz and triangle: exact scope and equality boundary

**CLASSICAL / proved here.** For a norm induced by a real or complex inner
product, expansion of `||x+y||^2` identifies the triangle inequality with
`Re<x,y> <= ||x||||y||`. Applying it to all scalar phase rotations of y
recovers `|<x,y>| <= ||x||||y||`. Conversely CS bounds the cross term and
gives the triangle inequality. This is equivalence of the inequalities
*for all pairs*, not identity of their equality cases for a given pair.

CS is also positivity of the two-vector Gram determinant:

    ||x||^2 ||y||^2 - |<x,y>|^2 >= 0.

In a real plane this is squared parallelogram area. Squaring forgets
orientation. For y=-x!=0, CS is an equality but triangle inequality is
strict. An arbitrary norm need not arise from an inner product: the l1
norm on R2 violates the parallelogram law. One may still write mixed-norm
inequalities for a separately chosen pairing, but that is a different claim.

## 2. Compactness, CLT, and an exact Rule 30 bridge

**CLASSICAL.** For iid real X_i of mean0 and variance1, the laws of
Z_n=(X_1+...+X_n)/sqrt(n) are tight by Chebyshev. A diagonal subsequence
of their distribution functions at rational points lies in `[0,1]^Q`.
Monotonicity supplies a distribution function and tightness prevents mass
escaping to infinity. This is the precise compact-product step related to
Tychonoff; bare product compactness is not the whole weak-convergence proof.

The additional identification is

    phi(t)=1-t^2/2+o(t^2),
    E exp(it Z_n)=phi(t/sqrt(n))^n -> exp(-t^2/2).

The finite-variance Taylor expansion and characteristic-function continuity
theorem identify the Gaussian. No higher moments or cumulants are assumed.
The standard Gaussian is the unique fixed point of normalized convolution
**at mean0 and variance1**, by iterating the same expansion. Unique fixed
point alone does not prove attraction: on the compact discrete set {0,1,2},
the continuous map 0->0, 1->2, 2->1 has one fixed point and a two-cycle.
Thus the old atlas sentence asserting that uniqueness alone makes every
cluster point a limit has been repaired.

**CITED repo theorem plus immediate corollary.** THM-3456 gives, for any
left-permutive finite-alphabet radius-one rule, the homeomorphism

    initial configuration -> (center trace, positive initial half-line).

Under uniform Bernoulli input this preserves product measure. For Rule30,
`f(l,c,r)=l XOR (c OR r)`, the center bits Y_t are therefore exactly iid fair
bits. Consequently their finite counts are Binomial(N,1/2), and

    N^(-1/2) sum_(t=0)^(N-1) (2Y_t-1) -> N(0,1) in distribution.

This is an exact random-input CLT, not a mixing conjecture. It does not
address the fixed single-one input: the prescribed inverse left boundary
is essential. Rule60 `l XOR c` has the same free law and a single-seed
center identically1. This is the same logical boundary as Haar-distributed
Collatz valuation words versus one marked positive integer. The respective
native maps and guards still differ; no cross-system convergence follows.

## 3. Reciprocity is an exchange law, not a four-class count

**CITED classical reciprocity.** Set epsilon(p)=(p-1)/2 mod2 for odd primes.
For distinct primes,

    chi(p,q)=(p/q),
    chi(p,q)=(-1)^(epsilon(p) epsilon(q)) chi(q,p).

On primes3 mod4 this is antisymmetry and hence a tournament. With an even
epsilon at either endpoint it is a symmetric pair. The sign is precisely
the usual graded interchange sign. This statement specifies the exchange
law, not a full super-vector-space bilinear structure.

The stronger known result is already in the incoming repo: CRT and
Dirichlet give every finite extension pattern on primes3 mod4, so this is
the countable random tournament. Thus every finite tournament occurs,
including the four classes at order4. This is not a theorem counting those
four classes, nor a frequency theorem for increasing-prime prefixes.
Primary sources checked: [Cameron, p.3](https://www.math.uni-bielefeld.de/~rehmann/ECM/cdrom/3ecm/pdfs/pant3/camer.pdf)
and [Jones, section6 remark4](https://arxiv.org/html/1912.02709v1#S6).

### An associative algebra and its precise loss of arithmetic

**PROVED elementary construction; no novelty claim.** For a finite ordered
set of distinct odd primes and a characteristic-not2 coefficient field,
take basis e_A indexed by subsets A. Define

    e_A e_B = 0                                      if A intersects B,
    e_A e_B = [product_(p in A,q in B) chi(p,q)] e_(A union B) otherwise.

The empty subset is the identity. If any two of A,B,C overlap both bracketings
of their product vanish. If they are pairwise disjoint both bracketings
multiply the same three cross-products. This proves associativity. Also

    e_A e_B = (-1)^(epsilon(A) epsilon(B)) e_B e_A,
    epsilon(A)=sum_(p in A)epsilon(p) mod2.

Define g(A)=product_(i<j in A)chi(p_i,p_j), and f_A=g(A)e_A. Then f_A is
the ordered product of its singleton generators. Its multiplication is
the standard sign for inversions of *odd* generators only. Therefore the
unmarked algebra is an exterior algebra on the odd generators tensored
with commuting square-zero variables for the even generators.

**Loss boundary:** if all four primes are3 mod4, every tournament sign
pattern yields the same unmarked exterior algebra on four generators.
It cannot distinguish the four tournament classes. Keeping the specified
subset basis e_A retains every edge through e_{p}e_{q}; forgetting that
basis erases the arithmetic signs. This is a concrete instance of why an
abstract algebra can be too coarse for a port-preserving recurrence.

## 4. What does survive among the three four-element sets?

**PROVED limited dictionary, not a native composition theorem.** For
a,b in F2 define the map on one bit

    F_ab(x)=(1+a+b)x+a.

These are all four maps of a two-element set. With J(x)=x+1,

    F_ba = J composed with F_ab composed with J.

| (a,b) | Return map | Assigned order-four class |
|---|---|---|
| (0,0) | identity | transitive TT4 |
| (1,1) | flip | strong class |
| (0,1) | constant0 | source above a directed triangle |
| (1,0) | constant1 | directed triangle above a sink |

Swapping the parity labels, conjugating the return map by J, and taking
tournament converse all fix two objects and exchange the other two.
The displayed assignment is an equivariant bijection for this C2 action.
It is a declared choice, not a canonical arithmetic classification.
It transports **converse symmetry only**. Tournament joins change size,
and no native operation realizing function composition on these four
isomorphism classes has been supplied. An arbitrary multiplication table
transported through the bijection would not supply such a realization.

This is more information than a shared numeral, but much less than an
equivalence of quadratic reciprocity, tournament enumeration and dynamics.

## 5. The area band, the actual Möbius bundle, and FTA

For a C1 real function f on [a,b], use the explicit map

    Phi(x,t)=(x,t f(x)),      0<=t<=1.

Then `Phi*(dx wedge dy)=f(x) dx wedge dt`. The Jacobian changes sign when
f does. At f=0 the vertical fibre collapses: this is degeneration of a
parametrized signed region, not a Möbius gluing. Its oriented area is
exactly integral_a^b f(x)dx. With F'=f, FTC gives F(b)-F(a).

A genuine Möbius bundle comes from a loop gamma along which a complex
function g never vanishes. Write w=wind(g composed with gamma,0). Continuing
log g around the loop adds 2pi*i*w. Continuing an m-th root multiplies it by
exp(2pi*i*w/m). Thus root monodromy is the quotient of integer winding by m.
The logarithm surface itself is biholomorphic to an ordinary copy of C:
use coordinate L and projection L->exp(L) in C*. The extra structure is the
covering projection and its deck translations, not a unique exotic topology.
For m=2, the **associated real sign line bundle** over the loop has gluing
(1,v)~(0,(-1)^w v). A bounded fibre gives a Möbius band exactly for odd w.
The square-root covering itself is not this band: over a circle it is a
circle for odd w, or two circles for even w. Its complex Riemann surface is
orientable. This corrects a type conflation in the older atlas.

**Exact FTC--FTA connection.** Along the parametrized loop choose a continuous
logarithm lift L. Ordinary calculus on its parameter interval gives

    integral_gamma p'(z)/p(z) dz = L(1)-L(0)=2pi*i*w.

For a degree-n polynomial choose a large circle on which its leading term
dominates. Interpolation to a_n z^n never crosses0 on that circle, so w=n.
If p had no zero inside, p/|p| would extend to the whole disk and its boundary
loop would have winding0. Contradiction. This proves existence of a root
without first assuming a factorization into roots. The spiral is a way to
store the lift L; it is not the topology of C, nor a unique possible topology
for an algebraically closed field. For example the algebraic numbers in C
form an algebraically closed, totally disconnected field in the subspace
topology. Any connected subset with two distinct points would have an
uncountable interval in a real coordinate projection, impossible here.

FTC above is stated over R with explicit regularity. Polynomial intermediate
value properties of a real-closed field alone do not justify arbitrary
analytic FTC. On the real algebraic numbers, the step at the missing cut
pi is locally constant and has derivative0 everywhere, yet changes between
3 and4. That is a hostile to silently replacing real completeness by real
closedness; it is not a counterexample to the ordinary real FTC.

## 6. Rule 30: where the flip is literal, and where it resets

THM-4047 uses `ell_r(t)=a_t(-t+r)` and proves

    ell_r(t+1)=ell_(r-2)(t)+ell_(r-1)(t)
                    +(1+ell_(r-1)(t))ell_r(t)        over F2.

For periodic lower-column drivers, composing one period gives x->Mx+C.
If any multiplier is0 this is a constant reset with exactly one fixed point.
If M=1,C=0 it is identity with two phases. If M=1,C=1 it is a flip with no
fixed point and a twofold return. These are precisely the four maps above.
Only the last two are invertible covering holonomy. Resets erase incoming
phase, so a theory of bundles or covering groups alone misses half the maps.

For the translation-only case x_(i+1)=x_i+b_i on a closed cycle, a section
exists iff sum b_i=0 in F2, with two sections then. Odd sum gives a flip;
the associated real line bundle is the Möbius case. This is an exact common
obstruction equation, not an identification of the whole spacetime diagram
with a surface. The physical center is `a_t(0)=ell_t(t)`, a moving diagonal;
fixed-depth periodicity and Haar CLT still do not control it.

## Next questions and verification

1. Which marked tournament boundary relations realize compositions of these
   four maps *natively*, without choosing fresh ports after each step?
2. Can the reset/holonomy distinction refine the moving-diagonal Rule30
   observer while retaining its actual inverse boundary? Existing fixed-depth
   results must not be repackaged as that missing uniform estimate.
3. For arithmetic transport, which subset-basis or ideal data survive a
   proposed gauge quotient? The algebra above supplies an all-size negative
   control: the unmarked exterior algebra forgets every tournament edge.

Run `python3 check.py` and `python3 -O check.py`; compare the saved receipts.
Universes: all1024 four-generator parity/sign systems; all64 labelled
four-vertex tournaments; all binary affine words through length6; complete
free Rule30/60 trace cones through depth6. These finite checks support the
formulas and loss examples; the all-size arguments are the proofs above.
No Collatz, Rule30 prize, LRC or black-core connectivity claim is made.
