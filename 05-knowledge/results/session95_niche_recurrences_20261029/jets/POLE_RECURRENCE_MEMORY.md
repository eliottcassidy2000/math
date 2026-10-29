# Two recurrence directions in a translated pole source

**Status: PROVED for the stated characteristic-zero family; exact finite controls.**
Source snapshot: `d3741fe6ed02f6072dbe6bfda2f3d04193d75ce0`.
No historical-priority claim is made for Prony reconstruction, Hankel rank,
finite differences, or rational summation. The new packet combines these
mechanisms into an exact iteration-memory law and an all-time bounded-jet
obstruction for the repository's two-flow family.

## 1. Inheritance and the underused bridge

The October 19 [metabelian invariant classification](../../metabelian_global_invariants_20261019/metabelian_global_invariants_and_logarithms.md)
uses canonical coordinates `z=1/x,w=-yx²`. Every relevant word has form

    Phi(z,w)=(z-a,w-F(z)),    F(z)=sum_t c_t/(z-t)^4.

For `a != 0`, a rational first integral, a rational normalized Hamiltonian
logarithm, and zero coefficient sum on each translation orbit `t+a Z` are
equivalent. That classification is inherited, not reproved as a new result.
Its finite-difference example hides global nonintegrability from a single
finite jet.

The [field-sensitive jet note](../../metabelian_field_sensitive_jets_20261019/field_sensitive_monomial_jet_lengths.md)
identifies normalized Taylor coefficients with finite weighted moments, while
[observable Hankel clocks](../../observable_hankel_clocks_20261020/OBSERVABLE_HANKEL_CLOCKS.md)
separates proved shift closure from a fitted recurrence. The older
[endpoint-Prony theorem, THM-2286](../../../../01-canon/theorems/THM-2286-endpoint-prony-lift-bank-and-sharp-owner-multiplier-landing.md)
uses the same finite exponential-sum mechanism for signed interval endpoints.
The present calculation measures what happens to that exact memory under
iteration, rather than treating a short list of moments as a faithful object.

## 2. General pole order and exact source iteration

Let K have characteristic zero, fix `a in K \ {0}`, a positive integer h,
and a finitely supported function `c:K -> K`. Combine equal sites first and
remove zero coefficients. Set

    F(z)=sum_t c_t/(z-t)^h,
    F_n(z)=sum_(j=0)^(n-1) F(z-ja),       n>=0.

The triangular map iterates exactly as

    Phi^n(z,w)=(z-na,w-F_n(z)).

Its merged coefficient at site u is

    c_u(n)=sum_(j=0)^(n-1) c_(u-ja).                 (1)

This is signed source addition. Counting source visits before merging gives
the wrong pole count when cancellation occurs. It is also essential to retain
the translation-orbit label, not just the total source weight.

Define normalized moments (including `0^0=1`)

    m_r(n)=[z^(-h-r)]F_n(z)/binom(h+r-1,r)
          =sum_u c_u(n) u^r,        r>=0.            (2)

The binomial divisor is nonzero in characteristic zero. All memory statements
below concern this normalized sequence, not the raw coefficients whose extra
polynomial factor changes their recurrence multiplicities.

## 3. Exact recurrence rank at a fixed time

Let `S_n={u:c_u(n)!=0}`, and `p_n=|S_n|`. The minimal monic homogeneous
recurrence that holds **for every r>=0** for `m_r(n)` is

    Q_n(E)m(n)=0,      Q_n(T)=product_(u in S_n)(T-u). (3)

For the zero sequence its degree is zero. For nonempty S_n the degree is
exactly p_n, even for signed or complex weights. Indeed any annihilator P
would give

    sum_u c_u(n) P(u) u^r=0 for all r>=0.

The first p_n equations form a nonsingular Vandermonde system, so every
`P(u)=0`. This proves minimality. Equivalently, the leading p_n-square Hankel
matrix has determinant

    det(m_(i+j)(n)) = product_u c_u(n)
                     product_(u<v)(v-u)^2 != 0.     (4)

An arbitrary ordering of the sites suffices in (4). No positivity is assumed.

**Zero-site boundary.** A surviving site zero contributes only to m_0.
It is counted in (3), whose validity begins at r=0. If one asks only for an
*eventual* recurrence, that one site disappears and its degree drops by one.
The infinite Hankel rank and the all-r recurrence convention used here retain
that initial information. The rational function F_n also retains the pole at
z=0, so no pole is actually removable.

A known bound `|S_n|<=B` makes the first `2B` moments a faithful representation
of the finite weighted source. Two candidates with those moments equal have
a difference supported at at most 2B sites; its square Vandermonde system
forces all difference weights to vanish. At exact size s, the first 2s moments
recover Q_n by the invertible Hankel system, then its roots and the weights.
This is uniqueness over K; an effective root algorithm needs a specified
presentation of K.

The uniform bound is sharp: at 2B distinct sites t_i, the weights
`1/product_(j!=i)(t_i-t_j)` have moments zero through degree 2B-2 and moment
2B-1 equal to one. Partitioning the sites into two B-element sets and moving
one side to the other gives distinct B-site sources agreeing in their first
2B-1 moments. This is the elementary Lagrange interpolation identity.

## 4. An exact affine law for recurrence memory

Partition the original support into translation orbits. In each occupied
orbit C, choose the least integer index so its source is written

    c_0,...,c_d,    c_0*c_d != 0,

at sites `t_C+i a`, with internal zero coefficients allowed. Let

    S_C=sum_(i=0)^d c_i,
    P_(C,j)=sum_(i=0)^j c_i,       0<=j<d.

For every `n>d`, convolution (1) consists of three disjoint blocks:

- prefix coefficients `P_(C,0),...,P_(C,d-1)`;
- `n-d` copies of S_C;
- suffix coefficients `S_C-P_(C,0),...,S_C-P_(C,d-1)`.

Therefore that orbit contributes exactly

    p_(C,n)=(n-d) 1_(S_C!=0)
            +sum_(j=0)^(d-1) [1_(P_(C,j)!=0)
                              +1_(S_C-P_(C,j)!=0)]. (5)

Distinct translation orbits never collide. For
`n>max_C d_C`, summing (5) gives the exact eventual law

    p_n = q n + b,
    q = #{C:S_C != 0},                              (6)

with b the displayed sum of boundary contributions minus
`sum_(S_C!=0)d_C`. The intercept can be negative; (6) is claimed only past the
explicit source-span threshold. The slope q is a nonnegative integer and is
not determined by the grand total `sum_C S_C`.

If every orbit has zero total, define its finitely supported primitive
`b_i=sum_(j<=i)c_j`. Then

    c_i(n)=b_i-b_(i-n).

For n>d the two boundary copies have disjoint support, giving
`p_(C,n)=2 #{j:P_(C,j)!=0}`. This explains bounded memory constructively.
Conversely each nonzero orbit total creates a nonzero plateau of length n-d,
so unbounded memory is unavoidable.

For h=4 the inherited integrability theorem now gives the exact equivalence:

    rational first integral / rational Hamiltonian logarithm
      iff bounded all-r moment-recurrence order along the iterates
      iff q=0.                                     (7)

For other h, the source growth theorem still holds; (7)'s Hamiltonian wording
is stated only for the inherited h=4 two-flow family.

### Small hostile pair

Take a=1. The source `delta_0-delta_1` has one zero-total orbit and p_n=2
for every n>=1. The source `2delta_0-2delta_(1/2)` has two nonzero-total orbits
and p_n=2n. Both have

    m_0(n)=0,    m_1(n)=-n  for every n.

Thus even the complete time histories of these first two moments miss the
difference between bounded and linearly growing exact recurrence memory.

## 5. The second recurrence direction: time instead of jet degree

The same source has a separate triangular recurrence in n:

    m_r(n+1)-m_r(n)
      =sum_(ell=0)^r binom(r,ell) (a n)^(r-ell) m_ell(1). (8)

This follows by expanding the new translated layer `(t+a n)^r` and summing
its original weights. Consequently m_r(n) is a polynomial in n of degree at
most r+1, and

    Delta_n^(r+2) m_r(n)=0.                         (9)

These finite-degree time recurrences coexist with arbitrarily large and
linearly growing all-r recurrence order in (3). They are different axes.
A fixed query may have short memory while the full evolving object does not.
If low original moments vanish, (8) preserves that vanishing at every time.

### An all-time finite-jet obstruction

For any R>=0 set `L=R+1`, a=1, and take

    t_j=j/(R+3),   c_j=(-1)^j binom(L,j),  0<=j<=L. (10)

Finite differences give `m_r(1)=0` for 0<=r<=R. Equation (8) then gives

    m_r(n)=0 for every n>=0 and 0<=r<=R.             (11)

All R+2 sites lie in distinct cosets of Z and have nonzero weights, so
`p_n=(R+2)n` for n>=1. For h=4 this is a nonintegrable exact two-flow word,
while the source-free translation has a polynomial first integral.
Nevertheless every iterate has the same bounded jet as that translation.
In the original coordinates,

    Phi^n(x,y)=(x/(1-nx),
               y(1-nx)^2+(1-nx)^2 x^(-2) F_n(1/x)).

Equation (11) makes the correction start at degree R+3, so the maps agree
through total degree R+2 for **every** n. Their first differing coefficient
is nonzero for n>=1: at r=R+1 equation (8) has constant nonzero increment,
and `m_(R+1)(n)=n*m_(R+1)(1)`.

Thus arbitrarily long observation time cannot rescue a fixed jet cutoff.
A support-size promise, a growing jet cutoff, or the actual orbit-labelled
source is needed. This is not an impossibility theorem for all decision
certificates: the finitely listed orbit totals themselves are short witnesses
for this specified input representation.

## 6. Relation to the recent cycle inventories and Frobenius work

The transfer is precise but limited. Finite exponential sums occur in the
old endpoint-Prony bank, in pole moments, and in finite linear observables.
Their Hankel rank measures the faithful linear shift representation of the
specified readout. It does not recover every structure that generated it.
The recent C9 versus C3×C3 hostile has equal fixed-point moment sequences but
different three-subset inventories: even all moments of *that* observable
cannot replace the missing cycle data. Here, by contrast, all normalized pole
moments do reconstruct the finite weighted pole source because the evaluation
sites themselves are the exponential bases. The observables are different.

The [scalar Frobenius template theorem](../../wronskian_frobenius_20261019/SCALAR_TEMPLATE_DEPTH_CLASSIFICATION.md)
shows that coefficient precision and spatial dilation must also be retained.
Nothing above authorizes reducing characteristic-zero ranks modulo p without
a new audit: distinct sites can collide, source weights can vanish, and the
normalizing binomial factors can become nonunits. For instance the two-site
source `delta_0-delta_p` has characteristic-zero rank two but all moments
vanish after reduction mod p. Its lost separation has valuation v_p(p)=1,
not evidence that the original sources coincide. The inherited Smith-Hankel
framework is the correct route when an integral global recurrence and a
chosen lattice are available. No Frobenius-compatible lift is asserted here.

## 7. A finite cyclic source carrier has a different exact clock

Here change the carrier explicitly: spatial sites are `Z/M`, coefficients are
also `Z/M`, and `M>=2`. This is a finite source vector, not its moment readout.
Let E cyclically translate the coordinates by one, let c be any source vector,
and set `C_n=sum_(j=0)^(n-1) E^j c`. Let d be the least positive spatial period
of c, so `d|M`, and let S_d be the sum of one d-coordinate period, interpreted
modulo M. Then the least positive time period of the entire sequence C_n is

    P = d * M / gcd(M,S_d).                         (12)

Use any integer representative of S_d in the gcd. For the zero vector the
spatial period is d=1, S_d=0, and (12) gives P=1, as required.

Proof: if P is a time period, differencing `C_(n+P)=C_n` gives
`E^(n+P)c=E^n c`, hence d divides P. One d-step sum of translates is the
constant vector with value S_d. For P=d*l, therefore, `C_P=l*S_d*1`, which
vanishes precisely when M divides l*S_d. This gives the displayed minimum.
Conversely that P has `E^P c=c` and `C_P=0`, so
`C_(n+P)=C_P+E^P C_n=C_n` for every n. This proves both directions.

A single unit source has d=M, S_d=1, and least time period M², even though the
spatial shift itself has order M. A constant source with value c has d=1 and
period `M/gcd(M,c)`. For example, constant value two modulo four has period
two; replacing its one-period sum by the grand total zero would lose the
clock. The d-period sum and the coefficient ring must remain explicit.

This gives a source-preserving finite recurrence `C_(n+d)=C_n+S_d*1` and exact
clock; passing to moment measurements may produce a strictly smaller quotient.
It neither contradicts characteristic-zero linear pole growth nor establishes
an arithmetic orbit period before the reduction and observer have been fixed.

## 8. Reproduction and scope

Run `python3 check_pole_recurrences.py` and `python3 -O check_pole_recurrences.py`.
The checker reports the exact check count and finite universes in its JSON receipt (including 162 signed source words of spans zero through four). It uses explicit exceptions rather than Python assertions. It checks
signed source convolution against (5), exact Hankel determinants and recovered
recurrence polynomials, bounded versus unbounded hostiles, the sharp finite
moment threshold, all-time bounded-jet controls, and every cyclic source vector for moduli two through five. It uses only the standard
library and rational arithmetic.

This is not an unrestricted A000568 formula, a tournament-class reconstruction
theorem, a Collatz convergence statement, or a new rational-integrability
classification beyond the inherited family. Its reusable contribution is the
exact source-labelled affine memory law, together with the separation between
time-query recurrences and full object recurrence memory.
