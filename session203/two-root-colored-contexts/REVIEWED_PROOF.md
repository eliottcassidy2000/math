# Two sensitive roots retain the repair exactly on ordered color blocks

Status: **PROOF CANDIDATE, pending independent audit.** The statement is an all-order identity and equality classification for two explicit physical operations. It is not an obstruction to all two-root contexts, and it proves no eventual odd-spectrum threshold.

## 1. Two physical promotions of the audited repairs

Let T be any tournament whose vertices have colors0 or1; T may be empty and need not be transitive. A context consists of a fixed core Q and two masks M0,M1. Retain T and Q disjointly. A body vertex of color i beats precisely the core vertices selected by Mi; all remaining cross edges point from the core into that body vertex. Masks and outgoing adjacency rows use vertices numbered from0.

Use either of these parent/child pairs:

- A parent core (2,24,3,37,13,23), masks(15,63); child core (2,8,67,101,15,87,19), masks(79,127).
- B parent core (32,5,1,39,47,70,31), masks(0,5); child core (0,5,129,39,47,199,159,27), masks(0,5).

In A, exactly old vertices4 and5 are sensitive to the body color. The child reverses the old arc1->4 and adds vertex6. In B, exactly old vertices0 and2 are sensitive. The child reverses0->5 and adds vertex7. Every body edge and every retained core/body incidence is unchanged. These are actual one-old-arc-reversal/one-birth operations, not differences of formal graphs.

The monochromatic cases inherit the previously audited transitive and arbitrary-module repairs. The point of promoting two roots is to allow a freely mixed color word on a common spine, beyond the one-sensitive-root obstruction. No uniformity or transitivity assumption is imposed on the colored body here.

Write P_A(T),C_A(T),P_B(T),C_B(T) for the four actual contexts. Put H(empty)=1.

**Theorem.** For each pair i in {A,B}, every finite colored tournament T satisfies

    H(C_i(T)) >= 2 H(P_i(T)) + H(T).               (1)

Equality holds if and only if EVERY color1 vertex beats EVERY color0 vertex. Empty color classes are allowed.

Consequently, if T is a transitive spine with colors in its increasing order, the desired +1 repair holds exactly for words1*0*. Every occurrence of a color0 vertex before a later color1 vertex makes the correction strictly larger than1. The two explicit promotions therefore do not turn the monochromatic repair into a free binary encoder.

This does not exclude other two-root contexts, changing the old spine edges, target-dependent choices, or using the strictly positive excess as an additional controlled parameter.

## 2. The colored path-cover algebra

For a colored tournament T, let x,y,z,w count its Hamilton paths whose first and last colors are, respectively,00,01,10,11. Each function is zero on the empty tournament. Let1 denote the function equal to1 only on the empty tournament and0 otherwise.

Products below are **subset convolution**, not multiplication of the scalar values on T:

    (fg)(T) = sum_(S subset V(T)) f(T[S]) g(T[V(T) minus S]).  (2)

This product is associative and commutative. For instance xy counts ordered spanning two-path covers whose first named path has endpoint colors00 and whose second named path has endpoint colors01. A monomial of degree r counts r labeled nonempty paths covering the body. Equal letters still give distinct labeled slots. Thus ordinary commutative polynomial notation encodes exact, nonnegative path-cover counts. The constant function1 in this algebra is the empty-set unit, not the numerical constant1 on every input.

In this notation the ordinary Hamilton-count function is

    H = 1+x+y+z+w.                                 (3)

Nonnegative functions have nonnegative convolution products.

## 3. A strictly positive cut defect

Adjoin an apex q beating every color0 vertex and losing to every color1 vertex. Define

    W(T)=H(T+q)-H(T).

Splitting Hamilton words according to whether q is first, last, or internal gives the exact polynomial

    W = y-z+xy+xw+y^2+yw.                          (4)

At the empty body both sides are0. At a nonempty body, the first and last cases contribute x+y and y+w; the internal case contributes (y+w)(x+y). Subtract (3) to obtain(4).

**Cut lemma.** W(T)>=0, with equality exactly when every color1 vertex beats every color0 vertex.

Here is a self-contained positive form. For distinct a,b in T, let Q_ab count ordered spanning nonempty two-path covers whose first path ends at a and whose second path starts at b. Let epsilon_v be1 if q beats v and0 otherwise. Then

    H(T+q)-H(T)
      = sum_(a->b in T) (Q_ab+Q_ba) epsilon_a(1-epsilon_b).   (5)

To verify(5), direct insertion gives

    H(T+q) = sum_v epsilon_v start_v
             + sum_v (1-epsilon_v) end_v
             + sum_(a!=b) (1-epsilon_a)epsilon_b Q_ab.

The constant term is H(T). For each vertex a, the path-flow identity

    start_a-end_a
      = sum_(a->b) Q_ab - sum_(b->a) Q_ba

follows by counting outgoing minus incoming used path edges at a. It identifies the linear and quadratic coefficients with(5).

Every capacity Q_ab+Q_ba is positive. Given a Hamilton path on T minus {a,b}, find the first vertex beaten by either a or b. Use that endpoint immediately before the remaining suffix and put the other endpoint after the preceding prefix, producing one of the two required spanning covers. If neither endpoint beats a path vertex, put the entire path before one endpoint and leave the other endpoint as a singleton second path. The empty remainder also gives two singleton paths. Thus at least one of the two covers exists.

In(5), epsilon_a(1-epsilon_b)=1 precisely for a directed edge from color0 to color1. This proves the cut lemma. It is the inherited weighted-cut mechanism, included here to make the positivity and equality boundary explicit.

## 4. Exact excess identities

Let

    D_i = H(C_i)-2H(P_i)-H,
    F = y+xy+y^2,
    J = y+y^2+yw.

The exact identities in the convolution algebra are

    D_A = (4+x)W + (x+z)F
          + y(4+3x+4y+3z+4w),                     (6)

    D_B = (21+42y+20w)W + 20(z+w)J + 22yH.        (7)

The displayed integer constants multiply the empty-set unit. Thus, for example, the first summand of(6) contributes exactly4W(T)+(xW)(T).

Equations(6)-(7) are finite-core combinatorial identities valid for every colored body, including nontransitive bodies of arbitrary order. Their complete coefficient certificate is supplied in NORMAL.json and recomputed by verify_colored_identity.py. They are not inferred from fitting graph counts at selected body sizes.

For clarity, the expanded excess coefficients provide a compact independent algebra check:

    D_A = 8y-4z
          +9xy-xz+8y^2+4yz+4xw+8yw
          +2x^2y+2xy^2+xyz+y^2z+x^2w+xyw;

    D_B = 43y-21z
          +43xy+85y^2+83yw+21xw-20zw
          +42xy^2+42y^3+20y^2z+82y^2w
          +62xyw+20yzw+40yw^2+20xw^2.

Substitution of(4) and the definitions of F,J,H into(6)-(7) gives these expansions exactly.

### How the full core certificate is obtained

Fix an order of the core vertices. Its n+1 linear gaps may be empty or occupied by a nonempty body path. An interior gap may be empty exactly when its neighboring core vertices form a forward edge; either end gap may always be empty. A body path with first color f and last color l may occupy a gap exactly when its left core endpoint, if any, beats color f and its right core endpoint, if any, loses to color l.

Assign to each gap its zero-one empty permission plus the sum of its permitted endpoint-color variables x,y,z,w. Multiply these gap polynomials and sum over every core order. Every Hamilton path in the resulting context has a unique core order, unique occupied gaps, and a unique labeled spanning body-path cover. Conversely every permitted cover reconstructs a Hamilton path. This is a bijection for arbitrary internal body edges, proving that the polynomial is the exact context count.

The verifier enumerates every order of the four fixed cores (orders6,7,7,8), performs this polynomial calculation with integer coefficients, subtracts2 times the parent tensor and(3), and verifies precisely(6)-(7). A separately implemented C++ certificate retains ordered endpoint-type words before commutative aggregation and gives the same coefficients.

## 5. Positivity and the complete equality case

Every factor in(6)-(7), except W before applying the cut lemma, has visibly nonnegative coefficients. The cut lemma makes W nonnegative on every induced colored subgraph, so all convolution terms are nonnegative. This proves(1).

If any edge points from color0 to color1, then W(T)>0. Equations(6)-(7) give D_A(T)>=4W(T)>0 and D_B(T)>=21W(T)>0, so equality fails.

Conversely, suppose all cross edges point from color1 to color0. The same holds in every induced subgraph. W is consequently zero on every induced subgraph. No Hamilton path can start in color0 and end in color1, so y is also zero on every induced subgraph. Every term in(6)-(7) contains W or y, directly or through F or J. All those terms vanish, proving equality.

On a transitive spine, this condition is exactly that no0 precedes a1, hence words1*0*. In A the color1 block is then a global initial module, and in B the color0 block is a global terminal module. This explains why the surviving mixed cases reduce to ordered products of the original homogeneous repairs rather than supplying new free digits.

## 6. Hostile controls, scope, and reproduction

The first unsorted spine has body rows(2,0), colors(0,1), and H(body)=1. In A, its actual parent and child counts are133 and279, giving excess12 over2H(parent)+1. In B, the counts are115 and295, giving excess64. Literal rows are retained in PROMOTED_TWO_ROOTS.json.

The main verifier independently evaluates the convolution formulas on75 exhaustive colored inputs through body order3 and24 fixed-seed larger inputs, then compares them with literal graph-subset DP for both physical contexts. It checks198 graph pairs, including108 strict cases. The core tensor computation, rather than those finite body tests, establishes the arbitrary-order identity.

    python3 verify_colored_identity.py > NORMAL.json
    python3 -O verify_colored_identity.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

The exact operations inherit their core rows from the audited Session203 transitive-repair packet. This result does not classify all two-sensitive-root gates or show any globally missing scalar count. It records exactly why these two natural shared-spine promotions fail to isolate a unit correction, and what nontransitive bodies still preserve the equality.
