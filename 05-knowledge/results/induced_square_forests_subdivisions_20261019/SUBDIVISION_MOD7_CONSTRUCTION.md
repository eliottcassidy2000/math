# Every finite graph has an induced square-sum two-subdivision at quartic height

Status: **PROVED**, independently [reviewed](review/INDEPENDENT_REVIEW.md); independent exact controls accompany this note. Date: 2026-10-19. No external priority claim.

## Statement and inherited interface

For every finite simple graph G with n >= 1 vertices and m edges, write K=n+2m for the number of vertices in its exact two-subdivision: replace each edge by a path of length three, with two new internal vertices private to that edge. There is an induced embedding of this graph in the positive-integer square-sum graph with maximum label O(K^4), with an absolute implied constant. All intended edges and all nonedges are preserved.

The inherited two-vertex bridge identity is

    v = r^2-B,  s = (A+B+1-r^2)/2,  w = s^2-A,
    B+v = r^2,  v+w = (s-1)^2,  w+A = s^2.

The new ingredient is simultaneous use of three residue classes modulo seven. These remove all within-family edges without the divisor-function bound. A common parameter window then removes every off-diagonal cross-family edge by the gap between consecutive squares.

## Explicit construction

Index the original vertices by i=0,...,n-1 and give them labels

    a_i = 14i+3,       H = 14n-11.

Index the edges arbitrarily by e=0,...,m-1 and orient each arbitrarily, writing its endpoint labels as A_e,B_e. Let R be the least integer congruent to 1 modulo 14 satisfying

    R >= 100(H+m+1).

Thus R < 100(H+m+1)+14. Put

    r_e = R+14e,
    v_e = r_e^2-B_e,
    s_e = (A_e+B_e+1-r_e^2)/2,
    w_e = s_e^2-A_e.

Each s_e is an integer, since A_e,B_e,r_e are odd. The new path for edge e is

    A_e -- w_e -- v_e -- B_e.

For m=0 there is nothing to add. The empty graph on zero vertices is trivial separately.

## Arithmetic ranges and injectivity

We have H <= R/100 and m <= R/100, and

    R <= r_e <= R+14m <= 1.14R,
    t_e := -s_e = (r_e^2-A_e-B_e-1)/2 > 0.

In particular t_e >= (R^2-2H-1)/2. Since R >= 400, this is greater than R^2/3. Consequently

    H < v_e < 2R^2 < w_e < 5R^4.

For the upper bound, t_e < r_e^2/2 < 2R^2, so w_e < 4R^4; the displayed 5 is deliberately loose. The lower bound follows from w_e > R^4/9-H > 2R^2.

For e > f, r_e^2-r_f^2 >= 28R. Therefore

    v_e-v_f >= 28R-H > 0,
    t_e-t_f >= 14R-H > 1.

Also t_e+t_f >= R^2-2H-1 > H. Hence

    w_e-w_f = (t_e-t_f)(t_e+t_f)-(A_e-A_f) > 0.

The V and W families are separately injective, and their disjoint size bands separate them from one another and from all original labels. Thus all K labels are distinct.

## Complete nonedge audit

The square residues modulo seven are 0,1,2,4. The three families satisfy

    original labels: 3 modulo 7,
    V labels:        5 modulo 7,
    W labels:        6 modulo 7.

Indeed r_e is 1 modulo 7, and s_e is 3 modulo 7. Sums of two distinct original labels, two V labels, or two W labels have respective residues 6,3,5 modulo seven. All these pairs are nonedges.

For an original label y != B_e,

    v_e+y = r_e^2+(y-B_e).

The deviation is nonzero and has absolute value less than H < 2r_e-1. Thus this sum lies strictly between the neighboring squares around r_e^2 and is nonsquare. The intended pair y=B_e is a square.

Likewise, for an original label y != A_e,

    w_e+y = s_e^2+(y-A_e)

has nonzero deviation of absolute value less than H < 2|s_e|-1, so is nonsquare. The intended pair y=A_e is a square.

It remains to inspect V-W pairs. The same-edge pair is a square by the bridge identity. For e != f,

    v_e+w_f = (s_f-1)^2+Delta_ef,
    Delta_ef = r_e^2-r_f^2+B_f-B_e.

The minimum nonzero squared-parameter gap is at least 28R, while |B_f-B_e| < H, so Delta_ef != 0. Moreover

    |Delta_ef| <= 28Rm+196m^2+H
                 <= 0.2996R^2+H
                  < R^2-2H.

The strict final inequality follows from H <= R/100 and R >= 400. Since s_f<0,

    2|s_f-1|-1 = r_f^2-A_f-B_f >= R^2-2H.

Thus the nonzero deviation is smaller than the lower adjacent-square gap at (s_f-1)^2, so this off-diagonal sum is nonsquare. This exhausts all pairs and proves inducedness.

## Uniform bound

Every label is below

    5[100(H+m+1)+14]^4
      <= 5[100(14n+m+1)+14]^4.

Since n,m <= K and n>=1, this is an absolute O(K^4) bound. This improves the generic common-window divisor argument's O_eta(K^(4+eta)) specifically for the free placement of the original vertices. It does not prescribe arbitrary original labels; for arbitrary fixed seeds the generic divisor-avoidance theorem remains useful.

No near-quadratic bound for two-subdivisions follows here, and no arbitrary graph is claimed to embed without subdivision. Subdivision does not preserve chromatic number.

## Independent exact controls

Run:

    python3 -B review/verify_mod7_subdivision.py
    python3 -B -O review/verify_mod7_subdivision.py

The controls inspect every pair of constructed labels, checking the equivalence between a square sum and an intended edge. They cover every labeled simple graph through order five; complete graphs of orders six through twenty-four; one hundred deterministic randomized edge lists, orientations and orders on at most forty original vertices; and the generic bridge against every label of several full interval seeds, allowing equal endpoints for triangle attachments. Checks are explicit exceptions and remain active under Python optimization. Finite tests support the implementation; the all-graph claim is proved above.
