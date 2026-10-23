# Powering, discriminant square classes, and the recovered Fibonacci forest operator

Status: PROVED polynomial identities and elementary applications; named transcendence, VOA and matrix-forest inputs remain CITED. This is a separate addendum; the initial Coulomb note is unchanged.

## 1. Powering descends through the reciprocal quotient

Write C(t)=13−6(t+t^−1), Δ(c)=(c−1)(c−25), and x=(13−c)/12. For every positive integer m, put

    P_m(c)=13−12T_m((13−c)/12),

where T_m is the Chebyshev polynomial defined by T_m((t+t^−1)/2)=(t^m+t^−m)/2. Then

    C(t^m)=P_m(C(t)),
    P_m(P_n(c))=P_(mn)(c).

In particular P_1(c)=c and P_2(c)=25−(13−c)^2/6. These are polynomial maps over Q, not branch-dependent choices of square roots.

The second-kind Chebyshev identity

    (t^m−t^−m)/(t−t^−1)=U_(m−1)((t+t^−1)/2)

and Δ(C(t))=36(t−t^−1)^2 give the polynomial identity

    Δ(P_m(c))=Δ(c) U_(m−1)((13−c)/12)^2.

Although the quotient expression is initially written away from t=±1, the polynomial identity extends to the endpoints by continuity or polynomial equality.

For rational c<1 we have x>1. Choose t>1; then the quotient defining U_(m−1)(x) is positive, hence nonzero. Therefore Δ(P_m(c)) and Δ(c) represent the same nonzero square class in Q×/(Q×)^2. Also t^m>1, so P_m(c)<1. The algebraic-versus-transcendental mirror-base dichotomy of the preceding note survives every positive integer powering operation.

Load-bearing boundary: this preservation is not claimed for every rational c. For c=13, Δ(c)=−144 is nonsquare, but P_2(13)=25 has zero discriminant because U_1(0)=0. Similarly m=0 collapses every input to the branch value1, so the positive-m hypothesis matters.

## 2. Golden-ratio orbit and exact integer recurrences

Let φ=(1+sqrt5)/2 and t=φ²=(3+sqrt5)/2, so C(t)=−5. Standard Fibonacci/Lucas initial conditions F0=0,F1=1 and L0=2,L1=1 give

    t^m+t^−m=L_(2m),
    t^m−t^−m=sqrt5 F_(2m).

Hence

    c_m=P_m(−5)=13−6L_(2m),
    Δ(c_m)=180 F_(2m)^2.

The first c_m, m≥1, are−5,−29,−95,−269,−725. Set c0=1; then

    c_(m+1)=3c_m−c_(m−1)−13.

The square class is always5, since180=36·5. Thus every c_m has quadratic irrational reciprocal parameters t^±m and a pair of transcendental positive mirror bases. This is an infinite exact family, not an extrapolation from its first terms.

Since c_m are negative integers, the simple Virasoro VOAs L(c_m,0) are all nonrational by the integer-central-charge corollary in the preceding note. This statement does not extend to arbitrary VOAs with those charges.

## 3. The exact recovered square-sum/forest connection

The prior proved note `05-knowledge/results/session32_square_monster_20261023/CAT_MAP_PATH_FORESTS.md` identifies the path-forest transfer matrix

    B=[[3,−1],[1,0]]

with the cat map A=[[2,1],[1,1]] by an explicit integral conjugacy. Its distinct observables are

    f_m(1)=det(I+L(P_m))=F_(2m),
    tr(A^m)=L_(2m),

where P_m here denotes the undirected path on m vertices; this graph notation is separate from the polynomial map P_m(c) above. To avoid confusion one may call the graph Path_m.

Therefore the present bridge is the explicit pair of identities

    c_m =13−6 tr(A^m),
    Δ(c_m)=180 det(I+L(Path_m))^2.

The trace supplies the central-charge coordinate; the rooted-forest count supplies its discriminant factor. They are distinct observables of the same transfer operator, not equal counts. For example m=4 gives trace47 and forest count21, hence c4=−269 and Δ(c4)=180·21².

This transports an actual invariant and arithmetic predicate: powering the transfer operator preserves the discriminant square class, which determines the algebraic/transcendental type of the corresponding mirror pair. It does **not** prove a VOA realization attached to the path, nor Hamiltonian existence for arbitrary square-sum graphs, nor a Monster correspondence. The inherited graph note's endpoint/port caveats remain intact.

## 4. Small integer classification and a simpler c=−5 rational VOA

The complete rational-parameter list for integer c<1 is

    c∈{0,−2,−7,−24}, with t>1 respectively3/2,2,3,6.

Indeed Δ(c)=(13−c)^2−144 is a square iff (13−c−d)(13−c+d)=144 for nonnegative integer d. Same-parity positive factor pairs, excluding the endpoint c=1, give exactly the four values. All other integer c<1 yield transcendental mirror pairs. Only c=0 among the four is a rational simple Virasoro model.

The earlier many-factor rational VOA at c=−5 can be replaced by the simpler

    L(c_(2,5),0) ⊗ L(c_(3,5),0),
    c_(2,5)=−22/5, c_(3,5)=−3/5.

Their sum is−5. Each is a minimal model and their tensor product is regular, hence rational. Dong–Li–Mason prove tensor closure in Proposition3.3 and regularity of minimal models in Theorem3.13: [primary paper](https://arxiv.org/pdf/q-alg/9508018). This simplifies the example without retracting the valid earlier construction.

For a verified accessible statement and proof of the transcendence input, see Richard Elman's author-hosted [Lectures on Abstract Algebra, §74](https://www.math.ucla.edu/~rse/algebra_book.pdf). The Gelfond–Schneider hypotheses here are precisely algebraic base t≠0,1 and algebraic irrational exponent1/(t−1), with the positive real branch.

Reproduce: `python3 session35-lie/powering_checks.py`. The exact checks verify the polynomial recurrences on a stated rational grid, the matrix/Fibonacci identities through m=30, and the hostile c=13 boundary. The all-m statements follow from the displayed symbolic identities, not finite extrapolation.
