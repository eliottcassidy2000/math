# Binary carries classify H-neutral native improvements

Status: **PROVED scoped**, extending the balanced two-apex sources and full edge-response formula of `apex_wedge_constructor_20261103/THEORY.md`, `PROFILE.md` and `CARRIER.md`. No inverse for arbitrary H is asserted.

## D1. Literal source and complete carry formula

For positive P,R let p be the ordinary leading-one binary expansion of P and s the reverse-complement of the leading-one binary expansion of R. Put

    w=p01s,  u=p10s.

Take the transitive base chain0->1->...->m-1. The apex a beats base vertex i iff w_i=1; b beats it iff u_i=1. Let T+ have a->b and T- have b->a. Both sources are strong, and the inherited scalar jump is

    H(T+)-H(T-)=2(ceil(R/2)-ceil(P/2)).

For every odd P, set R=P+1. The flip is H-neutral. Write k=bitlength(P) and j=v_2(P+1), the number of trailing ones. With F the ordered-bin path-cover polynomial, the complete response is

    F_(T+)-F_(T-) = -r^a (r^2-1)[(r+1)^(j-1)-(r-1)^(j-1)],       (D1.1)

where

    a=2k-j if P!=2^k-1,   a=k+1 if P=2^k-1.

In particular the full profile is unchanged **iff** j=1, or equivalently P=1 mod4. Scalar equality holds for every odd P. This gives a complete classification inside the specified balanced adjacent-chart family.

Proof. Put L_h=((r+1)^h+(r-1)^h)/2 and, for binary digits n_i of N,

    C_r(N)=sum_i n_i r^i L_(bitlength(N)-1-i).

The inherited balanced response is

    2r[r^k C_r(R)-r^bitlength(R) C_r(P)].

For a nonoverflow carry, unchanged prefix digits cancel. The remaining difference is

    C_r(P+1)-C_r(P)
      =r^(k-j-1)[L_j-r sum_(h=0)^(j-1) r^(j-1-h)L_h]
      =-r^(k-j-1)(r^2-1)[(r+1)^(j-1)-(r-1)^(j-1)]/2.

The second equality follows by telescoping the two geometric sums with roots r+1 and r-1. Multiplication by2r^(k+1) gives D1.1. In the overflow case C_r(2^k)=L_k and C_r(2^k-1)=sum_(h=0)^(k-1)r^(k-1-h)L_h; the same identity gives a=k+1. For j=1 the bracket vanishes; for j>=2 it is a nonzero polynomial with nonnegative coefficients.

## D2. Dominance survives every homogeneous module context

For j>=2 write g=F_(T-)-F_(T+), the positive of D1.1. Its coefficients in the falling-factorial basis(r)_k are nonnegative integers. To see this, use

    r(r^2-1)=(r)_3+3(r)_2,
    r(r)_k=(r)_(k+1)+k(r)_k.

Here a>=1 and the remaining bracket has nonnegative integer coefficients in ordinary powers of r. Repeated multiplication by r preserves the nonnegative falling-factorial cone. Consequently

    p_k(T-)>=p_k(T+) for every k,
    p_1(T-)=p_1(T+),    p_2(T-)>p_2(T+).

The strict last assertion follows from g(0)=g(1)=0 and g(2)>0. This proves more than pointwise positivity at a selected integer r.

In any fixed tournament context where the substituted source is a homogeneous module, a Hamilton path breaks its visits into k nonempty source paths. The number of ways to connect these visits through the outside vertices is a nonnegative coefficient multiplying the number of ordered k-path covers k!p_k. Therefore replacing T+ by T- weakly increases H in **every** such context. The same statement holds for path-cover observations by the corresponding visit decomposition. This is an order on responses for a specified operation language, not an assertion about arbitrary named-port attachments.

The triangle wrapper G strictly distinguishes the pair:

    H(G(T-))-H(G(T+))=g(1)+g(2)=g(2)>0.

The inherited operator for G preserves the nonnegative falling-factorial cone. Its coefficients satisfy

    p'_k=k(k+1)p_(k+1)+k(k+2)p_k+(2k-1)p_(k-1)+p_(k-2)

with absent indices zero and p_0=0 for nonempty sources. Thus every additional G iteration preserves strict scalar improvement after the first one.

## D3. Physical ports and exact native cost

The sources start with w_0=u_0=1. The inherited reset introduces one helper x beating only a, and uses the two actual Hamilton paths

    P=(a,b,0,1,...,m-1,x),
    Q=(a,x,b,0,1,...,m-1).

Apply native path complement along P, then along Q, then delete the helper. The endpoint is exactly T- on the original labels. The physical apexes a,b are retained, and the vertex clock is(n,n+1,n+1,n+1,n). This is a two-native-move program with one temporary birth and deletion; it does not mean that H is constant at the intermediate graphs.

The following sources are concrete:

| P,R | Order | Common H | F_(T-)-F_(T+) | Gain after G |
|---|---:|---:|---|---:|
| 1,2 | 7 | 89 | 0 | 0 |
| 3,4 | 9 | 497=7*71 | 2r^3(r^2-1) | 48 |
| 5,6 | 10 | 2389 | 0 | 0 |
| 7,8 | 11 | 2741 | 4r^5(r^2-1) | 384 |
| 11,12 | 12 | 19189 | 2r^6(r^2-1) | 384 |

The order-ten example is a full-profile-neutral physical reset, while the neighboring carry classes supply strict improvement. Both belong to one binary-carry law, with overflow explicitly retaining its extra bitlength clock.

`neutral_carries.py` checks all odd P through255 symbolically, reconstructs every response in the falling basis, and independently replays all128 native programs. Endpoint subset DP checks the five literal pairs in the table; complete path-cover counting checks P=1,3,5; literal G counts check P=1,3,7. The all-order proof is the carry identity and cone argument above.

## Next use

These moves give controlled variation inside an H fibre. They can repair a supplied source's future response while preserving both old apex names. The remaining constructive obligation is to select a pair that changes a desired rooted return or divisibility condition, then match the next source again. Homogeneous response dominance alone does not establish that rooted closure. The 81/27 and failed55->111 controls remain relevant precisely at this boundary.
