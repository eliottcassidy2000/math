# A digit formula for physical retained-root blocks

Exploratory positive reduction, Session198. The formulas below are proved from the physical O/Z grammar. They do not prove cofinite scalar coverage and are not a new obstruction theorem.

The physical input is the reviewed [three-letter source grammar, Theorem1 and equations(4)](https://github.com/eliottcassidy2000/math/blob/6f3f93a3915ef2a7eba9a611f67a1bb75bd392c9/session198/bilinear-descent/PROOF.md). Its Theorem1 proves completeness of the O/Z/J orbit of the two-vertex source. The C2 identities below are the fixed-left specialization of the reviewed [complete C2 bilinear theorem](https://github.com/eliottcassidy2000/math/blob/6f3f93a3915ef2a7eba9a611f67a1bb75bd392c9/session198/bilinear-descent/dependencies/BILINEAR_C2.md), equivalently Section4 of that source-grammar proof.

Let a nontrivial physical source have profile s=(D,E,U,d), and apply a chronological word of k>=1 retained-root letters O,Z. Write epsilon_i=1 for Z and 0 for O, i=0,...,k-1. Define

    A=sum epsilon_i 2^i, B=2^k, R=rev_k(A),
    p=A mod2, t=floor(A/2^(k-1)).

Here rev_k reverses exactly k binary digits, retaining leading zeroes. Every integer 0<=A<B specifies an actual word. Let g be the already proved one-apex scalar function:

    g(1)=1,
    g(2l)=g(l),
    g(2l-1)=g(l)+l-(l mod2) for l>=2.

## Exact block formula

The resulting actual source profile is

    D'=D+(U-d)A+d(g(A+1)-1)+(d-E)p,
    E'=dt,
    U'=BU-dR,
    d'=d.                                      (1)

The physical interpretation is restricted to nontrivial sources: Z at the singleton does not preserve the start-one condition. No relaxed integer-cone membership is assumed.

To prove it, let j_i be the integer represented by epsilon_0,...,epsilon_(i-1) in chronological order. The physical recurrences give U_i=2^i U-d j_i. For i>=1, E_i=d epsilon_(i-1). At a Z letter, D increases by U_i-E_i; at an O letter it is unchanged. Hence

    D'=D+UA-E epsilon_0-d C(A),
    C(A)=sum_(h<i) epsilon_h epsilon_i 2^(i-h-1)
         +sum_(i>=1) epsilon_(i-1)epsilon_i.

Put f(n)=g(n+1). The recurrence for g is equivalently

    f(2x+1)=f(x),
    f(2x)=f(x)+x+(x mod2), f(0)=1.

Expanding this over the bits gives

    f(A)-1=sum_(h<i)(1-epsilon_h)epsilon_i 2^(i-h-1)
           +sum_(i>=1)(1-epsilon_(i-1))epsilon_i.

Adding the expressions for C and f-1 yields A-p, because the contribution of each set bit i>=1 is (2^i-1)+1=2^i. Thus C=A-p-g(A+1)+1, proving D' in (1). The U,E,d formulas are immediate from their recurrences.

## Complement identity

For C=B-1-A, the two one-apex values satisfy

    g(A+1)-g(B-A)=A-R+t-p.                    (2)

Subtract the complemented-bit expansion from the preceding original expansion. Each pair h<i contributes (epsilon_i-epsilon_h)2^(i-h-1); their sum is A-R. The adjacent correction telescopes to t-p. This proves (2) for every k>=1 and 0<=A<B, including leading-zero cases.

Equation (2) cancels the quadratic bit interactions in a difference. It does not by itself construct a tournament whose Hamilton-path count equals that difference.

## Explicit fixed-partner arithmetic readouts

Fix the actual two-vertex source tau=(1,0,2,1). Its two physical C2 counts on a nontrivial source w are

    N(tau,w)=3D-2E+6U-2d,
    M(tau,w)=5D-3E+13U-2d.

Applying (1) to s gives the following two arithmetic values, each physically realized whenever the output source is nontrivial:

    N_block(s;k,A)
      =3D+3(U-d)A+6BU
       +d[3g(A+1)-6R-2t-5]+3(d-E)p,

    M_block(s;k,A)
      =5D+5(U-d)A+13BU
       +d[5g(A+1)-13R-3t-7]+5(d-E)p.          (3)

A zero-length retained-root block is the identity and is treated separately, rather than assigning a top bit to it.

Together with the genuine source operation

    J(D,E,U,d)=(D-E+U,D-E,D+U,D),

formula (1) gives an exact free-integer parameterization of the complete source language: start at tau and alternate arbitrary retained-root blocks with J. Every allowed parameter tuple is physical, and every nontrivial physical source admits such a tuple by splitting its O/Z/J word at J letters. Formula (3) therefore eliminates source membership constraints from the fixed-left C2 scalar problem. The remaining problem is arithmetic range coverage of these explicitly specified digit expressions, not mere residue coverage.

A sufficient unproved statement for the original coverage objective is that the two fixed-left ranges in (3), allowing all such physically generated prefix states s, contain every sufficiently large odd integer. A separate finite certificate would then handle the remaining targets using the known C1/C2/C3 witnesses. The reduction supplies no conductor and does not assume that the observed bounded conductor is the true one.
