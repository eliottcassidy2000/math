# Pinned inherited divisor-balance theorem excerpt

Source: 05-knowledge/results/trijection_reset_20261020/arithmetic.md, section1, snapshot 8f76b84d2960f9e14ad62a154b8533c77593702f.

For `N>=2`, let `D(N)={d:d|N,1<d<N}`. Define `F=|D|`, `S` as the number
of squarefree elements of `D`, and `U` as the number of prime elements of
`D`. Thus the prime set is contained in the squarefree set. `F=U+S` is an
equality of counts, not the literal disjoint union of those two subsets.

For `N=prod p_i^a_i` with `r` distinct primes, the inherited formulas are

    F=prod(a_i+1)-2,
    S=2^r-1-[N squarefree],
    U=r-[N prime].

They give the complete inherited classification

    F=U+S  iff  N=p, p^3, or p^2 q r,

with three distinct primes in the last shape. For completeness, a
nonsquarefree solution needs `prod(a_i+1)=2^r+r+1`. At `r=1,2,3`, this
forces respectively `(3)`, no profile, and `(2,1,1)`; at `r>=4`,
`prod(a_i+1)>=3*2^(r-1)>2^r+r+1`. Squarefree composites fail because
`F=S` and `U>0`; primes have all three counts zero.

The expression `{p^q r,p,p^3}` must retain its literal typing. If `p,r`
are distinct primes and `q>=2` is an integer exponent, then `p^q r` has
`(F,S,U)=(2q,3,2)`, so equality would require `2q=5`. It is never a
solution. At `q=1`, the squarefree product also fails. The recovered
three-shape classification uses `p^2 q r`, with `q` a separate prime
factor. These two strings describe different families.

