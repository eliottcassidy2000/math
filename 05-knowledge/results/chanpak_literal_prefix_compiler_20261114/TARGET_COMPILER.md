# Try the whole exact target line, not just its least arithmetic point

For a prime cofactor \(p\ge159^2\), choose the largest \(k\) with
\(A_k^2\le p\), where \((u_k,v_k)=(47,56)M^k\),
\(A_k=u_k+2v_k\), and \(M=D_2C_1\).
The terminal 7-dividing construction asks for
\[
                 u_k a+v_k b=p,\quad
                 a\text{ odd},\quad b\text{ even},\quad b\ge2a>0.
\]
The row is primitive odd/even, so invert \(u_k\bmod2v_k\).
All candidate first coordinates are in one progression of
step \(2v_k\), up to \(p/A_k\).
The [boundary analysis](PHYSICAL_BOUNDARY.md) gives at most
49 points on the whole line at this maximal depth. A short
alternative proof puts \(t=u_k/v_k\in(0,1]\). Since
\(A_{k+1}=5u_k+12v_k\) and \(p<A_{k+1}^2\), the interval's
length in multiples of \(2v_k\) is below
\[
    \frac{p}{2v_k A_k}<
    \frac{(5t+12)^2}{2(t+2)}\le\frac{289}{6}<49.
\]
The last rational expression increases on \([0,1]\).
There is no reason to choose just the least progression point.

At **each** point the compiler tests two actual source classes:

1. The complete terrace state recognizer, including length-0 and
   length-13 resets where appropriate. Its structural filters use
   the terminal dyadic slope band and the forced near-boundary
   predecessor; they are not arbitrary search cutoffs.
2. An [arbitrary-root bridge](ARBITRARY_ROOTED_BRIDGE.md), starting
   from the provided bank of strong cores. It factors by actual
   strong components and checks the integer \(R_+\) profile.

For the second class split any proposed rooted source \(T\) at the
strong component of the marked vertex \(v\), writing
\(T=A\to B\to C\) with \(v\in B\).
A path starting at \(v\) cannot visit \(A\). Its selected vertices
in the later \(C\) are an arbitrary path subset with the remainder
handled in the other slot. Thus
\[
 H(T)=H(A)H(B)H(C),\qquad
 R_+(T,v)=H(A)R_+(B,v)F_2(C).
\]
At a prime cofactor, \(\gcd(a,b)=1\), because any common factor
divides \(p\) and \(a<p\). Thus \(A\) consists only of singletons
and can be deleted. The marked \(B\) may be a singleton or
any strong core from the bank with a marked root. Divide out its
two responses \(H(B),R_+(B,v)\), then find an ordinal product
for the remaining two responses using strong components of \(C\).
Singletons of \(C\) are absorbed into the bridge's width.
The recursion decreases \(a\) on every nonsingleton component.
A complete strong inventory up to the count-dependent Busch bound
would decide the full bridge-source class; a small named bank
can build successful points without enumerating big tournaments.
Neither class need meet every target line.

For a successful point, the code outputs the literal prefix and
all arcs of the result after \(k\) copies of \(q=1,r=1\)
and one \(q=7,r=2\) seam. Two examples exercise different repairs:

| prime cofactor | actual source coordinates | repair | final total |
|---:|---|---|---:|
| 1,137,067 | \((1105,2634)\) | bypass the unphysical least \((469,3176)\) on the same line | 7,959,469 |
| 1,512,281 | \((23,4736)\) | change to an arbitrary-root core; this pair is impossible on the terrace | 10,585,967 |

The first full witness has 28 vertices; the second has 25.
Both prefix profiles were checked by direct subset path DP,
and the exact bridge transfer verifies the full totals.
The output adjacency masks are zero-based in both witnesses.

Run `python construct_target.py 1137067 1512281`.
For any other target the default bank may return `null`;
that is not a claim of nonattainment in the full spectrum.
