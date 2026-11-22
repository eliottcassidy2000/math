# Sharp Pell approximation with optimal E61 pruning on both spines

Call a primitive positive normal \((a,b)\) optimal if the sum \(ax+by\) attains the smallest possible residue support on **each** canonical E61 spine under the fixed-parity projection modulo \(2h\). These are the actual [projected Farey congruences](SQUARE-PRIMES.md). Let \(\theta=1/\sqrt2\).

**PROVED with finite phase inputs.** Among optimal directions converging to \(\theta\) with \(b\to\infty\), the one-sided lower limits are
\[
 \liminf_{a/b\to\theta^-} b^2(\theta-a/b)=\frac1{\sqrt2},\qquad
 \liminf_{a/b\to\theta^+} b^2(a/b-\theta)=\frac{31}{\sqrt2}.
\]
Both are attained as limits of explicit sequences of optimal normals. The least absolute defects \(|2a^2-b^2|\) on the two sides are 2 and 62, respectively. This calculates the exact cost of keeping the strongest modular sieve while approximating one specific irrational tangent. It does not prescribe its absolute weighted sum.

## From an odd square to the Collatz direction

Put \(c=b/2\). The pair \((a,b)=(7,10)\) satisfies \(a^2+1=2c^2\) with \(c=5\). Among positive odd \(c\), also imposing the shortened minus-Collatz step \((3c-1)/2=a\) gives \(c^2-6c+5=0\). Its only values are 1 and 5. At 5 the next two steps are \(7\mapsto10\mapsto5\). Thus the nontrivial intersection of this Pell conic with the actual first cycle step gives precisely the small optimal normal. At exponent 61 the prime 17 distinguishes the two signs of the Pell defect, even though their unconstrained real slopes approach the same limit.

## The local constraints and the absolute defect

**FINITE-EXACT inputs.** Intersect the minimizing projective directions from the [two finite tables](square-directions.py). At both 2 and 5 the only possible normal is \([1:0]\). At 17 the allowed \(b/a\) are the eight nonzero squares. The local test uses the shared physical parity choices before taking the minimum. The [replay](pell_directions.py) also checks, by cardinality, that the minimizing tuples on each spine are the full product of their local projections.

Thus \(a\) is odd, \(b\) is a multiple of ten, and \(5\nmid a\). Write \(D=2a^2-b^2\). It is a nonzero even integer. For \(a/b<\theta\), \(D=-2\) is realized by \((7,10)\). For the other side, \(D\equiv2,14\pmod{16}\) and \(D\equiv2,3\pmod5\). The only positive possibilities below 62 are 2 and 18. If \(D=18\), reduction modulo three forces \(3\mid a,b\), violating primitivity.

To exclude \(D=2\), take the integral determinant-one transformation
\[
 T(a,b)=(3a+2b,4a+3b).
\]
It preserves \(D\). Every nonnegative solution for \(D=2\) with odd \(a\) is on the forward orbit of \((1,0)\). Indeed, while \(b\ge4a/3\), the inverse \((3a-2b,3b-4a)\) stays nonnegative and decreases \(a\) strictly. Positivity of the first entry follows from \(b/a<\sqrt2\). Below that threshold \(2=D>2a^2/9\), leaving only \((1,0)\) with odd \(a\). The projective orbit of \((1,0)\) modulo 17 gives slopes \(0,7,\infty,10\). None is an allowed nonzero square. This rules out 2 for **every** candidate, at arbitrary size.

By contrast \(T^{14}(7,6)=(293465818991,415023341310)\) is primitive, optimal at both spines and has \(D=62\). The action of \(T\) on the projective lines over \(2,5,17,29,257\) has orders \(1,3,4,5,32\). Hence \(T^{480}\) preserves all local directions and generates infinitely many optimal normals from each displayed witness. Since
\[
 b^2\lvert a/b-\theta\rvert=\frac{\lvert D\rvert}{2(a/b+\theta)},
\]
the defect bounds and the two attaining sequences prove both stated lower limits. The factor 31 is the ratio \(62/2\) of the permitted Pell defects.

There are 96 optimal indices in a full period of 480 on the \((1,2)\) orbit, which contains \((7,10)\). On each of the \((7,6)\) and \((9,10)\) orbits of defect 62 there are 30. These are counts of indices in projective orbits. The actual positive lattice vectors grow under every nonzero power of \(T\).

~~~sh
python3 -B pell_directions.py
~~~

## Absolute search at both irrational slopes

**FINITE-EXACT.** For either spine, any E61 solution satisfies both
\[
 |y-\sqrt2\,x|>2^{36},\qquad |y-x/\sqrt2|>2^{36}.
\]
The search in fact excluded every pair in the following four complete intervals of the smaller normalized address, regardless of its slope. Each interval includes more than 27 billion consecutive values, with no conditions on the two binary tapes.

| Spine | \(\alpha\) | Inclusive normalized interval | Literal height tests |
|---:|:---:|---:|---:|
|4|\(1/\sqrt2\)|128,012,277,430,877–128,039,765,221,572|2,406,697|
|4|\(\sqrt2\)|90,514,323,952,759–90,541,811,743,455|4,385,390|
|8|\(1/\sqrt2\)|16,097,761,995,573–16,125,249,786,269|56,330|
|8|\(\sqrt2\)|11,378,811,175,249–11,406,298,965,945|102,508|

Here is why a tube around an irrational line admits a search in consecutive integer rows. Let \(Y_0=(N-k)/q\). The chart bounds \(0\le f(z)\le z\) and the strict inequalities \(t+h,c+h,r+e<q\) give \(xy\le Y_0<(x+1)(y+1)\) at a solution. Set \(s=\lfloor\sqrt{Y_0/\alpha}\rfloor\), \(W=2^{36}\). If \(|y-\alpha x|\le W\), then \(s-2W-3\le x\le s+2W+3\). Above that interval both \(x>s+1\) and \(y>\alpha(s+1)\), contradicting the product upper bound. Below it, both \(x+1\le s\) and \(y+1\le\alpha s\), contradicting the strict lower bound on \((x+1)(y+1)\). This uses only \(\alpha\ge1/2\).

Both values of \(s\) can be found by two integer square roots since \(\alpha^2\) is rational. The [native integer slab inverse](../chanpak_farey_two_three_clocks_20261122/NATIVE-SLABS.md) uses the smaller of \(X,Y\) when \(x=10X+i,\ y=2Y+j\). In these tubes \(5(y-j)>x-i\), checked already at the two extremal row bounds with the smallest permitted \(y\). Thus \(X\) is the smaller normalized address. The replay sends the whole interval in the chart to the native inverse across all eight branch assignments. This also clears solutions outside the specified tubes when they lie on the same rows.

~~~sh
c++ -std=c++20 -O3 -Wall -Wextra -Werror ../chanpak_farey_two_three_clocks_20261122/farey-native-slab.cpp -o /tmp/chanpak-pell-slabs
for n in 4 8; do
  for a in 1 2; do python3 -B pell_slabs.py /tmp/chanpak-pell-slabs "$n" "$a"; done
done
~~~
