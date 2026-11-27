# One order-28 frequency and exact modulus 48 for six-path counts

2026-11-27. This strengthens the [unlabelled switching histogram
result](../pants_unlabelled_six_path_histograms_20261127/README.md).
All counts here are directed paths on six distinct vertices, with
no ordering factor for the other vertices. A cut is identified with
its complement.

## 1. The one-bin test at order 28

The frequency of **8,474,100** paths suffices to identify the
skew conference tournament of order 28 up to switching, relabelling
and converse. The previous result used this bin only to resolve a
tie. Here that bin alone has 41 different frequencies across the
65 signed-congruence representatives, hence separates all 41 classes
up to converse. Five other bins also work: 8,474,004, 8,473,956,
8,473,908, 8,473,860 and 8,473,764.

The replay uses [Spence's published representatives](https://www.rangevoting.org/SpenceMat.html)
in web-list order and the signed-equivalence census justified in the
previous note. It checks both the direct signed row-sum score and
the quartic formula below at every one of the 63,882,000 fourth-shell
locations across the representatives. The previous score separation
at order 28 ensures the fourth shell includes all cuts at these bins.
The 65 frequencies are in [replay.out](replay.out).

## 2. Every switching class has exact modulus 48

Let \(A\) be a skew conference tournament matrix of order \(n\ge8\),
so \(4\mid n\). For a cut \(\varepsilon\), let \(P(\varepsilon)\) be
its raw six-path count and write
\[
S=\sum_i((A\varepsilon)_i^2-1)^2.
\]
The [fourth-moment identity and extrema](../pants_extremal_cuts_golay_20261127/README.md)
give \(P=\mathrm{const}-S/32\). Fix any root \(r\); normalize it
as a source and hold its cut bit fixed to give one representative
per cut. The polynomial in the other \(0/1\) flip variables has
degree at most four. Its coefficients in positive degrees are
the corresponding forward differences at the unflipped cut:

| Degree | Coefficient for the indicated flip set |
|---|---|
| 1 | \(-8n(n-2)(n-4)\) |
| 2 | \(48n(n-4)\) |
| 3, triple \(B\) | \(-192(n-1+K_B+Q_{rB})\) |
| 4, four-set \(U\) | \(384Q_U\) |

The first two follow from \(s_0,s_1,s_2\), the first three
shell scores. The third is the rooted shell formula, with
\(K_B=-3\) on cyclic triples and \(+1\) otherwise. The fourth
follows from the Walsh coefficient \(24Q_U\) of \(S\) by
replacing \(\varepsilon_i\) with \(1-2x_i\). The
[quartic residue rules](../pants_skew_conference_quartic_detection_20261127/README.md)
give \(Q_{rB}\equiv n-4\pmod8\) in the cyclic case and
\(Q_{rB}\equiv n\pmod8\) otherwise. Each \(Q_U\) is divisible
by four. Consequently every listed coefficient is a multiple
of 1536: in degree 1 use \(32\mid n(n-2)(n-4)\) and a factor
of three from one of its three factors. Hence all \(P\) are
congruent modulo 48.

**Four cuts about any one minimum edge certify sharpness.**
Let \(P_+(B),P_-(B)\) be the counts obtained by flipping a triple
\(B\) from the two ends of that edge. The third-shell formula
has the same baseline at both endpoints, and
\(P_\pm(B)=P_*+6(n-1+K_B\pm Q_{rB})\), for \(P_*\) independent
of \(B\). For a transitive triple \(T\) and a cyclic triple \(C\),
\[
 P_+(T)+P_-(T)-P_+(C)-P_-(C)=48.
\]
Both types exist in the normalized core, a regular tournament
on at least seven vertices: the absence of a directed triangle
would imply transitivity, and an outneighborhood of size at
least two supplies a transitive triple. Thus the gcd of
differences is already exactly 48 on the local distance-three
neighborhood of **each** minimum edge, whether or not those
neighborhoods overlap elsewhere. This is a signed sum of four
counts, not an assertion that every progression value occurs.
The same endpoint sum identifies whether a triple is cyclic,
and the endpoint difference gives its signed quartic; the two
readings separate the structural bit from the quartic correction.

The extremal formula at \(n=4m\) gives
\[
 P_{\min}=2m(m-1)(4m-1)(2m-1)(8m^2-16m+7).
\]
It vanishes modulo three. Modulo 16 the polynomial depends
only on \(m\bmod8\), and direct substitution gives
\(0,0,12,12,8,8,4,4\). Hence the precise residue law is
\[
 P(\varepsilon)\equiv12\lfloor n/8\rfloor\pmod{48}.
\]
Multiply both the residue and the modulus by \((n-6)!\)
for the ordered path-cover convention.

## 3. The local polynomial and separated shells

Fix an edge with direction \(r\), in its source-normalized gauge.
A flip set \(B\) off \(r\) has local distance \(k=|B|\). Using the
same rooted \(Q\) values at the two endpoints, the Newton expansion is
\[
S_\pm(B)=s_0+k(s_1-s_0)+\binom{k}{2}(s_2-2s_1+s_0)
-192\sum_{\substack{T\subseteq B\\|T|=3}}(n-1+K_T\pm Q_{rT})
+384\sum_{\substack{U\subseteq B\\|U|=4}}Q_U.
\]
Starting from the sink toggles exactly the terms with root \(r\)
in their Walsh support. At \(k=4\) the first three terms give
\(3s_0-8s_1+6s_2\).

Write \(L=n(n-2)\). In the \(k\)-th shell the largest row
magnitude is \(n-1-2k\), so \(S\ge[(n-2k)(n-2k-2)]^2\).
For any farther cut, the largest row magnitude \(M\) satisfies
\(M\le n-2k-3\). As in the previous note, for \(M\le n/2\) use
\(S\le(n^2/4-1)L\); otherwise use the convex cubic
\[
 g(M)=(M^2-1)^2+((n-M)^2-1)(L-M^2+1).
\]
At \(k=4\) the previous proof handles \(n=28\). At \(n\ge32\),
putting \(x=n-32\), the fourth-shell lower bound exceeds
\(g(n/2)\) by \(3x^4/4+121x^3/2+1685x^2+17406x+33984\)
and exceeds \(g(n-11)\) by \(8x^3+528x^2+9216x+22784\).
Thus the fourth shell is separated at every \(n\ge28\).

More generally put \(z=n-2k-2\). The exact difference of the
lower shell bound and the far-end bound is
\[
 [z(z+2)]^2-g(z-1)
 =8z^3-16(k+1)^2(k+2)(n-k-2).
\]
In particular all shells through \(k\) are separated if
\(n\ge\max(8(k+1),3(k+1)\sqrt{k+2})\). Under this condition
\(z\ge3n/4\), so the lower shell bound beats \(g(n/2)\) and
\(z^3\ge27n^3/64>2(k+1)^2(k+2)(n-k-2)\). The condition
also holds for every smaller radius. Asymptotically, for any
\(c<2^{-1/3}\), taking \(k=\lfloor c n^{2/3}\rfloor\) makes
the displayed difference \(8(1-2c^3+o(1))n^3>0\) while the
midpoint bound is again smaller. Both comparisons hold for
every smaller radius too: the shell bound increases as the
radius decreases, and \((k+1)^2(k+2)(n-k-2)\) increases with
\(k\) throughout \(0\le k=o(n)\). These radii are below half
the minimum edge separation, so each counted shell cut has a
unique minimum edge. The constant is sufficient, not a
claim of the largest possible radius.

## Replay

    curl -L https://www.rangevoting.org/SpenceMat.html -o SpenceMat.html
    g++ -O3 -std=c++20 replay.cpp -o replay
    ./replay SpenceMat.html

On x86 with POPCNT, adding -mpopcnt speeds up the direct check.
