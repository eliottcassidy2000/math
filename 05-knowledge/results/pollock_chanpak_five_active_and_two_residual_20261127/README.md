# Five active indices and the two-residual Pollock shell inverse

This extends the [four-active exact bridge](../pollock_chanpak_first_non_cubic_overlaps_20261127/README.md) and the [shared-three shell barrier](../pollock_chanpak_small_edit_rigidity_20261127/README.md). Write \(F(x)=(x^3-x)/6\) with \(x\ge2\), so each summand is strictly positive.

## 1. The first fully active short-guard matches

Allow all five figurate indices to change, while preserving their sum and sum of squares. Require the two index multisets to be disjoint. The table gives all matches with the 1010-leading short-guard chart through lower height \(500001\); the raw count is \(H=2N+1\).

| Lower indices | Upper indices | \(N\to N+1\) | Payload \(v\) |
|---|---|---:|---:|
| 11,29,45,51,83 | 18,18,49,52,82 | 136844 → 136845 | 167 |
| 15,23,41,59,81 | 18,18,49,52,82 | 136844 → 136845 | 167 |
| 2,36,47,56,80 | 3,31,53,55,79 | 139647 → 139648 | 168 |

The first two rows have common moments \((219,12477)\). Their tournament move changes 273689 to 273691 on 18 vertices. The third changes 279295 to 279297, also on 18 vertices, with common moments \((221,13045)\). A transitive prefix with apex bits \(2v\) in binary is followed by the six-vertex fork, one guard and a final vertex. In the core every cross-block arc goes forward; the apex beats fork slots 0, 1, 4 and the guard, and loses to the final vertex. Fork masks 177 and 181 change only arc (0,3). The first bit is 1, so the apex reaches the whole core through the first vertex; every core vertex can reach the apex through the final vertex. These are literal strong tournament moves.

For completeness, the established height formula \(1662v-58q(v)+427\), bounded below by \(1604v+427\), gives 31 eligible lower heights. The replay solves the figurate inverse at each lower height and its successor by sorting the five indices, storing all two-sums, and searching the remaining three indices with monotone bounds. After matching the moments, it checks disjointness. It also independently counts paths in the actual tournaments by splitting at the apex. The smallest fully active figurate step without the guard constraint is already \(124\to125\): \((3,3,3,7,7)\to(2,4,4,5,8)\).

At the larger height 546598033 there are exactly two payloads, 329528 and 333037, both of order 29. They give the same exact endpoints 546598033 and 546598035. Their outdegree multisets differ at each endpoint even after allowing converse. On the Pollock side the unique equal-moment pair among all 59326 lower and 27104 upper positive five-term representations is the disjoint pair
\[
(141,485,496,879,897)\to(194,308,652,872,872).
\]
The moments are \((2898,2078372)\). Thus both a two-to-one and a one-to-two inverse occur in the same comparison, at different heights.

**Parity filter.** For a unit step let \(W(Z)\) be upper minus lower index generating polynomial. As in the four-active proof, \(W(-1)\equiv8\pmod{16}\); now \(|W(-1)|\le10\), so the number of even indices changes by exactly four in absolute value. One list is homogeneous and the other has four indices of the opposite parity and one exceptional index. If \(S\) is the common first moment, the exceptional index is \(S+2\pmod4\). To see this, use \(F(j)\) odd iff \(j\equiv2\pmod4\). If the homogeneous list is even, its total mod two is \(S/2\pmod2\), whereas only the exceptional even index contributes on the mixed side. If the homogeneous list is odd, the mixed list's four even indices must sum to two mod four. Both cases give the filter. It is necessary, not sufficient.

## 2. Exact divisor test for two changing residual indices

Write \(N=F(k+1)+r\), \(0\le r<\binom{k+2}{2}\). Lower the anchor to \(F(k)\), retain two residual terms, and replace a positive pair \((x,y)\) by \((x+a,y-a)\) with the same index sum. Since \(F\) is strictly convex on the positive indices, sorting the two pairs lets us assume \(x\ge y\) and \(a>0\): to increase the sum of their tetrahedra at the same index sum, the larger index must increase. If \(S=x+y\) and \(D=x-y+a\), then
\[
F(x+a)+F(y-a)-F(x)-F(y)=aSD/2.
\]
In particular the following is an **exact inverse** for this constrained move. Enumerate factorizations \(aSD=k(k+1)\) such that
\[
1\le a\le D,\qquad S\ge D+a+4,\qquad S+D-a\equiv0\pmod2 .
\]
Recover \(x=(S+D-a)/2\), \(y=(S-D+a)/2\), and test \(r-F(x)-F(y)\) for a sum of two strictly positive tetrahedral numbers. Each factorization passing these tests supplies both representations; every such representation appears. This needs no search over tournament or residue states.

**Prime-factor obstruction.** Every index in a residual representation is less than \(2(k+1)^{2/3}\): its tetrahedral value is below \((k+1)^2\) and \(F(j)\ge j^3/8\) for \(j\ge2\). Hence \(a,D<S<4(k+1)^{2/3}\). A prime divisor of \(k(k+1)\) at least \(4(k+1)^{2/3}\) excludes the conserved-sum move for the entire shell. In particular all sufficiently large prime \(k\), or prime \(k+1\), are excluded. This is not an obstruction to changing the residual index sum, and is not a Pollock nonrepresentation criterion.

## 3. An infinite family of positive shell repairs

For every integer \(t\ge11\) take the following factorizations:

| \(t\) | \(k\) | \(S\) | \(a\) | \(D\) |
|---|---|---|---|---|
| even | \(t^3\) | \(2(t^2-t+1)\) | \(t^2/2\) | \(t(t+1)\) |
| \(1\bmod4\) | \(t^3\) | \(2t^2\) | \(t(t+1)/2\) | \(t^2-t+1\) |
| \(3\bmod4\) | \(t^3-1\) | \(2t^2\) | \(t(t-1)/2\) | \(t^2+t+1\) |

Factor \(t^3\pm1\) to get \(aSD=k(k+1)\). Each row satisfies the parity condition, \(D\ge a\), and \(S-D-a\ge4\). Also \(S\le2t^2\) and \(0\le D-a=x-y\le2t^2/3\) (substitute the three quadratic expressions and use \(t\ge11\)). The lower pair costs
\[
F(x)+F(y)=\frac{S^3+3S(x-y)^2-4S}{24}\le4t^6/9.
\]
Let the two retained indices \(u,v\) range independently from 2 through \(\lfloor t^2/2\rfloor\). They cost at most \(t^6/24\), giving \(r\le35t^6/72<t^6/2\le\binom{k+2}{2}\). So both \((k+1,x,y,u,v)\) and \((k,x+a,y-a,u,v)\) are positive representations of the same \(N\) in the stated shell.

For \(R=F(u)+F(v)\), \(u+v\) divides \(6R\), and fixing \(u+v\) fixes \(uv\), hence fixes the unordered pair. The multiplicity of each retained sum is at most \(\tau(6R)\). Out of \(\Theta(t^4)\) unordered pairs, this gives \(\gg_\varepsilon k^{4/3-\varepsilon}\) **different exact remainders** in every selected shell. For \(t=11\), \(k=1330\) and \((160,82)\to(215,27)\) is the moved pair; any two retained indices from 2 through 60 work. This proves an infinite positive construction, not a full-shell coverage result.

Run the two JavaScript replays with Node. The large two-plus-three search can take a minute. In combination the positive family and the prime-factor obstruction explain why keeping the index sum is productive for finding examples but too restrictive for a uniform shell proof.
