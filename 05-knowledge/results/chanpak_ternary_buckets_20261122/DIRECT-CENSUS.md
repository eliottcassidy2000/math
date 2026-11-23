# All direct tripling forks on cores of at most eight vertices

**FINITE-EXACT, with an unbounded prefix inverse.** Use the right-flank, one-apex tournament of the [moving-fork search](PHYSICAL-FORKS.md). Choose any ordinary core of order 1 through 8, its apex mask, one mode with positive slope factor, and **any** positive odd integer \(w\). A presentation qualifies if three fixed individual arc reversals at the center give heights \(H-2,H+2,H/3\) for every \(s\ge1\). Each reversed arc lies either inside the ordinary core or between a core vertex and the apex. The center and source have the same vertex set and the same right address; the inverse move needs no padding or relabelling.

Under these quantifiers the only center curves are
\[
                 H=63\cdot 2^s-39,\qquad H=126\cdot 2^s-39 .
\]
Thus there is one source ray \(21\cdot2^s-13\), together with its shift by one place. Exactly two of the 62 distinct center curves found with \(w\le31\) survive the direct-source requirement. Allowing arbitrary odd prefixes adds no curve. All successful presentations have \(w=1\), right address \(2^s-1\), and a core-to-core tripling edge. At \(s=1\) the first fork is \(29\to87\to85,89\).

Here are the three presentations among the supplied ordinary isomorphism representatives. A row integer encodes its outgoing neighbors; the vertex labels start at zero. Different presentations may be isomorphic as complete unmarked tournaments.

| Center coefficient of \(2^s\) | Ordinary adjacency rows | Mask | Arc to \(H-2\) | Arc to \(H+2\) | Arc to \(H/3\) |
|---:|---|---:|---|---|---|
|63|0, 1, 11, 19, 71, 31, 47|11|(2,3)|(4,5)|(4,6)|
|126|0, 1, 3, 23, 39, 143, 63, 95|22|(3,4)|(5,6)|(5,7)|
|63|0, 1, 11, 19, 71, 31, 47, 127|11|(2,3)|(4,5)|(4,6)|

## How an integer intercept replaces the unbounded search

For an ordinary core of height \(h\), write its full right-flank height as \(K+cy+e[y\text{ odd}]+hf(y)\), with \(f\) the [binary read](../chanpak_decimal_kernel_20261122/REVERSED-PHASE.md). The two modes give \(H=a w 2^s+b-hd\) with the following data.

| Address \(y\) | Slope factor \(a\) | Fixed intercept \(b\) | Variable correction \(d\) |
|---|---:|---:|---:|
|\(w2^s\)|\(c+h\)|\(K\)|\(w-1-f(w)\)|
|\(w2^s-1\)|\(c\)|\(K+e-c\)|\(-f(w-1)\)|

Take \((a',b',h')\) from the center after a candidate fixed reversal. A move to offset \(j\in\{-2,2\}\) has to satisfy \(a'=a\) and \((h'-h)d=b'-b-j\). A move to a third of the height has to satisfy \(3a'=a\) and \((3h'-h)d=3b'-b\). These conditions are equivalent to the identities for all \(s\), by equality of the two affine coefficients. Intersect the solutions from one arc of each type. Each equation in \(d\) has all integers, one integer, or no integer as its solution. A zero slope factor can occur only in the second mode with the apex beating no core vertices. Then the center has height \(h(1+f(w-1))\), and the original 62-curve catalogue explicitly restricted to positive slopes. **All conclusions here concern growing centers \(a>0\).**

For every fixed \(d\) the remaining word inverse is finite. Put \(w=2v+1\). In the first mode \(d=2v-f(v)\ge v\), since \(0\le f(v)\le v\). In the second mode \(-d=f(v)+v+[v\text{ odd}]\ge v\). Hence \(w\le2|d|+1\) in either case, and \(d\) must be even. Enumerate this explicit interval and check the exact read. The all-integers case would accept all odd prefixes; no such triple arose. The scan covers all tag assignments on the [ordinary representatives](../../../session212/cycle-levels/MANIFEST.json) and each of their \(\binom{n}{2}+n\) eligible arcs. It checks the arcs on the actual labelled center, so two moves on different center presentations cannot combine to form a fork. This proof of completeness allows the three choices of arc to depend on the core, mask, mode and prefix.

The tagged presentation counts before screening at orders one through eight are \(2,4,16,64,384,3584,58368,1761280\). The replay below prints every successful presentation. The [direct-fork replay](direct_fork.py) gives a complementary vertex-row check.

~~~sh
c++ -O3 -std=c++17 -Wall -Wextra -Werror direct_census.cpp -o /tmp/chanpak-direct-census
/tmp/chanpak-direct-census ../../../session212/cycle-levels 8
~~~

## The boundary of the apex obstruction

A single apex arc *can* triple a moving height. On the two-vertex ordinary core \((0,1)\), use the even address \(y=2^s\) and switch the mask between 0 and 2. With mask 0 the core has \(h=1,c=0,K=1\) and height \(2^s+1\). With mask 2 it has \(h=1,c=2,K=3\) and height \(3(2^s+1)\). The obstruction in the census is the joint demand for the other two physical exits from that same center. The census leaves open moves involving vertices of the growing flank, larger ordinary cores, and coincidences at a single specified value of \(s\).
