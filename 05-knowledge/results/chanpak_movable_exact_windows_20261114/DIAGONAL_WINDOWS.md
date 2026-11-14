# Exact bilateral windows by balancing two source coordinates

This is a physical tournament construction of several exact moving windows, not
an absolute-target theorem. At a fixed slope and dyadic scale the center still
moves quadratically. Relative/locally residue-complete control alone does not
place a window over an arbitrary odd integer.

## Bilateral composition

Let (B) have at least one vertex and mark a subset (U) beaten by a new apex
(v); all other vertices of (B) beat (v). Prepend a transitive wing (L)
and append a transitive wing (R) to (B); old cross-arcs go forward. Label
each wing in its transitive order by whether (v) beats that vertex. On (L)
use the ordinary MSB-first binary expansion of (x\ge1). On (R) use the
reversed complemented binary expansion of (y\ge1).

Write (h=H(B)), and (f=\sum_{D\subseteq B}H(D)H(B\setminus D)), with
(H(\varnothing)=1). A complementary pair means a directed Hamiltonian path
on each of two complementary subsets, with order of the subsets retained.
(P) counts complementary pairs whose second path starts in (U), and (Q)
counts those whose first path ends outside (U). Write (b) for the number
of whole-core paths starting in (U), and (e) for the number ending outside
(U). Define (S=h-e,T=h-b). The core puncture count (K) is
(b+e+I), where (I) counts pairs with both paths nonempty that can be joined
via (v). Finally (E(U)) counts pairs with first end and second start both
in (U). Then (P=b+I+E(U)), hence (K-P+S=h-E(U)).

Let (p(x)=x\bmod2) and let (g(x)) be the one-apex count on just the left
wing. It satisfies (g(0)=1), (g(2z+1)=g(z)), and
(g(2z)=g(z)+z+p(z)) for (z\ge1). Appending a terminal one forces that
new vertex last, leaving the count unchanged. Appending a terminal zero
allows (z+1) ways to place it immediately before the apex and an end-vertex
copy of all old paths except the (1-p(z)) paths ending in the apex.
The number of nonempty subsets of the left wing whose first vertex is beaten
by (v) is its binary value (x); the right-wing dual gives (y).

The full count is
$$
H=fxy+Qx+Py+h(g(x)+g(y)-2)+K+Tp(x)+Sp(y).
$$
To see this, delete (v) from a path, leaving prefix and suffix. If both
core pieces are nonempty, each complementary pair contributes
((y+a)(x+d)), with (a) the indicator that its first path ends outside
(U), and (d) the indicator that its second path starts in (U).
If the entire core is in the suffix, the contribution is
(hxy+by+h(g(x)-1+p(x))+b(1-p(x))): the last two terms cover absent right
prefix with respectively present and absent left suffix. If the entire
core is in the prefix, the reflected contribution is
(hxy+ex+h(g(y)-1+p(y))+e(1-p(y))). Summing yields the formula.
All these full tournaments are strong: the first ordinary vertex is beaten
by (v), the last ordinary vertex beats (v), and every core vertex is
between the wings.

## Tilt and offset

Take (x=MLu,y=MRu-1), where (M\ge2) is a power of two and (L,R,u\ge1).
Put (W=RP+LQ). The recurrences give (g(x)=g(Lu)+(M-1)Lu+p(Lu)) and
(g(y)=g(Ru-1)). Thus all cores with common ((h,f,W)) have common
height
$$
G=fM^2LRu^2+Mu(W-fL)
+h[(M-1)Lu+g(Lu)+g(Ru-1)+p(Lu)-2]
$$
and offsets (h-E(U)). Finding a same-height window means finding one
linear level (W) on which the quadratic end/start count (E) takes
consecutive step-two values.

* (L=R=1): ((h,f,W)=(5,64,35)) and also ((5,64,83)), each with offsets
  (-7,-5,-3,-1). At (M=2,u=1) these give respectively 201–207 and
  297–303, odd values only.
* (L=1,R=2): ((h,f,W)=(9,672,981)) gives offsets
  (-59,-57,-55,-53,-51). Here (G=1344M^2u^2+309Mu+
  9[(M-1)u+g(u)+g(u-1)+p(u)-2]). Five consecutive odd values fill every
  odd residue modulo ten with span eight. An order-7 core bank with
  (L=1,R=4), ((h,f,W)=(9,336,645)), produces the identical heights and
  offsets. Increasing the core by absorbing a leading right-wing zero bit
  explains the redundancy: this vertex is a new core sink outside (U).

A separate three-odd construction uses only two transitive wings and no
ordinary core. There the count is (xy+g(x)+g(y)-1+p(x)p(y)), obtained by
partitioning paths according to whether the prefix has a right vertex and
the suffix has a left vertex. If both events occur, the contribution is
(xy); if exactly one holds, the respective contributions are
(g(x)-1+p(x)) and (g(y)-1+p(y)); if neither holds the contribution
is ((1-p(x))(1-p(y))). At (A=8u), choose (x=A+a,y=A-1-a) for
(a=0,1,2). The product decreases by (a(a+1)); the sum of the two (g)'s
changes by (0,0,2), by the binary recurrences. Exact outputs are therefore
(64u^2-u+g(u)+g(u-1)+p(u)-1-2a).

## Reproduction

```sh
python replay_diagonal.py
g++ -O3 scan_diagonal.cpp -o /tmp/anchor-scan
g++ -O3 scan_tilt.cpp -o /tmp/anchor-tilt
python canonical_cores.py > /tmp/anchor-cores
/tmp/anchor-scan < /tmp/anchor-cores
/tmp/anchor-tilt 1 2 < /tmp/anchor-cores
/tmp/anchor-tilt 1 4 < /tmp/anchor-cores
```

The scanners generate every order-8 tournament at least once from the order-7
canonical bank. The diagonal scan has 58,900 inputs including duplicates,
714,132 profiles, and maximum contiguous window four; the slope-two scan
has 1,115,310 profiles and maximum five. The profile calculation uses subset
Hamiltonian DP; the separate Python replay builds full literal tournaments
and checks their Hamiltonian counts directly at selected small (M,u).
No census of target magnitudes is used as evidence for eventual coverage.

## Moving a phase without changing the quadratic baseline

At slope two let (X=Mu), and use (x=X+q,y=2(X-q)-1) with integer
(q) and (M>2(|q|+1)). Because (g(2t+1)=g(t)), the wing contribution
is (g(X+q)+g(X-q-1)). For (qge0) the length-(\log_2 M) low words of
these two inputs are complementary; the high words are (u,u-1).
At each low-bit step exactly one zero triggers a contribution with
linear (u)-coefficient (M/2^{i+1}). The high-word parity enters only
at the final step, and enters on the (u) side just as in the baseline
because the high bit of (q) is zero. Therefore the difference
(\delta_q=g(X+q)+g(X-q-1)-g(X)-g(X-1)) does not depend on (u).
The other signed phases follow by (q\mapsto -q-1), which exchanges the
two arguments. The full difference from the same core at (q=0) is
$$
-fq(2q+1)+q(Q-2P)+h\delta_q+T p(q).
$$
In particular (\delta_{-1}=\delta_1=0). These shifts keep the same
absolute baseline while retaining a different finite offset selector.
An exact order-through-eight scan including phases (-1,0,1) finds no
longer window than five at slope two. This rules out only the fixed
profile and these phases at those orders. Reproduce with
`g++ -O3 scan_tilt_phase.cpp -o /tmp/anchor-phase; /tmp/anchor-phase 1 2 < /tmp/anchor-cores`.
