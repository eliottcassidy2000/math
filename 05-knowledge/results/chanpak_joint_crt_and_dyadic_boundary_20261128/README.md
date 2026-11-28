# Moving-root descent: exact CRT gluing and the integer side of the boundary

2026-11-28. This develops the [moving-root Chan–Pak grammar](../chanpak_moving_root_and_module_interfaces_20261128) with two ideas from the repo's Hadwiger–Nelson and planar Jacobian work. It proves full reachability **simultaneously modulo every odd modulus**, not just prime-by-prime. It also gives a sharp positive-integer inverse for fixed-root stretches and an infinite family of hostile ports that pass every odd residue test but have no positive predecessor. These are statements about the specified tournament grammar, not proofs of Chan–Pak, Hadwiger–Nelson, or the planar Jacobian conjecture.

## 1. Three charts, one common wall

The state (h,e,u,g) records the Hamiltonian-path count, the number of paths ending at the root, the ordered split count used in the preceding note, and the root-deletion count. From the singleton (1,1,1,1), the actual tournament operations are
\[
P:(h,0,2u,g),\quad
K:(h-e+u,g,2u-g,g),\quad
M:(h-e+u,h-e,h+u,h).
\]
Here P appends a sink, K reverses just the root contact and keeps the root, and M reverses it and moves the root to the new vertex. Use \(\mathcal P:e=0,\ \mathcal K:e=g,\ \mathcal M:u=h+g-e\) for the three destination charts. In all charts take h,u,g as coordinates. Direct differentiation gives

| Source | Append P | Append K | Append M |
|---|---:|---:|---:|
| P chart | 2 | 2 | 0 |
| K chart | 2 | 2 | 1 |
| M chart | 2 | 0 | 1 |

The two singular transitions land on the *same* wall \(\mathcal K\cap\mathcal M=\{h=u,e=g\}\); both have rank two with that entire rational image. They are P followed by M and M followed by K, in chronological notation. The other restrictions have determinant 1 or 2. In particular, prime 2 is the only possible lattice-index obstruction for a nonsingular transition. This calculation does not say those transitions map the positive state region onto itself.

The analogy with the [planar delta collision-wall calculation](../../../01-canon/theorems/THM-4141-delta-d-collision-wall-boundary-monodromy-exclusion.md) and the [two-term wall](../../../01-canon/theorems/THM-4143-two-term-collision-wall-critical-boundary-monodromy-exclusion.md) is method, not implication: normalize and label the live charts, compute the collapsed boundary separately, then test whether the formal inverse lands on an allowable branch. Here everything is linear, and on a fixed-g slice of \(\mathcal K\) the active maps below are genuinely planar affine maps. No Keller-map or monodromy assertion is transferred.

## 2. Exact simultaneous odd-modulus theorem

**Theorem.** For every odd integer m≥3, the states reached by finite words in P,K,M modulo m are exactly \(\mathcal P\cup\mathcal K\cup\mathcal M\) over \(\mathbb Z/m\mathbb Z\). There is a constructive witness of length \(O(m^3)\). The size of the union is \(3m^3-3m^2+m\).

Inclusion in the union follows from the destination charts; the seed is already on \(\mathcal K\). It suffices to fill \(\mathcal K\), since P and M restricted to that chart have determinants 2 and 1. Write its coordinates as \(x=u-g,\ y=h-u,\ g\). The following calculations work over any odd modulus. On \(\mathcal K\), K acts by \((x,y)\mapsto(2x,y)\). As in the previous prime calculation, follow PK with \(K^{-2}\), or PPK with \(K^{-3}\), to get
\[
E:(x,y)\mapsto(x+g/2,y-x),\qquad
F:(x,y)\mapsto(x+3g/4,y-3x-2g).
\]
Negative powers are shorthand for positive iterations in the finite ring. This shorthand is never applied to the positive integer inverse.

### Vertical corrections, including powers of 3

The commutator of E and \(KEK^{-1}\) is a pure vertical translation by \(\pm3g/4\); call it C. Let J be F followed by \(E^{-r}\), with \(r=3/2\bmod m\). The x-shifts cancel exactly, so J is a vertical translation by a value dependent on x and g. Modulo 3 that value is g: F's y-increment is g there and \(E^r\) is the identity there. Thus on a component \(q=p^a\) where g is a unit, C corrects y when \(p\ne3\), and J corrects y when \(p=3\). Their translations have period dividing m. CRT on their exponents turns either correction on in the desired prime-power factors and off in all other factors. If g is a unit on the whole ring, E first sets x arbitrarily, and these two corrections then set y arbitrarily.

From the singleton's g=1 slice we can reach all such unit slices: MK applied to a K-state with h=t lands at \(g'=t,x'=2u,y'=0\). One entrance is enough when t is a unit.

### A single supplying slice handles the nonunit factors together

Take a target (x,y,g) with \(\gcd(x,y,g,m)=1\). Start on a supplying K slice \(g_{\rm old}=t\), where t will be a unit. Set h=g and choose its u freely. Appending MPK gives the target-g slice with initial coordinates
\[
x_0=2g+4u,\qquad y_0=-x_0/4+g/2-t.
\]
Now apply \(E^j\). To arrive at x, take \(x_0=x-jg/2\). Then the resulting y is
\[
-x/4+g/2-t-jx+gj/8+g j(j+1)/4 .
\]
Choose t and j coherently by CRT:

* At a prime power where g is a unit, take t=1 and j=0 on that component and correct y afterward by C or J.
* Where g is a nonunit but x is a unit, take t=1. The displayed function of j is a permutation modulo that prime power: the difference at j and k is (j−k) times a factor congruent to −x modulo p. Choose its unique preimage of the desired y.
* Where x and g are both nonunits, y must be a unit. Take \(t=-y-x/4+g/2\) there and \(j=0\). This t is a unit on that component, and the output y is exact.

These choices give **one** unit t modulo m, **one** initial x_0, and **one** exponent j. The g=t slice is fully reachable, so select the required u there and execute the same word for all components. Correct y only on the unit-g factors. No separate incompatible choice of last edit or motion is made at different primes.

For a nonprimitive target put \(\delta=\gcd(x,y,g,m)\). We can reach the scaled seed \(\delta(1,1,1,1)\bmod m\): on g=1 choose h=δ, u=0 and apply MK. Use the word for the divided primitive target modulo \(m/\delta\) from that scaled seed. Homogeneity of the three integer matrices transports it to the desired target modulo m. This also handles the zero target.

Finally, \(\operatorname{ord}_m(2)\le m\) makes each use of E and F cost \(O(m)\) actual edits. Forming one C or J costs \(O(m^2)\), and the CRT powers of C and J can each be selected once modulo m, for \(O(m^3)\) total edits. Only a fixed number of unit-slice and scaling calls is needed. Each pair of charts intersects in \(m^2\) points and the triple intersection has m points, giving the stated size. The replay uses compressed words of actual positive moves; it does not claim this is a practical way to realize a specified integer h.

## 3. The dyadic boundary does select a positive-integer branch

Stay at the same root from one K-state to the next, using \(P^rK\) for r≥0. Both u and g at a K-state are odd for a genuine word (K gives \(u'=2u-g\), and deletion counts are odd). The literal update is
\[
u'=2^{r+1}u-g,\qquad
h'=\begin{cases}h+u-g,&r=0,\\h+2^ru,&r>0.\end{cases}
\]
Consequently the exact integer \(u'+g\) determines the one possible value \(r=v_2(u'+g)-1\), the predecessor \(u=(u'+g)/2^{r+1}\), and the predecessor h. In the nontrivial region the gap \(u-g\) at least doubles with each block, so an exact inverse down a fixed-root stretch uses only logarithmically many nontrivial blocks. The only zero-growth K block is r=0,u=g; it fixes the state and can be deleted from the word. This peels runs of actual moves, not rational choices of inverse sheet.

There is a stronger one-step exclusion off the singular wall. Suppose a proposed positive state is on \(\mathcal K\), has \(h\ne u\), and \((u+g)/2\) is odd. Its last letter must be K. A precursor on the M chart would land on the wall and is excluded. A literal P-precursor has even u and is excluded. Thus the K-precursor is forced:
\[
(h_\mathrm{old},e_\mathrm{old},u_\mathrm{old},g_\mathrm{old})
=\left(h-\frac{u+g}{2}+g,\ g,\ \frac{u+g}{2},\ g\right).
\]
For example, for every positive odd h, the proposed port \((h,1,2h+3,1)\) satisfies the simple guards \(h\ge g\ge e\), \(u\ge h\) and all odd congruence tests from the theorem, yet it is not reached from the singleton. Its unique allowable predecessor would have height −1. This is an obstruction to a **port state**, not to the scalar h. It explains why proving unrestricted simultaneous odd-modulus coverage is not enough for Chan–Pak.

## 4. Hadwiger–Nelson: keep the support test, not just the inverse data

The repo's [five-interior extension theorem](../paired_budget_selectors_matrix_20261028/FIVE-INTERIORS.md) first forms the exact boundary quotient and then checks whether its formal critical graph has a *literal geometric realization*: the \(K_4\vee C_5\) quotient fails on allowed contacts, not because the abstract quotient is uncolorable at some unrelated parameter. Here the three chart labels and their common rank-drop wall are the boundary quotient; the parity and positivity of the backward register state are the literal support test. The family above has genuine modular lifts at every odd modulus but fails that support test at a unique forced branch.

The [unit-distance two-port ear work](../two_port_unit_distance_20261103) makes a second useful distinction: full reconstruction of boundary weights from positive ear probes does not create a forcing gadget. In the Chan–Pak chart, the positive CRT words reconstruct every odd residue vector, but a specified absolute target also needs compatible positive integer precursors. Conversely, the two-dimensional Jacobian on each live slice and the explicit CRT calculation explain **which** compatibility has already been settled: odd-prime gluing, not positive descent.

Proposed next inverse test: choose two different incoming faces with the same scalar h, retain their exact dyadic signatures, and seek a bound that guarantees a positive precursor for at least one of them on all large targets. Terminal h alone, and the three face labels with odd residues alone, both omit constraints proven necessary here.

## Replay

Run python3 replay.py. It constructs compressed positive P/K/M words for all 881 K-chart points across moduli 3, 5 and 9, and for 100 seeded points at each of fifteen further odd moduli, including 3-power and mixed composite cases. It checks their literal modular matrix images, the append-to-P and append-to-M images, the dyadic block inverse, and the hostile negative predecessors. The proof is uniform; these checks are not its basis.
