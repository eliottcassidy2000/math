# A diagonal norm carrier and the two different mod-8 gates

**Proved.** Write (h_A=H(A)), (q_A=P_2(A)+2h_A) for two ordered optional path slots, $x_A=q_A-h_A$ and $y_A=h_A$. For nonempty tournaments (A,B), orient the three blocks $A\to B\to c\to A$, with (c) a singleton. Then
$$
F(A,B)=x_Ax_B+2y_Ay_B.
$$
Indeed, deleting (c) breaks each path into its ordered before- and after-runs. Each run is a path in (A), in (B), or in $A\Rightarrow B$ using a portion of each. The allowed block words, deleting (c) from the display only after the type has been recorded, are (cAB,ABc,BcA,ABcA,BcAB,ABcAB). The resulting sum is (3h_Ah_B+P_2(A)h_B+h_AP_2(B)+P_2(A)P_2(B)). Every resulting tournament is strong (the directed triangle of modules provides paths between and within modules). Converse leaves both coordinates invariant. Thus diagonal and converse-diagonal joins give (D(A)=x_A^2+2y_A^2).

Every nonempty induced tournament has odd (H); the empty one has (H=1). Pairing complementary subsets in $q_A=\sum_{S\subseteq V}H(A[S])H(A[V\setminus S])$ gives $q_A\equiv2^{|A|}\pmod{4}$, even for a singleton. Also (P_2(A)\ge (|A|-1)H(A)) by cutting each Hamiltonian path at each internal gap: concatenation recovers the path and the lengths recover the gap. Hence $x_A\ge |A|y_A>0$. For order two, the exact pair is ((x,y)=(3,1)), so $x\ge3y$ at every order at least two.

## The prime restriction is stronger than a residue class

Actual diagonal outputs are always **3 mod 8**, not merely 1 or 3: both $x_A$ and $y_A$ are odd. General positive norms in $\mathbb Z[\sqrt{-2}]$ that are odd are 1 or 3 modulo 8; the indefinite form $x^2-2y^2$ instead gives 1 or 7. The zero equation $x^2=2y^2$ has no nonzero rational solution.

If a prime $p\equiv5,7\pmod{8}$ divides $x^2+2y^2$, then it divides both (x,y) by the supplementary laws of quadratic reciprocity. Peeling common powers shows the valuation at every such prime is even. In particular (7,21,63) all fail this test: each has v_7=1. Thus the all-spectrum exceptions cannot be characterized by this norm because 63 is attained outside the diagonal family. Conversely, 27=5²+2·1² is a norm but is not a diagonal output: its only representation in positive odd coordinates is (5,1), which would give q_A=6. The congruence above rules that out at orders at least two, and the singleton has q_A=2. The residue-class gap is likewise sharp as an obstruction, not a characterization: a singleton diagonal has 3, TT2 has 11, TT3 has 51, and the directed triangle paired with a singleton yields 15 off diagonal. TT2 paired with TT3 gives 23 (prime 23 is inert), showing the restriction does not survive recontact.

Prime diagonal outputs have **relative density zero even among primes 3 mod 8.** Here is the structural reason, invoking the classical angular prime-ideal theorem for fixed sectors of an imaginary quadratic field. The ring $\mathbb Z[\sqrt{-2}]$ is norm-Euclidean (round the two real coefficients independently; norm of the error is at most (3/4)); a split prime has, up to signs and conjugation, only one representation as a norm. Its positive angle is $\arctan(\sqrt2\,y_A/x_A)\le\arctan(\sqrt2/|A|)$. For any fixed (K), all but finitely many diagonal prime outputs must have $|A|\ge K$, so lie in a sector of angular width (O(1/K)). The fixed-sector theorem gives upper relative density (O(1/K)) among primes (and hence among the positive-density class 3 mod 8), then let $K\to\infty$; no shrinking-sector uniformity is assumed. Primary source for sectorial distribution in imaginary quadratic fields: [Coleman, The distribution of points at which binary quadratic forms are prime](https://doi.org/10.1112/plms/s3-61.3.433). For elementary support counting, (D(A)\le X) and $|A|\ge K$ imply $x_A\le\sqrt X, y_A\le\sqrt X/K$; hence the whole diagonal family has $O(X/\sqrt{\log X})$ distinct outputs by taking $K\asymp\sqrt{\log X}$ and bounding all smaller orders by labeled tournament counts.

## The physical diagonal family is locally complete on its one dyadic sheet

For every \(L\ge3\), actual diagonal outputs hit *every* residue \(3\bmod8\) modulo \(2^L\). Choose a source \(A\) which is an ordinal sum of \(k\) directed triangles and \(s\ge L\) singletons. Hamilton counts and two-optional-slot counts multiply under ordinal sum, so \(H(A)=3^k\) and \(q_A=12^k2^s\). Consequently
\[
 D(A)=9^k\big((2^{2k+s}-1)^2+2\big)\equiv3\cdot9^k\pmod{2^L}.
\]
The powers of 9 fill precisely \(1+8\mathbb Z\) modulo \(2^L\), by \(v_2(9^{2^j}-1)=j+3\), obtained by repeated squaring. This is a local-tournament realization that genuinely matches the Collatz multiplier gate after fixing the parity of the exponent of 3. It coexists with the relative prime-density-zero theorem: for \(k>0\) the explicit local filler has a factor 9. Local saturation cannot replace absolute placement or the exact \(+2\) switch.

## Recontact is a Gram operation, not norm multiplication

The exact identity
$$
D(A)D(B)-F(A,B)^2=2(x_Ay_B-x_By_A)^2
$$
says what the off-diagonal operation preserves. It is the real part of $(x_A+y_A\sqrt{-2})(x_B-y_B\sqrt{-2})$; the determinant is its missing imaginary coefficient up to sign. Norm multiplicativity is a property of the *pair*, not of the real part. Holding the diagonal shell (D(A)=D(B)=N) gives
$$
N-F(A,B)=\tfrac12[(x_A-x_B)^2+2(y_A-y_B)^2].
$$
For distinct profiles of orders at least two the loss is at least **38**. Put $a=|x_A-x_B|/2$, $b=|y_A-y_B|/2$. Both positive, and the two coordinate differences have opposite signs. Equal norms give $a(x_A+x_B)=2b(y_A+y_B)$, so $a\le2b/3$. The common congruence $x+y=q\equiv0\pmod{4}$ gives $a\equiv b\pmod{2}$. For (b=1,2) no positive (a) is allowed; the first possible pair is ((1,3)), costing $2(a^2+2b^2)=38$. This is a lower bound, not a claim that physical profiles attain it. The singleton shell is unique. A prime shell is unique in profile by norm-Euclidean factorization. In particular no diagonal-shell sheet gives an exact +2 by comparing its diagonal to a mixed recontact. This says nothing about the gap between two independent mixed recontacts, or between unequal shells.

## A fixed switch cannot carry exact +2 to arbitrarily large hosts

Fix alternative source profiles \(A,A'\), let \(a=(x_{A'}-x_A)/2\), \(b=(y_{A'}-y_A)/2\), and vary the same host \(B\) in both joins. The half-step is \(a x_B+2b y_B\). Exact step \(\pm2\) requires \(a\ne0\) odd, since \(x_B,y_B\) are odd. If both source orders are at least two then \(a\equiv b\pmod2\), so \(b\) is odd as well. For host order at least two in that case \(a,b\) have opposite signs, since adding same-sign terms gives absolute value at least five. In every case with \(|B|\ge 2\),
\[
 |a|\,|B|\le |a|\,x_B/y_B \le 2|b|+1/y_B\le 2|b|+1.
\]
Thus a fixed pair only supplies exact unit half-steps for a bounded range of host orders, including when the two source norms do not agree. Any unbounded single-source-switch construction in this join has to change the source pair or the reader: the difference of its Hamilton counts must grow at least linearly with host order, relative to the change in the path-cover coordinate. This filters the proposed repair much more strongly than mod 8 alone.

## The relation to Collatz: a common projector, different transfer

For $k\ge3$, the subgroup generated by 3 in $(\mathbb Z/2^k)^\times$ is precisely the inverse image of ${1,3}\bmod8$: indeed $v_2(3^{2^j}-1)=j+2$ for $j\ge1$ by squaring, so the group has size $2^{k-2}$, and all its elements are 1 or 3 modulo 8. The same subgroup is the kernel of the split character $\chi_{-8}$. This is an equality of a *homogeneous multiplier subgroup* with a *quadratic norm gate*. It is not an invariant of the inhomogeneous dynamics. For any exact division depth $d\ge3$, the odd-step relation $3n+\varepsilon=2^dm$ fixes $n\equiv-3\varepsilon\pmod{8}$, i.e. (5) for plus and (3) for minus. Nevertheless **all** odd classes occur as $m\bmod8$, even at fixed (d): choose positive $m\equiv\varepsilon 2^{-d}\pmod{3}$ in the desired mod-8 class, then $n=(2^dm-\varepsilon)/3$ is a positive odd integer for (m) sufficiently large. Thus the plus/minus clocks distinguish source residue in deep edges, but do not enforce an angular/norm obstruction on destinations.

Repository leads: [October 4 recursive dependency kernel](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_recursive_dependency_kernel_20261004.md) already records the dyadic subgroup; [THM-3455’s Berggren sieve](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum.md) uses the same quadratic sign, independently of tournament gluing. The synthesis principle is to retain the lost coordinate: determinant/sign under off-diagonal gluing, and affine carry under odd-edge dynamics. The residue projector is useful for pruning candidate families, but neither theorem produces an unrestricted absolute target selector or a Collatz convergence claim.
