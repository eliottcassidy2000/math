# Two scalable exact sevenfold squares and a joint ratio/Jacobian selector

**PROVED** for the formulas at every nonempty filler, the simultaneous modular selector, and the fixed-modulus scale estimate. **FINITE-EXACT** for the sample prime counts. Compare the [one-apex sparsity theorem](ONE-APEX-INVERSE.md), the [two-port Jacobian](UNRESTRICTED-JACOBIAN.md), and [prior two-apex run signatures](../chanpak_hadwiger_planar_folds_20261122/TWO-APEX.md).

## Two modules, two marked reversals, exact at every width

Use ordered spine positions with forward orientation and two apices \(a,b\) with \(b\to a\). A spine tag 0,1,2,3 is the bitmask of apices beating that position (1 for \(a\), 2 for \(b\)). The symbol \(M_j\) denotes **any nonempty tournament** inserted as a uniform block with external tag \(j\); all earlier spine positions beat it and it beats all later spine positions. Let \(D=H(M)\), \(Q=Q_2(M)=\sum_{S\subseteq V(M)}H(M[S])H(M[V(M)\setminus S])\), including empty slots. The two frames and their \(00,10,01,11\) counts are:

| Frame | Ordered tags | Mark 1; mark 2 | Height table |
|---|---|---|---|
| A | \(2,M_0,3,0,3,3,3\) | \(a\) to penultimate position; \(a\) to last | \(6Q+5D,\ 20Q+17D,\ 34Q+29D,\ 42Q+35D\) |
| B | \(1,0,M_1,3,2,3,3\) | \(a\) to first fixed 3 after block; \(b\) to last | \(15Q+D,\ 21Q-3D,\ 75Q+23D,\ 105Q+7D\) |

Both upper heights are exactly seven times the lower. Both orders are strictly increasing since \(Q\ge2D>0\). Both upper tournaments are strong: the spine marches from first to last through the module (each module vertex can enter from its predecessor block and exit to its successor block). At A the last position reaches \(a\), then the first 3 after the block, then the following 0, then \(b\) and the first position. At B the last reaches \(b\to a\), then the first position. These routes also place both apices in the spine's SCC. The lower last position is a sink, so the upper result does not depend on the lower being strong on the entire displayed vertex set.

**Run proof and coefficient certificate.** A Hamilton path reads the apices in one of two orders, with three spine runs. Assign each of the six fixed spine positions to a run. In a run, these positions must be increasing, and any used part of the module occupies a single contiguous block in the sorted order of external positions. Choose which runs receive nonempty module parts. The first run after \(a\) in A, or after \(b\) in B, cannot meet the module: that apex loses to the module and to all fixed positions before the module, so it could not start the run. No three occupied runs are possible. Test the boundary arcs against one placeholder in each occupied run, including the direct apex arc if the middle is empty. One occupied run contributes \(D\); two specified occupied runs contribute \(Q-2D\) by taking all ordered nonempty two-slot covers. This also explains why replacing the transitive block by an arbitrary tournament changes only the two readings, not the coefficients.

For each corner, the counts of accepted skeletons with one and two occupied runs are:

| Frame / apex order | 00 | 10 | 01 | 11 |
|---|---:|---:|---:|---:|
| A, \(a,b\) | (3,1) | (9,3) | (15,5) | (21,7) |
| A, \(b,a\) | (14,5) | (48,17) | (82,29) | (98,35) |
| B, \(a,b\) | (13,6) | (18,9) | (73,24) | (98,39) |
| B, \(b,a\) | (18,9) | (21,12) | (100,51) | (119,66) |

[Replay](jacobian-squares-check.js) enumerates the constant-size skeletons and independently evaluates literal tournaments with transitive and nontransitive modules by subset DP. In particular, with \(M=TT_t,\ x=2^t\), we have \(D=1,Q=x\), giving the two stated exponential rails of exact squares. Some prime starts are \(101,197\) in A and \(241\) in B. The nontransitive module with masks \((30,4,8,18,6)\), \(D=5,Q=64\), gives A's exact square \((409,1365,2321,2863)\); masks \((18,4,9,19,6)\), \(D=11,Q=100\), give B's \((1511,2067,7753,10577)\). These two primes lie off the respective transitive-filler rails. This is not a claim of infinitely many prime starts, or of minimal edit distance for those examples.

## Exact joint steering of the ratio and curvature

Let \(h=6Q+5D\) in A. The half-mixed derivative \(d=(h_{11}+h_{00}-h_{10}-h_{01})/2\) is \(-3k\), with \(k=Q+D\). The linear change 
\[
 (Q,D)\longmapsto(h,k),\qquad Q=h-5k,\qquad D=6k-h
\]
is unimodular. Thus the [physical joint slot controller](../chanpak_fork_curvature_20261123/ALL-PRIME-WORDS.md#complete-two-reading-image-including-powers-of-three) prescribes **both** \(h\) and \(k\) modulo any fixed odd modulus, with no constraint away from 3. At powers of 3 the only mod-three condition is \(h\equiv k\pmod3\) for even filler order or \(k\equiv0\pmod3\) for odd filler order. All selected physical squares still satisfy the **exact integer** equality \(h_{11}=7h\). In particular, any desired height residue mod 7 can coexist with a nonsingular gain plane there, by choosing \(k\not\equiv0\pmod7\). This contrasts with tuning the residue of an approximate sevenfold height. In B, \(d=6(2Q-D)\); for dyadic \(Q=2^t,D=1\), its Jacobian is nonsingular mod 7 iff \(t\not\equiv2\bmod3\), while A is nonsingular mod 7 at every dyadic width.

**Both middle readings can avoid a whole fixed wheel.** Let \(m=9\prod_{5\le\ell\le y}\ell\) and fix a unit \(r\bmod m\) for the residue of \(h\). Use even filler order. At each \(\ell\ge7\), put \(D=(r-6Q)/5\). Keeping \(D,Q,k\) and both middle readings nonzero forbids at most five choices of \(Q\): the middle numerators are \(17r-2Q\) and \(29r-4Q\), and \(5k=r-Q\). At 5 take \(Q=r\) and exclude \(D=0,-r\), which makes the middle values \(2D,4(r+D)\) nonzero and also protects \(k\). At 3 take \(Q\equiv2r\bmod3\), \(D=(r-6Q)/5\bmod9\); then \(D\equiv Q\bmod3\), \(k\) and the first middle are units, and the second middle is divisible by 3. Of the three lifts to 9, discard the unique one making the second middle divisible by 9. The literal square then has both \(h\) and its first middle \(y\)-rough, and its second middle divided by **exactly 3** is \(y\)-rough. (Only the odd part of the wheel is imposed on \(Q\), which is even.) At primes \(\ell\le y\) other than 3, the gain plane is nonsingular too. The three height forms \(h,h_{10},h_{01}\) have pair determinants of magnitude 2,4,2, so all three pair gcds equal the common filler content \(\gcd(D,Q)\); their heights are odd and an odd prime cannot divide two independent forms after common content is removed. Division here normalizes arithmetic readings; it does not construct the divided counts.

**Absolute scale, fixed wheel.** Put \(\theta=1-\tfrac12\log_3(3/2)\). Any permitted joint residue pair \((h,k)\) admits a literal exact-ratio square with \(|h-X|=O_m(X^\theta)\) at every sufficiently large real scale \(X\). Indeed the [native height selector](../chanpak_heegner_absolute_targets_20261123/BALANCED-WINDOWS.md#native-factor-seven-fork-height-a-second-base-32-window) approximates \(280Q+7D\) at scale \((280/6)X\), with the inferred \((Q,D)\) residues and with \(3^n\asymp X\) for its chain length \(n\). Its proof bounds \(D=O(2^n)=O(X^{\log_3 2})\). As \(6Q+5D=(6/280)(280Q+7D)+(5-42/280)D\) and \(\log_3 2<\theta\), the same window applies. Taking two nearby target scales on opposite sides with a sufficiently large constant gives a straddle, in the same exact cofactor residue class. Unlike the earlier approximate fork alone, each such cofactor is already the actual lower corner of its own literal \(h\to7h\) square. This still does not say \(h=p\) for a prescribed prime \(p\).

**Absolute residual to solve.** For a prescribed lower integer \(P\), A reduces exact admission to finding a single module with slot pair \((Q,D)=(P-5k,6k-P)\), for an integer \(k\). Since \(D\ge1,\ Q\ge2D\), the interval is \((P+1)/6\le k\le3P/17\). The constraints are necessary; modular surjectivity in that interval is not a proof of integer realization. This is the point at which a stronger path-slot inverse, not a real local Jacobian, is needed.

## Subsequent absolute capacity barrier

The [glued one-apex profile theorem](MODULE-CAPACITY.md) bounds the integer image even after arbitrary ordinal joins and a fixed additional bank. It also gives explicit six-vertex modules outside the glued profile set and proves the first scalar and prime escapes of the unaugmented class.
