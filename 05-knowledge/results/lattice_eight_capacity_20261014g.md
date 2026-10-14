# Marked lattice congruences and a 57-prime extension bound for 30

**Status:** PROVED elementary block/capacity consequences; CITED van der Blij input for the lattice filters; FINITE-EXACT rational and integer certificates. No new Giuga number, higher-index example, or Collatz convergence theorem is asserted. The new 57-prime lower bound in §4 does not require the lattice theorem.

## 1. Inheritance, domain, and the useful new obligation

The closest proved mechanism is the strict index increase and reciprocal-capacity bound in [Giuga block lattice](giuga_block_lattice_20261014f.md), which gives at least 56 new primes for a proper Giuga multiple of 30. The retained sidecar is the **complete labelled old-prime set**, not just its reciprocal sum. The canonical hostile is an apparently admissible reciprocal-maximizing block that fails the old local congruences. Here those congruences exclude the entire 56-prime case.

For squarefree $N>1$, set

$$S(N)=\sum_{p\mid N}\frac1p,\qquad k_\mu(N)=S(N)+\frac\mu N,\qquad \mu\in\{-1,1\}.$$

The lattice domain requires $k_\mu(N)\in\mathbb Z$. The Giuga sheet is $\mu=-1$; a prime has index zero, but is not a composite Giuga number. The code keeps these distinctions. The signed index is not defined as an integer merely because a source is a Mersenne number or a Collatz endpoint.

Live objects in this packet are: the source-owned prime set; two characteristic vectors; the same-source block product; exact reciprocal capacity; and the bounded residual CRT obligation. Their maps preserve necessary extension predicates, while dropping labels or replacing the entire remaining product by a partial product loses those predicates.

## 2. Two characteristic coordinates and same-parity block guards

The inherited star matrix is

$$B=\begin{pmatrix}k&\mathbf1^T\\\mathbf1&\operatorname{diag}(p_1,\ldots,p_r)\end{pmatrix}.$$

Its Schur complement is $\mu/N$, so its determinant is $\mu$ and signature is $r+\mu$. Van der Blij's characteristic-vector congruence is the cited input. Mavrotheris, Theorem 1 (p.2; characteristic-norm calculation p.8), obtains

$$r+\mu\equiv k\pmod8\quad(N\text{ odd}),\qquad
N\equiv4\mu(r-k)+4-2\mu\pmod{32}\quad(N\text{ even}). \tag{1}$$

Sources: [Mavrotheris, *Characteristic-Vector Congruences*](https://zenodo.org/records/22085830), [van der Blij's original note](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/vdb.pdf). This is inherited lattice input, not a new theorem about signature. The present script checks its exact characteristic coordinates independently on its stated finite universe.

For the implementation, put $\epsilon_i=p_i\bmod2$, $t=\sum\epsilon_i/p_i$, and

$$c_0=\mu N(k-t),\qquad c_i=(\epsilon_i-c_0)/p_i.$$

Then $Bc=(k,\epsilon_1,\ldots,\epsilon_r)^T\equiv\operatorname{diag}(B)\pmod2$. Unimodularity makes $c$ integral. Direct multiplication gives norm $k$ in the odd case and $k+(\mu N+2)/4$ in the even case. In particular, for $N=30$, $\mu=-1$, the exact vector is $(-14,7,5,3)$, its norm is $-6$, and the signature is $2$.

Now suppose $NQ$ is a same-sign, same-parity squarefree extension, with $Q$ a product of $u$ distinct new odd primes and index jump $m=k_\mu(NQ)-k_\mu(N)$. Subtract (1). For odd $N$,

$$u\equiv m\pmod8. \tag{2}$$

For even $N=2n$ with $n$ odd,

$$N(Q-1)\equiv4\mu(u-m)\pmod{32},\qquad
Q\equiv1+2\mu n^{-1}(u-m)\pmod{16}. \tag{3}$$

This is a block-product guard, not a step-eight recurrence. It does not apply to a parity-changing extension.

A second convenient characteristic vector has centre zero, coefficient one on every odd-prime leaf, and coefficient $\delta=(k+r-1)\bmod2$ on the leaf 2. Its norm is $\sum_{p\text{ odd}}p+2\delta$. Writing $\delta'=(\delta+m+u)\bmod2$, subtraction gives

$$\sum_{q\mid Q}\frac{q-1}{2}\equiv\delta-\delta'\pmod4. \tag{4}$$

For the hypothetical 56-prime extension of 30 below, $m=1$, (3) gives $Q\equiv15\pmod{16}$ and (4) gives half-leaf sum $1\pmod4$. The familiar old-prime Giuga guards also give $Q\equiv1\pmod{30}$, hence

$$Q\equiv31\pmod{240}. \tag{5}$$

These are necessary constraints. Neither (3), (4), nor their conjunction makes the target index integral or checks every new prime.

## 3. Exact capacity with labels and residue constraints

Let $NQ$ be a proper Giuga multiple of a Giuga source $N$. The inherited exact identity is

$$m=\sum_{q\mid Q}\frac1q+\frac{Q-1}{NQ},\qquad
0<m<\frac1N+\sum_{q\mid Q}\frac1q. \tag{6}$$

Thus $m$ is a positive integer. For a fixed number $u$ of added primes, let $\ell_1<\cdots<\ell_u$ be the smallest primes not dividing $N$. Define

$$U_u=\frac1N+\sum_{i=1}^u\frac1{\ell_i}.$$

If $U_u\le m$, the block is impossible. If $m-U_{u-1}>0$, its largest prime obeys the strict bound

$$q_{\max}<\frac1{m-U_{u-1}}. \tag{7}$$

Indeed, the other $u-1$ primes contribute at most the smallest-prime sum. For any proposed omitted prime $p$, replacing the eligible universe by one that also excludes $p$ gives another exact upper bound. If that bound is at most $m$, $p$ is mandatory. Both tests retain strictness at equality.

For a **finite** pool and prescribed product residue modulo $D$, `bounded_capacity` computes the exact maximum reciprocal sum. Group primes by residue. For fixed class counts, the best choice consists of the smallest primes in each class. Dynamic programming over the class counts and the product residue therefore exhausts every possible subset, without enumerating its orderings. The packet retains the universe and the selected maximizer; it is not an unrestricted prime-existence certificate. Independent exhaustive-subset controls check the optimizer on small pools, including nonunit classes.

For 30 and $u=56$, the product-constrained maximizer under (5) is the ordinary first 56 eligible primes with 271 replaced by 317. Its relaxed upper bound is approximately $1.0018146886824166>1$. Thus this statistic alone does **not** exclude 56. This maximizer also fails (4), demonstrating that the product coordinate by itself can discard another characteristic constraint. We do not need to solve that further optimization to prove the next result.

## 4. A proper Giuga multiple of 30 needs at least 57 new primes

**Theorem.** If $30Q$ is Giuga, $Q>1$, and $\gcd(30,Q)=1$, then $Q$ has at least **57 distinct prime factors**. In particular such a multiple has at least 60 prime factors in total. No existence assertion is included.

**Proof.** Squarefreeness is part of the Giuga condition. The inherited capacity argument gives

$$U_{55}<1<U_{56}<2.$$

Monotonicity excludes every positive block size at most 55. If the size is 56, (6) forces $m=1$, hence target index 2. All the following comparisons use exact fractions; decimal displays only aid reading.

From (7),

$$q_{\max}<1/(1-U_{55})<794,$$

so every added prime is at most 793. Omitting any one of the 36 primes from 7 through 167 gives capacity at most one: the best possible omission replaces that prime in the first 56 by 281. Therefore **all** these 36 labelled primes are mandatory. Omitting 173 already gives a bound greater than one, so the mandatory list here is exactly the list certified by this particular test.

Let

$$M=30\prod_{7\le p\le167\atop p\text{ prime}}p
=962947420735983927056946215901134429196419130606213075415963491270.$$

Write $30Q=MB$. Exactly 20 new primes remain. They avoid every prime of $M$ and are at most 793. Consequently

$$B\le B_{\max}
=1294399569600457862304602788432784367142015736802901139771,$$

the product of the 20 largest eligible primes, namely 653 through 787 as listed in the JSON receipt.

For **each** prime $p\mid M$, the final Giuga condition is

$$\frac Mp B\equiv1\pmod p.$$

Set $C=\sum_{p\mid M}M/p$. Because $M$ is squarefree, $C$ is a unit modulo $M$. Combining all these old-prime conditions gives

$$B\equiv C^{-1}\pmod M.$$

Its least nonnegative residue is

$$b_0=241154391290894464663467115463743362761848622423921060203644724609.$$

The exact comparison is $B_{\max}<b_0<M$. No positive integer $B\le B_{\max}$ can belong to this residue class. This excludes the entire 56-prime case. $\square$

The final step is a simultaneous old-prime obstruction. Merely testing $Q\equiv1\pmod{30}$ loses the 36 forced labels and misses it. The proof does not rely on the optional lattice restrictions (3)–(5), a primality heuristic, a chosen search cutoff, or Korselt normality. In fact an even number cannot be Carmichael; adding Korselt to this example would be an unrelated immediate obstruction.

The program independently checks the same giant residue through `adaptive_giuga_completion_20261014g.pending(core(M), -1)`, retaining the entire remaining-product variable $B$. This is a direct consumer of the [adaptive fixed-count compiler](adaptive_giuga_completion_20261014g.md). It is not permission to apply a remaining-product guard to an unfinished selected prefix.

At 57 added primes, $1-U_{56}<0$. This particular largest-prime bound supplies no cap, and the present argument stops. The theorem does not classify 57-prime blocks or claim that 57 is attainable.

There is nevertheless a typed **OPEN** next obligation. The same omission test forces the 24 new primes 7 through 103. Their enlarged core is
$M_{57}=23984823528925228172706521638692258396210$; the target index is still 2 and 33 primes remain. Its positive integer deficit is
$A=2M_{57}-\sum_{p\mid M_{57}}M_{57}/p=4259058771137487325644873173720211915023$.
For the smallest remaining prime $q$, (6) in partial-core form gives $A/M_{57}<33/q$, hence $103<q\le185$. The whole remaining product must be congruent to
$6351005307643534964051443105592355986743$ modulo $M_{57}$. `next57_obligation()` retains the original source 30 and these exact data. No 57-prime completion tree has been run; a finite first branch is not a completed exclusion.

## 5. Source and inference boundaries

The lattice constructor rejects squarefree sources with nonintegral signed index. For example, $15=2^4-1$ has $k_{-1}(15)=7/15$, so no integral unimodular star is exported. The Mersenne prime 31 has index zero, but the actual odd Collatz chain $31\to47\to71\to107\to161$ leaves that integral-index domain at $161=7\cdot23$. The nonsquarefree Mersenne source 63 is also rejected. Thus the lattice marker is neither a universal Mersenne coordinate nor a preserved Collatz certificate.

Incoming [THM-4651, Giuga sum at $p^2$ and the parity divide](../../01-canon/theorems/THM-4651-giuga-sum-at-p-squared-is-the-wilson-quotient-and-the-korselt-giuga-parity-divide.md) was read as a pre-audit idea source; no step above depends on its audit status. The 56-to-57 improvement uses the independently replayed capacity/CRT certificate. Count 8 in (2), modulus 8 for characteristic norms, modulus 32 in (1), and the prime-count bound 57 are different objects, not interchangeable clocks.

## 6. Exact controls and reproduction

The script checks signed squarefree sources through 2500, independent rational determinants and characteristic parity, known larger source vectors, actual same-sign PPN extensions, independent exhaustive small-pool subset optimization, the 56-prime forced labels and strict cap, all 39 old-prime CRT guards, and malformed typed receipts. Normal and optimized execution use always-active gates, not `assert`.

Files: [script](../../04-computation/experiments/lattice_eight_capacity_20261014g.py), [JSON receipt](lattice_eight_capacity_20261014g.json). Reproduce from the repository root:

```powershell
python -B -X utf8 04-computation/experiments/lattice_eight_capacity_20261014g.py
python -B -O -X utf8 04-computation/experiments/lattice_eight_capacity_20261014g.py
```

Normal/optimized/saved output agree: **5,143 checks**, LF SHA256
`1137b17b278b3474c617880d04895b1bc5442bd169b094bcf4333b064d90cd26`.
