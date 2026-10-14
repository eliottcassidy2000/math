# Adaptive prime completion with a retained index and pending guards

2026-10-14. **PROVED:** a complete terminating mathematical search for each fixed squarefree core, sign and added-prime count, optionally with a fixed integer index and the full Korselt target. **FINITE-EXACT:** the catalogs below and their independent controls. No bound on an unknown successful added-prime count is proved; the classical Giuga conjecture remains open.

The useful change is to retain an unfinished reciprocal deficit and the conditions on the *remaining product*. A proper prefix need not already be a signed integer. Requiring that would discard real completions.

## 1. Inheritance and the precise target

The [signed one-prime packet](giuga_korselt_defects_20261014f.md) gives a complete one-prime classifier and distinguishes old-prime CRT conditions from new-prime obligations. The [two-prime/lattice packet](giuga_block_lattice_20261014f.md), §2, removes a prime cutoff by an exact divisor identity. Its source-specific reciprocal budget, §4, also excludes small blocks above an already Giuga core. We extend that fixed-count method recursively, rather than demand a signed intermediate at each step.

The closest hostile is the opposed pair 30/561: signed cofactor conditions and Korselt periods are separate. The sharper procedural hostile is

\[
6\longmapsto6\cdot7\cdot41=1722.
\]

At the initial core, the remaining product must be $5\bmod6$. The first chosen prime $7$ fails that test, while $7\cdot41$ passes. A filter applied prematurely to a proper prefix loses this valid Giuga completion.

Our concept board is **retained core / integer index / reciprocal deficit / remaining-product phase / all-prime authentication / decreasing remaining count**. The anchor is completion; the niche is the combined pending CRT; the wildcard is the relation to a finite proof search with an unbounded outer parameter. The relevant [META-PATTERNS](../../00-navigation/META-PATTERNS.md) are retaining the source under a quotient, and distinguishing local support from bounded coverage.

Fix a squarefree integer $N>1$, an added count $t\ge0$, and a sign $\mu\in\{-1,1\}$. Seek ordered distinct primes

\[
q_1<\cdots<q_t,\qquad q_i\nmid N,
\]

such that, for $T=N\prod_iq_i$,

\[
\kappa_\mu(T)=\sum_{p\mid T}\frac1p+\frac\mu T=k\in\mathbb Z.
\tag{1}
\]

The negative sign gives Giuga numbers when $T$ is composite; positive sign and $k=1$ gives primary pseudoperfect numbers. We do not assume index one in the general algorithm. In the older cofactor convention $p\mid T/p-\alpha$, the sign is $\alpha=-\mu$.

The optional joint target additionally requires $T$ composite and

\[
p-1\mid T-1\quad\hbox{for every }p\mid T.
\tag{2}
\]

This is the full Korselt test, including all newly adjoined primes. For $\mu=-1$ the composite joint target is exactly the classical Giuga-counterexample condition. Prime targets of index zero may occur in the signed-only zero-count interface; the joint interface explicitly rejects them.

## 2. Finite indices and exact recursive state

Write $C(M)=\sum_{p\mid M}M/p$ and $S(M)=C(M)/M$ for a squarefree partial product $M$. Suppose $s\ge1$ primes remain and the last chosen prime is $b$ (initially $b=1$). For a fixed proposed target index $k$, retain

\[
A=M(k-S(M))=kM-C(M).
\]

If $R$ is the product of the remaining primes, (1) is exactly

\[
\boxed{\frac A M=\sum_{q\mid R}\frac1q+\frac\mu{MR}.}
\tag{3}
\]

The right side is strictly positive. For the negative sign, its negative term is smaller than $1/q$ for even the first remaining prime, since $M\ge2$ and $R\ge q$. Thus **every unfinished successful node has $A>0$**. At a finished node $s=0$, the condition is instead $A=\mu$; positivity must not be imposed there. The two boundaries are tested separately.

Let $\ell_1<\cdots<\ell_s$ be the least $s$ primes above $b$ not dividing $M$. The right side of (3) is at most

\[
U(M,b,s,\mu)=\sum_{i=1}^s\frac1{\ell_i}
                  +\frac\mu{M\prod_i\ell_i}.
\tag{4}
\]

Indeed every ordered candidate satisfies $q_i\ge\ell_i$, and the expression is strictly decreasing in every coordinate. Its partial derivative in a positive real coordinate $q_i$ has the sign of

\[
-\left(1+\frac\mu{M\prod_{j\ne i}q_j}\right)<0.
\]

Initially this gives the finite index interval

\[
S(N)<k\le S(N)+U(N,1,t,\mu)\qquad(t\ge1).
\tag{5}
\]

The interface may enumerate these integers or select one supplied exact integer $k$. If it lies outside (5), the answer is empty. For $t=0$ the interface directly checks $C(N)+\mu=kN$, with no strict inequality from (5).

Now let $q$ be the next remaining prime. Multiplying (3) by $Mq$ yields the safe branch bounds

\[
\begin{array}{ll}
\mu=-1:& Aq<sM,\quad q\le\left\lfloor\dfrac{sM-1}{A}\right\rfloor,\\[2mm]
\mu=+1:& Aq\le sM+1,\quad q\le\left\lfloor\dfrac{sM+1}{A}\right\rfloor.
\end{array}
\tag{6}
\]

For the positive case use $q\sum_i1/q_i\le s$ and $q/R\le1$. For the negative case the subtracted term makes the inequality strict. These bounds are derived at the current node; they are not chosen prime cutoffs.

Selecting a native prime $q>b$, $q\nmid M$, updates the entire state by

\[
(M,A,s,b)\longmapsto(Mq,Aq-M,s-1,q).
\tag{7}
\]

This follows either by subtracting $1/q$ from the gap or from $C(Mq)=qC(M)+M$. The same target index is retained. If another prime remains, a nonpositive updated deficit is impossible and is rejected. The deficit itself need not decrease; the decreasing integer is the **number of remaining factors**.

### 2.1 Exact final consumers

For one remaining prime, solve

\[
Aq=M+\mu.
\tag{8}
\]

It gives at most one rational candidate. Authenticate integrality, primality, $q>b$ and $q\nmid M$.

For two remaining primes $q<r$, the inherited identity becomes

\[
\boxed{(Aq-M)(Ar-M)=M^2+\mu A.}
\tag{9}
\]

The factors are positive: for example $Aq-M=(Mq+\mu)/r>0$. Enumerate positive divisor pairs $u<v$ of the right side, require $A\mid u+M,v+M$, and recover $q=(u+M)/A$, $r=(v+M)/A$. Retain the ordering and all native prime guards. A nonpositive right side has no candidate.

**Completeness and termination.** Every solution has one ordered prime address and one index in (5). Equation (6) includes its first choice, (7) retains the exact remaining equation, and (8)/(9) include its terminal choices. Conversely an authenticated leaf satisfies (1), and the code independently checks every local condition $p\mid T/p+\mu$ as well as the exact integer index. Each nonterminal node has finitely many children and decreases $s$; the terminal divisor lists are finite. Thus the mathematical algorithm terminates for every fixed input, with no prime-height cutoff.

Trial division is sufficient for the primality and factorization steps. No practical time, memory or polynomial complexity bound is claimed. A resource failure in a large implementation run is not an empty catalog. Trying every added count is only an outer semidecision search: fixed-count termination supplies no universal stopping count.

## 3. Pending old-prime conditions and the joint target

For squarefree $M$, $C(M)$ is a unit modulo $M$: at a prime $p\mid M$, all terms of $C(M)$ except $M/p$ vanish, and that last term is a unit. The signed conditions for the primes *already in $M$* are therefore exactly

\[
R\equiv-\mu C(M)^{-1}\pmod M.
\tag{10}
\]

This congruence constrains the whole remaining product $R$. It is not an instruction to require a newly selected prime or an unfinished prefix to satisfy the final signed criterion. Equation (3) implies (10), so the deficit state already retains it; `pending` makes it explicit and `marked_trace` checks it at every stage of an accepted completion.

For the joint target, put $\lambda(M)=\operatorname{lcm}_{p\mid M}(p-1)$. The old Korselt conditions add

\[
MR\equiv1\pmod{\lambda(M)}.
\tag{11}
\]

They are impossible when $\gcd(M,\lambda(M))>1$. Otherwise

\[
R\equiv M^{-1}\pmod{\lambda(M)},
\]

and this combines with (10) to one explicit class modulo $M\lambda(M)$, because the two moduli are coprime. The code may safely prune an impossible old Korselt state. Every final accepted joint leaf nevertheless checks (2) for *all* primes, including the last two selected by (9).

The pending congruence does not encode the target index, the remaining factor count, primality, or the new-prime periods. Those coordinates are retained separately. At count two, $561\cdot2\cdot59=66198$ is a positive signed-only control whose full Korselt audit fails. At core six, testing $7$ rather than $7\cdot41$ in (10) is the canonical premature-filter hostile.

The API's joint result is compared against the signed-only result followed by a separate full Korselt test throughout the declared small universe. This verifies the pruning direction without turning an old-prime guard into a full completion certificate.

## 4. Complete three/four-prime catalogs for the requested cores

These are complete for the specified core, sign and exact added count, with no bound on the added primes. All possible indices from (5) are included.

| Core | Added primes | Giuga sign $\mu=-1$ | Positive sign $\mu=+1$ |
|---:|---:|---|---|
|30|3|none|none|
|30|4|none|none|
|561|3|none|none|
|561|4|none|$(2,101,149,3109)$, index 1|
|1105|3|none|none|
|1105|4|none|none|
|1729|3|none|none|
|1729|4|none|none|

Every joint signed-plus-Korselt catalog in this table is empty. For core 30, (5) contains no integer even before branching. This agrees with, but is much weaker than, the inherited minimum 56-new-prime budget for a proper Giuga multiple of 30. The [marked lattice/capacity companion, §4](lattice_eight_capacity_20261014g.md) sharpens that minimum to 57: capacity forces an enlarged old-prime core, and its pending residue exceeds the maximum remaining product. That stronger result consumes (10); it is not needed to prove these catalogs.

The positive control is

\[
561\cdot2\cdot101\cdot149\cdot3109=52\,495\,396\,602.
\]

This is a known primary pseudoperfect number, explicitly listed in [Sondow–Tsukerman, §4, Definition 4 and the following list](https://arxiv.org/html/1401.0322). The packet does not claim a new number or a new general Egyptian-fraction method. Its contribution is the retained-core fixed-count classification and an executable adaptive consumer with complete native guards.

The marked successful branch is:

| Chosen new primes | Partial product $M$ | Remaining count | Deficit $A$ | Remaining product |
|---|---:|---:|---:|---:|
|none|561|4|290|93574682|
|2|1122|3|19|46787341|
|2,101|113322|2|797|463241|
|2,101,149|16884978|1|5431|3109|
|2,101,149,3109|52495396602|0|1|1|

The first four partial products are not positive-sign integral-index integers. Rejecting them for that reason would lose the completion. At the terminal, the signed local residues vanish at all seven primes, while the Korselt residue at $p=3$ is $1\bmod2$. This is a PPN completion of a Carmichael source, not preservation of the Carmichael property.

The complete catalogs were independently reproduced by a second rational-gap recursion. It subtracts $1/q$ as a `Fraction`, uses a coarser initial index bound and the common loose bound $q\le(s+1/M)/(k-S(M))$, and goes down to the one-prime equation. It uses neither the two-prime factorization nor the main capacity pruning nor the integer-deficit update. Its 2,852 visited nodes give exactly the same sixteen signed-only catalogs.

## 5. Interfaces, controls and the next obligation

[Script](../../04-computation/experiments/adaptive_giuga_completion_20261014g.py), [saved JSON](adaptive_giuga_completion_20261014g.json).

`compile_completion(core,count,mu,joint=False,target_index=None)` authenticates the inherited exact `Core`, computes the finite index interval, and returns the complete ordered catalog for those inputs. `Completion` retains the original source, sign, exact index, added prime tuple, target and joint-mode marker. `audit_completion` independently verifies all of them. `pending` and `accepts_pending` expose the remaining-product guard; `marked_trace` exposes the successful recursive state without claiming each unfinished product is already terminal.

All public source/receipt consumers reject Boolean and floating aliases before equality or arithmetic. Warm-cache factorization is reached only through the inherited exact-integer validator. Zero-count terminals and the signed terminal value $A=-1$ have explicit controls.

The finite control universe consists of:

* all 49 squarefree cores in $2\ldots80$, both signs and added counts $0\ldots3$;
* 8,464 independent tuples from the nine-prime bank through 23, compared through every individual local signed condition, with a separate all-prime Korselt comparison;
* agreement with the inherited complete one- and two-prime compilers;
* the sixteen requested signed catalogs, their sixteen joint counterparts, and the independent rational-gap tree;
* wrong core, sign, index, phase, target, duplicate/nonprime/order and exact-type hostiles; the successful unfinished-prefix and unsuccessful signed-to-Korselt controls above.

The finite prime bank is an independent test universe only. It does not limit the production compiler or the requested catalogs.

```powershell
python -B -X utf8 04-computation/experiments/adaptive_giuga_completion_20261014g.py
python -B -O -X utf8 04-computation/experiments/adaptive_giuga_completion_20261014g.py
```

Both runs agree with saved JSON: **17,743 checks**, LF SHA256 `8af10c34eb5b4bb1de656ae2bac9c4a4dea77ce0266f268573287b0f38d559f5`.

The next global obligation is explicit. To turn adaptive block search into a proof that a fixed core has no joint completion of any size, one needs a valid bound on the successful added count or an invariant excluding the entire remaining tree. Repeatedly finding empty finite catalogs does not supply that bound. Conversely any accepted joint leaf would be a completely authenticated arithmetic witness, without relying on a prime-density or CRT-existence heuristic.
