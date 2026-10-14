# The actual Mersenne child: lattice exclusion and a descent that does not pay its source

**PROVED scoped:** an exponent-only necessary Carmichael condition; a stronger
conditional Giuga index constraint and an order-group capacity bound for odd
Mersenne exponents; complete rejection
of one fixed two-twos head transplant; and an exact infinite exponent phase on
which a later local descent remains above the original source.
**FINITE-EXACT:** the stated modular readers and controls.
The current giant Mersenne ROOT obligation and universal Collatz remain **OPEN**.

Artifacts: [checker](../../04-computation/experiments/mersenne_lattice_child_guards_20261014g.py)
and [saved output](mersenne_lattice_child_guards_20261014g.json).
No giant source integer is expanded, and no ROOT word is discovered or asserted.

## 1. Inheritance and the unchanged source

The current minimum of the inherited 1,949-label Mersenne fan is

\[
E_*=99\,708\,991\,741,\qquad M_E=2^E-1.
\]

The [clocked-completion packet](mersenne_clocked_completion_20261014c.md)
retains this source, the complete absence of a contracting ordinary odd
predecessor through 64 odd steps, and the finite-head clock obstruction.
The corrected [THM-4605, Mersenne exits by pair-chain absorption](../../01-canon/theorems/THM-4605-mersenne-exits-by-pair-chain-absorption-the-t23-fan-escapes-at-depth-1.9e7.md)
must be read with its pre-ROOT, zero-defect scope. It does not supply a terminal
proof for this minimum. We do not repeat the saturated inverse search.

The closest positive mechanism is the actual eight-step first descent on
`8911+32768t` from the [successor-guard packet](collatz_successor_guards_20261014f.md).
The canonical hostile is its empty direct Mersenne intersection. The least-used
sidecar here is the original source: a rule can pay a later checkpoint without
paying the initial integer. Our concept board is unchanged exponent / native
word / affine carry / original-source payment / full arithmetic eligibility /
supplied ROOT terminal.

## 2. An exponent-only Carmichael exclusion

For every integer $E\ge2$,

\[
\operatorname{ord}_{2^E-1}(2)=E.
\]

Indeed $2^E\equiv1$, whereas $0<a<E$ gives $0<2^a-1<2^E-1$.
If $M_E$ is Carmichael, its base-two Fermat condition therefore forces

\[
\boxed{E\mid2^E-2.}\tag{1}
\]

The same necessary condition holds when $M_E$ is prime. It is not sufficient
for either conclusion: $E=11$ passes, but $M_{11}=2047=23\cdot89$ fails
Korselt because $88\nmid2046$.

For the unchanged exponent the exact small-integer calculation is

\[
E_*=47779\cdot2086879,\qquad
2^{E_*}\bmod E_*=13\,454\,057\,925\ne2.
\]

Thus **$M_{E_*}$ is not Carmichael**, and cannot enter a full joint
Giuga--Korselt route. This does not decide its Collatz trajectory.

There is also an all-height exponent exclusion. If an odd prime $q\mid E$ and
$d=\operatorname{ord}_q(2)$, (1) forces $E\equiv1\pmod d$.
The checker proves $d=15926$ for $q=47779$. Consequently every

\[
E\equiv27\,377\,367\pmod{760\,928\,354}
\tag{2}
\]

is excluded: this class is zero modulo $q$, but is $573$, rather than one,
modulo $d$. It contains $E_*$. The API authenticates the prime, its exact
order, the source divisibility and the incompatible phase; passing another
phase yields no positive arithmetic certificate.

## 3. Odd-exponent joint sources would need index divisible by four

Suppose $E\ge3$ is odd and $N=M_E$ is both Giuga and Carmichael, with distinct
prime factors $p_1,\ldots,p_r$. Write its positive Giuga index as

\[
\sum_{p\mid N}\frac Np=\kappa N+1.
\tag{3}
\]

Korselt and $v_2(N-1)=1$ imply $p\equiv3\pmod4$ for every factor.
Also $\operatorname{ord}_p(2)\mid E$ is odd. An odd divisor of $p-1$ divides
$(p-1)/2$, so $(2/p)=1$. The supplementary quadratic-reciprocity formula
$(2/p)=(-1)^{(p^2-1)/8}$ now gives

\[
\boxed{p\equiv7\pmod8\quad\hbox{for every }p\mid N.}\tag{4}
\]

Since $N\equiv7\pmod8$, the number $r$ is odd; every cofactor $N/p$ is one
modulo eight. Reducing (3) gives $r\equiv1-\kappa\pmod8$.
Independently, the inherited odd Giuga star-lattice constraint is
$r\equiv\kappa+1\pmod8$. It follows that

\[
\boxed{4\mid\kappa,\qquad\kappa\ge4.}\tag{5}
\]

For completeness, the star has central diagonal $\kappa$, prime diagonals
$p_i$, and one on every central spoke. Its Schur complement is $-1/N$,
so it is integral unimodular of determinant $-1$ and signature $r-1$.
The central basis vector is characteristic and has norm $\kappa$.
The only non-elementary input is the **CITED**
[van der Blij signature congruence](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/vdb.pdf),
used exactly as in [the inherited star-lattice construction](giuga_block_lattice_20261014f.md#5-a-literal-graph-shaped-lattice-and-its-eightfold-constraint).
That construction credits
[Mavrotheris's characteristic-vector application](https://zenodo.org/records/22085830).

Neither (4) nor (5) asserts that such a joint source exists, or that every
Mersenne source is excluded. The current source is excluded separately by (1).
No Collatz terminal follows from the signature congruence. The checker verifies
the residue implications and prime-order controls, not a nonexistent positive
joint example.

### An independent exclusion from the exponent's divisor count

There is a further source-preserving bridge that needs neither a factorization
of $M_E$ nor the Korselt hypothesis. For odd $E\ge3$, group the distinct prime
divisors of $M_E$
by $d=\operatorname{ord}_p(2)>1$, where $d\mid E$. Every such $d$ is
odd and at least three. Since $d\mid p-1$ and $p$ is odd, $p\ge2d+1$.
If a group contains $r_d$ primes, their product divides $2^d-1$, giving

\[
r_d\log_2(2d+1)<d,\qquad
S_d:=\sum_{\operatorname{ord}_p(2)=d}\frac1p
\le\frac{r_d}{2d+1}
<\frac1{2\log_2(2d+1)}<\frac15.
\]

The last strict inequality follows from $2d+1\ge7$ and $7^2>2^5$.
Empty groups contribute zero. Summing over the $\tau(E)-1$ nontrivial
divisors proves the exact bound

\[
\boxed{S(M_E):=\sum_{p\mid M_E}\frac1p<\frac{\tau(E)-1}{5}.}\tag{9}
\]

Thus **$\tau(E)\le6$ excludes Giuga entirely**, since a positive integral
Giuga index requires $S=\kappa+1/M_E>1$.
For a joint source, (5) requires $S>4$, so (9) forces
**$\tau(E)\ge22$**. These are necessary conditions with unbounded prime
factors; they do not assert existence beyond either threshold.

The fully authenticated exponent factorization above has $\tau(E_*)=4$.
Consequently $S(M_{E_*})<3/5$: **the actual source is not Giuga independently
of its failure of Korselt**. Because $E_*$ is composite, $M_{E_*}$ is composite,
so $S-1/M_{E_*}>0$. Both signed indices satisfy

\[
0<S\pm1/M_{E_*}<\frac35+\frac17=\frac{26}{35}<1.
\]

Neither signed star construction has an integral central index at this source.
The checker consumes a full prime-power factorization of the *unchanged
exponent*, rejects partial or altered factorizations, and verifies (9) by an
independent small support factorization for odd exponents from 3 through 31.
For the minus sign, positivity is explicitly restricted to composite $M_E$;
a prime $M_E$ instead has minus index zero.

## 4. Why the 8911 word cannot be inserted into the two-twos template

Retain the actual word

\[
u=(1,1,1,3,1,1,1,5),\quad
F_u(x)=\frac{6561x+9425}{16384}.
\]

Its native odd source cell is $8911\pmod{32768}$ and it maps
$8911+32768t$ to $3569+13122t$. All seven proper prefix slopes exceed one,
while the full word pays $5342+19646t>0$. No positive Mersenne lies in that
source cell: for $E\ge15$ its residue is 32767; for smaller $E$, the equation
$2^E=8912$ has no solution. These are inherited facts, retained here as guards.

After a Mersenne's two-twos prefix, the D3/D4 child clearing in
[the signed-ladder packet](collatz_complement_ladders_20261007c.md)
has marked relation $X=81Y+10$. A partner $v$ with signed sibling gap $g$
would have length 12, cost $14-2g$, and carry

\[
B_v=9425+10\cdot6561+\frac{16384}{3}(1-4^{-g}).\tag{6}
\]

Positive valuations require $14-2g\ge12$, hence $g\le1$.
Positive carry excludes $g\le-2$. The remaining gaps $-1,0,1$ have respective
carries $58651,75035,79131$, all less than $3^{11}$, the first unavoidable
carry contribution of a twelve-letter word. Hence **no signed gap works in
this entire fixed-head template**. The complete carry decoder independently
rejects those three candidates. This is not a claim about other heads,
longer continuations or changed anchors.

## 5. The rule does occur later, but its payment is local

Write

\[
X_E=U^{E-1}(M_E)=2\cdot3^{E-1}-1.
\]

A bounded exploratory read with 200,000 dyadic bits and at most 65,536
post-ascent odd steps located the 8911 cell. The production checker retains
only the fixed occurrence at depth $d=18843$: it reads exactly $d+8$
letters with **37,431 bits**, rather than searching until a hit occurs.
Let $w$ be the first $d$ letters and $A=\sum w=37416$. Exact modular division
proves

\[
Z_E=U^d(X_E)\equiv8911\pmod{32768}
\quad(E=E_*),\qquad
U^8(Z_E)=F_u(Z_E)<Z_E.
\tag{7}
\]

The prefix is reproducibly specified by the source formula, depth, bit budget
and its ASCII comma-separated SHA256 in the saved JSON. A second path composes
its complete affine carrier and checks the final oddness congruence. Every
endpoint residue query keeps the full denominator $Q$ in its working modulus
$Qm$ before division. The giant integer is never constructed.

This finite read has an exact infinite-family scope. The total head $wu$ has
cost $A'=37430$ and a single native source residue modulo $2^{A'+1}$.
Substituting $X_E$ reduces its guard to one class of $3^{E-1}$ modulo $2^{A'}$.
Since the order of three is $2^{A'-2}$, it is equivalent to

\[
\boxed{E\equiv E_*\pmod{2^{37428}}.}\tag{8}
\]

Every positive exponent in (8) is at least $E_*$ and obeys
$e=E-1\ge2A'$. For each head prefix of length $j$ and cost $a\le A'$, the
full affine slope from the original $M_E$ is

\[
\frac{3^{e+j}}{2^{e+a}}
=\left(\frac32\right)^e\frac{3^j}{2^a}>1,
\]

because $(3/2)^{2A'}=(9/4)^{A'}>2^{A'}$. The affine carry is positive.
The initial ones run also grows strictly. Therefore there is no earlier ROOT,
and **every retained checkpoint and the endpoint of (7) exceed $M_E$**.
The eight-step suffix is a genuine first descent relative to $Z_E$ but is
unpaid relative to the immutable original source. Its child still has no
supplied ROOT certificate here.

The arithmetic exclusion (2) and this native phase coexist: their intersection
is

\[
E=E_*+380464177\cdot2^{37428}t,\qquad t\ge0.
\]

Thus one can simultaneously authenticate a local descent, an original-source
payment failure and a full-joint-arithmetic exclusion. None is a substitute
for the missing original-source paid exit or a supplied terminal proof.

## 6. Reproduction and remaining obligation

Run from the repository root:

```text
python -B 04-computation/experiments/mersenne_lattice_child_guards_20261014g.py
python -B -O 04-computation/experiments/mersenne_lattice_child_guards_20261014g.py
```

The explicit universes are literal order controls $2\le E\le20$; odd primes
below 2000; index/count residues 1 through 32; 64 literal members of the 8911
ray; small formal readers $3\le E\le24$; and the one fixed giant-source
modular prefix. Invalid exact types, insufficient precision, wrong exponent
phases, unsupported stage tags and ROOT padding are rejected. Assertions are
not used for validation. Normal and optimized output must agree exactly.

The next genuine completion must retain the same source and supply either a
paid child with its own grounded proof or an independent ROOT word. The present
packet gives new arithmetic restrictions and a source-authenticated stopping
reason for the proposed local-rule repair. It does not turn a finite quotient,
the lattice index, or a later local descent into global Collatz coverage.

Frozen normal and optimized output: **2,320 checks**, LF SHA256
`aca6765bdde0e61d96127cdbd4ff69be75c4354d9c351333e671cb8c25af0d3c`.
An independent peer proof/API audit and optimized replay passed with the same
count and hash, including the order-group bounds and original-source payment
boundary. No correction remains.
