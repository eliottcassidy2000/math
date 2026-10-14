# Fermat primes, Wieferich delay, and exact source-residue lifting

**Status:** PROVED elementary identities and complete phase compiler;
FINITE-EXACT implementation controls; CITED classical arithmetic and rotation
groups. No completeness claim for either prime list, no Collatz convergence
claim, and no new graph-colouring consequence.

The useful distinction is between two defects. The order of 2 modulo a prime
can already miss many nonzero residues. Independently, its first nonzero
prime-adic defect can delay the appearance of new digits. Fermat primes have
the first defect, except at 3 and 5, but never the second. Wieferich primes
have the second. Keeping both yields an exact compiler for prescribed source
residues of odd inverse edges.

## 1. Inheritance and the actual connection

- [THM-4631, level factorisation](../../01-canon/theorems/THM-4631-level-factorisation-of-the-compressed-twisted-transfer-operator-of-the-px-plus-1-maps-fast-modes-are-the-doubling-cycles-on-the-units-mod-p-i.md)
  already uses the cycles of multiplication by 2 on prime-power units. We
  rederive only that elementary permutation layer, not its full operator.
- [The corrected functional-graph note, Theorem R_a and audit remark](collatz_functional_uniqueness_20261001.md)
  distinguishes primitive-root separation from non-Wieferich lifting.
  MISTAKE-553 repaired their conflation; the broader rigidity remark is only
  sketched there and is not used here.
- [The Fermat/Platonic note, sections 1, 2, 4](gilbreath_fermat_platonic_20260926.md)
  already separates the binary kernels, constructible polygons and five
  Schläfli types. [THM-871, Fermat-rung rigidity](../../01-canon/theorems/THM-871-fermat-rung-rigidity.md)
  concerns tournament multipliers, not transitivity of doubling.
- Cheapest hostile: 17 is a Fermat prime and is non-Wieferich, but 2 has
  order 8, not 16, modulo 17. The endpoint 3 has no inverse edge under either
  odd map $17n+1$ or $17n-1$.
- Least-used retained datum: the first **nonzero** defect
  $(2^d-1)/p^s\pmod p$, together with the target residue and exponent phase.

The concept board is: binary prime shape; order and lifting delay; native
source residue; pointed orbit count; Platonic group type. The common move is
to retain a phase before passing from a count to an individual source. The
packet adds an executable source-bound compiler and its precise failure
branches. The underlying lifting and non-Wieferich arithmetic are classical.

## 2. The complete split/grow tower

Let $p$ be an odd prime and put

$$
d=\operatorname{ord}_p(2),\qquad s=v_p(2^d-1)\ge1,
\qquad u=(2^d-1)/p^s\not\equiv0\pmod p.
$$

For every $k\ge1$,

$$
L_k=\operatorname{ord}_{p^k}(2)=d p^{\max(0,k-s)},\qquad
C_k=\frac{(p-1)p^{k-1}}{L_k}
=\frac{p-1}{d}p^{\min(k-1,s-1)}.                 \tag{1}
$$

**Proof.** Binomial expansion for odd $p$ gives
$v_p((1+p^su)^t-1)=s+v_p(t)$ for every positive integer $t$:
an exponent prime to $p$ preserves the valuation, and each multiplication of
the exponent by $p$ raises it by exactly one. Any exponent returning to 1
modulo $p^k$ must first be a multiple of $d$. The least permitted multiplier
then proves (1). Multiplication by 2 is a permutation of the units, and every
orbit has length $L_k$, giving $C_k$.

Under reduction from level $k+1$ to $k$, each cycle has $pL_k$ points above
it. For $k<s$ these form **$p$ separate cycles of unchanged length**; for
$k\ge s$ they form **one cycle of length $pL_k$**. Reduction maps each such
cycle onto the original cycle. Thus the two independent quantities are
the initial number $(p-1)/d$ of cycles and the delay $s-1$.

Since $(p-1)/d$ is prime to $p$,
$v_p(2^{p-1}-1)=s$. Consequently the base-2 Wieferich condition is exactly
$s\ge2$. This does not assert that 2 generates the units modulo $p$.

The following rows are independently checked by modular exponentiation,
including minimality of the orders and nonvanishing of the next defect.

| $p$ | $d$ | $s$ | $u\pmod p$ | $C_1$ | $C_2$ | $L_3$ |
|---:|---:|---:|---:|---:|---:|---:|
|1093|364|2|891|3|3279|397852|
|3511|1755|2|1781|2|7022|6161805|

These two primes illustrate a one-level delay. Equation (1) also covers any
higher delay without assuming such a prime exists.

## 3. Fermat primes have no vertical delay

Suppose $p=2^{2^m}+1$ is prime, with $m\ge0$. Then

$$
d=2^{m+1},\qquad 2^d-1=p(p-2),\qquad s=1.       \tag{2}
$$

Indeed, $2^{2^m}\equiv-1\pmod p$, so its order is exactly $2^{m+1}$:
it divides this power of two but does not divide its half. The literal
factorization in (2) proves that the defect is nonzero modulo $p$.
Therefore every Fermat prime, including any additional one, is non-Wieferich.

This is a special case of a result explicitly recorded by Crandall, Dilcher
and Pomerance: if the binary ones of a prime are equally spaced, it is not
Wieferich. For completeness, if

$$p=1+2^h+\cdots+2^{(t-1)h}\quad(h\ge1,\ t\ge2)$$

is prime, then $2^{ht}-1=p(2^h-1)$. A proper divisor of $ht$ is at most
$ht/2$, whereas $p>2^{h(t-1)}\ge2^{ht/2}$; thus no such divisor can be
the order. Hence $d=ht$, and $0<2^h-1<p$ gives $s=1$.
This includes both Fermat and Mersenne prime shapes. It is a proved
exclusion, not a reason to identify the two prime families.

For a Fermat prime, the Fermat quotient has the explicit residue

$$
q_p(2)=\frac{2^{p-1}-1}{p}
\equiv\frac{p-1}{d}(p-2)
\equiv\frac{2}{d}=2^{-m}\pmod p.                \tag{3}
$$

The inverse powers in (3) are residues modulo $p$, not real fractions.
The number of doubling cycles at **every** prime-power level is
$2^{2^m-m-1}$. Only $m=0,1$, hence $p=3,5$, give one cycle.

| $m$ | verified prime $p$ | order modulo $p$ | cycles at every level | $q_p(2)\pmod p$ |
|---:|---:|---:|---:|---:|
|0|3|2|1|1|
|1|5|4|1|3|
|2|17|8|2|13|
|3|257|16|16|225|
|4|65537|32|2048|61441|

The supplied five numbers are exactly these five verified entries. The
next Fermat number is $2^{32}+1=641\cdot6700417$. Neither that composite
nor this table proves that there are no further Fermat primes.
In fact the binomial congruence (3), interpreted as the indicated integer
quotient, holds for every Fermat **number**, including that composite:
$d\mid p-1$ and $2^d=1+p(p-2)$ suffice. Thus satisfying this congruence
cannot certify primality. The prime-dependent APIs still reject composites.

The companion [pointed tournament-path packet](tournament_path_defects_20261014e.md)
constructs an action on a particular set of path presentations whose free
orbit count is $2q_p(2)$. Where that companion's hypotheses apply, (3) reads
this **actual** count modulo $p$ as $2^{1-m}$, never zero for a Fermat prime.
It is not an action by automorphisms of the underlying tournament.

## 4. The phase compiler and its source-owned inverse edges

Suppose $2^a\equiv z\pmod{p^k}$, with $0\le a<L_k$ and $z$ a canonical
unit residue. A requested next target $Z\equiv z\pmod{p^k}$ has:

- **one permitted target digit** if $k<s$: it must equal
  $2^a\pmod{p^{k+1}}$. The exponent period remains $L_k$; every other
  target digit is impossible.
- **a unique exponent-digit lift for every target digit** if $k\ge s$.
  Put $v=2^a\pmod{p^{k+1}}$ and
  $c=(2^{L_k}-1)/p^k\pmod p\ne0$. Then the new exponent is

$$
a+jL_k,\qquad
j\equiv\frac{Z-v}{p^k}(vc)^{-1}\pmod p,\quad0\le j<p.       \tag{4}
$$

The first assertion follows from the unchanged order. For the second,
$2^{a+jL_k}\equiv v(1+jcp^k)\pmod{p^{k+1}}$ proves both existence and
uniqueness. The implementation finds the initial phase by an exact search
through $d$ residues, then applies (4). Thus it is complete, with no
discrete-log search over the full prime-power group. Failure at either the
initial orbit or a delayed digit is a certificate of nonmembership.

Here is a literal generalized-Collatz consumer. Fix an odd positive endpoint
$q$ with $p\nmid q$, a sign $b\in\{-1,1\}$ and a prescribed source
residue $c_0\pmod{p^{k-1}}$. Then

$$
n=\frac{q2^a-b}{p},\qquad
n\equiv c_0\pmod{p^{k-1}}
\quad\Longleftrightarrow\quad
2^a\equiv(b+pc_0)q^{-1}\pmod{p^k}.               \tag{5}
$$

For $a\ge1$, integral $n$ is odd. The explicit cut $q2^a>p+b$ gives
$n>1$, and then $pn+b=q2^a$ certifies the actual odd edge with valuation
exactly $a$ and endpoint exactly $q$. The cut excludes a padded edge from
1. For $p=3,b=1$ this is an ordinary Collatz inverse edge; for other
parameters it is the specified different map, not the classical map.
Nothing supplies a later terminal proof for $q$.

`compile_port_guard` returns the marked prime, depth, endpoint, sign,
source class and exponent phase. `make_port` and `audit_port` retain and
recompute the exact source, not merely a congruence witness. In particular:

- At $p=17,q=3$, neither sign admits even a first-level phase: $3$ is not
  in $\langle2\rangle$, which already contains $-1$.
- At $p=3,q=b=1$, the class $n\equiv1\pmod9$ requires
  $a\equiv2\pmod{18}$. The formal $a=2$ gives $n=1$ and is rejected;
  $a=20$ gives the valid source $349525$.
- At $p=1093$ or $3511$, $q=b=1$, the source class $n\equiv1\pmod p$
  is impossible, although $n\equiv0\pmod p$ is possible.

### The retained source address is a scaled isometry

Fix a reachable endpoint/sign and one sufficiently large permitted exponent
$a_*$ so that the first source exceeds 1. All native exponents in this
phase at or above $a_*$ have the form $a_*+dt$, $t\ge0$. Write

$$n_t=\frac{q2^{a_*+dt}-b}{p}.$$

For distinct nonnegative integers $t,u$, the same valuation proof gives

$$v_p(n_t-n_u)=s-1+v_p(t-u).                     \tag{6}$$

Thus all these sources occupy **one** class modulo $p^{s-1}$, but the
normalized coordinate $(n_t-n_0)/p^{s-1}$ is an exact isometry of finite
residue addresses: modulo every $p^r$ it is a bijection of the $p^r$
parameter classes. The valuation proves injectivity and well-definedness;
equal finite cardinalities give surjectivity. No assertion that this
geometric sequence contains every positive integer is made.

For $q=b=1$, every positive native exponent is $a=dt$, and

$$v_p\!\left(\frac{2^{dt}-1}{p}\right)=s-1+v_p(t).            \tag{7}$$

Wieferich therefore means every endpoint-1 inverse source is divisible by
$p$, with deeper strata prescribed by $t$. At a Fermat prime $s=1$, every
finite source residue is reachable, despite the possible initial
endpoint-orbit obstruction. Its first positive source is $p-2$;
at $p=3$ this is 1 and the strict port consumer excludes it.

Equation (6) connects the order defect to the repository's source-parameter
isometries. The preserved object is an actual native edge with its residue
and exponent address. The lost data in a bare cycle count are which endpoint
and which source class. Restoring those coordinates repairs the local
compiler; it does not prove a global trajectory reaches 1.

## 5. What the Platonic comparison does and does not preserve

The Schläfli inequality $(r-2)(t-2)<4$, $r,t\ge3$, has exactly the five
solutions $(3,3),(3,4),(4,3),(3,5),(5,3)$. Duality pairs the cube with
the octahedron and the dodecahedron with the icosahedron. Consequently
there are three abstract rotation-group types: $A_4,S_4,A_5$, of orders
12,24,60. These inherited classical identifications are checked here only
through the finite permutation groups and their element-order lists.

The earlier group connection $\operatorname{PSL}_2(3)\cong A_4$ and
$\operatorname{PSL}_2(5)\cong A_5$ does not continue along Fermat primes.
For every prime $p$, the unipotent matrix
$\left(\begin{smallmatrix}1&1\\0&1\end{smallmatrix}\right)$ has order
$p$ in $\operatorname{PSL}_2(p)$. At any Fermat prime $p\ge17$, none of
12,24,60 is divisible by $p$. Lagrange's theorem therefore forbids an
injective group map into a Platonic rotation group. This is an exact
obstruction to a group-preserving extension, not a prohibition on arbitrary
maps or on unrelated geometric constructions.

## 6. Sources, verification and next obligation

**CITED primary sources, not new results here.** Crandall, Dilcher and
Pomerance, [A Search for Wieferich and Wilson Primes](https://math.dartmouth.edu/~carlp/PDF/paper111.pdf),
*Math. Comp.* 66 (1997), 433–449: definition and the two examples at p433;
equally spaced binary-one exclusion, including Fermat and Mersenne primes,
at p434, section 1. The elementary proof is reproduced above. Dorais and
Klyve's [2011 authored article record](https://cs.uwaterloo.ca/journals/JIS/VOL14/Klyve/klyve3.html)
reports a search with no additional Wieferich prime below $6.7\cdot10^{15}$.
That historical bounded result is not an assertion that no others exist.
Euler's [translated original paper](https://arxiv.org/abs/math/0501118)
records the factor 641 of the first composite Fermat number. Our small
integer factorization and primality checks are independent.

Run:

```text
python -B -X utf8 04-computation/experiments/fermat_wieferich_order_20261014e.py
python -B -O -X utf8 04-computation/experiments/fermat_wieferich_order_20261014e.py
```

Declared controls: all odd primes below 200, plus 257,1093,3511,65537,
at levels 1–4 by independent minimal-order certificates; direct cycle
walks for primes 3 through 19 at levels 1–3; complete target phases for
$p=3,5,7,11,13$ through level 3; every delayed and next live digit for
1093 and 3511; source-address valuation and normalized digit-bijection
controls; all five displayed Fermat primes; prime binary repunits with
spacing 1–8, length 2–12 and highest bit at most 20; 1194 literal native
port replays in the declared small source-class universe; exact malformed
type, source, phase and ROOT-boundary hostiles; and the three permutation
group element-order lists. No floating arithmetic or external prime table.

Normal and optimized runs agree on **29,435 checks** and the saved LF JSON,
SHA256 `f3d9535e2545b0ba421f239cfeb7233b982a76a8342e9fe4f732087450793201`.
This is the finite verification scope; the all-parameter claims above are
supported by their proofs, not extrapolated from these runs.

The next meaningful question is source-specific: does a retained local
phase combine with an authenticated terminal proof and a decreasing
source-owned budget? The compiler answers the phase question completely.
Neither non-Wieferich lifting nor an orbit-count congruence answers the
terminal or budget question. A proposed bridge must pass all three.
