# Nineteen residue rays with completed Collatz routes of sharp depth four

2026-10-04. **PROVED:** every residue modulo every power of 19 has infinitely
many positive odd representatives with a first-hit route to 1 of at most four
accelerated steps. Four is the smallest possible uniform bound. The explicit
bank below has natural density zero. **FINITE-EXACT:** the bounded controls
listed in section 6. **OPEN:** convergence of arbitrary supplied integers.
No historical novelty is claimed for inverse Collatz constructions or
prime-adic lifting.

The output is a new, explicitly defined integer and its certificate.
Prescribing the residue of an existing integer does not certify that integer.

## 1. Inheritance and the question being answered

The closest proved mechanism is the first-hit parent-pointer codec in
[inverse rays and ternary addresses](inverse_ray_ternary_addresses_20261004.md),
especially its mod-9 inverse guards and root-self-return exclusion.
The earlier [word-completion theorem](arithmetic_braids2_20260917_inverse_completion.md)
already constructs complete routes in arbitrary binary/ternary classes;
its section 3 explicitly separates residue support, source identity, and
ordinary size. The current [boundary compiler](collatz_boundary_compiler_20261004.md)
instead retains a supplied source and a magnitude guard. Neither direction
can be silently substituted for the other.

Canonical hostile: an inverse source divisible by 3 has a completed forward
route but cannot be used as the target of another inverse step. Corrected
near miss: a residue-dense bank need not cover the integers. Least-used
sidecar: the selected inverse exponent, together with the first-hit parent
certificate. The live board is root / inverse leaf / residue disk / digit
parameter / exact odd rank / ordinary height.

The [golden phase tower](difference_families_20261004.md), equations 14–16,
has 19 child cycles per parent above denominator 76. The
[cubic depth tower](cyclotomic_depth_towers_20261004.md), section 1, has a
specified multiplicative orbit missing an index-19 factor in dimension 18.
Our route construction uses the exact base-2 datum

\[
 \operatorname{ord}_{19}(2)=18,\qquad
 2^{18}-1=19\cdot13797,\qquad13797\equiv3\pmod {19}.
\]

This is the same factorization behind that finite-field index, whereas the
golden tower uses its separate matrix identity \(G^{18}-I=76G^9\).
The objects and maps differ: our children are digits of a parameter for an
actual completed inverse route. No golden-phase or finite-field orbit is
declared to be that route.

## 2. Eighteen easy residue disks and one exceptional disk

Write \(U(n)=(3n+1)/2^{v_2(3n+1)}\) on positive odd integers. For \(t\ge0\)
use these families:

1. For each even \(k\in\{4,6,\ldots,20\}\),
   \[
     n_t=\frac{2^{k+18t}-1}{3},\qquad U(n_t)=1.
   \]
2. For each odd \(k\in\{1,3,\ldots,17\}\),
   \[
     n_t=\frac{5\cdot2^{k+18t}-1}{3},\qquad
     n_t\longrightarrow5\longrightarrow1.
   \]

All divisions are integral, all sources are odd and greater than 1, and the
displayed exponents are exact because the targets are odd. Using \(k=20\)
instead of \(k=2\) removes the root self-return without losing a residue.
Modulo 19, the values \(3n_t+1\) in the first group are precisely the nine
nonzero squares. In the second group they are \(10\) times those squares;
10 is a nonsquare, so this gives the other nine units. Thus the eighteen
families occupy distinct source residues and omit exactly
\[
 r=6=-1/3\pmod {19}.
\]

The exceptional family changes an internal exponent:

\[
 \boxed{\quad
 n_t=\frac{53\cdot2^{8+18t}-5}{9},
 \qquad \text{valuation word }(1,\,7+18t,\,5,\,4).
 \quad}                                                     \tag{1}
\]

To verify its guards, put \(u_t=(53\cdot2^{7+18t}-1)/3\).
Since \(2^{18}\equiv1\pmod9\), \(u_t\equiv2\pmod3\).
Consequently \(n_t=(2u_t-1)/3\) is an odd integer. The exact route is
\[
 n_t\longrightarrow u_t\longrightarrow53\longrightarrow5\longrightarrow1.
\]
All nodes before 1 exceed 1. The first instance is
\[
 1507\longrightarrow2261\longrightarrow53\longrightarrow5\longrightarrow1.
\]
Modulo 19, \(u_t=0\) and \(n_t=6\). In particular, (1) repairs the missing
residue without trying to extend a ternary leaf.

## 3. Why four odd steps are necessary

**PROVED obstruction.** No positive odd integer congruent to 6 modulo 19
has first-hit odd rank at most three.

If \(n\equiv6\pmod {19}\), its odd target \(u=U(n)\) is divisible by 19.
We show every completed source of odd rank at most two that is divisible
by 19 is also divisible by 3. Such a source cannot be an inverse target:
the identity \(3n+1=2^k u\) would fail modulo 3.

At rank one, \(u=(2^a-1)/3\), with \(a\) even. Divisibility by 19 implies
\(18\mid a\), hence \(u\equiv0\pmod3\).

At rank two write
\[
 v=(2^a-1)/3,\qquad u=(2^b v-1)/3.
\]
The intermediate target \(v\) must be a ternary unit, so
\(a\bmod18\in\{2,4,8,10,14,16\}\).
Both \(v\bmod9\) and \(v\bmod19\) are determined by this residue:
the first uses \(2^a\bmod27\), whose period is 18.
For each of the six cases, \(19\mid u\) uniquely fixes \(b\bmod18\).
The complete table is:

| \(a\bmod18\) | \(b\bmod18\), with 18 for 0 | Integral second inverse step? | \(u\bmod3\) |
|---:|---:|:---:|---:|
|2|18|yes|0|
|4|2|no|—|
|8|10|yes|0|
|10|9|yes|0|
|14|15|no|—|
|16|11|yes|0|

This is a finite proof of all exponent heights, not a bounded-height
experiment. The class \(a=2\) is retained because \(a=20,38,\ldots\) gives
nonroot predecessors. Each integral case is a ternary leaf. The rank-zero
target 1 is not divisible by 19. This exhausts the possibilities and proves
the obstruction. Equation (1) attains rank four. Therefore four is the
sharp uniform residue-coverage bound at every precision \(19^a\), \(a\ge1\).

## 4. A digit compiler that never expands the source

Each ray is \(n_t=(A2^{18t}-B)/D\), with \(19\nmid AD\).
For distinct nonnegative parameters \(s,t\), elementary binomial lifting
at the odd prime 19 gives
\[
 v_{19}(n_t-n_s)=1+v_{19}(t-s).                         \tag{2}
\]
Indeed \(v_{19}(2^{18h}-1)=1+v_{19}(h)\), initially from the displayed
factorization and then by taking nineteenth powers and prime-to-19 powers.
The other factors in the difference are 19-adic units.

More precisely, for \(a\ge1\) and \(d\in\{0,\ldots,18\}\),
\[
 n_{t+d19^{a-1}}\equiv n_t+c_r d19^a\pmod {19^{a+1}},
 \quad
 c_r=\begin{cases}3r+1&r\ne6,\\7&r=6.\end{cases}          \tag{3}
\]
The coefficient is independent of \(t\) and is nonzero modulo 19.
For ordinary outer rays it is
\((A/3)\cdot3=A=3r+1\pmod {19}\).
For (1) it is \((256\cdot53/9)\cdot3=7\pmod {19}\).

Given a canonical address \(z\bmod19^a\), choose its unique ray
\(r=z\bmod19\). Start \(t=0\). Once the first \(j\) source digits match, set
\[
 d=\frac{z-n_t}{19^j}c_r^{-1}\pmod {19},\qquad
 t\gets t+d19^{j-1}.
\]
In this calculation \(n_t\) need only be read modulo \(19^{j+1}\).
This retains all previous digits and fixes the next one. It returns the
unique \(0\le t_0<19^{a-1}\) with the prescribed address. All representatives
on that ray are exactly
\[
 t=t_0+h19^{a-1},\qquad h\ge0.                         \tag{4}
\]
They are distinct, positive, and carry the same fixed-length first-hit
grammar. This proves the all-height surjectivity and infinitude.

The script's compile_address(address, precision, lift=0) returns the ray,
parameter and inherited Certificate. Its modular reader checks the result
independently by walking the parent pointers. The source is never expanded.
The canonical parameter takes \(O(a)\) bits; the certificate has at most
four nodes. It uses \(a-1\) digit corrections, with modular powers at growing
precision, rather than a search through \(19^{a-1}\) parameters.

For the canonical choice, ordinary first-hit rank is at most
\(18\cdot19^{a-1}+5\). This exponential bound concerns the trajectory clock,
not the size of its symbolic description. Our precision-80 control has a
336-bit parameter and a source with more than \(10^{100}\) binary digits.

Equations (2)–(3) also extend each ray continuously from \(\mathbb Z_{19}\)
onto its residue disk \(r+19\mathbb Z_{19}\). A nonterminating 19-adic
parameter is not an ordinary nonnegative exponent. An infinite branch of
the address tree therefore does not produce a fixed positive integer with
an infinite or completed ordinary route; the actual certificates above
use finite nonnegative parameters.

## 5. Full residue support, zero density, and the selector boundary

For each fixed ray,
\[
 \#\{t\ge0:n_t\le X\}=\frac{\log_2 X}{18}+O(1).
\]
The nineteen rays are disjoint because their residues modulo 19 differ.
Their union therefore contains
\[
 \frac{19}{18}\log_2X+O(1)
\]
sources below \(X\), and has natural density zero. Inside a fixed address
at precision \(a\), equation (4) gives only
\(\log_2X/(18\cdot19^{a-1})+O(1)\) sources.

The map's source is (residue address, lift quotient). Its target is a
new positive odd integer with an explicit completed route. It preserves
the requested residue and exact route semantics. It does not preserve an
independently supplied integer or its magnitude. Retaining the parameter,
ray and parent certificate recovers the constructed source, but it cannot
repair a false equality to the supplied source.

For example, asking this bank for the residue of 19 modulo 19 returns
87381, which reaches 1 in one odd step. This says nothing by itself about
the integer 19. An adaptive source selector may use these rays as certified
targets only after proving its actual endpoint equals one of them. Agreement
at arbitrarily chosen finite modular precision is not that equality proof.

## 6. Reproduction and independent controls

Run from the repository root:

    python 04-computation/experiments/mod19_route_lifts_20261004.py
    python -O 04-computation/experiments/mod19_route_lifts_20261004.py

The [script](../../04-computation/experiments/mod19_route_lifts_20261004.py)
and [saved output](mod19_route_lifts_20261004.out) check:

- all 7,239 addresses modulo \(19,19^2,19^3\), including compatibility of
  parameter prefixes and the independent inherited certificate reader;
- all nineteen rays at parameters 0–24, with 475 literal forward replays,
  exact valuation words, first-hit root checks and codec round trips;
- all 5,700 pairwise differences in those parameter samples;
- all nineteen child digits in all nineteen rays at eight precision levels,
  giving 2,888 exact digit checks;
- three lift quotients at every base residue, the complete six-case
  obstruction table, one symbolic precision-80 case, and eleven
  type/domain/root/integrality hostiles.

The inverse-leaf obstruction and exceptional-ray derivative were also
independently checked by the sibling arithmetic lane. The all-height
conclusions follow from the proofs above, not extrapolation from these
finite universes.
