# A precise aliquot three-cycle frontier, with Fermat and tournament connections

**2026-10-19 update:** Hypotheses A3 and D3 below are now proved for all exponents in [the dyadic three-cycle exclusion](dyadic_aliquot_three_cycle_exclusion_20261019/DYADIC_THREE_CYCLE_EXCLUSION.md). The earlier status and finite census below are retained as the historical frontier. This does not exclude general aliquot three-cycles.

**PROVED scoped:** repeated dyadic cofactors cannot support a primitive aliquot three-cycle, and a resultant bounds the third exponent once the other two are fixed. **FINITE-EXACT:** no remaining candidate with second-smallest exponent at most 128, with the third exponent unrestricted. **OPEN:** the conjectured exclusion for all dyadic-times-prime triples, and the general existence of an aliquot three-cycle. No global Collatz or aliquot termination claim follows.

Artifacts: [checker](../../04-computation/experiments/aliquot_dyadic_frontier_20261014e.py), [output](aliquot_dyadic_frontier_20261014e.json), [general core compiler](aliquot_guarded_cycles_20261014e.md), [tournament path encodings](tournament_path_defects_20261014e.md), and [Fermat/Wieferich order lifts](fermat_wieferich_order_20261014e.md).

## 1. Inheritance and the live target

The closest mechanism is the multiplicative divisor-sum formula retained in [the aliquot/Lehmer packet](collatz_aliquot_lehmer_five_20260927.md). Its [later audit](collatz_sticky_20260927_audit.md) limits the claims about sampled drift and driver persistence: they do not prove orbit fate. The analogous current Collatz mechanism is the [source-owned affine word and clock classifier](mersenne_clocked_completion_20261014c.md). Both require an algebraic candidate followed by native guards.

Our anchor is a sharply delimited part of the unknown aliquot three-cycle problem. The niche is the actual path-encoding relation to Fermat quotients. The wildcard is the contrast between five regular solids and five known Fermat primes. The board is **factorization template / ordered affine word / integrality / native prime guard / orbit quotient / surviving parameter**.

The least-used coordinate is the cofactor left after removing one prime that appears exactly once. Fixing that coordinate makes the apparently unbounded search for prime values linear. The canonical hostile below is the apparent cycle `46,26,16`: its final edge works only after pretending that 1 is prime.

The incoming repair-algebra audit ([MISTAKE-606](../../01-canon/MISTAKES.md) and the repaired [THM-4640, repair algebra at minus one](../../01-canon/theorems/THM-4640-the-repair-algebra-at-minus-one-collisions-are-increasing-two-three-representations-ladders-and-gersonides-carries.md)) is relevant methodologically: uniform switch coverage needed the actual root switch and its corrected quantifiers. The incoming Hadwiger--Nelson [hinge/mono-set experiment](hn_type_m_hinge_rotations_and_forced_radii_20261014.md) likewise distinguishes a forced radius from enough cross-edge constraints. Here a rational cyclic solution is not enough; the prime and coprimality conditions are the retained constraints.

## 2. From abundance labels to an exact cyclic compiler

Write `s(n)=sigma(n)-n`. Deficient, perfect and abundant mean respectively `s(n)<n`, `s(n)=n` and `s(n)>n`. These are local directions in the functional graph, not three possible eventual orbit fates.

If `c>=2` and `p` is prime with `gcd(c,p)=1`, then

\[
s(cp)=(\sigma(c)-c)p+\sigma(c),\qquad
s(cp)-cp=p(\sigma(c)-2c)+\sigma(c).                 \tag{1}
\]

Thus a deficient core has an exact transition threshold `p=sigma(c)/(2c-sigma(c))`; perfect or abundant cores produce abundant prime extensions. For an ordered core word `(c_0,...,c_(l-1))`, requiring `n_i=c_i p_i` and `s(n_i)=n_(i+1)` gives

\[
p_{i+1}=\frac{(\sigma(c_i)-c_i)p_i+\sigma(c_i)}{c_{i+1}}.
\]

The [companion](aliquot_guarded_cycles_20261014e.md) composes these affine maps. It gives at most one rational candidate per fixed core word and checks the remaining arithmetic guards. This is an exact classifier with no cutoff on the prime variables, not a claim that a finite core bank covers every integer.

We specialize to `c_i=2^(a_i)`, `a_i>=1`, with the odd prime occurring to exponent one. Formula (1) becomes

\[
s(cp)=(c-1)p+2c-1.                                \tag{2}
\]

For a three-cycle with ordered cores `(c,d,e)`, let the prime variables be `(p,q,r)`. Direct elimination of the three equations gives

\[
p=\frac{N(e)}{D(e)},\qquad D(e)=Le+K,
\]

\[
\begin{aligned}
L&=c+d-1,&K&=(c-1)(d-1),\\
N(e)&=Ae^2+Be+C,&A&=2d,\\
B&=2cd-2c+2d^2-3d+1,\\
C&=-2cd+2c-2d^2+2d-1.
\end{aligned}                                                    \tag{3}
\]

All denominators here are positive. Once `p` is known, (2) determines `q` and `r`. An integer first candidate is necessary; primality and three distinct source integers must still be checked. The checker compares (3) with independent rational affine composition on the source variables.

**Integrality propagates in this dyadic class.** The unreduced common denominator of all three cyclically rotated candidates is

\[
D=cde-(c-1)(d-1)(e-1)=cd+ce+de-c-d-e+1,
\]

which is odd. The core compiler writes each candidate as `B_i/D` with `B_i` a positive integer. If one candidate is integral, its edge equation makes the next candidate times its dyadic core integral. Coprimality of that core with `D` forces the next candidate to be integral too. All three are then positive, and their edge equations force oddness. Thus one integrality test decides all three; this predicate is invariant under cyclic rotation. Primality and distinctness remain separate conditions. This simplification relies on the common odd denominator and need not hold for general cores.

## 3. All-exponent theorem: a repeated exponent is impossible

**Theorem.** A primitive aliquot three-cycle entirely of the form `2^a p`, with `a>=1` and `p` an odd prime, must have three pairwise distinct exponents of 2.

Rotate a hypothetical cycle with two equal cores into `(C,C,D)`, where `C,D` are powers of two at least 2. For this proof temporarily allow positive odd integers as the prime variables. The first equation implies

\[
p=Ck-1,\qquad q=(C-1)k+1,
\]

where `k` is a positive **even** integer, since `q` is odd. Eliminating the remaining variables gives

\[
k=\frac{4CD-3C+2D^2-3D+2}{C^2+2CD-2C-D+1}.        \tag{4}
\]

If `D<C`, subtracting twice the denominator from the numerator gives `(D-C)(2(D+C)-1)<0`. Hence `0<k<2`, impossible. If `D=C`, (4) gives `k=2`, and all three variables are `2C-1`: only a fixed point is possible.

Suppose `D>C`. The third variable is

\[
r=\frac{4C^3+2C^2D-5C^2-4CD+3C+2D-1}
        {C^2+2CD-2C-D+1}.                         \tag{5}
\]

Subtracting `(2C-1)` times its denominator gives `(C-D)(2C^2-1)<0`. Subtracting `(C-2)` times the denominator gives `3C^3-C^2+CD-2C+1>0`. Therefore

\[
C-2<r<2C-1.
\]

Because `C|D`, the numerator and denominator in (5) are respectively `-1` and `1` modulo `C`. An integral `r` must satisfy `r=-1 mod C`, leaving exactly `r=C-1`. Substitution forces

\[
D=\frac{C^2(3C-2)}{C-1}.
\]

Its integrality requires `C-1|1`, so `C=2`, `D=16`. This gives `(p,q,r)=(23,13,1)` and the formal source cycle `(46,26,16)`. The final variable is not prime: actually `s(16)=15`, not 46. Thus this last case is nonnative. This proves the theorem and classifies the sole nonconstant positive-odd-variable exception with repeated cores.

The failed implication is concrete: multiplicativity used `sigma(p)=p+1`; at `p=1`, the true value is `sigma(1)=1`. Retaining the factorization guard repairs the cyclic compiler.

## 4. The third exponent has an exact finite bound

Define

\[
R=A K^2-B K L+C L^2.
\]

The polynomial identity

\[
L^2N(e)-R=(Le+K)(ALe+BL-AK)                       \tag{6}
\]

shows that `D(e)|N(e)` implies `D(e)|R`.

**The retained resultant never vanishes for `c,d>=2`.** Put `u=K/L`, so `0<u<d-1`. Since

\[
B=(d-1)(2c+2d-1)>2d(d-1),\qquad C<0,
\]

we have `R/L^2=2du^2-Bu+C<0`. Consequently every integral first candidate obeys

\[
\boxed{e\le\left\lfloor\frac{-R-K}{L}\right\rfloor.}           \tag{7}
\]

No primality assumption or floating-point estimate enters this bound. It covers every larger third exponent by divisibility, not by extrapolating a search.

For three distinct exponents, place the largest cofactor last. Given the two smaller exponents `a<b`, both cyclic orientations are represented by ordered pairs `(c,d)=(2^a,2^b)` and `(2^b,2^a)`. We tested every `1<=a<b<=128` and every larger power of two allowed by (7). There were **3,156,331 oriented core words**. None had even an integral first prime variable. The largest exponent requiring an explicit test was 512; every still larger exponent is excluded by (7) for the same lower pair.

Combining this finite exact calculation with section 3 gives the following scoped conclusion:

> Any primitive aliquot three-cycle entirely of the form `2^a p` would have three distinct exponents, and its second-smallest exponent would be at least **129**.

This does not exclude general aliquot three-cycles, whose odd parts can contain several primes or repeated odd prime powers. It is a parameter-space exclusion, not an arbitrary-height census of all aliquot inputs.

## 5. Two precise hypotheses and their failure tests

**OPEN working hypothesis D3.** No primitive aliquot three-cycle lies entirely in the class `{2^a p : a>=1, p odd prime}`. The repeated-exponent part is proved above; the surviving case has three distinct exponents and second-smallest exponent at least 129. Novelty relative to the full literature is not asserted.

**OPEN stronger arithmetic hypothesis A3.** For pairwise distinct powers of two `c,d,e>=2`, the denominator `D(e)` in (3) never divides `N(e)`. This would exclude a cyclic solution even before the prime guards. The stated computation supports A3 in its exact universe; it is not a proof for all parameters.

The next decisive test is therefore a divisibility mechanism for A3 or an explicit integral exception. By the propagation lemma, an exception to A3 would already supply three positive odd integer variables. It would still have to pass primality and source-distinctness before it could challenge D3. This two-stage target prevents a failed stronger hypothesis from being mistaken for a new aliquot cycle. Raising the bound alone is useful only if it reveals structure in the remaining resultant divisors.

The general core compiler also identifies a different constructive possibility: assemble different cores whose individual direction labels fail to predict the whole cycle, then solve the coupled word exactly. This is the same operational lesson as the [two-prime repair of 561](carmichael_successor_blocks_20261014b.md). The actual maps and guards differ; no arithmetic identification is claimed.

## 6. The graph and five-object connections that survive exact maps

The [tournament companion](tournament_path_defects_20261014e.md) retains a useful actual family. Reverse only the extreme edge in a transitive `n`-vertex tournament. Its Hamiltonian-path count is `1+2^(n-2)`, with an explicit subset encoding. At orders `3,4,6,10,18` these counts are precisely `3,5,17,257,65537`. This explains the user's corrected Fermat list by a construction, not just by the fact that it has five entries.

For an odd prime `p`, cyclically rotating the `p` subset coordinates gives exactly `(2^p-2)/p=2q_p(2)` free orbits. Thus base-two Wieferich divisibility is precisely divisibility by `p` of this marked-path orbit count. This action is on encoded path presentations; the companion proves it is not generally a vertex automorphism. The inherited exclusions `H=7,21` concern a different, global Hamiltonian spectrum question. Its old claim that these were *exactly* the missing values has been repaired to the established exclusion, with completeness still open.

The [order-lift companion](fermat_wieferich_order_20261014e.md) explains the other side. A Fermat prime `p=2^(2^m)+1` has `ord_p(2)=2^(m+1)` and

\[
2^{\operatorname{ord}_p(2)}-1=p(p-2).
\]

Its prime-power lifting defect is exactly one, so it is not Wieferich. The known examples 1093 and 3511 instead have lifting defect two. Retaining that defect gives a source-residue compiler for generalized signed Collatz inverse ports. The Fermat assertion is classical and is attributed in the companion; the reusable artifact is the guarded lift.

There are five Platonic solids, but duality groups them into three rotation-group types `A4,S4,A5`. The companion retains the spherical curvature condition and constructible-polygon connection. It does not assert a bijection from the five known Fermat primes to the solids or classify future Fermat primes by their number.

Finally, every aliquot cycle satisfies `sigma(n_i)=n_i+n_(i+1)`. For odd length `l`,

\[
n_i=\frac12\sum_{j=0}^{l-1}(-1)^j\sigma(n_{i+j});
\]

for even length there is instead an alternating-sum compatibility condition and a free alternating coordinate. This is the elementary `I+P` distinction for a cyclic permutation matrix. It treats lengths 3 and 5 identically at the linear level, so it cannot explain the absence of a known three-cycle by itself. The factorization guard is where the arithmetic enters. The companion also recovers Pollack's stronger gcd restrictions for odd-length cycles and supplies a period-four control: the historically famous lengths 1,2,5,28 do not exhaust the known periods.

## 7. Reproduction and scope

```text
python -B 04-computation/experiments/aliquot_dyadic_frontier_20261014e.py
python -B -O 04-computation/experiments/aliquot_dyadic_frontier_20261014e.py
```

Both runs reproduce **3,178,849 explicit checks**, with LF-normalized JSON SHA256
`1c52c86b0f13c6ae70304b9ec9047be5ceb5b6bcf5c67cdd008317a265061f64`.
Independent affine and cyclic-integrality controls use all exponent triples in `1..8`; repeated-core controls use every pair in `1..64`; the main distinct-core universe and its unbounded tail proof are stated in section 4. Actual divisor enumeration verifies the ghost failure and fixed-point controls. Invalid cofactor types are rejected, and no optional assertion can disable a check.

Independent audits by the aliquot and order-lifting lanes both passed the repeated-core proof, the resultant sign and tail bound, the two-orientation coverage, and optimized-versus-saved reproduction. The three companion packets were also independently peer-audited; the cycle-length and exponent-tail wording corrections were applied before this checkpoint.
