# A recursive common-future kernel with an exact guard and a rooted exit

**Status: PROVED** for the identities, guard equivalences, substitution, and
constructed completed families below. **FINITE-EXACT** for the stated script
universes. This does not prove that arbitrary positive odd inputs enter this
kernel, nor that an arbitrary exit dependency reaches 1.

The companion experiment is
[collatz_recursive_dependency_kernel_20261004.py](../../04-computation/experiments/collatz_recursive_dependency_kernel_20261004.py);
its deterministic output is
[collatz_recursive_dependency_kernel_20261004.out](collatz_recursive_dependency_kernel_20261004.out).

## 1. Inheritance and the useful four-vertex object

The starting family is the later-checkpoint quarter-child switch in
[collatz_checkpoint_reroute_20261004.md](collatz_checkpoint_reroute_20261004.md).
For the accelerated odd Collatz map
\[
 U(n)=\frac{3n+1}{2^{v_2(3n+1)}},
\]
write
\[
 v=(1,2,1,1,1,2),\qquad
 F_v(n)=\frac{729n+925}{256}.
\]
On precisely the binary source cylinder
\[
 n=155+2048t,\quad t\ge0,
\]
the six valuations \(v\) are exact and
\[
 x=F_v(n)=445+5832t=4h+1,\qquad
 h=111+1458t=\frac{729n+669}{1024}<n.
\]
Here \(x\equiv5\pmod8\). If \(a=v_2(3x+1)\ge3\), then
\[
 U(x)=U(h),\qquad v_2(3h+1)=a-2.
\]
Thus four vertices \(n,x,h,U(x)\) carry an actual six-edge macro, two actual
odd edges to the common future, and a proof dependency \(n\Rightarrow h\).
The dependency is useful once a rooted certificate for \(h\) is supplied.
It is not itself an actual Collatz edge. There is no intrinsic observable for
all six unordered pairs, so these four vertices are a typed proof diagram,
not a tournament.

The earlier first-common-future cutting and quarter-child switch supply this
kernel. The new contribution here is its exact repeat guard, its compilation
into a repeated actual word, and a reusable certificate record. The general
idea of valuation fuel and the obstruction to replacing it by a finite
unmarked state machine are inherited from
[collatz_guarded_pumping_memory_20261004.md](collatz_guarded_pumping_memory_20261004.md).
The ROOT codec and grounded graph interface are inherited from
[inverse_ray_ternary_addresses_20261004.md](inverse_ray_ternary_addresses_20261004.md)
and [fair_frontier_extension_20261004.md](fair_frontier_extension_20261004.md).

## 2. One finite kernel, an unbounded arithmetic register

Define the dependency map and its rational fixed point by
\[
 H(n)=\frac{729n+669}{1024},\qquad
 \alpha=\frac{669}{295},\qquad
 \Delta(n)=295n-669.
\]
Since \(295\cdot155-669=22\cdot2048\), and 295 is odd,
\[
 n\equiv155\pmod{2048}
 \quad\Longleftrightarrow\quad
 v_2(\Delta(n))\ge11.
\]
No integer has \(\Delta(n)=0\), so its valuation always exists. Moreover
\[
 \Delta(H(n))=\frac{729}{1024}\Delta(n).
\]

**Exact repetition theorem.** For every positive odd integer \(n\) and every
integer \(m\ge1\), the kernel can be applied successively \(m\) times if and
only if
\[
 \boxed{v_2(295n-669)\ge10m+1.}
\]
Each application is on its original binary cylinder, produces a positive odd
integer, and strictly decreases that current source. Its terminal value is
\[
 c=H^m(n)=
 \frac{669+729^m(295n-669)/1024^m}{295}.
\]

For necessity, each application consumes exactly ten units of binary
valuation and requires eleven at its input. For sufficiency, the source
congruence gives \(n\ge155\) and \(H(n)\ge111\); the fuel identity inductively
keeps every required intermediate input on the same cylinder. The strict
inequality follows from
\[
 n-H(n)=\frac{295n-669}{1024}>0.
\]
Consequently the exact available repeat count is
\[
 \left\lfloor\frac{v_2(295n-669)-1}{10}\right\rfloor.
\]
The cylinders for at least \(m\) repetitions are nested single residue classes
modulo \(2^{10m+1}\), with density \(2^{-10m}\) relative to the odd integers.
The final oddness bit is essential: valuation exactly \(10m\) does not suffice.

This dependency chain never reaches ROOT by itself: every kernel output is at
least 111. At the end of its available fuel it needs another rule or a supplied
rooted certificate. Repetition does not ground its own proof.

For height and storage accounting,
\[
 c<n<\left(\frac{1024}{729}\right)^m c,\qquad
 \operatorname{bits}(c)\le\operatorname{bits}(n)
 \le\operatorname{bits}(c)+\lceil m/2\rceil.
\]
The first inequality follows by writing
\(n=\alpha+(1024/729)^m(c-\alpha)\); the bit bound uses
\(1024/729<\sqrt2\). These bounds retain the actual integer sizes.

## 3. The actual repeated word and the conjugacy square

Set
\[
 W=(3,2,1,1,1,2),\qquad
 F_W(x)=\frac{729x+2971}{1024}.
\]
Exact affine composition gives
\[
 \boxed{F_W\circ F_v=F_v\circ H.}
\]
Thus the dependency iteration becomes an actual word iteration after the
marked head \(v\). The two fixed points differ:
\[
 \alpha=669/295,\qquad F_v(\alpha)=2971/295.
\]
The source guard and the head are part of this conjugacy; forgetting them
does not create a legal integer trajectory.

Suppose a supplied terminal first-hit word for \(c=H^m(n)>1\) is
\((b,\tau)\). One quarter-child splice replaces it by
\((v,b+2,\tau)\). Repeating the splice gives the exact first-hit word
\[
 \boxed{v\,W^{m-1}\,(b+2)\,\tau.}
\]
This follows inductively because replacing the leading 1 of the previously
substituted \(v\) by 3 produces \(W\). Each splice adds six odd edges, ten
halvings, and sixteen ordinary Collatz edges.

There is no early-root padding. Every nonempty prefix of \(v\) has coefficient
greater than 1 and positive carry, so its source states exceed the current
kernel input, which is at least 155. The last splice edge reaches exactly
the supplied child's first successor; the supplied first-hit suffix then
finishes the route. Induction preserves this property at every splice.

**Cheap type hostile.** The dependency \(H\) has diagonal data
\((3^6,2^{10})\), identical to \(W\), but carry 669 instead of 2971.
It is not any positive valuation word. The ordered carry decoder from
[geometry_collatz_drift_carriers_20261004.md](geometry_collatz_drift_carriers_20261004.md)
rejects it. Concretely, \(H(155)=111\), while \(U^6(155)=445\).
Identical multiplier, trace, or odd-edge/cost counts do not repair this type
error.

## 4. An explicitly completed family at every repeat depth

For each \(m\ge1\), choose a positive integer \(s\ge2\) satisfying
\[
 4^s\equiv2302\cdot295^{-1}\pmod{3^{6m+1}},
 \qquad b=2s,\qquad c=\frac{2^b-1}{3}.
\]
The right side is a principal unit modulo 3. The order of 4 modulo
\(3^d\) is \(3^{d-1}\): the binomial identity inductively gives
\(v_3(4^{3^j}-1)=j+1\), and multiplication by a ternary unit preserves that
valuation. Hence the exponent \(s\) is one class modulo \(3^{6m}\).
Taking further positive periods gives infinitely many choices at every \(m\).
The script finds this class by lifting one ternary digit at a time.

Now
\[
 295c-669=\frac{295\,2^b-2302}{3}
\]
is divisible by \(729^m\). Define
\[
 \boxed{n=\frac{669+1024^m(295c-669)/729^m}{295}.}
\]
This is an integer: \(1024\equiv729\pmod{295}\), and both are units modulo
295, so its numerator vanishes modulo 295. It is positive and \(n>c\),
because \(c\ge5>\alpha\). It is odd because its numerator is odd. Also
\[
 \Delta(n)=\frac{1024^m}{729^m}\Delta(c),\qquad
 v_2(\Delta(n))=10m+1.
\]
The equality holds because \(b\ge4\) gives
\(v_2((295\,2^b-2302)/3)=1\). Thus every constructed source has exactly
\(m\) available kernel repetitions, followed by its explicit one-edge
terminal route; its final exit is not hidden in further repetitions.
The repetition theorem therefore applies and gives \(H^m(n)=c\).
The terminal \(c\) has the one-edge first-hit word \((b)\), so the source
has the completed first-hit word
\[
 vW^{m-1}(b+2).
\]
Its odd rank is \(6m+1\), its halving cost is \(10m+b\), and its ordinary
rank is \(16m+b+1\). This proves an infinite constructed completed family at
every depth, without an orbit-search premise. It is an explicit use of
inherited inverse-fibre/valuation lifting, not a claim of new universal
coverage.

Two literal controls are:

| Repeats \(m\) | Terminal exponent \(b\) | Source bits | Terminal bits | Odd / ordinary rank |
|---:|---:|---:|---:|---:|
| 1 | 888 | 887 | 887 | 7 / 905 |
| 2 | 786750 | 786750 | 786749 | 13 / 786783 |

Higher controls use exact modular readers without expanding those source
integers. The reader enlarges its requested modulus by the full denominator
before division; thus it works even when that modulus shares factors with
3 or 295. An independent reader reverses each actual word letter. Both also
agree with the inherited ROOT codec at powers of 2 and 3.

### A positive power-of-three intersection at every depth

**PROVED guard coverage, with no terminal-home conclusion.** The same
repeat cylinders also contain infinitely many powers of three at every
depth. This connects the construction to the subgroup observation in
[HYP-9175, powers-of-three-shadow-only-cycle-points-1-or-3-mod-8-resisting-share-constant](../hypotheses/HYP-9175-powers-of-three-shadow-only-cycle-points-1-or-3-mod-8-resisting-share-constant.md),
Proposition E (incoming commit 6961974fe). That file has overall OPEN status;
only its proved modular fact motivates this connection, and the short proof
below is independent of its conjectural coverage statements.

For \(k\ge3\), the subgroup generated by 3 modulo \(2^k\) is exactly the
units congruent to 1 or 3 modulo 8. Indeed every power has that form, and
\[
 v_2(3^{2^j}-1)=j+2\quad(j\ge1)
\]
follows from \(3^2-1=8\) and successive squaring. Odd exponents have
\(v_2(3^e-1)=1\); factoring an arbitrary exponent into its odd part and a
power of 2 then shows the order is exactly \(2^{k-2}\). This equals the
number of residues in the two displayed classes, proving equality of the
subgroup with that set.

Since \(\alpha=669/295\equiv3\pmod8\), there is a unique exponent
\[
 e_m\pmod{2^{10m-1}},\qquad
 3^{e_m}\equiv\alpha\pmod{2^{10m+1}}.
\]
It is odd. Consequently every source
\[
 n=3^{\,e_m+t\,2^{10m-1}},\qquad t\ge0,
\]
admits at least \(m\) decreasing kernel replacements. The exponent classes
are nested, with \(e_{m+1}\equiv e_m\pmod{2^{10m-1}}\); each depth increment
retains ten more binary exponent digits. At \(m=1\) the class is the inherited
\(e_1=483\pmod{512}\). The binary-lifting API computes higher refinements
without constructing the source integers.

This is a different parameter family from the explicitly completed plans
above: its repeat guard is paid, but no rooted proof for its terminal
dependency has been supplied. It neither certifies all powers of three nor
asserts that these power-of-three sources also lie in the completed family.

## 5. What the compression preserves and costs

The new RepeatCertificate(source, repetitions, terminal, terminal_word)
stores a supplied source, one repeat count, its exact terminal identity, and
the supplied terminal first-hit proof. Its verifier checks the sharp fuel
guard, the closed-form endpoint, and the terminal first-hit word. It does not
observe or replay the \(6m\) added odd edges. The optional word/ROOT-codec
exporters replay the expanded source route and preserve identity and ranks.

The macro changes neither the grounded-closure requirement nor the existing
controller. It supplies a compact checked proof which that controller can
consume. A missing terminal proof remains missing. The standalone test helper
that observes terminal routes is outside this verifier; its observations
are counted separately.

Any positive repeat count up to the available fuel is legal. A controller
should retain earlier certified exits as alternatives: the largest available
count is not a theorem that the final dependency is the best place to search.
The record takes a chosen count and checks its supplied terminal identity.

The expanded word adds \(6m\) letters; the shared templates \(v,W\) and a binary
repeat count describe those letters without duplicating them. Binary
composition builds \(F_W^{m-1}\) with \(O(\log m)\) matrix products and shared
squares. This is a structural compression statement, not a logarithmic
bit-complexity or runtime claim. Summaries, source integers, the terminal
exponent, and terminal proof still require their actual bits. Expanded
certificate export deliberately pays the expansion cost.

The familiar reciprocal/trace diagram is retained rather than rediscovered:
for \(P=729,Q=1024\), let
\[
 \lambda=P/Q,\qquad J=\lambda+\lambda^{-1},\qquad
 X=(P+Q)/\sqrt{PQ}.
\]
Then \(X^2=J+2\), and doubling a word gives
\(\lambda\mapsto\lambda^2\) and \(J\mapsto J^2-2\).
This is the mechanism already proved in
[quadratic_escape_rank_atlas_20261004.md](quadratic_escape_rank_atlas_20261004.md),
section 6. It is suitable for repeated block summaries. The carry and source
guard, however, must accompany that trace: \(H\) and \(W\) are an explicit
same-trace pair with different mathematical roles.

## 6. Reproduction and bounded controls

Run from the repository root:

    python -B 04-computation/experiments/collatz_recursive_dependency_kernel_20261004.py
    python -B -O 04-computation/experiments/collatz_recursive_dependency_kernel_20261004.py

The finite universe is:

- 128 prefix cases: \(m=1,\ldots,32\), with lifts \(0,1,7,31\) of the exact
  binary guard; closed-form endpoints, every local guard/decrease, the actual
  substituted prefix, fuel consumption and integer bit bounds are checked.
- Four larger prefix-only controls \(m=64,128,256,512\); these do not claim a
  terminal route. Linear and binary word composition agree for \(m=1,\ldots,32\).
- Twelve completed supplied-terminal cases, one canonical binary representative
  at each \(m=1,\ldots,12\). Finding their terminal proofs costs 1,622 explicitly
  counted odd observations. Their added 468 odd letters are subsequently
  represented by twelve counts and shared templates. The 1,622 terminal
  letters remain.
- The two literal completed-family controls in the table.
- Forty-eight symbolic completed plans: \(m=1,\ldots,16\), exponent-family
  lifts \(0,1,7\), five moduli each, giving 240 independent residue comparisons
  plus the inherited binary/ternary codec checks.
- Eight power-of-three exponent refinements \(m=1,\ldots,8\), each tested
  at five exponent lifts \(0,1,2,7,31\), giving 40 modular guard checks.
  Exact clock orders and retention of previous exponent digits are checked.
  None of these controls expands its source or supplies a terminal proof.
- Explicit malformed-source/count/terminal/phase controls; ROOT-padded words,
  corrupted terminal identity, insufficient fuel, and an exceeded expansion
  cap are rejected. Another 128 sharp controls omit exactly the final
  oddness bit. The non-word dependency \(H\) is rejected by the actual-word
  decoder.

All checks use explicit exceptions and execute under optimization. There are
no import-time experiments, hidden file writes, new ROOT seeds, or modifications
to the existing controller. The next useful operation is to let a controller
retain this typed repeat node alongside its literal graph, expanding only at
the boundary of an interface that actually requires the full word.

**Local API correction.** An audit found that the initial carrier decoder
accepted the float-alias tuple (729.0, 1024, 2971), because Python arithmetic
and equality treated its first coordinate as the corresponding integer.
The decoder now requires a tuple of three exact positive integers before
decoding, and the hostile input is retained as a regression. The mathematical
carrier theorem and the guarded certificate verifier were unaffected.
