# Multiplier transport, arithmetic cycle spectra, and moving Fermat certificates

Session: `collatz_fermat_reframes_20261029`. Started from `92514b16e2f6`.

**Status.** PROVED below: elementary cycle transport and cone separation (C1),
the forcing-divisor relation and exact recurrence (C2), and the one-letter
failure of a denominator-only carrier (C3). These refine classical cycle-word
formulas; no literature-first claim is made. FINITE-EXACT: the stated census
and hostile tests in [results.json](results.json), reproduced by [verify.py](verify.py).
[FERMAT.md](FERMAT.md) proves a restriction on fixed algebraic certificates and
gives a constructive projection onto the odd-order unit subgroup. Collatz
convergence, exclusion of all nontrivial integer Collatz cycles, and
compositeness of every Fermat number beyond the fifth known prime remain OPEN.
This is a research packet with explicit proofs and same-agent exact checks,
not independently refereed canon or a Lean formalization.

## Inheritance and live concept board

The closest proved mechanisms are the classical parity-word cycle formula and
the cycle product identity: THM-4634
[`collatz-parity-number-dictionary-irrationality-exponents-of-parity-vectors-and-jarnik-strata`](../../../01-canon/theorems/THM-4634-collatz-parity-number-dictionary-irrationality-exponents-of-parity-vectors-and-jarnik-strata.md)
and THM-4633
[`collatz-cycles-as-rational-approximations-of-log2-3-with-the-admissible-minimum-budget`](../../../01-canon/theorems/THM-4633-collatz-cycles-as-rational-approximations-of-log2-3-with-the-admissible-minimum-budget.md).
The old [weighted-mediant law, section 3.2](../collatz_pentagon_fixed_points_argument_styles_20260930.md)
already supplies concatenation and warns that splitting an integral cycle need
not preserve integrality. The rational cycle `19/101,79/101,169/101` is already
in [the September squarefree-guard calculation](../collatz_guards_20260921_squarefree.md).
It is inherited, not newly discovered here.

Canonical hostile examples: the positive `5x+1` cycle `13,33,83`; two different
cycles with the same denominator and clocks; the tournament pair with the
same 67 labelled cyclic triples but automorphism orders 81/27. Corrected near
misses: MISTAKE-649 (generalized parity conjugacy/ergodicity is an old barrier,
not a convergence proof), MISTAKE-605 (rational parity and integer-cycle
exclusion are distinct obligations), and MISTAKE-673 (no assumed
squarefreeness of Fermat numbers). Least-used relevant sidecar recovered here:
the *numerator/forcing* behind the older rational-cycle catalogue.

The portfolio is **Anchor:** integer-selective Collatz transport; **Niche:**
Möbius inversion over forcing divisors; **Wildcard:** cyclotomic rigidity of
Fermat certificate recipes. The live board has five objects: marked valuation
words, reduced denominators, the real positive cone, future attachments, and
odd-order unit certificates. The results below compare each new representation
against these five. We do not manufacture a tournament from arithmetic ties.

## C1. A bijection that preserves the clock and exposes the arithmetic loss

For odd `a >= 3`, use the odd-step map

`U_a(x) = (a*x+1)/2^v`, where `v=v_2(a*x+1)`.

Here `x` is rational with odd numerator and denominator. We exclude the
terminal input `a*x+1=0`; it never occurs in the cycles below. For a nonempty
word `w=(v_0,...,v_(k-1))` of positive integers, write

`s_j = v_0+...+v_(j-1)`, `H=s_k`,

`C_a(w) = sum_(j=0)^(k-1) a^(k-1-j)*2^(s_j)`,

`D_a(w)=2^H-a^k`, `x_a(w)=C_a(w)/D_a(w)`.

**Cycle-word fact (classical).** Composition of the indicated affine branches
is `(a^k*x+C_a)/2^H`. Its unique fixed point is `x_a(w)`. Cyclic rotations
give the other marked points, and their valuations are *exactly* the letters
of `w`. Indeed `C_a` and `D_a` are odd and

`a*C_a(w)+D_a(w)=2^(v_0)*C_a(rotate(w))`.

After division by the odd `D_a`, the next point has valuation zero, proving
exactness. A primitive word has exact odd period `k`; a repeated word does
not. Thus matching the same primitive word gives a bijection between marked
rational periodic orbits for any two odd multipliers. It also respects cyclic
re-marking. Its unrestricted two-adic parent is classical; see
[Bernstein–Lagarias, The 3x+1 Conjugacy Map](https://websites.umich.edu/~lagarias/doc/bernstein.pdf).

**Cone-separation theorem.** Let `b>a` be odd. If the word has a positive
`U_b` cycle, its `U_a` cycle is positive and

`min_i x_(a,i) < 1/(b-a)`.

Conversely, if the `U_a` cycle consists of positive integers, then its
same-word `U_b` cycle is negative.

*Proof.* Multiplication around a positive cycle gives the classical identity

`2^H = product_i (a+1/x_(a,i))`.

Positivity at multiplier `b` is exactly `2^H>b^k`. If every `a`-point were at
least `1/(b-a)`, every factor in the displayed product would be at most `b`,
a contradiction. If every `a`-point is a positive integer, then
`2^H <= (a+1)^k < b^k`, proving the converse sign claim. The weak equality
in the first bound requires all `a`-points to be 1. The strict inequalities
at multiplier `b` hold because odd `b>a` implies `b>=a+2`. No probabilistic
assumption is used. QED.

In particular **every positive `5x+1` cycle transports to a rational
`3x+1` cycle with some value below `1/2`**. This provides an actual arithmetic
separation from the hostile map. It does not delete that map's cycles by a
purported conjugacy of the positive integers.

| Valuation word | Positive `5x+1` cycle | Same-word rational `3x+1` cycle | Integer forcing after scaling |
|---|---|---|---|
| `(1,1,5)` | `13,33,83` | `19/101,79/101,169/101` | `3x+101`: `19,79,169` |
| `(1,3,3)` | `17,43,27` | `31/101,97/101,49/101` | `3x+101`: `31,97,49` |
| `(1,4)` | `1,3` | `5/23,19/23` | `3x+23`: `5,19` |

The positive integer Collatz fixed point 1, word `(2)`, instead goes to the
negative `5x+1` fixed point -1. Any hypothetical positive integer Collatz
cycle would likewise go to a negative cycle. Excluding positive `5x+1`
cycles therefore does **not** exclude unknown positive `3x+1` cycles.

The exact retained arithmetic quantity is

`d_a(w)=abs(2^H-a^k)/gcd(C_a(w),2^H-a^k)`.

This reduced denominator is rotation-invariant and coprime to `2a`. Positivity
and `d_a=1` characterize positive integer cycles. Scaling by `d_a` takes the
rational orbit to an integer orbit of `a*x+d_a`, with all clocks and phases
unchanged. Changing the forcing from `d_a` to 1 is precisely the missing
arithmetic operation; it cannot be omitted from the representation.

### Incoming Session 100: exclude all positive hostile cycle basins

While this packet was being checked, the independently audited
[Session 100 selective target](../session100_selective_collatz_targets_20261029/bijections/VALUATION_CONJUGACY_AND_INTEGER_SIDECAR.md)
arrived on `origin/main`. Its exact all-source bijections are inherited here:
for an infinite positive valuation word `w`,

`Phi_a(w)=-sum_(j>=0) 2^(s_j)/a^(j+1)` in the 2-adic integers.

The least positive odd representatives of its finite 3-source cylinders
are nondecreasing. They remain bounded exactly when `Phi_3(w)` is a positive
integer. Define `W_3^+` by this arithmetic boundedness condition and
`Y=Phi_5(W_3^+)`. There are exact clock-preserving bijections

`positive odd Collatz sources <-> W_3^+ <-> Y`.

Session 100 already excludes the named `13,33,83` basin by its transported
denominator 101. **C1 extends the same exclusion to every positive rational
`5x+1` cycle and every finite prehistory of its itinerary.** This includes
any presently unknown positive integer 5-cycle, without enumerating it.

*Proof of the extension.* C1 forces the transported 3-cycle to have a value
strictly between 0 and `1/2`, hence reduced denominator `q>1`. Denominators
are invariant around this cycle. Under an inverse 3-branch a reduced `p/q`
becomes `(2^v*p-q)/(3q)`; its numerator is coprime to `q`, so the reduced
denominator remains divisible by `q`. Every finite prehistory is therefore
nonintegral. No such word belongs to `W_3^+`. QED.

This supplies the requested target without the positive 5-cycle obstruction,
with a noncircular arithmetic membership condition. Its intended terminal
state is **-1**, the image of Collatz's 1. Proving that all of `Y` reaches -1
is still equivalent to Collatz convergence. The new exclusion handles all
positive hostile cycle basins, not arbitrary aperiodic words or negative
5-cycles. Source arithmetic has been made explicit, not automatically solved.

Session 100 also proves that a rational coordinate change preserving the
valuation-1 and valuation-2 affine branches cannot change the multiplier.
That is compatible with this construction: `Phi_5 o Phi_3^(-1)` is a
2-adic map, not a rational coordinate change on the real line.

## C2. All-order recurrence and forcing-divisor inversion

Let `P_a(d;k,H)` count primitive cyclic valuation words with `k` letters,
sum `H`, positive fixed point and reduced denominator exactly `d`. Let
`Z_(a,b)(k,H)` count positive integer cycles of exact odd period `k` and total
valuation `H` for `U_(a,b)(n)=(an+b)/2^v`, with odd `b>0`.

**Theorem.** For every such `a,b,k,H`,

`Z_(a,b)(k,H) = sum_(d|b) P_a(d;k,H)`,

`P_a(d;k,H) = sum_(e|d) mu(d/e)*Z_(a,e)(k,H)`.

*Proof.* The rational fixed point with forcing `b` is `b*x_a(w)`, integral
exactly when `d_a(w)|b`. Odd scaling preserves each valuation and exact
period. Cyclic word orbits correspond bijectively to cycles, so this has no
unaccounted automorphism factor. Möbius inversion proves the second formula.
Even `b` is outside the stated odd-state universe. QED.

Here is an exact recurrence for these quantities without listing cycles.
For a chosen target `(k,H)`, if `D=2^H-a^k<=0`, set the marked count `M=0`.
Otherwise use residues modulo the *target's* positive `D`. Initialize
`B_(0,0)(0)=1`, all other states zero. For `j<k`, append each `v>=1` that
leaves room for the remaining positive letters:

`B_(j+1,h+v)(a*r+2^h mod D) += B_(j,h)(r)`.

Then

`M_(a,b)(k,H)=sum_(r: b*r=0 mod D) B_(k,H)(r)`,

`Z_(a,b)(k,H)=(1/k)*sum_(e|gcd(k,H)) mu(e)*M_(a,b)(k/e,H/e)`.

The first recurrence compiles the numerator `C`. For the last identity,
every word is a repetition of a unique primitive word. Repetition preserves
the rational fixed point and integrality. Möbius inversion removes repeats;
a primitive cyclic orbit has exactly `k` marked words. Each smaller marked
count uses its own `D`, not the original target modulus. These are finite,
all-order identities, **not** a claimed polynomial-time counting algorithm.

At `(k,H)=(3,7)`, all five primitive word necklaces for `a=3` have denominator
101. Consequently `Z_(3,b)(3,7)=5` if `101|b`, and zero otherwise. The five
cycles are recovered by forcing 101; they do not become native forcing-1
cycles. This is a precise use of Möbius inversion, not an analogy based only
on the signs of the Möbius and Liouville functions.

## C3. Future attachments require more than the denominator

Concatenation has the inherited exact law

`C_a(uv)=a^(len(v))*C_a(u)+2^(H_u)*C_a(v)`, `H_uv=H_u+H_v`.

The pair `u=(1,1,5)`, `u'=(1,3,3)` has the same `(k,H,d)=(3,7,101)` at
multiplier 3, but `C(u)=19` and `C(u')=31`. Attach the **same** letter `(1)`:

`x_3(u(1))=185/175=37/35`,

`x_3(u'(1))=221/175`.

Thus a one-letter future test separates denominators 35 and 175. The first
failed implication is "same reduced denominator and clock means same future
integrality data." The strongest survivor is the forcing-divisor formula
for *closed* cycles. A sufficient composition carrier is `(a,k,H,C)` with its
marked source, together with the valuation word if internal phases are needed.
No claim of minimality under all later operations is made. The minimal
possible attachment length, one, already refutes the compressed alternative.

This is structurally the same information-loss test as the 81/27 tournament
control, not a bijection between these two mathematical categories. Re-running
the inherited labelled test confirms 67 common cyclic triples, automorphism
orders 81 and 27, and two-run response ranks 24 and 22 modulo 3.

## Finite-prefix cost and what remains open

The [older repeated-block fuel result](../session53_g27_collatz_cooper_20261024/collatz27-repeated-block/REPORT.md)
already distinguishes arbitrarily long realizable prefixes from a single
infinite integer orbit. For the word `(1,1,5)`, repeating it `L` times requires

`n = 19*101^(-1) mod 2^(7L+1)`.

Every finite cylinder contains positive integers. For its least representative
`r_L`, define `101*r_L=19+c_L*2^(7L+1)`. The carry satisfies
`c_(L+1)=2^(-7)*c_L mod 101`; its clock has order 100. The full source/carry,
not just this finite clock, determines how long the block can repeat. For
large `L`, `c_L>=1`, so `r_L >= (2^(7L+1)+19)/101`. This is a growing
source-height cost, not a finite universal cutoff. The checker realizes all
20 prefixes tested, including ones well beyond an arbitrary short search.

More generally, for any positive integer `U_a` segment,

`2^H <= (a+1)^k*n_0/n_k <= (a+1)^k*n_0`.

So `L` consecutive copies of a word require
`n_0 >= (2^H/(a+1)^k)^L` whenever that ratio exceeds 1. A proof for a single
repeated word does not cover an aperiodic orbit that continually changes words.
This is where the proposed convergence strategy still needs a new theorem.

## Exact audit, sources, and next decisive tests

The census contains all 21,699 positive words with `k<=5,H<=20`, reduced to
4,788 primitive rotation orbits. It checks the rational steps and signs
exactly for multipliers 3,5,7; 13,224 positive-pair cone comparisons; 588
independent residue-DP/denominator-spectrum comparisons; forcing Möbius
inversion; 20 repeated-word source cylinders; the two canonical attachment
and tournament counterexamples; and the Fermat tests specified separately.
The native cycles in this finite box are precisely the trivial 3-cycle and
the three listed 5-cycles. This is a box statement, not an exhaustive theorem
over all clocks. Gates remain active under `python -O`.

Normal and optimized runs agree byte-for-byte: 240,402 explicit local gates,
plus the invoked inherited tournament checker. The saved JSON SHA256 is
`0a95ea87b8b224105fc331e412d6e229dc24724c7145d3f0f0596e5b67dd6fac`.

```text
python -B 05-knowledge/results/collatz_fermat_reframes_20261029/verify.py --output 05-knowledge/results/collatz_fermat_reframes_20261029/results.json
python -B -O 05-knowledge/results/collatz_fermat_reframes_20261029/verify.py --output <temporary-output.json>
```

The recent [Session 99 pole bridges](../session99_collatz_pole_bridges_20261029/README.md)
make the same distinction needed here: their all-multiplier pole construction
retains the `5x+1` cycle and their pole-order Collatz readout is not a descent
proof. C1 adds a real-cone selector, C2 retains arithmetic forcing, and C3
shows why reducing the carrier too early destroys future tests. The older
mediant census already defeats a uniform "distance from integrality" gap;
we do not revive that failed approach.

Next decisive Collatz test: find an exact positive-integer selector that
controls **changes between** primitive blocks, rather than another bound for
one repeated block. Any claimed word surgery must preserve forcing 1 and
positive integrality and supply a well-founded decreasing quantity. Test it
on the inherited mediant near misses `2-1/242461` and on the one-letter
denominator collision before scaling up. Next decisive Fermat test: quantify
how certificate order, coefficient height, or chart depth grows, as detailed
in [FERMAT.md](FERMAT.md); neither finite examples nor witness density supplies
existence for every `n>=5`.
