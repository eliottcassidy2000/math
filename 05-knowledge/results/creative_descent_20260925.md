# Affine ports, an explicit residual, and counter-sized root certificates

**Status: PROVED elementary schemas and reduction rules; FINITE-EXACT controls.**
The global positive Collatz root-certificate obligation remains **OPEN**.
This is a compilation and specialization of inherited results, with no
priority claim. No convergence assumption, density-to-all inference, or
external cycle-classification theorem is used.

## Inheritance and the live board

The closest proved mechanisms are the ordered affine carry in
[THM-4469-mahler-bridge-adjacent-block-pairs.md](../../01-canon/theorems/THM-4469-mahler-bridge-adjacent-block-pairs.md),
the exact first-reset boundary in
[the valves note, sections 3–4](collatz_guards_20260921_valves.md), and
[the inverse-completion theorem](arithmetic_braids2_20260917_inverse_completion.md).
The last already completes any prescribed finite halving word to any
admissible target. The root family below is its particularly transparent
one-rise specialization, with an explicit exponent phase and integer rank.

The canonical hostiles are arbitrarily long growing all-one prefixes and
the positive `3n-1` cycles containing 5 and 17. The corrected near miss is
an unguarded existential affine certificate: allowing its word counters
to vary independently of the actual orbit merely restates descent. See
[the compression audit, S15](collatz_mod6_20260922_compression_and_lean_audit.md).
Bounded-lookahead obstructions are proved in
[THM-4471-kawasaki-fixed-point-collatz-proof-refuted.md](../../01-canon/theorems/THM-4471-kawasaki-fixed-point-collatz-proof-refuted.md)
and [the energy note](collatz_blueprint_20260921_energy.md).
The least-used sidecar is an ordinary integer height attached to a
symbolic port, together with its exact source cylinder.

| Live concept | Exact data retained | Test or boundary |
|---|---|---|
| Affine port | parity word, sign, residue, ordered carry | independent forward expansion |
| First reset | `n+sigma=2^a u`, actual zero-run length | a zero run is forced, not chosen |
| Adaptive cover | guarded exit maps and residual words | exact disjoint cylinder cover |
| Counter grammar | arbitrarily large rise and fall counters | root certificates beyond every fixed cutoff |
| Well-founded rank | ordinary positive integer at port boundaries | every reduction strictly decreases it |
| Reflected sign | carry sign and residue assignment | 5 and 17 remain hostile cycles |

Anchor: a reusable proof-edge interface with an explicit residual. Niche:
closed-form certified ancestors with unbounded first growth. Wildcard:
the same compiler on the reflected sign exposes what a structural argument
still fails to distinguish.

Throughout use the shortcut map, for positive integers and `sigma=+1,-1`,

```text
C_sigma(n) = n/2                    if n is even,
             (3n+sigma)/2           if n is odd.
```

For this shortcut convention the plus terminal cycle is `1<->2`;
4 enters it and is not itself recurrent. A universal terminal-core
certificate must therefore accept 1 and 2 separately or identify them
with the root bug. The explicit families below do literally visit 4.

## 1. A port is a cylinder with its actual exit map

Let `w=(e_0,...,e_(L-1))` be a binary parity word. Define

```text
p = 3^(sum e_j),          q = 2^L,
R = sum_j e_j 2^j 3^(sum_(i>j) e_i),
r = -sigma R p^(-1) mod q,          0 <= r < q,
s = (p r + sigma R)/q.
```

The empty word has `p=q=1`, `R=r=s=0`. The exact port is

```text
n = r+qk       -->       C_sigma^L(n) = s+pk,       k>=0.       (P1)
```

Only positive inputs are used for a reduction. This is an iff source
guard: the first L parities equal w exactly when `n=r mod q`.
The affine identity follows by updating `R <- 3^e R+e 2^j` at step j.
For the converse, the displayed congruence forces the first parity
(`R mod2` is the first bit); after the first actual step the remaining
congruence is the suffix congruence. Induction gives the entire word.
There is exactly one source residue since p is odd.

If `q>p`, define the explicit integer threshold

```text
K = max(0, floor((s-r)/(q-p))+1).                           (P2)
```

Then, among the nonnegative lifts, the exit is strictly smaller than its
source **iff** `k>=K`. This follows directly from
`s+pk < r+qk`. The finitely many lifts `k<K` are retained as a fringe;
they cannot be silently accepted. When `p>q` no descent is claimed by
this compiler. On the minus sign some such ports can descend or return,
so this is a sufficient bank, not a complete descent classifier.

The source is an actual finite orbit segment; the target is its guarded
affine port. The map preserves source membership and the exact endpoint.
Discarding w loses the order of intermediate states; discarding sigma
loses the reflected carry; discarding k loses the size threshold. The
necessary sidecar is therefore `(w,sigma,r,s,K)` or equivalent exact data.
The cheap decisive test is direct expansion of both a small lift and a
lift just across the threshold.

Composing ports requires the exit of one to satisfy the source guard of
the next. At accepted boundaries the rank is the positive integer n.
Strict decrease of this rank forbids an infinite sequence of accepted
reductions; intermediate forward values are allowed to rise. Termination
of the compiler can end in a residual input, so this alone is not a root
theorem.

## 2. An adaptive finite cover leaves a constructible residual

Start with the empty word, split a word by appending 0 and 1, and accept
it at the first prefix for which `3^(number of ones)<2^(length)`. Stop
splitting unaccepted nodes at a chosen integer depth D. This gives a
prefix-free finite set of ports `A_D` and a residual language

```text
W_D = {w in {0,1}^D : 3^(ones in w[0:j]) >= 2^j
                         for every j=1,...,D}.
R_(D,sigma) = {r_sigma(w)+2^D k : w in W_D, k>=0} intersect Z_(>0).
                                                               (P3)
```

The actual uncertified set for this bank is `R_(D,sigma)` together with
the positive fringe of accepted ports. These are explicit arithmetic
progressions plus a finite set, not an unspecified exceptional set.
The split construction partitions all residue classes modulo `2^D`.
For any fixed positive input, repeatedly applying accepted ports is a
finite algorithm: it either reaches a supplied terminal certificate or
halts at this represented residual/fringe. All generated edges are real
Collatz segments.

**FINITE-EXACT at D=14:** each sign has 142 accepted ports and 734
residual classes modulo 16384. The positive fringe is `{1}` on plus and
empty on minus. The residual classes differ between signs; for example
the minus residual contains 1, 5, and 17. The identical counts are
expected because the stopping language depends only on word slopes.
This recovers the inherited 734-class control in
[the minus-sheet note](collatz_mod6_20260922_minus_sheet_positive_control.md).
No claim identifies slope stopping with actual stopping for every input.

The script constructs the full sets by the stated algorithm and checks
them independently by native forward iteration of every one of the
16384 input residues. The formulas (P1)–(P3), rather than only the first
few printed residues, specify the reproducible residual.

## 3. A counter grammar compiles an infinite descent bank

For plus, choose `a>=1` and positive odd u with
`3^a u=1 mod2^a`. Starting at `n=2^a u-1`, the first a odd steps reach
`3^a u-1`; at least a even steps then give

```text
y=(3^a u-1)/2^a,
n-y = ((4^a-3^a)u-2^a+1)/2^a.                             (P4)
```

The exact parity prefix is `1^a0^a`. For `a>=2`,
`4^a-3^a>2^a-1`, so every such source strictly descends. At `a=1`
only `n=1` gives equality; the other lifts descend. The inequality starts
at a=2 and follows inductively since
`(4^(a+1)-3^(a+1)-2^(a+1)+1)
 =3(4^a-3^a-2^a+1)+4^a+2^a-2`.

Equivalently this is the single residue
`n=2^a(3^(-a) mod2^a)-1 mod4^a`, with threshold zero for `a>=2`
and threshold one for a=1. These are infinite arithmetic progressions.
They are a simple sufficient subcase of the stronger inherited
first-reset boundary in the valves note, not a new first-reset theorem.
The unbounded a makes this a counter grammar rather than one fixed
lookahead table. The actual first zero run may exceed a, which is harmless
because this certificate stops after its first a zeros.

## 4. Exact root families and reusable target ports

Let v be any positive odd target with `3` not dividing v, let `a>=2`,
and choose an integer `b>=a` satisfying

```text
2^b v = -sigma mod3^a.
u=(2^b v+sigma)/3^a,
N=2^a u-sigma.                                             (P5)
```

Then u is positive and odd. For `0<=j<=a`, the first j odd steps give

```text
C_sigma^j(N)=3^j 2^(a-j)u-sigma.
```

For j<a these numbers are odd. At j=a the result is `2^b v`, so the
next b steps are even and reach the odd target v. Thus `1^a0^b` is an
exact maximal rise/fall passport. Moreover

```text
3^a(N-v)=(2^(a+b)-3^a)v+sigma(2^a-3^a) > 0.               (P6)
```

For minus both terms are positive. For plus, the lower bound obtained
from `b>=a` and `v>=1` is `4^a+2^a-2*3^a>0`. Its value at a=2 is
2 and its next value is three times the preceding value plus
`4^a-2^a>0`. Hence this passport strictly decreases the integer rank.
It can be prepended to any supplied certificate for v.

**The exponent phase is exact.** The order of 2 modulo `3^a` is
`2*3^(a-1)`. An elementary proof starts with
`v3(4^h-1)=1+v3(h)` for h>0: for h not divisible by3, factor
`4^h-1=(4-1)(1+4+...+4^(h-1))`, whose second factor is h mod3;
cubing a number `1+3^j z`, with `3` not dividing z and j>=1,
increases that valuation exactly by one. An odd exponent of 2 is
`-1 mod3`, so the claimed order follows. There are exactly that many
units modulo `3^a`, so powers of 2 exhaust them. Consequently (P5)
has one and only one b phase modulo `2*3^(a-1)`. Adding periods obtains
infinitely many b satisfying `b>=a`.

For the target v=1 this gives complete, explicit root families:

```text
plus:  b=3^(a-1)(2t+1),     N=2^a(2^b+1)/3^a-1,
minus: b=2*3^(a-1)(t+1),    N=2^a(2^b-1)/3^a+1,
                             a>=2, t>=0.                 (P7)
```

For every positive b, the plus numerator is divisible by `3^a` **iff**
`b=3^(a-1) mod 2*3^(a-1)`; the minus numerator is divisible by
`3^a` **iff** b is a multiple of `2*3^(a-1)`. For plus, the element
`2^(3^(a-1))` is the nonidentity element of order two and therefore is
-1: the only square roots of 1 modulo an odd prime power are +/-1,
as follows by factoring `(x-1)(x+1)`. This proves both iff statements,
not merely sufficiency of the displayed examples.

Every (P7) source reaches 4 after `a+b-2` shortcut steps and 1 after
`a+b` steps. The first three plus examples at t=0 are
`N=3,151,26512143`, with `(a,b)=(2,3),(3,9),(4,27)`.
On minus the corresponding sources are
`29,77673,3558399705576689`, with b equal to 6,18,54.
For each fixed a the sources grow strictly with t, giving a genuine
infinite family of complete root certificates. The certificates can
store the counters and the power identity, without listing all halvings.

**These certificates exceed every fixed cutoff.** Given D, choose
`a>max(D,1)` in either family. The first D bits are all 1 and every
one of those steps increases the source, so it lies in the residual
`R_(D,sigma)` of (P3). Nevertheless its counter passport proves the
complete route to 4. For D=14 the script checks a=15 on both signs;
the plus certificate has `b=4782969` and a source with 4782961 bits.
It verifies the 15 odd steps and the terminal power identity exactly,
without traversing the millions of halving steps.

This resolves explicitly selected residual points of arbitrarily large
lookahead. It does not resolve every residual point. The inverse
construction changes the starting integer; it must not be used to alter
the future of a given input. That distinction is already central in the
inherited inverse-completion theorem.

## 5. The remaining obligation stays visible

The minus cycle through 5 has shortcut word `110`, with
`p=9,q=8,R=5`; the cycle through 17 has word `11110111000`, with
`p=2187,q=2048,R=2363`. Their exact equalities are

```text
9*5-5=8*5,
2187*17-2363=2048*17.
```

They coexist with all the minus root families in (P7). Therefore
existence of arbitrarily large, explicitly root-certified ancestors,
or of a large collection of guarded descent ports, does not prove
universality. Sign-blind counts cannot remove this obstruction.

A concrete open target is now: enlarge the accepted port/schema bank
until every positive plus input is either terminal or has a verified
port to a smaller certified input. It is also legitimate for another
proof object to provide a direct root certificate for a residual input.
Any such claim must use the input's actual word and carry. The represented
residual is available for further investigation; no density estimate is
substituted for its emptiness, and no real-valued non-well-founded rank
is used.

## 6. Exact reproduction and scope of the controls

```text
python 04-computation/experiments/creative_descent_20260925.py
python -O 04-computation/experiments/creative_descent_20260925.py
```

Both runs agree with [the frozen output](creative_descent_20260925.out).
The script uses integer arithmetic and active checks even under `-O`.

- Both signs, depth 14: all 16384 native residues, all accepted/residual
  port carries by an independent sum, and lifts k=0,1,2,17 (3504 checks
  per sign). Threshold fringes and cover mass are checked exactly.
- Balanced grammar: a=1..64 and k=0,1,3,100, including the equality n=1.
- Complete root families: a=2..9 and t=0..3 on both signs, plus the
  a=15 witnesses beyond the cutoff. Only compressed terminal powers
  replace iteration of the long halving tails.
- Integrality iff: a=1..8 and every b in two complete positive exponent
  periods, both signs, 26240 comparisons including failures; every
  unit in each period is also checked to occur exactly once.
- Generic target port: every positive odd v<=101 not divisible by3,
  a=2..6, both signs (340 cases); unique phase and direct odd-step
  endpoint verified.
- Hostiles: growing Mersenne prefixes through length 64 and both
  reflected cycles, derived directly from their forward orbits.

The infinite statements follow from the elementary proofs above, not
from extrapolating these finite universes. The independent sibling
[finite-tape transducer](creative_transducer_20260925.md) supplies a
different local implementation of the same signed rise/fall passport.
An independent read-only audit checked (P5)–(P7), the order and integrality
iff proofs, strict integer rank, the residual inclusion beyond every
fixed cutoff, and an optimized replay: PASS. `Port` is an internal
constructor used here with legal binary words and signs; it is not
advertised as a validator for arbitrary malformed external records.
