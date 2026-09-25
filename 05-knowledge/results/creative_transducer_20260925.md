# A finite-word carry machine with signed run passports

**Date:** 2026-09-25. **Status:** PROVED elementary transducer and
verification statements; INHERITED run/cylinder mechanisms and descent
corollaries; FINITE-EXACT controls. Repeated outer Collatz iteration is
OPEN. No novelty claim for binary multiplication, run acceleration, or
the known one-odd-step root family.

## 1. Inheritance and the chosen object

The closest proved mechanism is
[THM-4473, Collatz digit chains, rotation, and repunits](../../01-canon/theorems/THM-4473-collatz-digit-chains-rotation-repunits.md):
the initial plus odd-run length is v2(n+1), and its all-ones boundary
shadows the 2-adic fixed point -1. The stronger exact first-reset
calculus, residue cylinders, and spatial descent density already occur
in [guards and valves, sections 3--4](collatz_guards_20260921_valves.md).
The previous [root decoder](decoder_mahler_catalan_20260925.md) retains
the actual parity word and ordered carry. These are dependencies, not
new results being renamed here.

Canonical hostile: a finite input word can expand even though its
individual digit scan terminates. Corrected near miss: finite digit
information and a favorable density do not imply every-orbit descent;
see the 2026-09-23 refreshed-digit correction in MISTAKES. Least-used
sidecar here: the finite end marker separating an ordinary integer
from an infinite 2-adic digit tail.

Anchor: a sound structure-to-arithmetic map. Niche: one small machine
for both signs. Wildcard: compress arbitrarily long, correctly guarded
runs into independently expandable proof records.

| Live concept | Exact operation | Boundary / decisive control |
|---|---|---|
| Finite binary tape | read least significant digit first | end marker is part of the object |
| Three carry states | multiply suffix by3 and add1 or2 | all states remain in0,1,2 |
| Halving | delete an initial zero | intrinsic digit operation |
| Run passport | replace maximal1^a0^b by an endpoint | both oddness guards must remain |
| Compiled cylinder | choose source class before orbit discovery | local descent is not universal coverage |
| Signed controls | change initial carry2 to1 | minus cycles5 and17 must survive |

## 2. Six local transition rules and three terminal rules

Let C_sigma be the shortcut map on positive integers:

```
C_sigma(n)=n/2                     when n is even,
           (3n+sigma)/2             when n is odd,
sigma in {+1,-1}.
```

Represent n by its finite binary word, **least significant bit first**,
followed by an end marker `#`. Thus n=6 is `011#`. Canonical positive
words end in1 immediately before `#`; the empty word represents zero
only inside the multiplier implementation.

The multiplier machine M_c computes `3x+c` on the remaining word,
with carry c in `{0,1,2}`. Its entire transition table is:

| Carry | Read0: emitted bit, new carry | Read1: emitted bit, new carry |
|---|---|---|
| 0 | 0,0 | 1,1 |
| 1 | 1,0 | 0,2 |
| 2 | 0,1 | 1,2 |

In one formula, reading b emits `(3b+c) mod2` and replaces c by
`floor((3b+c)/2)`. At the end marker, flush:

```
[0]# -> #,       [1]# -> 1#,       [2]# -> 01#.
```

These terminal outputs are also written least significant bit first.
High zero digits, if noncanonical inputs were permitted internally,
can then be removed without changing the value.

**PROVED arithmetic invariant.** After j input digits, with emitted
prefix d, unread suffix value r, and current carry c_j,

```
3x+c_0 = value(d) + 2^j(3r+c_j).                     (D1)
```

One substitution of `3b+c=2c'+d` proves preservation. The terminal
rules prove that the complete output is exactly3x+c_0.

Every scan has a unique normal form: each transition consumes one
unread digit, and the terminal rule emits at most two more digits.
The deterministic local rules therefore terminate on each finite tape.
This is a genuine finite-control normalization, not an assumed
termination of the outer dynamical system.

### The signed Collatz wrapper

```
C_sigma(0w#) = w#;
C_+(1w#)    = M_2(w#);
C_-(1w#)    = M_1(w#).                               (D2)
```

Indeed n=1+2x gives `(3n+1)/2=3x+2` and `(3n-1)/2=3x+1`.
Thus halving is literally deletion, while the odd step is a certified
finite local carry scan. No tournament relation is needed. The sign
is an actual control input; deleting it would conflate different maps.

For the smallest finite word, `C_+(1#)=01#`, whereas
`C_-(1#)=1#`. This correctly gives the plus1<->2 cycle and the minus
fixed point1. An infinite all-ones word has no terminal marker; on
the plus side carry2 loops on input1 and emits1 forever. That word
represents -1 in the 2-adics and is outside the positive finite-tape
domain. The end marker and sign separate these cases.

## 3. Signed run passports: exact compression with independent expansion

Let n be positive odd, excluding n=1 when sigma=-1. Define

```
a=v2(n+sigma),             u=(n+sigma)/2^a,
z=3^a u-sigma,
b=v2(z),                  m=z/2^b.                  (D3)
```

Then a,b>=1 and u,m are positive odd. The record `(a,u,b,m)` is a
**run passport** from n to m. Its typed verifier checks only

```
n+sigma=2^a u,             3^a u-sigma=2^b m,
a,b>=1,                   u,m positive odd.         (D4)
```

**PROVED:** these checks are equivalent to the stated maximal parity
block `1^a0^b`, with endpoint m. During the first a steps,

```
n_j=3^j 2^(a-j)u-sigma,        0<=j<=a.
```

For j<a this is positive odd; n_a=z is positive even. The next b
steps divide z by2, with even inputs until the final odd m. The
oddness of u and m proves maximality at both ends. This is the
signed integration of the inherited plus run/valve identity.

For a supplied source, the stored proof can omit redundant u and m
and recompute them from a,b; the implementation retains them for
readability. A record compresses a structured run, not arbitrary
integer data. The separate expander replays ordinary integer steps
and checks each claimed parity, so it does not merely repeat the
closed-form verifier.

The known cycles give exact hostile passports:

```
plus:  1 --(a=1,u=1,b=1)--> 1;
minus: 1 is fixed, with no finite odd-run length;
minus: 5 --(a=2,u=1,b=1)--> 5;
minus:17 --(a=4,u=1,b=1)-->41
         --(a=3,u=5,b=3)-->17.                       (D5)
```

Any claimed termination rule for this machine must accommodate these
sign-specific facts. The constant number of carry states does not
remove the unbounded tape or run lengths.

## 4. A compiler, not just an after-the-fact verifier

For a,b>=1, choose the unique odd residue

```
u_0 = 3^(-a)(sigma+2^b) mod2^(b+1).
```

For every integer t>=0, set

```
u=u_0+2^(b+1)t,             n=2^a u-sigma.            (D6)
```

Then (D3) has exactly the requested a and b: u is odd, and
`3^a u-sigma=2^b mod2^(b+1)`. Thus the machine compiles an infinite
arithmetic family with a guarded maximal run **before following an
orbit**. The plus version is inherited from valves; the common signed
implementation supplies explicit proof records for either sign.

An intentionally simple inherited descent sub-bank asks for b>=a.
Its sources are the single residue classes

```
n = 2^a(sigma*3^(-a) mod2^a)-sigma mod4^a.            (D7)
```

They admit the actual block `1^a0^a`, followed perhaps by more
halvings. For sigma=+1 their first rows are
`1 mod4`, `3 mod16`, `23 mod64`, `15 mod256`.
For sigma=-1 they are the negated residue classes.

The drop after 2a steps follows directly from

```
4^a n - [3^a n+sigma(3^a-2^a)]
  = 2^a [u(4^a-3^a)-sigma(2^a-1)] > 0,              (D8)
```

except for plus n=1. For the plus sign the minimum u=1 gives
`F_a=4^a-3^a-2^a+1`; F_1=0 and
`F_(a+1)-4F_a=3^a+2^(a+1)-3>0`. For the minus sign all terms in
the bracket are positive. This is only a sufficient local bank.

The classes (D7) are disjoint, because they have different a values,
and have densities4^(-a). Their union has natural density1/3 among
all positive integers, or2/3 among odd integers. The omitted a>A
tail lies inside one residue class modulo2^(A+1), so the density
passage is justified. Including immediate even halving gives5/6
of all starts, ignoring finite terminal exceptions. This **does not
improve** the stronger first-reset density in valves, sections3--4,
and is not a frequency statement along any orbit. The same simple
density holds on the minus side despite its known nontrivial cycles.

### Generated root-certificate control

The familiar one-odd-step family

```
n_k=(4^k-1)/3,                k>=2,
```

has passport a=1, b=2k-1, m=1. Its expansion reaches shortcut root4
after exactly2k-2 steps. This is a genuine generated infinite family,
but is inherited from the inverse sibling ladder rather than a new
coverage theorem. The implementation also checks k=4096: an8191-bit
input receives an8190-step root4 certificate from the short run counts.
It is a constructive large control for the transducer/passport interface.

## 5. A precise boundary for local termination rules

For every k>=1 and either sign, set

```
n=2^(k+1)-sigma.
```

Then n is a positive nonterminal odd integer and

```
C_sigma(n)=3*2^k-sigma = n mod2^k.                   (D9)
```

Hence no rank depending only on a **fixed number of low binary digits**
can strictly decrease at every nonterminal Collatz step: (D9) forces
equal rank on one such edge. This does not exclude ranks with an
unbounded height/register, delayed descent rules, automata with richer
acceptance, or a successful global proof. It identifies the exact
information that this bounded-prefix proposal would lose.

The terminating scan and the outer rewrite have different ranks.
Unread input length decreases during one scan, but an odd outer step
can emit a longer word. The plus root family above shows useful
normalization; the minus cycles show that finite local rules alone
do not establish global termination.

## 6. Connection contract and reproduction

```
source: a finite positive binary word, sign, and optional run passport;
target: the next word / compressed arithmetic path;
map: (D2), or equations(D4) plus independent expansion;
preserved: every actual shortcut edge, parity guard, source and endpoint;
lost on compression: intermediate words and clock length;
restoration sidecar: sign plus run counts; source fixes all intermediates;
decisive test: expand using integer arithmetic and compare every parity.
```

Run the [checker](../../04-computation/experiments/creative_transducer_20260925.py):

```
python 04-computation/experiments/creative_transducer_20260925.py
python -O 04-computation/experiments/creative_transducer_20260925.py
```

The [output](creative_transducer_20260925.out) records all sources0..4096
with each carry, all positive sources1..16384 with both signs, every
odd source1..16383 with independent maximal-run expansion, all864
compiled cylinder controls with a,b1..12 and t0,1,3, plus signed cycles,
malformed/sign-transplanted records, bounded-prefix hostiles, and the
large root family. No random samples, floating point, or removable
assertions. Positive finite words and the minus fixed point are handled
explicitly.

The remaining obligation is not the local arithmetic decoder: it is
to produce a finite chain of sound records ending at the root for
every requested positive source. These records are suitable edges
for a common-future proof DAG, and their compiler supplies infinite
verified subfamilies. Universal reachability/termination remains OPEN.
