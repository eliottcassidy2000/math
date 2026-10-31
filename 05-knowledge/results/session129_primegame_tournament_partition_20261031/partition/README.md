# The proposed PRIMEGAME partition: a genuine height grading, not a finite quotient

**Session129. PROVED structural statements; FINITE-EXACT saturated census.**
No claim of literature priority is made. All instructions retain their original
first-enabled order. The most useful positive result is that the proposed
partition admits a unique primitive additive height character. Its strongest
limitation is infinite continuation memory, not a mismatch of small counts.
A further positive surprise is an exactly classified infinite family of genuine
six-cycles, outside the one-control-token prime-generating domain. In fact
every periodic orbit on all positive integers has least period2 or6.

## 1. Carrier, exact order, and factorization

For prime p write e_p=ν_p(N). Primes outside
`(2,3,5,7,11,13,17,19,23,29)` are inert spectator factors. Set data counters
`(a,b,c,d)=(e_2,e_3,e_5,e_7)` and call the other six coordinates control
coordinates. This terminology does **not** restrict the theorem to one token.

| Letter | Fraction | Exponent increment | Proposed block |
|---|---|---|---|
| A | 17/91 | −7−13+17 | X1 |
| B | 78/85 | +2+3−5+13−17 | X3 |
| C | 19/51 | −3−17+19 | X2 |
| D | 23/38 | −2−19+23 | X3 |
| E | 29/33 | −3−11+29 | X2 |
| F | 77/29 | +7+11−29 | X2 |
| G | 95/23 | +5+19−23 | X3 |
| H | 77/19 | +7+11−19 | X2 |
| I | 1/17 | −17 | X1 |
| J | 11/13 | +11−13 | X0 |
| K | 13/11 | −11+13 | X0 |
| L | 15/2 | −2+3+5 | X3 |
| M | 1/7 | −7 | X1 |
| N | 55/1 | +5+11 | X1 |

Here `+p` means add one to e_p, not add p to the integer. Thus the partition is
`X0=JK, X1=AIMN, X2=CEFH, X3=BDGL`. Letters are always executed in the table's
order, not block order. The next letter is the first whose denominator divides
the current positive integer. No halt exists because the last denominator is1.
The classical products `AB=6/35`, `EF=7/3`, `DG=5/2` are inherited from Conway,
not new discoveries of this packet. Product identities alone discard priority.

## 2. Unique additive height character — the positive bridge

Define

```
T = e_2 + e_5 + e_23,
U = e_11 + e_13 + e_17,
V = e_19 + e_23 + e_29,
Φ = 2e_7 + 2U + 3V − 4T.
```

**Theorem 1.** On every positive integer, including multi-control states:

* T increases by1 at N and is unchanged by every other instruction;
* Φ changes by0 on X0, by−2 on X1, by+1 on X2, and by0 on X3;
* up to scalar, Φ is the unique additive valuation character supported on
  the ten program primes which is zero
  on X0 and X3 and constant separately on X1 and X2.

The primitive coefficient vector is
`(-4,0,-4,2,2,2,2,3,-1,3)` in the displayed prime order. For any finite
executed word w, or indeed any formal ordered product,

```
Φ(end)−Φ(start) = #X2(w) − 2#X1(w).
```

In particular every return word satisfies **#X2=2#X1**. This is an exact
structural relation between the proposed blocks, not a count-based analogy.

**Proof and uniqueness.** Every increment follows by substitution. For
uniqueness let the unknown character weights be w_p and block increments q,r
on X1,X2. The zero increments J,K,D,G,L,B force

```
w_11=w_13=w_17=t,  w_2=w_5=s,  w_3=0,
w_23=s+w_19.
```

The X1 equations M,I,N give `w_7=t=−q` and `s=2q`. The X2 equations C,E give
`w_19=w_29=r−q`, while F gives `−q−r=r`, hence `q=−2r`.
These determine precisely r times the stated vector; H and A introduce no
new restriction. The primitive orientation r=1 yields Φ. This uniqueness is
relative to the stated block conditions, not a classification of all possible
partitions of the program. If spectator-prime weights are allowed, uniqueness
is instead modulo arbitrary spectator valuation invariants.

The physical nonnegative height

```
Ψ = Φ+4T = 2e_7+2U+3V ≥ 0
```

satisfies `ΔΨ=#X2−2#X1+4#N`. Thus every actual prefix obeys
`2#X1−#X2−4#N ≤ Ψ(start)`. On an N-free interval the walk has horizontal,
+1, and −2 steps. Formal nonnegative excursions of this step set are a useful
ternary-walk model, but the extra counter guards and priority remain essential;
this packet does not identify that formal language with executable return words.

The six horizontal instructions have a further intrinsic split: J,K are
exactly the horizontal instructions leaving all four data counters unchanged.
The other four are X3. Thus the proposed 2+4+4+4 partition is recovered from
four lowering moves, four raising moves, and six horizontal moves split into
two pure-control exchanges and four data-changing moves.

## 3. Complete nonnegative circulation cone

Let z_A,...,z_N count occurrences in any formal product with net exponent
increment zero. The following is an **integer-semigroup classification**, not
just a rational nullspace computation. There are unique nonnegative integers
`g=z_D, j=z_J, t=z_C, m=z_M` such that

```
(z_A,...,z_N) =
(g+2t+m, g+t+m, t, g, g+t+2m, g+t+2m, g, t,
 0, j, t+j, t+m, m, 0).
```

Equivalently the nonnegative kernel is the free commutative monoid on

```
g: A B D E F G,
j: J K,
t: A² B C E F H K L,
m: A B E² F² L M.
```

These are multisets, not executable word orders. In particular **I and N
never occur in a return word**, and `z_L=z_C+z_M`. The block counts also obey
`#X3−#X1=2g`, besides `#X2=2#X1`.

**Proof.** Prime23 balance gives G=D; prime2 gives B=D+L; prime5 then gives
N=0. Prime29 gives F=E; prime19 gives H=C; prime11 gives K=C+J; prime13 gives A=B+C; and prime17 then gives I=0. Prime3 gives E=D+2L−C; prime7 gives L=C+M.
Substitute g,j,t,m. Conversely substitution makes all ten increments zero.
Since g,j,t,m are actual nonnegative coordinates, the characterization is
integral and exact. It also proves full increment-matrix rank10 and kernel
rank4 without a numerical rank assumption.

## 4. Complete priority-cycle classification: only periods2 and6

Let `k=U+V` be the total number of control tokens. It changes only at N (+1)
and I (−1).

**Theorem 2.** Every genuine periodic orbit is of one of the following types.

1. A two-cycle
   `S·2^a5^c11^h ↔ S·2^a5^c11^(h−1)13`, where a,c≥0, h≥1, and S contains only
   spectator primes. Its word is JK up to cyclic rotation.
2. A cycle of length6g, g≥1, using each of A,B,D,E,F,G exactly g times and no
   other instruction. It has at least three control tokens, with at least one
   token in each of `{13,17}`, `{11,29}`, `{19,23}` throughout.

Theorem3 below sharpens case2 to g=1, completing the global classification.
Already this theorem implies that all periodic orbits with at most two
control tokens are exactly the two-cycles in case1. In particular it fully
classifies periodic orbits on the invariant zero-or-one-control-token domain
containing seed2.

**Proof.** The circulation theorem excludes I,N. If k=0, only L,M,N can fire;
excluding N leaves strict decreases of a or d, so no cycle exists. Thus k>0.
At least one control prime always enables one of I,J,K,F,G,H before L,M;
hence L=M=0 in the cycle. The kernel equations imply C=H=0 and leave only
`g·ABDEFG+j·JK`. If g>0, the positive number of19/23 tokens is constant
because C,H are absent. One of G,H is always enabled, precluding J,K, so j=0.
Without J,K the three displayed control pairs cannot exchange tokens, and
each must be nonempty to perform its pair of instructions. If g=0, only J,K
fire. Their priority requires all other control coordinates zero, b=d=0,
and J drains13 until it is0; K then creates one13, immediately removed by J.
This gives exactly case1.

**Positive sharp witness.** The three-token state4301 has the legal cycle

```
4301 --G→ 17765 --B→ 16302 --D→ 9867
     --E→ 8671 --F→ 23023 --A→ 4301.
```

It realizes the classical product identity AB·DG·EF=1 with the correct
first-enabled order `GBDEFA`. Its block word is `X3 X3 X3 X2 X2 X1`.
The one-token cycle theorem therefore cannot be extended to arbitrary states.

**Theorem 3 (all remaining periodic orbits).** Every non-JK periodic orbit
has least period6 and is exactly a constant integer rescaling of the displayed cycle by

```
X = S·2^a11^r13^v19^w23^s,
a,r,v,w,s ≥ 0,  a·w=0,
```

where S is a spectator factor. The unique state just before G is4301X.
All such scales really give first-enabled six-cycles. There are no other
six-cycle word orders except cyclic rotations of GBDEFA. The three-token
lower bound and period6 lower bound for non-JK cycles are both attained.

**Why longer cycles cannot occur.** Work in case2 of Theorem2. Put

```
x=e_17,  u=e_13+e_17>0,
y=e_29,  v=e_11+e_29>0,
z=e_23,  w=e_19+e_23>0.
```

The totals u,v,w remain constant. Choose any state just before a G in the
cycle; such a state exists because g>0. F is disabled, so y=0. E is disabled
and e_11=v>0, so b=0. D is disabled, so a(w−z)=0. Also z>0, and A is disabled.

If x=0, then e_13=u>0, so A being disabled forces d=0. From this state onward
A,B,E,F remain disabled: only D,G can fire before H. D conserves a+z and
G reduces a+z by1. Each D strictly reduces a and neither D nor G increases
a. Thus only finitely many
D/G steps can occur, after which a=z=0 and H fires. This contradicts H=0
throughout a cycle. Therefore x>0.

Now B being disabled at the G-state forces c=0. If d>0, A being disabled
forces x=u. G then B then A are forced: G supplies5 for B; B creates13 for A.
After those three steps x=u, b=1,c=0,e_13=0, so C is first enabled, again
contradicting case2. Hence d=0. After the forced G,B pair, if the original
x were at least2, C would likewise be first enabled. Consequently x=1.

We have derived exactly y=b=c=d=0, x=1, e_11=v≥1, z≥1, and a(w−z)=0.
These are precisely the GBDEFA guards proved next. They give an immediate
six-step return to the chosen G-state. A periodic orbit containing that state
therefore has least period6, not6g with g>1.

**Exact guards.** Write the proposed G-state as
`S·2^a11^u13^v17·19^w23^z`, with u,z≥1, v,w,a≥0 and aw=0. Its other supported
exponents b,c,d,e_29 are0. Initially G beats I because23 is present; B and C
are disabled by c=b=0; D is disabled by aw=0. G creates one5 and one19,
so B fires; B consumes the unique17 and creates one2, one3 and one13,
so D fires; D restores a and w, so aw=0 permits E; E removes that3 and
produces29, so F fires; F produces7, so A fires and restores the start.
Conversely a GBDEFA word forces exactly these conditions: its B needs17;
initial G then forces b=c=0; the later D forces the original e_17=1 because
otherwise C would preempt it; B's created13 forces initial d=0; initial F
must be disabled so e_29=0; E requires u≥1; and initial G requires z≥1 and
aw=0. Thus this family is exact, not just sufficient.

**Independent complete six-order certification.** Even without the preceding longer-cycle exclusion, the kernel and Theorem2
imply that a primitive
six-cycle must use each of A,B,D,E,F,G once. In such a word every prime
coordinate is consumed at most once. Consequently all prefix enabling tests
are unchanged by capping each initial exponent at2: a coordinate initially
at least2 never becomes0. The checker exhausts all3^10 capped states and
finds only the six rotations of GBDEFA; it independently checks the exact
G-state guard family on every capped state. This finite saturation is a proof
for the entire unbounded state space, not a bounded experimental inference.
The base cycle's gcd is1, and it has six distinct states, so the scaling
parametrization is unique and preserves least period6.

## 5. What the four blocks do and do not quotient

### Exact possible-next-block graphs

The table gives all possible next block labels after a currently chosen label.
The one-hot domain allows zero or one control token and arbitrary data counters.

| Current | One-hot possible next | All positive integers |
|---|---|---|
| X0 | X0,X1,X2 | X0,X1,X2 |
| X1 | all four | all four |
| X2 | X0,X2,X3 | all four |
| X3 | all four | all four |

This is a nondeterministic abstraction, not a deterministic quotient and not
a tournament. There are loops and many bidirectional pairs. On the one-hot
domain the only strict pairwise directions are X1→X2 and X3→X0; the other
pairs have both directions. An arbitrary choice in each bidirectional pair
would discard genuine transitions.

The census is complete because the first two instructions depend only on
`min(e_p,2)`: denominators are squarefree and one step changes a coordinate
by at most1. The checker exhausts3^10 states, and independently compares
exact integer divisibility against exponent-vector updates. For the one-hot
subdomain its complete universe is7·3^4=567 states. Smallest witnesses in
these universal capped carriers are retained in results.json.

The globally forbidden X0→X3 transition also has a direct proof: when J or K
fires, control17,19,23,29 are absent (each has an earlier enabled denominator).
J or K leaves a11/13 token, so the next instruction is J,K,A, or E, never X3.
The extra one-hot prohibition X2→X1 fails with multiple markers:
`247=13·19 --H→1001 --A→187`.

### The reciprocal pair can be contracted — but only with a guard

Treat J,K as silent steps. Compare states with identical a,b,c,d and identical
spectator factors, each with **exactly one control token**: one at11 versus
one at13, with all other control coordinates zero. Record the first non-X0
instruction's block and resulting full state, or record infinite silent
divergence. This comparison does not allow additional control tokens.

**Theorem 4.** These two weak observations agree exactly when **b·d=0**.
If b>0,d>0, control11 emits X2 via E while control13 emits X1 via A. The
minimal witness is231 versus273. If b>0,d=0, both reach the same E transition,
with a preceding J on the13 side; if b=0,d>0, both reach the same A transition,
with a preceding K on the11 side. If b=d=0, both diverge silently forever.
After a non-silent exit the states coincide, so this is a future-preserving
weak contraction, not just one-step agreement. On the interior b,d>0 the
phase bit is indispensable. Step timing is deliberately not preserved.

### No finite autonomous observer preserves the block sequence

**Theorem 5.** There is no finite deterministic quotient of all states,
or even of the seed2 reachable orbit, that preserves all future block labels.

For the all-state assertion, state2^a emits exactly a consecutive X3 labels
(the L instructions), followed by X1 (N). Distinct a are separated by their
first exit position, giving infinitely many continuation-distinguishable
states. This is not a cardinality objection: an infinite state space can have
a finite quotient, but these particular future observations distinguish it.

For seed2, the inherited divisor-loop proof reaches, for every candidate
n≥2, a successful-divisor state `2^n7^(d−1)` immediately after I. It then
executes exactly n L instructions before M or N. These reachable states are
again distinguished by their first non-X3 output. Equivalently, the block
word has arbitrarily long X3 runs and infinitely many exits, so cannot be
ultimately periodic. The finite checker repeats the inherited first-ten
prime-output trace as a control; the all-time claim uses the inherited
symbolic divisor-loop proof, not that finite trace.

A four-node **label graph with unbounded counter/phase sidecars** is compatible
with these results. A four- or ten-state autonomous deterministic simulation
preserving those observations is not. No tournament-node identification or
n5 determination is claimed by this packet; such a bridge needs an additional
explicit map and its own preservation theorem.

## 6. Reproduction, evidence, and inheritance

From any working directory, Python3 standard library only:

```
python -B /path/to/check_partition.py > /tmp/partition-normal.json
python -B -O /path/to/check_partition.py > /tmp/partition-optimized.json
cmp /tmp/partition-normal.json /tmp/partition-optimized.json
cmp /tmp/partition-normal.json /path/to/results.json
```

All gates are explicit exceptions, active under optimization. The checker
imports no inherited code or data and requires no network. results.json is
compact; it retains the exact graphs, witnesses, matrix ranks, cone generators,
six-cycle census and source-control trace. SHA256SUMS hashes the retained
packet, excluding itself. verification.json records normal/optimized and
relocated isolated-mode replay.

The closest inherited mechanisms are Session126's exact ordered exponent
machine, four-counter/seven-state one-hot table, and divisor-loop theorem.
The hostile inherited examples are priority at91, the common seed7/17/55
basin, and the failure of a ten-node observer to preserve arbitrary controls.
The corrected near miss is identifying finite controller labels with
unbounded register states. The least-used sidecar pursued here is the full
nonnegative increment kernel and the multi-token domain, which reveals the
six-cycles rather than erasing them.

Source pin: `sources.json` records the SHA-256 of the exact inherited report
read for this work. Classical attribution: J.H. Conway, *FRACTRAN: A Simple
Universal Programming Language for Arithmetic* (1987), printedpp4,9–10,
[author chapter scan](https://gwern.net/doc/cs/computable/1987-conway.pdf).
The arithmetic/control partition results here are supplied with proofs and
self-contained computation; the only imported all-time statement is the
explicit Session126 divisor-loop route used in the seed2 observer corollary.
