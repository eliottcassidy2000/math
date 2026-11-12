# New clock rank beyond the old stable bank

**PROVED by full finite skeleton certificates; exact family scope.**
This positive construction is separate from the old fixed-bank union in
[PROOF.md](PROOF.md), whose generating function is unchanged.

The old bank can have six stable vertices, but after canonical unselected
transitive-twin compression it has at most **two genuine independent length
clocks**. That rank does not persist at larger orders.

## 1. A first three-clock strong core

Take outgoing rows

    R9=(158,156,8,16,36,335,159,60,223), O={0,6}, P={1,2,8}.

It is strong, its marked automorphism group is trivial, and it has no
unselected transitive twins. The ports induce a transitive tournament. The
complete compressed-run defect tensor is

    a_(1,1,1)=1, all other a_r=0.

Hence delta=1 for **every** independent positive triple of port lengths.
This gives exactly `binom(n-7,2)` marked classes at order n>=9 and series
`z^9/(1-z)^3`. The parameterization is injective because the twin quotient
is canonical and marked rigid. The base has H=143, W=142 and U=2004, checked
independently by the audit's permutation/split calculation.

The complete inherited bank through8 has no three-clock case, so order9 is
the first occurrence of **three jointly free canonical unselected twin
lengths in a strong marked unit seed**, relative to that inherited exhaustive
lower-order classification. This is not a new exhaustive order9 census.

## 2. Four clocks at order ten

Adjoin vertex9 to R9, with its outgoing old-vertex mask148={2,4,7}. The result is

    R10=(670,668,8,528,36,847,671,60,735,148),
    O={0,6}, P={1,2,8,9}.

Again the graph is strong and unselected-twin-reduced, the marked quotient
is rigid, P is transitive, and its **entire** tensor is
`a_(1,1,1,1)=1`, zero elsewhere. Thus it gives exactly `binom(n-7,3)` marked
classes at order n>=10 and series `z^10/(1-z)^4`. No minimal-order assertion
for four clocks is made: the full order9 seed universe was not enumerated.

These two families and the old fixed-bank union are pairwise disjoint as
marked isomorphism families, because their canonical twin quotients have
orders9,10, and at most8 respectively. Their exact combined series is

    C(z) + z^9/(1-z)^3 + z^10/(1-z)^4.

By Pascal's identity its coefficient at every n>=8 is exactly
`19n-49+binom(n-6,3)` (123,145,170 at n9,n10,n11).
In particular its coefficients grow cubically. This does not say all seeds
have polynomial growth, nor that these families exhaust large orders. Each
finite union still obeys the inherited physical height obstruction; the
new seed-tail contributions are at most L^4+L^5 below X in the notation of
PROOF.md, section6.

## 3. The immediate clock-extension rule stalls

A tempting extrapolation is to keep the same six exterior vertices and add
one new stable type at every step. We test exactly that restricted move from
R10. A new vertex can be inserted in one of five positions in the existing
four-port transitive order, and has64 possible orientations to the six fixed
exterior vertices. These320 templates exhaust all choices for which the
five candidate ports remain transitive and the original graph/marks stay
unchanged. Exclude nonstrong outputs and any output whose new port merges
with an old unselected twin block, since that is not a new canonical clock.

The complete calculation in `clocks.json` finds no surviving unit seed that
passes even all the necessary one-port length2 tests. Therefore no template
in this specified list yields a new jointly stable five-clock chart. This is
a finite exact obstruction to this immediate move, not a bound on global
clock rank, not a restriction on two-vertex replacements, and not an
exhaustive search over larger cores.

A global source or sink can add a cheap clock, but would destroy strongness;
those candidates are explicitly excluded. The useful next question is whether
one can change an exterior signature, add a coupled pair, or move the marks
to obtain an iterable rank-raising operation. That remains OPEN here.

## Certificate

`clocks.py` recomputes both full tensors, strong connectivity, twin reduction,
marked rigidity, literal all-ports-doubled counts (order12 and14), and each
single-port length3 control. It exhausts the320 continuation templates. Its
normal and optimized receipts are identical. Positive-length universality
comes from the proved acyclic finite skeleton bound, not those literal spot
checks. The original signature hunt is not asserted exhaustive outside the
precisely enumerated continuation templates.
