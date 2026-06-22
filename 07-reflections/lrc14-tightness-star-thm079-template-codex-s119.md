# LRC14 tightness star in the THM-079 template

The owner's proposed template is the right shape:

```text
THM-079: reduce H=21 to one strong atom, then force the atom impossible.
LRC14:   reduce to one bounded/top-balanced atom, then force the atom loose.
```

The exact S119 audit makes the boundary visible.  Both known tight rows,

```text
AP = {1,...,13}
GW = {1,...,11,13,24}
```

have `M=1/14`, and their full argmax sets are the same six denominator-14
points:

```text
1/14, 3/14, 5/14, 9/14, 11/14, 13/14.
```

The binding pairs are also the same: `{1,13}`, `{5,9}`, and `{3,11}` at the
paired unit points.  The denominator-14 grid obstruction is exact: for those
six units `k`, `14|vk` if and only if `14|v`.  Hence any row with no multiple
of `14` has all six denominator-14 survivor points and cannot be a strict
counterexample.

That proves the "apex-7 floor" part cleanly.  It also explains why AP and GW
are non-covering tight boundary rows: they do not contain a multiple of `14`.
The q-covering rows are on the other side of the wall.  In the exact finite
window `[1,18]`, every row satisfying the necessary q-covering condition
(a multiple of every `q=2..14`) has slack; the minimum is `M=1/12`, and no row
is tight or below threshold.

Incoming THM-568 sharpens the state after this audit.  It proves the
apex-denominator lemma: if a tight optimum is `t=a/D`, then `14|D`, the
binding runners have sum divisible by `D`, and in fact `D=14*gcd(S)`.  Thus a
primitive tight row is forced to optimize at denominator `14`.

That means the structural half of the star target is now theorem-level.  The
remaining branch is no longer "prove all tight rows are AP/GW" in the broad
sense; it is the localized statement:

```text
if S contains multiples of 14, then S is not tight.
```

The 14-free case is done by THM-523 plus THM-568: a 14-free tight row has a
denominator-14 survivor and therefore optimizes at the apex.  The case with at
most six multiples of 14 is also done by the 14-free core's `1/13` margin plus
the comb-teeth union bound.  What remains is the `>=7` multiples-of-14 branch,
an apex-localized second-moment/equidistribution problem.

The useful correction remains that the equality theorem alone is not the proof.

```text
M(S)=1/14 => AP/GW
```

does not by itself exclude an `M(S)<1/14` bounded over-cover.  To make the
THM-079 template a proof, the atom step must include one of:

```text
strict over-cover -> tight boundary by compression/compactness;
strict over-cover -> forbidden K3/H7 packet by state lift;
or an exact bounded-core finite atlas with positive minimum margin.
```

This is where HYP-2908 enters.  The K3 slogan is not generic binary data.
Loose digraphs can realize Hamiltonian-path value `7`, and `K3` subgraphs can
occur inside larger conflict graphs.  The forbidden object is a connected
tournament/OCF packet with `I(.,2)=7`.  So the LRC packet vertices should not be
raw runners.  They should be cover arcs, exact-period packets, measured
sector-state words, support-six relation packets, or whatever quotient keeps
the LRC over-cover predicate and the two-state activity.

So the next sharp target is:

```text
bounded apex-7 over-cover
  -> either AP/GW boundary plus non-covering contradiction,
     or exact slack,
     or tournament-conflict-realizable K3 packet.
```

That is the real LRC14 Moon step.

After THM-568, a more concrete formulation is:

```text
S = R union M14,  |M14| >= 7,
R has at most six speeds and is 14-free,
R has an interval with margin >1/14 by LRC(<=13).

Show the danger combs from the multiples of 14 cannot cover that interval.
```

That is the same Moon step, but now in exact arithmetic rather than broad
tight-locus language.
