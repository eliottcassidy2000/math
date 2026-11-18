# A physical integer decoder, with a separate absolute-count register

**Status:** all-order statements below are for the specified rooted and
metagraph constructions. None claims Chan–Pak completeness. Controls are
literal path counts. The inherited ingredients are the
[retained-root binary carry](PAIRED_AXIOMS.md) and
[two-line rim exchange](RIM_LIFT.md). The hostile control is a scalar collision
with different continuations, and the positive control is its converse in
spirit: different totals with the same continuation.

## 1. Which register remembers the word?

Fix a tournament rooted at a. Let h be its Hamiltonian path count, e its
number of paths ending at a, and d=H(T−a), with H(empty)=1. Let t be the
weighted start read: sum, over subsets S containing a, of the ways to
start a path of T[S] at a times H(T−S). In particular t≥h and d is odd.

Append a sink for letter 0. For letter 1 reverse just the arc from a to
the new sink. Keep a as root. Both operations are physical, and the
known transfer is
```
0: (h,e,t,d) -> (h,0,2t,d).
1: (h,e,t,d) -> (h+t-e,d,2t-d,d).
```
For 1, the old concatenations lose e, the uses of the new back edge
contribute t, and the paths ending at a use the new vertex just before
a, contributing d. In the weighted start sum, a new sink can be in
either part of the split. The back-edge edit removes the contributions
where the part beginning at a would otherwise consist of just a and the
new sink (loss d). The root deletion merely gains a sink.

Because d is odd, the last bit is exactly t' mod 2, and
t=(t'+d·bit)/2. After at least one preceding bit, e=d(t mod 2), so
the entire profile is reversible until the known seed. In fact
t_k=2^k t_0−d·V, where V reads the bits in chronological order,
most significant first. At a fixed depth the t register runs through
one exact arithmetic progression, with one distinct physical word
per term. If t_0>d, t increases at every step, so t alone determines
its predecessor chain back to the known seed even without the length.
Depth is not side information.

This is the precise [ANS](https://arxiv.org/abs/1311.2540) connection:
one integer retains the word, the next symbol is selected by its
residue, and the quotient returns the prior state. This is a binary
specialization, not a claim of nonuniform compression efficiency.
**The Hamiltonian register itself is not such a code.** For the old
physical seed (29,11,31,13), words 1010011 and 1111111 both give h=2317,
but their t registers are 2889 and 2317. Letter 1 then gives 5193 and
4621. These are literal rooted tournaments at both stages, including
simultaneously at sevenfold scale.

## 2. Inverting only the requested absolute count

Suppose the seed has t_0≥2d. Fix a word length k and take a prefix
of length p=ceil(k/2); let q=k−p. Read the q suffix bits least
significant first as W=Σ b_i 2^i. At prefix state (h,e,t,d) the final
Hamiltonian count is exactly
```
h+tW−dK(W)−e b_0,
K(W) = Σ (2^(j−i−1) + [j=i+1]), over i<j with b_i=b_j=1.
```
Each preceding bit 1 changes the later t by −d·2^(j−i−1), and
adjacent ones pay the same end debit d again. The empty suffix has
W=K=b_0=0. For q≥1, K≤2^q−2 by setting all bits to one and summing
the geometric series. Here e≤d once a prefix bit has been taken
and t≥(t_0−d)2^p+d≥d(2^q+1). Thus the entire debit is less than t.

Consequently **each prefix has at most one candidate suffix at the
prescribed absolute target N**:
```
W = ceil((N−h)/t),
accept iff 0≤W<2^q and dK(W)+e b_0 = h+tW−N.
```
This is a disjoint integer-bin test *inside each prefix*. Checking the
last equality is essential; the bin does not guarantee the sum lands
on N. The spelling of the word is returned, not just a membership bit.
The first formula reads the suffix in the opposite bit orientation
from the t-register formula.

Trailing 0s change no count and can be dropped. If a nonempty word
ends in 1, its last step has h'≥t and t'=2t−d≤2h'−d. Hence
2^k(t_0−d)+d≤2N−d. Testing the seed and lengths up to this exact
bound gives a **total decision algorithm for this infinite physical
language**, with O(sqrt(1+N/(t_0−d)) log(2+N)) arithmetic operations
for a fixed seed, using a direct recurrence for each prefix and
suffix. This is a bound in arithmetic operations, not in the bit
length of the input. A negative return excludes the rooted
language, not all tournaments. In particular the decoder rejects
30011 and 30013 for the (29,11,31,13) seed despite other physical
constructions hitting them. It accepts the earlier 33749 and
33751 targets with explicit words.

There is an all-order counting obstruction to trying to repair this
by adding prepared descendants of the same seed. If M is the maximum
possible depth through X, there are at most 1+Σ_(k=1)^M 2^(k−1)=2^M
different totals, hence at most 2(X−d)/(t_0−d) when M exists.
The seed (29,11,31,13) can cover at most a limsup fraction 2/9 of the
odd targets. Prepared descendants do not enlarge the original
language at all. More generally a finite bank has limsup fraction
at most 4 Σ 1/(t_0−d), for seeds with t_0>d. This bound doesn't
obstruct sufficiently different sources or a genuinely new move.

## 3. A physical +2 switch invisible to the decoder's future reads

For every suffix, its increment in h is a function of (e,t,d) only.
Conversely probes 1, 01, 11 have increments t−e, 2t,
and 3t−e−2d. These recover every coordinate. Thus (e,t,d) is the
**exact continuation equivalence** for this alphabet when comparing
increments relative to the absolute count. It is not claimed to be
the interface for arbitrary attachments.

Here is a switch that varies h without changing this equivalence
class. Use outgoing masks (126,108,56,112,98,64,4), root 3,
and reverse arc 2→5. Literal profiles are (23,2,50,9) and
(25,2,50,9). Every shared word keeps exactly the same single-arc
difference, continuation profile, and +2 count difference. For
example the common word 1000101101 gives the twin primes 32561,
32563. This is an ordinary edge edit, not a claim that it is one
merged metagraph line. To invert the switch, test the same suffix
decoder on N and N−2. To generalize, a physically supplied interval
of absolute shifts with the other three coordinates fixed remains
an interval of exactly that width at the output; the shifts do not
get scaled by the word. The next proof obligation is to *select*
such a palette at the right center, or build several independent
switches with different scales without exposing them to the later
reader. Perfect prediction of the future increments by itself does
not give the palette.

## 4. An exact quarter-step decoder using actual metagraph lines

An independent two-template decoder works on some signed
quarter-steps y=(3x+s)/4. The destination itself chooses the
template: y≡2 mod 3 selects s=−1 and y≡1 selects s=+1. Accept the
branch only on the exact power-of-two test below. Then the
quotient gives w, the signed equation gives x, and the template
constructs a physical route.

| Sign | Source x | Destination y | Core masks D | Contact mask I | Reversed edge | Core path P |
|---|---|---|---|---|---|---|
| −1 | 168·2^w−17 | 126·2^w−13 | (56,17,3,6,12,30) | 25 | 0→5 | (0,5,2,1,4,3) |
| +1 | 172·2^w−79 | 129·2^w−59 | (36,1,34,7,15,26) | 29 | 2→5 | (1,0,2,5,4,3) |

For w≥1 wrap the core with a transitive block X and apex b with
b→X→D and I→b→D−I. A path starting at b contributes H(D).
Otherwise the part just before b ends with a nonempty path of D
on a support S and ending in I. Choose the subset of X before b,
with the rest immediately after b. All 2^w choices work except
that with no X after b, the first vertex of a nonempty remainder
of D must not be in I. Thus H(w)=G·2^w+H(D)−Q, where G is the
weighted end read over I and Q is the number of nonempty ordered
two-path covers with the first ending and the second starting
in I. This proves the two linear forms from the finite core
coefficients in the replay.

For w≥2, take the source frame consisting of the first vertex
of X, the rest of X except its second vertex, P, b, and the
second vertex of X. Reverse the frame edges; up to converse this
is a merged metagraph line. The reversed frame closes using the
unchanged first-to-second edge of X. In the resulting cycle omit
the reversed core edge to get a second frame. Reversing it cancels
the entire first change except for that core edge and the
closing edge. The latter just swaps adjacent transitive twins,
so the result is the displayed target up to relabeling.
At w=2 the two routes are 655→787→491 and 609→833→457.
Width one has separate literal frame certificates (319→255→239
and 265→199) in the replay. Overshoot is allowed here.

Any source frame has at most two X-runs because the only reentry
from outside is through b. Reversing the frame creates at least
ceil((w−4)/3) vertex-disjoint triangles wholly within X, by taking
disjoint triples of consecutive vertices within its runs. At the
target every triangle meeting X must spend b and one core vertex.
The target packing number is at most
max(floor((n+1)/3),1+floor((n−1)/3)) for core order n.
Thus for w≥11 these representative pairs require **exactly two**
merged lines; this is not a statement about other representatives
of either count. For the older beta-59 family with n=8 the same
argument gives w≥14, improving its sufficient threshold of 32.

## 5. Prequential evaluation without confusing compression and coverage

Following [Dawid](https://rss.onlinelibrary.wiley.com/doi/10.2307/2981683)
and the [Hutter Prize](https://prize.hutter1.net/), compare exact
certificate generators on a predetermined stream of target height
blocks. At the start of each block freeze the template bank and
state interface. Score the bank description, the first-use cost
of any new templates, and the bits of independently checkable
literal certificates used as escapes. Targets themselves can be
public side information. Report unresolved targets separately;
assigning a prior probability to a move is not a certificate.
A bank that compresses an exceptional ray extremely well can still
miss almost every height. A useful structural metric is the number
of exact physical returns predicted *before* a block opens, not the
number for which a graph can be fitted after seeing the answers.

The count-space analogue of arithmetic coding is to predict a
contiguous interval only when every integer in it has a literal
constructor. The ceiling decoder above is weaker but honest: it
produces one candidate exact correction rather than assuming its
entire bin is legal. A successful extension would choose a new
physical prefix from the failed debit equation, then return to a
known legal continuation state. That would close the remaining gap
between reversible word memory and a universal selector for counts.

Replay: `python3 ans_replay.py` and `python3 fixtures.py`.
