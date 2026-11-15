# The reversed carry is the real phase; select it from the target

For a chosen sum and an actual address interval, the algorithm below
chooses a nearest feasible carry word in a certified first-order model.
It can also reject the entire interval. It does not assert that every
target has a successful interval.

## A carry address equals the reversed checksum phase

Fix s=x+y>0, r with s<2^r, and carries c_i into bit i from the
right, so c_0=c_r=0. Put
```
rho(n) = sum_{i>=0} n_i 2^(-i-1),
C = sum_{i=1}^r c_i 2^(-i-1).
```
The column equation x_i+y_i=s_i+2c_(i+1)-c_i gives exactly
`rho(x)+rho(y)=rho(s)+3C` after summing with weights 2^(-i-1).
No independence assumption enters.

The recurrence for p is a sum, over zero bits i, of
`floor(n/2^(i+1))+n_(i+1)`. Replacing the floor by its argument
causes error at most one per position, including the correction from
the adjacent bit, and the remaining geometric tail is less than
one. Thus |p(n)−n(1−rho(n))|<=r+1 when 0<=n<=s. If t=x−s/2,
the literal tournament transfer is
```
F(x,y) = B+Q(s²/4−t²)+(T+S)s/2+(T−S)t
       +45s(1−rho(s)/2) − (135s/2) C
       −45t(rho(x)−rho(y)) + epsilon,
|epsilon| <= 90(r+1).
```
The main digit-dependent part is linear in the *reversed carry*,
not an independently selectable checksum. The remaining difference
of reversed addresses is weighted by t, small near the balanced
direction.

## Deterministic selection in a prescribed physical address interval

Take integers 0<=L<=U<=s, put t0=(L+U−s)/2, w=(U−L)/2, and let
D be the part of the formula with t=t0, C=0 and the last line
omitted. Every x in the interval, with y=s−x, satisfies
```
|F(x,s−x) − [D−(135s/2)C]| <= E,
E = Q(2|t0|w+w²)+|T−S|w+45(|t0|+w)+90(r+1).
```
So choose the allowed carry nearest C*=(D−N)/(135s/2).
This is an exact optimization of the carry, using O(r)
bounded-state transitions after the threshold is read into place-value
bits. Addresses and arithmetic have O(r)-bit size (or O(log N+r) for
an arbitrarily distant target).

Scan columns i=0,...,r−1. Keep c_i, comparisons of the current
lower-bit prefixes of x against L and U, and the comparison of
c_1...c_i to the corresponding prefix of the threshold carry.
There are at most 2·3³=54 states per layer. Choose the column
bits b,b' and next carry with b+b'+c_i=s_i+2c_(i+1).
If s_i differs from c_i, the next carry must stay c_i and one
unit may go to either tape. If s_i equals c_i, the next carry
can be either digit but b,b' are forced for each choice.
A new differing address bit overrides old comparisons because
it is more significant; the first differing carry bit takes
precedence because the carry is reversed.

For the best predecessor of C*, require the carry numerator
2^(r+1)C at most floor(2^(r+1)C*). Mark live suffix states
backward, imposing c_r=0 and both address bounds. At each
layer take the largest next carry bit among all live successors
of the currently retained states; retain all states for that
choice. A successor uses the ceiling and smallest bit. Handle
thresholds outside the representable range by rejecting the
empty side and clipping the other. The winning predecessor
or successor recovers one full pair by stored transitions,
with a fixed tie break. This selects the target's carry rather
than computing a menu of every attainable checksum.

Let R be the resulting model distance. The physical pair is
within R+E of the target. If R>E, the **entire prescribed
address interval** has no exact witness, because no feasible
carry model is closer and every physical discrepancy is
bounded by E. This is a quantitative boundary test, not a
conclusion that any particular missing target is absent from
the full spectrum. It keeps the scale envelope and the
discrete carry constraint as separate mathematical data.

## Replay

Run `replay_carry_selector.py` to check the carry identity and
the nearest-carry states against brute-force physical addresses.
