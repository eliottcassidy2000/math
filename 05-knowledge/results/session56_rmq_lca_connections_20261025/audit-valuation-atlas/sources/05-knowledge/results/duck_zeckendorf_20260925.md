# Fibonacci colors survive distinct decompositions; diagonal carries need one bit

**Status: PROVED elementary normal-form, charge, carry, and action statements;
FINITE-EXACT controls.** This note specifies a precise three-color convention
compatible with Fibonacci recursion. No current repository theorem explicitly
named “three-color Zeckendorf decomposition” was located in the targeted
inheritance search, so the convention is stated here rather than attributed to
an unidentified result. No Collatz, square-sum, or graceful-tree equivalence is
claimed.

## Inheritance and scope

The closest established mechanism is Zeckendorf's unique normal form, using
weights `1,2,3,5,8,...`, with coefficients zero or one and no adjacent occupied
positions. The primary account
[Zeckendorf, A Generalized Fibonacci Numeration (1972)](https://www.fq.math.ca/Scanned/10-4/zeckendorf.pdf)
also describes generalized Fibonacci sequences determined by two initial
values. A modern primary paper records the standard uniqueness convention:
[Kologlu--Kopp--Miller--Wang (2010)](https://arxiv.org/abs/1008.3204).
The proofs needed here are included below; no novelty claim is made for
Fibonacci numeration, its companion matrix, or its Beatty coordinates.

The precise geometric inheritance is
[THM-3339, Fibonacci three-ray Berggren transplant](../../01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md).
It proves a period-three primitive-content split, three exact golden
Berggren rays, and the need to retain a four-state affine owner/current
coordinate. Section 5 below gives the exact interface to our charge; it
does not identify two different carry cocycles. The relevant repository
correction is explicit in the warning on
[the historical summand/Zeckendorf reflection](../../07-reflections/summand-graph-fermat-zeckendorf.md)
and the repaired
[arithmetic-braids summand note](arithmetic_braids_20260917_summand.md): apparent
Fibonacci resets of the unrestricted distinct-summand depth were a cutoff
artifact. That history cannot justify a coloring theorem. The current
[ternary clock](ternary_berggren_20260925.md) supplies a separate exact
three-adic residue action, not a Fibonacci representation.

The canonical hostile is `2=1+1`: repetitions and distinct decompositions are
different universes. The corrected near miss is assuming that two distinct
integer summands have disjoint Fibonacci supports. The least-used sidecar is
their shared Fibonacci support, compressed below to an exact carry and, for
colors, one bit.

| Live concept | Precise object | Boundary |
|---|---|---|
| Canonical representation | binary word without adjacent ones | unique for every integer |
| Relaxed distinct representation | arbitrary finite binary word | multiple words, same two-coordinate charge |
| Repeated representation | finite nonnegative multiplicities | charge needs a carry correction |
| Diagonal color | XOR of the colors of two summands | overlap can change it |
| Ternary action | a cycle on three nonzero binary vectors | an action, not addition modulo three |
| Sum/difference graphs | actual integer edge labels | a color does not certify label admissibility |

Anchor: resolve whether decomposition ambiguity breaks the color. Niche:
find the smallest carry coordinate that repairs the diagonal. Wildcard:
identify the exact order-three action shared with ternary residue words and
the ten-state unordered-pair object.

## 1. Three representation conventions

Use standard Fibonacci numbers `F_0=0,F_1=1` and define

    f_i=F_(i+2), i>=0: 1,2,3,5,8,... .                   (Z1)

There is only one weight 1 in this alphabet. Define:

1. Canonical: `n=sum e_i f_i`, `e_i in {0,1}`, with `e_i e_(i+1)=0`.
2. Relaxed distinct: the same formula with arbitrary finite binary digits.
3. Repeated: coefficients are arbitrary nonnegative integers of finite support.

Canonical representations are unique for Fibonacci and non-Fibonacci
integers alike. Relaxed distinct representations need not be unique even
for Fibonacci numbers: `3=3=1+2`. Repeated representations are larger again.

For completeness, greedy subtraction constructs the canonical form. If
`f_j<=n<f_(j+1)=f_j+f_(j-1)`, its remainder is below `f_(j-1)`, so it skips
the neighboring weight. To see uniqueness, the greatest sum of a
nonadjacent subset of `f_0,...,f_(j-1)` is `f_j-1`, by induction. The largest
differing selected weight therefore cannot be canceled by the smaller
weights of the other canonical form.

There is also a terminating, confluent rewrite on relaxed distinct words:

    011 -> 100                                             (Z2)

where the printed positions run **high to low**, with weights
`f_(i+2),f_(i+1),f_i` from left to right; equivalently replace the two
lower occupied positions by the empty higher position. Each rewrite
preserves n and decreases the number of
occupied positions by one. Whenever adjacent ones exist, the top two of
their maximal run have an empty position above, so some rewrite applies.
Every terminal word is canonical. Uniqueness then implies every legal
rewrite order has the same terminal word. Thus the graph of all distinct
representations of a fixed n is connected through its normal form.

## 2. A two-coordinate invariant of every distinct decomposition

Give each atom an integer charge

    v_0=(1,0), v_1=(0,1), v_(i+2)=v_i+v_(i+1).          (Z3)

Its second coordinate is `F_i`; write the first as `a_i`, so
`f_i=a_i+2F_i`. For a finite distinct representation define

    Q=sum e_i v_i=(A,B),          n=A+2B.

By (Z2)--(Z3), Q is unchanged under every legal rewrite, hence is independent
of the choice of relaxed distinct representation. It has an exact closed
form. Put

    phi=(1+sqrt5)/2, alpha=1/phi^2=(3-sqrt5)/2,
    b(n)=floor(alpha(n+1)).

Then

    Q(n)=(n-2b(n),b(n)).                                  (Z4)

Here is a proof independent of rewrite confluence. Let `r=-1/phi`. The
recurrence and the two initial values give

    alpha f_i-F_i=alpha r^i.

For any finite binary digits, the sum of selected powers `r^i` is strictly
between the sum of all negative powers and the sum of all positive powers:

    -1 < sum e_i r^i < phi.

The bounds follow from the two geometric series and are strict for finite
support. Multiplying by alpha gives

    -alpha < alpha n-B < 1-alpha,
    0 < alpha(n+1)-B < 1.

Thus B is exactly the floor in (Z4). This proof applies to every distinct
representation, including consecutive terms and the empty representation
of zero. No floating-point arithmetic is needed in the experiment:

    b(n)=(3m-isqrt(5m^2)-1)//2,           m=n+1.           (Z5)

Since `sqrt(5m^2)` is irrational, its integer floor brackets the expression
strictly; considering the parity of `3m-isqrt(5m^2)` gives (Z5).

More generally, additive atom weights in any abelian group preserve every
distinct-representation fibre **iff** they satisfy the Fibonacci recurrence.
Necessity follows from the representations `f_(i+2)` and `f_i+f_(i+1)`;
sufficiency follows from (Z2). Therefore two initial group values describe
the entire family of such invariants. This is a two-register statement,
not an additive homomorphism on all integers: distinct representations are
not closed under ordinary addition when their supports overlap.

## 3. The precise three-color convention

Reduce (Z3)--(Z4) modulo two. Let

    R=(1,0), G=(0,1), B=(1,1) in V=F_2^2,
    c(n)=Q(n) mod2.                                      (Z6)

The atom colors are `R,G,B,R,G,B,...` because their sums use XOR:
`R+G=B`, `G+B=R`, `B+R=G`. Every relaxed distinct Fibonacci decomposition
has the same XOR color. The presence of multiple such decompositions does
not break this invariant.

There is a fourth, neutral state `(0,0)`. For example `6=5+1` has color
`R+R=0`. Thus this is three nonzero atom colors in a four-state charge
space, not a proper three-coloring of every integer. Treating the neutral
state as absent discards actual arithmetic.

This convention is also not additive coloring by `Z/3Z`. A three-periodic
sequence of `Z/3Z` weights satisfying the Fibonacci recurrence must be
identically zero: if the period is `(a,b,c)`, then `c=a+b`, `a=b+c`,
`b=c+a`, hence `2a=2b=0`. In characteristic three that forces a=b=c=0.
The nonzero cyclic XOR colors work precisely because their ambient group
has characteristic two.

## 4. Repeated atoms and the exact diagonal defect

For repeated coefficients `d_i`, let `Q_raw=sum d_i v_i=(A_raw,B_raw)`
and `n=sum d_i f_i`. Define the integer defect

    D=B_raw-b(n).

Then exactly

    Q_raw-Q(n)=(-2D,D).                                  (Z7)

Consequently only the G coordinate can change modulo two, and the single
bit `D mod2` repairs it. The minimal hostile is

    2: singleton f_1 has color G,
       repeated f_0+f_0 has color R+R=0.

If an additive atom invariant is required to survive **all** repeated
representations, it must additionally satisfy `w_1=2w_0`. Together with
the Fibonacci recurrence this gives `w_i=f_i w_0`. Conversely those weights
clearly survive every representation. Over V this loses the independent
G channel and cannot retain the three distinct atom colors.

For a pair of integer summands a,b, write

    delta(a,b)=b(a)+b(b)-b(a+b).                          (Z8)

The repeated representation obtained by concatenating their canonical
supports satisfies

    Q(a)+Q(b)-Q(a+b)=(-2delta,delta),
    c(a) XOR c(b) XOR (0,delta mod2)=c(a+b).             (Z9)

The exact two-summand carry has only three possible values:

    delta(a,b) in {-1,0,1}.                              (Z10)

Indeed, if u and v are the fractional parts of `alpha(a+1)` and
`alpha(b+1)`, then `delta=-floor(u+v-alpha)`. Since
`-alpha<u+v-alpha<2-alpha`, the floor is -1,0,or1.
All three cases occur: `(a,b)=(1,1),(1,2),(2,2)` give -1,0,1 respectively.

This pinpoints the diagonal failure. If a+b=n, the naive cell color
`c(a) XOR c(b)` can change along the diagonal. The first strict-summand
split is already n=5:

    1+4: R XOR G=B,         2+3: G XOR B=R=c(5).         (Z11)

The canonical representation of 4 is `3+1`, so `1+4` repeats the atom 1.
Its defect is -1. In contrast `2+3` uses disjoint atom supports, and its
defect is zero. Distinct *integer summands* therefore do not guarantee a
distinct combined *Fibonacci representation*.

If the canonical supports of a and b are disjoint, their union is a
relaxed distinct representation, so delta=0. Overlap is not itself a
complete classifier; retaining the exact carry makes the statement iff.
Along any fixed diagonal the first color coordinate is fixed by n mod2,
so at most two naive colors occur; the carry bit makes all cells agree
with c(n).

Regrouping has no residual path ambiguity once the carry is retained:

    delta(a,b)+delta(a+b,c)
      =delta(b,c)+delta(a,b+c).                          (Z12)

Both sides telescope to `b(a)+b(b)+b(c)-b(a+b+c)`. Thus the accumulated
correction is independent of the addition tree. The uncorrected statistic
was not a well-defined arithmetic sum color; the repaired statistic is.

## 5. The shared order-three action, and the distinct carries

The Fibonacci companion matrix on V is

    M(a,b)=(b,a+b),     M^3=I.                           (Z13)

It fixes neutral zero and cycles R->G->B->R. On actual integers, shifting
every canonical Fibonacci atom up one index gives the exact injection

    S(n)=2n-b(n),       c(S(n))=M c(n).                  (Z14)

Its integer values grow; the order-three color quotient does not make
`S^3(n)=n`. This separates the action on colors from arithmetic evolution.

There is now a precise link to the recovered geometric theorem. In
THM-3339, the golden Euclid pair is `u_k=(F_k,F_(k+1))`; our atom charge
is `v_i=u_(i-1)` (with the two small boundary indices interpreted by
(Z3)). Thus M is exactly its Fibonacci matrix reduced modulo two. The
three colors record odd/even, even/odd, and odd/odd parameter parity.
For i>=3, color B is exactly the content-two case of the raw Euclid
triple; normalization divides by two and changes the leg chamber.

THM-3339 equations (18)--(19) place the normalized triples on the true
geometric rays `(BA)^r`, `A(BA)^r`, and `C(BC)^r`, in that theorem's
branch convention. It explicitly does not produce a cyclic action
permuting the three Berggren children. Those golden rays differ from
the parabolic inverse-fibre sampling ray in the current ternary-clock
note. At non-Fibonacci labels, Q(n) need not be a coprime golden pair:
`Q(6)=(2,2)` is not primitive, and `Q(7)=(1,3)` has golden norm 5,
not plus or minus one. Three-color information cannot select the
golden geometric locus by itself.

The carry correction also selects a frame. Modulo two, the vector
`(-2,1)` is G; the integer weights 1 and 2 distinguish this channel.
If all colors are rotated by M, the correction channel must rotate
from G to B as well. Keeping G fixed breaks even the diagonal1+4
identity. In fact `S(a)+S(b)-S(a+b)=-delta(a,b)`, so the numeric
Fibonacci shift is not an automorphism of integer addition. This is
compatible with THM-3339's moving-owner warning: a quotient channel
action does not preserve all framed current operations. Its branch
owner cocycle and our overlap cocycle are different objects; neither
is identified with the other without an additional map.

There is an explicit, legitimate map from a ternary residue digit to a
nonzero color: `0->R,1->G,2->B`. It intertwines digit addition by one
modulo three with M. Applied letter by letter to least-significant-first
ternary words, it is an isomorphism of rooted address trees. Increment
acts by M on the first color, carrying to the next position exactly when
the old color was B. This is the full finite ternary odometer, with its
carry, in a three-color alphabet.

The connection to the proved residue clock is explicit. Fix
`c=2^(k0-1)` and let `H(j)=c(4^j-1)/3`. The inherited clock theorem
proves that `H` is a bijection modulo every `3^d`, compatible with
reduction, and `H(j+1)=4H(j)+c`. For a height residue h, first compute
its unique phase `j=H^(-1)(h) mod3^d`, then color its d ternary digits.
This composite decoder intertwines `h->4h+c` with the colored odometer
just described. It retains the entire phase word, not a single
Fibonacci charge. This is a concrete commuting square between the
residue-tree action and the color-word action.

It is an isomorphism of actions, not of additive groups: `0 mod3` maps
to R, while `R XOR R=0` is neutral and outside the digit alphabet. Nor is
the single Fibonacci charge c(n) itself this ternary address. It fails to
factor through any fixed residue quotient `n mod3^d`, since n and
`n+3^d` have opposite parity and hence different first coordinates in
c(n). A finite residue-tree node needs the ordinary integer/carry sidecar
if it is also to carry the Fibonacci color.

Thus the exact bridge has three levels: the common cyclic action M; an
equivariant encoding of ternary address words; and the separately computed
Fibonacci charge of their ordinary integer labels. The three-adic carry
of the odometer and the overlap carry delta are different operations.

## 6. The ten-state pair object and its color boundary

The root synthesis identifies an exact ten-state object:
`Sym^2(V)`, unordered pairs with repetitions allowed. Under M it has one
fixed point `{0,0}` and three free orbits of size three. Its seven-state
part consists of the six distinct pairs together with `{0,0}`; the remaining
three states are `{R,R},{G,G},{B,B}`. This is exactly the repeated-atom
diagonal channel that distinguishes distinct from repeated objects above.

The pair-XOR observable sends the four repeated pairs to zero, while each
nonzero charge occurs on two pairs. The abstract identification with the
ten proper divisors of `p^2 q r` is given in the companion prime/root
notes: prime labels are assigned virtual colors, not their numerical
Zeckendorf colors. The two assignments are generally different. An
equivariant ten-letter action therefore does not automatically preserve
integer sums, differences, divisibility, or numerical Fibonacci carries.

Even on this ten-state set there are distinct equivariant observables.
Assigning a lifted-prime letter its prime's virtual color differs from
pair-XOR, which makes its repeated pair neutral. Channel data must identify
which observable is being used; equality of the underlying group action
does not identify their charge-counting formulas.

## 7. What this transports to sums and differences

The useful source is a labeled pair of integers with its canonical support
or exact charge; the target is an additive edge label. Formula (Z9)
computes the color of a sum with one carry bit. For b>a and d=b-a, the
same formula gives

    c(a) XOR c(b) XOR (0,delta(a,d) mod2)=c(d).           (Z15)

So both sum and difference labels have an exact local decoder using the
same carry mechanism. This is a concrete shared arithmetic structure for
square-sum edges and graceful difference labels. The full charge retains
the integer because n=A+2B; its color quotient loses that information.

The quotient does not recognize square targets: the square/nonsquare
pairs `1/5`, `4/2`, `9/3`, and `16/6` have the same colors R,G,B,0
respectively. It also does not certify that a tree's n-1 differences are
the distinct integers `1,...,n-1`. Hamiltonian existence and graceful
labeling therefore require additional global constraints. The precise
repair here is local representation consistency, not an implication
between those conjectures.

## 8. Exact reproduction

    python 04-computation/experiments/duck_zeckendorf_20260925.py
    python -O 04-computation/experiments/duck_zeckendorf_20260925.py

The standard-library script retains all checks under `-O`; see
[the frozen output](duck_zeckendorf_20260925.out). Its universes are:

- All integers 0..100000, checking the floor formula against independently
  summed greedy charges and checking the exact Fibonacci shift.
- All 65,536 distinct subsets of the first 16 weights: 4,180 integer
  values and all 131,072 legal one-step rewrites preserve charge.
- All rewrite orders from each of the 4,096 subsets of the first 12
  weights terminate at the same canonical form.
- All `4^7=16384` repeated representations using multiplicities 0..3 on
  the first seven weights pass the one-bit correction, including the
  repeated-atom hostile.
- Every strict sum pair on diagonals n<=512: 65,280 pairs, all exact
  carry values, all 26,106 disjoint-support controls, and the first split
  at n=5. All `65^3` triples in 0..64 pass the cocycle identity.
- Every residue in ternary depths 1..7 passes the color-word odometer;
  each scalar residue quotient fails as predicted. The explicit
  height/phase/color commuting square passes on1,452 cases with
  `k0=2..5` and depths1..5. Square-filter and
  characteristic-two/three hostiles are checked explicitly.
- Atom indices3..49 match the inherited golden norm and odd/odd
  content-two seam; an explicit color rotation tests that the carry
  channel must rotate with its frame.

The universal claims are proved above. These finite controls do not turn
the color correspondence into a square-sum, graceful-tree, or Collatz proof.

Independent proof audit checked the binary-subset floor formula, both
parities of its integer-square-root implementation, the carry range and
cocycle, repeated-atom repair, and the first strict diagonal split, with
no mathematical findings. The root review corrected the printed rewrite
direction convention; (Z2) is explicitly written high to low.
