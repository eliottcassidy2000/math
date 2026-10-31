# Session129: actual triangular/Fibonacci bridges for PRIMEGAME

## Findings, with the mathematical types kept separate

1. Every PRIMEGAME **integer state** has an injective harmonic/Fibonacci-recurrence lift whose third harmonic edge is twice a triangular number. It conjugates every enabled PRIMEGAME step exactly.
2. Exponent parity gives an honest tournament observer after a declared assignment of ten prime registers to the ten **edge slots** of K5. Passing to isomorphism and converse gives ten nodes, but the observer is not autonomous. The least counterexample in the stated gauge is `1,4`.
3. The canonical harmonic lift of an **instruction fraction** gives a Fibonacci quadratic norm. It proves uniform triangular identities, and there is one genuine golden-unit pairing of the fourteen instructions: `11/13` with `15/2`. The proposed reciprocal pair `11/13,13/11` is not that pairing.
4. None of this makes the proposed `2+4+4+4` instruction blocks determine ten merged nodes. No new adjacency correspondence, CA law, or finite-state PRIMEGAME quotient is asserted.

All statements below are proved algebraically or by exhaustive finite enumeration, not extrapolated from numerology. The code uses integer/Fraction arithmetic and explicit exception gates, identically under ordinary Python and `python -O`.

## 1. Two uses of the same harmonic lift

For positive coprime integers `a,b`, compactify the positive fraction by

`a/b -> a/(a+b)`.

The THM-3509 seed is `(x,y)=(b,a)`, giving the recurrence window

`W=(b,a,a+b,b+2a)`

and harmonic triangle

`H=(u,v,z)=(ab,b(a+b),a(a+b))`.

Indeed `vz=u(v+z)`, equivalently `1/u=1/v+1/z`. Coprimality gives the exact decoder

`(gcd(u,v),gcd(u,z),gcd(v,z))=(b,a,a+b)`.

Thus this is injective, and genuinely preserves the original reduced fraction. The word *Fibonacci* here means the recurrence `w2=w0+w1, w3=w1+w2`, not necessarily consecutive terms of the distinguished Fibonacci sequence.

There are two differently typed inputs:

- An instruction `a/b` gives the window above. There are fourteen of these.
- A state `N=N/1` gives the special U-spine window `(1,N,N+1,2N+1)` and face

  `H(N)=(N,N+1,N(N+1))=(N,N+1,2T_N)`, where `T_N=N(N+1)/2`.

In the second case the gcd decoder recovers `(1,N,N+1)`, hence N, hence its full prime factorization. This is not a quotient and does not discard unbounded exponent information.

## 2. An exact triangular dynamical conjugacy

Let the selected PRIMEGAME rule be `a/b`. Its original divisibility guard is `b|N`, including the first-enabled priority rule. Put `M=aN/b`. Then the transported operation is

`(N,N+1,2T_N) -> (M,M+1,2T_M)`.

This is an exact commuting diagram because `H` is injective and the same priority is retained. An explicit triangular-coordinate update is

`T_M = (a^2/b^2) T_N + a(b-a)N/(2b^2)`.

Proof: substitute `M=aN/b` into `M(M+1)/2` and use `N^2=2T_N-N`. Dropping the linear N term would be false. If one stores only the triangular number, N is still recovered on its image by

`N=(sqrt(8T_N+1)-1)/2`.

Thus no approximation or irrational arithmetic is required; the square root is integral on the image. All guards can be recovered, but they have not become local or finite-state.

For example, the first state update is `2 -> 15`, and the lifted triangular update is `T_2=3 -> T_15=120`. The numerical label 55 should not be confused with its index: a state `N=55` lifts to `T_55=1540`, not `T_10=55`.

## 3. The exact Fibonacci boundary

The antipodal difference in the harmonic K4 construction is

`kappa=e03-e12=x^2+xy-y^2=:Q(x,y)`.

For the U-spine state lift, `Q(1,N)=1+N-N^2`. The equation `|Q(1,N)|=1` with positive integral N has exactly the solutions `N=1,2`:

- `Q=1` gives `N(N-1)=0`, hence N=1;
- `Q=-1` gives `(N-2)(N+1)=0`, hence N=2.

These produce `(1,1,2,3)` and `(1,2,3,5)`. Already the first PRIMEGAME step `2 -> 15` changes kappa from `-1` to `-209`. Therefore PRIMEGAME preserves the generalized recurrence-window carrier, not the unit-Cassini/golden ray. Ranking the four window entries or the six edge products gives only a transitive tournament away from the N=1 tie; this construction does not produce the ten n=5 classes.

## 4. A uniform triangular-to-Fibonacci-norm bridge

Let `phi=(1+sqrt(5))/2`. The algebraic norm satisfies

`Norm(x+y phi)=x^2+xy-y^2=Q(x,y)`.

For an eraser instruction `1/q`, its compactified harmonic seed is `(q,1)`, so

`kappa(1/q)=q^2+q-1=2T_q-1`.

This is a uniform theorem, not an isolated numeric observation. It has an exact inverse gate: a positive integer k occurs as this norm for a positive integer q iff `4k+5` is an odd square at least 9; then

`q=(sqrt(4k+5)-1)/2`.

For the printed instructions,

- `kappa(1/7)=55=F_10=T_10`, which is also the restart numerator `55=5*11`;
- `kappa(1/17)=305`, so the two erasers do not both produce listed restart factors;
- `kappa(13/11)=95`, a numerator elsewhere in the rack.

These are verified specializations of the fixed map. They do not classify Fibonacci norms or imply a hidden PRIMEGAME transition. In particular, occurrence of 305 elsewhere as a cubic continued-fraction coefficient supplies no map from this quadratic norm to that cubic process.

## 5. The unique golden-unit pairing among the fourteen instructions

Multiplication by phi and Galois conjugation act on primitive integer seeds by

`F(x,y)=(y,x+y)`, `C(x,y)=(x+y,-y)`.

Direct expansion proves `Q(Fv)=-Q(v)` and `Q(Cv)=Q(v)`. Both maps are unimodular, so primitivity is preserved, and absolute norm is invariant under the group they generate. Intermediate coordinates need not be positive. The same action may be read projectively on rational ratios; a primitive seed is determined up to simultaneous sign.

The fourteen instruction norms, in the given order, are

`9539, 7771, 3209, 1789, 1205, -2855, -6311, -4105, 305, 191, 95, -191, 55, -2969`.

There are thirteen distinct absolute values. The only repeated value is 191, from `11/13` and `15/2`. Necessary equal norm is sufficient for this pair, because

`F^3 C(x,y)=(x-y,2x-y)`

sends `(13,11)` to `(2,15)`. In fraction coordinates this is the involution

`r -> (2-r)/(1-r)`,

sending `11/13` to `15/2`. Its seed square is `-I`, and it reverses Q. Thus the fourteen printed instructions meet **exactly thirteen** Fibonacci-unit/Galois orbits.

This is a true algebraic relationship not asserted to be a program symmetry. The rack is not closed even under F: `1/7` is sent to `8/1`, which is absent. The reciprocal pair `11/13,13/11` has absolute norms 191 and 95, and therefore belongs to different such orbits. Consequently this golden action cannot explain the proposed `2+4+4+4` partition.

## 6. The n=5 triangular carrier that actually produces ten nodes

The number of edge slots in an n-vertex tournament is `C(n,2)=T_(n-1)`. A fixed directed Hamiltonian path forces n-1 orientations and leaves

`C(n,2)-(n-1)=C(n-1,2)=T_(n-2)`

independent off-path slots. These are different carriers from the quotient set of tournament types.

Independent exhaustive enumeration gives:

| n | edge slots | fixed-path free slots | labelled tournaments | isomorphism classes | self-converse classes | classes modulo converse |
|---|---:|---:|---:|---:|---:|---:|
|3|3|1|8|2|2|2|
|4|6|3|64|4|2|3|
|5|10|6|1024|12|8|10|

For the last column, converse acts as an involution on isomorphism classes, so the count is `(classes+self_converse)/2`. At n=5 this is `(12+8)/2=10`. At n=3 the transitive tournament and directed triangle are both self-converse, so the count is 2, not 1. Thus `1,3,10` mixes carriers unless a different quotient convention is supplied.

This enumeration checks every bitmask and every vertex permutation. The n=5 merged Hamiltonian-path counts, in canonical-representative order, are

`1,3,5,3,9,9,11,15,13,15`.

Neither these ten classes nor their H values follow from identifying their cardinality with the number of prime registers.

## 7. A genuine prime-exponent-to-tournament observer

Declare the following edge gauge: assign primes in increasing order `(2,3,5,7,11,13,17,19,23,29)` to the lexicographically ordered pairs

`01,02,03,04,12,13,14,23,24,34`.

For a supported state `N=prod p_i^e_i`, orient edge i forward exactly when `e_i` is odd. Its bitmask is

`B(N)=sum (e_i mod2)2^i`.

Then pass to `O(N)=[B(N)]` modulo vertex relabeling and global edge reversal. The map is surjective onto all ten merged classes: every labelled mask is realized by a squarefree supported integer. This is an explicit factorization-to-node map on **states**.

For rule a/b, let `alpha,beta` be numerator/denominator exponent vectors. Whenever this rule is selected,

`B(aN/b)=B(N) XOR parity(alpha-beta)`.

The printed parts are squarefree, so this toggle is exactly their union of prime supports (numerator and denominator are disjoint). The law remains valid with non-squarefree parts if parity vectors, rather than set unions, are used. The checker verifies it for all 1024 binary-exponent inputs and also checks the full enabled integer/exponent update and harmonic transport.

This is not a canonical register-to-node bijection. Ten registers label the ten **edges**, and the ten node classes arise only after using all 1024 configurations and taking a quotient. Indeed, in this gauge, the ten prime singleton states occupy only five merged classes. A bijection from register names to edges is declared data; distinct edge assignments generally are not induced by vertex permutations. Under prime renaming, the observer is transported by transporting its gauge as well.

## 8. Smallest hostile: the observer has no autonomous update

Take the two smallest states with equal exponent parity:

`N=1` and `N=4=2^2`.

Both have mask 0 and the transitive merged class. But priority gives

`1 -> 55=5*11` (rule14), `4 -> 30=2*3*5` (rule12).

Their next masks are 20 and 7. In the declared gauge, mask20 has canonical merged representative4 and H=5, whereas mask7 is transitive, canonical0 and H=1. Thus there is no function f on these ten observer nodes satisfying `O(F(N))=f(O(N))` for every N. The failure already occurs before taking the quotient because 20 differs from7.

The witness is smallest when ordered first by the larger initial state and then the smaller: below4, states1,2,3 have distinct parities; the only equal-merged pair among them is1,2, whose next states55 and15 have the same merged class. The checker searches every supported integer up to100 and pins the resulting minimum `(1,4)`.

The conditional XOR identity is still exact. It fails to close because parity discards divisibility and priority information, not because factorization ceases to be additive. This observer is a useful visualization/measurement interface; retaining full exponents supplies a faithful extension.

## 9. Numeric arithmetic versus prime-renaming invariance

Consistently relabeling the ten primes conjugates the FRACTRAN program exactly. The checker transports the entire 1024-state binary box to primes `(31,37,41,43,47,53,59,61,67,71)` and verifies the selected instruction and successor exponent vector.

The literal restart numerator55 becomes `41*47=1927`, which is neither triangular nor Fibonacci. The eraser prime7 becomes43, whose quadratic norm is `43^2+43-1=1891`, not1927. Therefore the numeric norm hits are real algebraic identities of the printed constants, not invariants of the exponent program. The exact state harmonic lift and parity observer remain meaningful, but their numerical coordinates or chosen gauge must be transported appropriately. No claim that arbitrary register renaming preserves the harmonic numeric values is made.

## 10. Provenance and reproduction

Only targeted `git show` reads were used; source paths and Git blob identities are recorded in SOURCE_PINS.json. Main inherited inputs are THM-3509, Session126's exact PRIMEGAME report and independent n=5 checker, and the fixed-Hamiltonian-path triangular carrier reflection. The new checker imports none of them.

Run from any directory:

```sh
python /path/to/check_figurate.py > normal.json
python -O /path/to/check_figurate.py > optimized.json
cmp normal.json optimized.json
```

Finite checks supplement the symbolic proofs in Sections1–5,7–9; the small tournament classification itself is finite exhaustive proof. All gates use exceptions, not removable assertions. No external publication or repository mutation is performed by this artifact.

`verify_release.py` also relocates the standalone checker into a fresh temporary directory, requires byte-identical ordinary/optimized output, and rejects six deliberately damaged variants: wrong Cassini cross term, wrong triangular index, presence instead of parity, replacing15/2 by15/14, dropping converse, and breaking the recurrence window. The fraction substitution is a hostile mutation, not a faithful implementation of the other published PRIMEGAME list. See REPLAY.json.
