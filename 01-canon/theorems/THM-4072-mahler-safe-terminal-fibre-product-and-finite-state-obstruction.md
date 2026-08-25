---
id: THM-4072
title: "Mahler safe-terminal fibre product and finite-state obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.  The strict Mahler
  carry condition is the reset-Buechi part of an exact countable follower
  graph, ordinary positive stabilization is a native-bit co-Buechi condition
  plus a positivity-seen reachability marker, and the coordinate change is
  a rooted-tree automorphism with an explicit unbounded integer sidecar.  Their
  exact fibre product gives an iff characterization of candidate Z-number
  codes, but every finite local terminal-prefix test is vacuous.  The rooted
  safe pair-prefix language is nonregular and the coordinate change has no
  finite synchronous transducer.
  None of these statements proves that the candidate language is nonempty or
  calls that unknown language nonsofic; Mahler's Z-number problem remains open.
source: codex-frontier-synthesis-creative-20260825c / Mahler fibre-product lane
audit: >
  PASS.  A carry-first replay exhausts all 131,070 nonempty binary words
  through depth 16, checks the suffix inequalities against the follower graph,
  proves every level map is a bijection by both affine and orbit identities,
  recovers safe counts 1,2,3,5,8,12,18,27,40,60,90,134,201,302,452,678,1018,
  sees exactly m+1 follower states at every depth m, and performs 2,752,470
  exact radix-clock descent gates.  A genuinely independent native-first
  replay enumerates ordinary binary prefixes, reconstructs carries by the
  ceiling orbit, tests safety lexicographically, independently recovers the
  same counts, separates all 136 pairs of depth-16 follower states, and
  performs 545,289 clock gates.  Both replays check the greedy equality ray,
  (100)^infinity, A=1, (01)^infinity, and the denominator-19 changing-start
  tower at horizons 18,36,54,72.  They use exact integers and rational numbers
  only, contain no Python assert nodes or floating constants, and normal and
  optimized runs byte-match their frozen outputs.
depends_on:
  - THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization
  - THM-2231-relation-carry-completion-and-exact-radix-clock
  - THM-2352-q-adic-prefix-residue-collision-spectrum
  - THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation
script: 04-computation/mahler_safe_terminal_fibre_product_thm4072.py
output: 05-knowledge/results/mahler_safe_terminal_fibre_product_thm4072.out
script_sha256: 6e07652a114cd259a8574c5ed9357c681bd6015e0c72e090a41b8089d4420a97
output_sha256: 3e6299922e2e04059ccbd7b464d06e26736cdee669ee4308a2872d6ce7c8c11f
independent_audit_script: 04-computation/mahler_safe_terminal_fibre_product_thm4072_independent_audit.py
independent_audit_output: 05-knowledge/results/mahler_safe_terminal_fibre_product_thm4072_independent_audit.out
independent_audit_script_sha256: 818ab917ef4e8e0268bf4821ec529f35efa21c09c67cedcbe2a89e5e45d2b25f
independent_audit_output_sha256: 019ed652a342108da1f31eb7f04bbdf88953fbc283cfdab0905052351ea891cc
hash_basis: raw working-tree bytes (LF)
---

# THM-4072 -- the Mahler safe-terminal product is exact but not finite-state

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.**  This theorem
joins the two coordinates separated by THM-2228.  The join is lossless and
gives an exact candidate-code theorem.  It does not decide whether any such
code exists.

## 1. The strict safe language is a reset-Buechi language

Use the notation of THM-2228 and THM-3848.  For

```text
c=c_0c_1... in {0,1}^N
```

put

```text
Y_i(c)=sum_(j>=0)c_(i+j)(2/3)^(j+1),
K={c:Y_i(c)<=1 for every i},
S={c:Y_i(c)< 1 for every i}.                           (1)
```

Let `d=d_1d_2...` be the greedy expansion of one,

```text
z_0=1,
d_k=floor((3/2)z_(k-1)),
z_k=(3/2)z_(k-1)-d_k,
1=sum_(k>=1)d_k(2/3)^k.                                (2)
```

Define the countable labeled follower graph `Q` with vertices
`q=0,1,2,...`, initial vertex zero, and the following transitions:

```text
q --d_(q+1)--> q+1;
q --0--> 0             when d_(q+1)=1.                 (3)
```

The second kind is called a **reset**.  When `d_(q+1)=0`, an input `1` has no
outgoing edge and is rejected.

> **Follower/Buechi theorem.**  A finite word labels a path from zero in `Q`
> if and only if it belongs to the strict truncated safe language `P_m` of
> THM-3848.  An infinite word labels a nonrejecting path if and only if it is
> in `K`.  On such a path,
>
> ```text
> c in S  iff  reset edges occur infinitely often.     (4)
> ```

The finite assertion is exactly the unique first-disagreement renewal
decomposition in THM-3848: matching `d` advances the state, a first downward
disagreement can occur only at a `1` of `d` and returns to a fresh safe suffix,
and a first upward disagreement is forbidden.  Passing to the inverse limit
gives the assertion for `K`.

If a path has only finitely many resets, then after its last reset it starts at
state zero and follows `d` forever; the no-reset case follows `d` from the
start.  Hence it lies in `K\S`, whose exact form in THM-3848 is

```text
K\S={c:sigma^n c=d for some n>=0}.                     (5)
```

Conversely, suppose a nonrejecting word has a suffix equal to `d`.  The
lexicographic maximum of the follower set at state `q` is `sigma^q d`.  For
`q>0`, shift invariance gives `sigma^q d<=d`, and equality would make `d`
eventually periodic.  THM-3848 proves it is not, so `sigma^q d<d`; consequently
the word `d` itself cannot be read from state `q>0`.  The equality suffix must
therefore begin at state zero and has no later reset.  This proves (4),
including the easy-to-miss reverse direction.

Every state `0,...,m` occurs at depth `m`: state `n` is reached by

```text
0^(m-n)d_1...d_n.                                      (6)
```

The states have distinct follower sets, with respective lexicographic maxima
`sigma^n d`.  Thus the countable state is mathematical information, not a
presentation artifact.

### The oriented strict/equality ledger

For a suffix beginning at `i` in a word of length `m`, set `h=m-i` and

```text
C_(i,h)=sum_(0<=j<h)c_(i+j)2^j3^(h-1-j).               (7)
```

The actual Mahler fractional phase `x_i=Y_i/2` of any infinite completion lies
in the sharp oriented interval

```text
x_i in [C_(i,h)/3^h, C_(i,h)/3^h+(2/3)^h].             (8)
```

The finite language uses only

```text
2C_(i,h)<3^h,                                          (9)
```

whereas an infinite strict suffix has some open-tail certificate

```text
2C_(i,h)+2^(h+1)<3^h.                                 (10)
```

Equation (10) is necessary and sufficient after allowing `h` to grow, by
THM-2228.  The equality ray `d` satisfies `Y_0(d)=1`, hence has the oriented
phase `x_0=1/2`.  Thus `Q` retains the strict/equality boundary through its
unbounded follower state and reset acceptance, while (8) is the quantitative
phase sidecar.  Neither may be replaced by a centered distance.

## 2. Carry coordinates and native binary coordinates are rooted-tree isomorphic

For a length-`m` carry word put

```text
C_m=sum_(0<=j<m)c_j2^j3^(m-1-j),
r_m=-3^(-m)C_m mod 2^m,       0<=r_m<2^m,              (11)
u_m=(3^m r_m+C_m)/2^m.                                 (12)
```

Set `r_0=u_0=0`.  Compatibility gives a unique digit `e_m in {0,1}` such that

```text
r_(m+1)=r_m+e_m2^m.                                    (13)
```

Then the exact carry-driven update is

```text
e_m=c_m xor (u_m mod 2),
u_(m+1)=(3u_m+c_m+3^(m+1)e_m)/2.                       (14)
```

Equivalently, when the native digit is supplied first,

```text
c_m=(u_m+e_m) mod 2.                                   (15)
```

Indeed, `C_(m+1)=3C_m+c_m2^m`.  Substitute (12)--(13) into the corresponding
identity at level `m+1`.  Divisibility by `2^(m+1)` says precisely

```text
e_m=u_m+c_m mod 2,
```

and division by two gives (14).  Also, THM-2228's affine orbit identity says

```text
u_m=T^m(r_m),              T(a)=(3a+(a mod 2))/2.      (16)
```

So `u_m` is both the algebraic quotient in (12) and the native integer reached
after reading the carry prefix.

Let

```text
theta_m:{0,1}^m -> {0,1}^m,
theta_m(c_0...c_(m-1))=e_0...e_(m-1).                  (17)
```

THM-2228 proves that `c -> r_m` bijects carry words with residues modulo
`2^m`; (13) is the ordinary binary expansion of that residue.  Hence every
`theta_m` is a bijection.  The maps commute with last-letter truncation, so
their inverse limit

```text
theta:{0,1}^N -> {0,1}^N                               (18)
```

is an automorphism of the rooted binary tree and a homeomorphism of its path
space.  It is **not** asserted to commute with the shift.  Its required online
sidecar is the generally unbounded integer `u_m`.

## 3. The native plateau clock and the exact fibre product

For a full native word `e`, the `1` positions are precisely the active
positions of THM-2352: `r_m` changes there and is constant across intervening
zero plateaux.  If an ordinary completion of the level-`m` cylinder is

```text
A=r_m+2^m k,                    k>=0,                  (19)
```

then its next native digit is `e_m=k mod 2` and its next remaining quotient is

```text
k'=(k-e_m)/2=floor(k/2).                               (20)
```

With the THM-2231 radix clock

```text
h_m=bit_length(k),              bit_length(0)=0,       (21)
```

one has `h_(m+1)=h_m-1` whenever `h_m>0`.  Thus ordinary stabilization is
exactly the co-Buechi condition that `e_m=1` only finitely often.  Stabilization
to a **positive** integer is

```text
E_+={e:e_m=1 finitely often but at least once}.         (22)
```

Define the synchronized graph

```text
G={(c,e):e=theta(c)}.                                  (23)
```

Then the exact carry/native code language for positive Mahler Z-numbers is

```text
Z_code
 ={(c,e) in G:
     the Q-path never rejects,
     reset edges occur infinitely often,
     e_m=1 only finitely often,
     e_m=1 at least once}
 ={(c,e) in G:c in S and e in E_+}.                    (24)
```

This is an **if and only if**.  The first two gates are exactly `c in S` by
Section 1.  The last two are exactly positive ordinary stabilization of
`Phi(c)` by (13) and (22).  THM-2228 then gives the unique number

```text
A=sum_(m>=0)e_m2^m,
alpha=A+Y_0(c)/2.                                      (25)
```

Conversely, every positive Mahler Z-number supplies these same coordinates
and satisfies all four gates.  The words "Buechi" and "co-Buechi" in (24)
describe the two acceptance predicates; they do not turn the unbounded graph
`G` into a finite automaton.

## 4. Finite local terminal tests prune nothing

The distinction between a local relaxation and an actual candidate projection
is essential.  Define the **naive level product**

```text
F_m={(w,e):w in P_m, e=theta_m(w),
             e is the prefix of some word in E_+}.     (26)
```

Every binary word `e` is a prefix of a word in `E_+`: append only zeros if a
`1` has already occurred, and otherwise append one `1` and then only zeros.
Therefore

```text
F_m={(w,theta_m(w)):w in P_m},
|F_m|=|P_m|=a_m,                                       (27)
```

and both projections in (27) are bijective onto their images.  In particular,
the local terminal-prefix condition removes **zero** safe nodes at every
finite level.  The first counts are

```text
a_0,...,a_16=
1,2,3,5,8,12,18,27,40,60,90,134,201,302,452,678,1018. (28)
```

Because the `theta_m` are compatible,

```text
inverse_limit F_m={(c,theta(c)):c in K}.                (29)
```

This is not (24): finite approximation closes `S` to `K` and closes `E_+` to
the full binary path space.  The equality pair `(d,theta(d))` and the
nonterminal safe pair over `(100)^infinity` are explicit points retained by
(29) but excluded from (24).  Moreover, `F_m` is **not** claimed to equal the
set of length-`m` prefixes of (24); that latter set could be empty, because the
existence of a Z-number is open.

The same obstruction is visible in the clock.  For a fixed safe prefix `w`
with residue `r_m`, every integer (19) has that carry prefix.  If `r_m>0`, one
may realize clock height zero with `k=0`, and every height `h>=1` by choosing
`2^(h-1)<=k<2^h`.  If `r_m=0`, positivity removes only `h=0`; every `h>=1`
still occurs.  These are terminal completions of the finite cylinder, not
claims of infinite safety.  Thus no finite follower/residue state bounds the
remaining terminal clock.

## 5. Two finite-state impossibility theorems

Let the finite rooted pair-prefix language be

```text
L_pair={((c_0,e_0),...,(c_(m-1),e_(m-1))):
        m>=0, c_0...c_(m-1) in P_m,
        e_0...e_(m-1)=theta_m(c_0...c_(m-1))}.          (30)
```

> **Nonregular rooted graph.**  `L_pair` is not a regular finite-word
> language.

If it were regular, its letterwise projection onto the carry coordinate would
be regular.  That projection is exactly `union_m P_m`, because every word in
`P_m` has the strict completion `w0^infinity`.  THM-3848 proves that the closed
shift `K` is nonsofic by its pairwise distinct follower sets, equivalently that
this finite language is nonregular.  This contradiction proves the claim.
The type is intentionally a **rooted pair-prefix language**: since `theta` is
not shift-commuting, no pair subshift is being declared nonsofic.

> **No finite synchronous coordinate changer.**  No finite-state synchronous
> letter-to-letter transducer realizes `theta` on every infinite carry word.

Suppose one did.  The inverse image of the omega-regular language `E_+` under
such a transducer would be omega-regular, by the standard product and
projection closure.  But

```text
theta^(-1)(E_+)=I_+,                                  (31)
```

where `I_+` is THM-2228's set of positive-integer carry itineraries.  That
theorem proves both that `I_+` is nonempty and that it contains no ultimately
periodic word.  Every nonempty omega-regular language contains an ultimately
periodic word, a contradiction.  Hence some unbounded coordinate memory, such
as `u_m`, is intrinsic.

Neither statement applies the word "nonsofic" to the unknown intersection
(24).  In particular, intersecting two individually explicit infinitary gates
does not establish either existence or a shift presentation for their
intersection.

## 6. Mandatory hostile ledger

The following five controls exercise different failed implications.

1. **Greedy equality boundary.**  Every finite prefix of `d` is in `P_m`, and
   appending `0^infinity` makes a strict member of `S`.  The full ray `d` has
   no resets and has `Y_0(d)=1`, so it is in `K\S`.  This separates finite
   strict truncation from the infinite equality endpoint.

2. **Safe nonterminal period.**  For `c=(100)^infinity`,

   ```text
   (Y_0,Y_1,Y_2)=(18/19,8/19,12/19),
   (x_0,x_1,x_2)=(9/19,4/19,6/19),
   Phi(c)=-9/19.                                       (32)
   ```

   There is a reset every third letter, so `c in S`.  Its native binary word
   is the pure least-significant-bit-first period

   ```text
   e=(101100001010011110)^infinity.                    (33)
   ```

   Indeed, the block value is
   `9(2^18-1)/19=124173`, so (33) represents
   `124173/(1-2^18)=-9/19` in `Z_2`.  Its cyclic zero plateau has length at
   most four, but it never terminates.  Strict safety therefore does not imply
   terminality, even with bounded plateaux.

3. **Terminal unsafe integer.**  The native word of `A=1` is `1000...`, while
   its carry itinerary starts `1011`.  The first three letters are safe, but
   the follower graph rejects the fourth, and their truncated tail is

   ```text
   94/81>1.                                             (34)
   ```

   Terminality therefore does not imply safety.

4. **Centered/oriented trap.**  For `c=(01)^infinity`,

   ```text
   Phi(c)=-2/5,
   e=(0110)^infinity,
   (Y_even,Y_odd)=(4/5,6/5).                           (35)
   ```

   The follower graph first rejects at length six; the offending length-five
   suffix `10101` has truncated value `266/243`.  The centered value `2/5`
   therefore hides the wrong oriented half-circle and does not repair either
   safety or terminality.

5. **Denominator-19 changing-start tower.**  For every positive multiple
   `m` of 18, put

   ```text
   A_m=9(2^m-1)/19,
   alpha_m=A_m+9/19=9*2^m/19.                          (36)
   ```

   The ordinary native digits of `A_m` agree with (33) in coordinates
   `0,...,m-1` and are zero thereafter.  Its carry itinerary agrees with
   `(100)^infinity` in the first `m` coordinates.  Thus the terminal and
   nonterminal objects have the same joint pair prefix of depth `m`.  The word
   "coordinates" is load-bearing: `A_m` has bit length `m-1`, not `m`.

   For `0<=n<=m`, the phase numerator modulo 19 follows the cycle
   `9 -> 4 -> 6 -> 9` under multiplication by `3/2=11 mod 19`.  Since
   `3^m=1 mod 38`, all these phases are in the lower half-circle, the
   index-`m` phase is `9/19`, and

   ```text
   fractional_part(alpha_m(3/2)^(m+1))=27/38>1/2.     (37)
   ```

   Hence each changing start passes the entire finite product test at its own
   horizon and fails immediately afterward.  This is a finite-horizon
   obstruction, not one infinite candidate.

## 7. Map, loss, sidecars, and exact test

```text
source:       a carry path c and its safe follower state q_m;
target:       the native binary path e=theta(c);
map:          equations (13)--(15), a compatible rooted-tree automorphism;
preserved:    exact 2-adic point Phi(c), every finite cylinder, and depth;
not preserved by e alone:
              oriented real suffix phase, strict/equality boundary, and
              safe follower state;
not preserved by q_m alone:
              native active positions, ordinary termination, and clock;
required sidecars:
              q_m plus the oriented interval (8), and u_m/r_m plus h_m;
cheapest decisive finite tests:
              rejection/upward disagreement, open-tail certificate (10),
              a native 1, or completion of a displayed clock;
irreducibly infinitary gates:
              infinitely many resets and eventually no native 1s.          (38)
```

The primary exact replay is

```text
python3 -B 04-computation/mahler_safe_terminal_fibre_product_thm4072.py
python3 -B -O 04-computation/mahler_safe_terminal_fibre_product_thm4072.py
```

and the independent native-first replay is

```text
python3 -B 04-computation/mahler_safe_terminal_fibre_product_thm4072_independent_audit.py
python3 -B -O 04-computation/mahler_safe_terminal_fibre_product_thm4072_independent_audit.py
```

Both outputs byte-match the frozen artifacts named in the front matter.

## 8. Dependency and scope audit

- **THM-2228** supplies the carry homeomorphism, affine orbit identity,
  strict-tail/positive-stabilization iff theorem, nonemptiness and absence of
  ultimately periodic words in `I_+`, and the original hostile mechanisms.
- **THM-3848** supplies the exact renewal decomposition, `K=closure(S)`, the
  equality boundary `K\S`, nonperiodicity of `d`, distinct follower sets,
  nonsoficity of `K`, and the count recurrence.
- **THM-2352** identifies native active positions and zero plateaux and warns
  that arbitrarily late active positions survive every finite cylinder.
- **THM-2231** supplies the typed exact radix clock and the reason that an end
  marker cannot be uniformly bounded.
- The new content here is the explicit transition (14), the compatible
  carry/native rooted-tree automorphism, the reset-Buechi proof including its
  reverse direction, the exact four-gate fibre product (24), the no-pruning
  inverse limit (29), the two finite-state impossibility results, and the joint
  denominator-19 hostile.

The proved consequence is an exact representation-and-obstruction theorem.
It neither constructs nor excludes a Mahler Z-number.
