# Denominator continuations: marked affine state and prime-adic cancellation balls

**Status:** elementary PROVED identities below; bounded FINITE-EXACT checks in
`checks.json`, independently implemented rational affine check in
`check_continuations.py`. These are research results, not formalized canon.
No claim of novelty over the inherited affine-word formulas. No Collatz
convergence or nontrivial-integer-cycle exclusion conclusion.

## Findings first

* The inherited words `(1,1,5)` and `(1,3,3)` have clocks `(k,H)=(3,7)`,
  denominator 101 and different numerator carries 19 and 31. Append `(1)`
  **and reclose the new word**: the values become `37/35` and `221/175`.
* The stronger symmetry obstruction already occurs **within one cycle**:
  rotations `(1,3,3),(3,3,1),(3,1,3)` all have denominator 101, but append-one
  denominators are `175,175,7`. An unmarked necklace is not an attachment state.
* Across all 15 marked positive compositions of 7 into 3 parts, denominator
  101 is constant, while append-one denominators have multiplicities
  `{7:1,25:1,35:5,175:8}`. This is the complete stated clock fibre.
* Every suffix denominator is a prime-by-prime capped distance from the
  retained numerator to a suffix-dependent residue centre. The old denominator
  records just one different cancellation test. This explains exactly why it
  cannot control future reclosure.

## Inheritance and scope

Closest mechanism: C1--C3 of
`05-knowledge/results/collatz_fermat_reframes_20261029/README.md` at
`e871dd4513c7dcaa75cb29bb9847aafac35a795c` (read by `git show`; sparse checkout
omits this path). C3 already gives the 35/175 witness and sufficient
`(b,k,H,C)` state; the concatenation law is explicitly inherited there from
`collatz_pentagon_fixed_points_argument_styles_20260930.md`, section 3.2.
Session100's
`bijections/VALUATION_CONJUGACY_AND_INTEGER_SIDECAR.md` supplies the inverse
prehistory distinction. New elaborations here are the explicit finite-menu
residue modulus, cancellation-ball profile, precision comparison, and
same-necklace attachment hostile. Literature priority is not asserted.

Concept board: affine monoid / closure denominator / cyclic conjugacy /
prime-power cancellation / attachment response. Anchor: denominator
continuation. Niche: p-adic cancellation profiles. Wildcard: quotient-versus-
attachment congruence. Canonical hostile: 35/175; corrected near miss:
confusing closed-cycle invariance with reclosure; least-used sidecar: marked
position. No artificial tournament is imposed on arithmetic ties.

## 1. Exact affine monoid carrier

Fix odd integer `b>=3`. For a nonempty word `u=(a_0,...,a_(k-1))` with positive
integer letters, set `H=sum a_i`, `s_j=sum_(i<j) a_i`, and

```
A_u=b^k, B_u=2^H,
C_u=sum_(j=0)^(k-1) b^(k-1-j)*2^s_j,
Delta_u=B_u-A_u,
M_u=[[A_u,C_u],[0,B_u]].
```

The branch composition is `(A_u*x+C_u)/B_u`, with unique fixed point
`x_u=C_u/Delta_u`. The denominator cannot vanish since `B_u` is even and
`A_u` odd. For the empty suffix use `(A,B,C)=(1,1,0)`; its own closure is
undefined, but composing it with a nonempty word is permitted.

For every suffix `v`, chronological concatenation satisfies

```
M_(uv)=M_v M_u,
A_(uv)=A_v A_u, B_(uv)=B_v B_u,
C_(uv)=A_v C_u+B_u C_v,
Delta_(uv)=B_u B_v-A_u A_v.
```

This follows by multiplying triangular matrices, or composing the affine
maps. A letter `a` updates `(A,B,C)` to `(bA,2^a B,bC+B)`.
Thus the triple is an exact sufficient state for **all** future suffixes,
not merely a closed-cycle denominator. Its coordinates retain chronological
order through the noncommuting top-right entry. For example, two single
letters `a,c` have carries `b+2^a` and `b+2^c` after swapping, while their
clocks agree. Equal determinant or equal diagonal entries discard precisely
this order-sensitive data.

### Exact carry is a lossless code, not genuine information compression

For fixed `b,k,H`, the numerator carry `C` uniquely determines the entire
marked word. For `k>1`, split off the first letter `a_0`:

```
C_b(u)-b^(k-1)=2^a_0 C_b(tail).
```

The tail is nonempty, so its carry is positive and odd. Hence the difference
is positive and `a_0=v_2(C-b^(k-1))`. Divide the difference by `2^a_0`, reduce
`k` by one and `H` by `a_0`, and recurse. At length one the carry is always 1
and the last letter is the remaining `H`. This proves uniqueness and gives
an exact decoder. The zero-difference case occurs only at the terminal
single-letter stage of valid states, where the procedure does not take a
valuation. A validity decoder also rejects nonpositive remainders, zero
letters, terminal carry other than 1, or nonpositive final `H`.

Thus `(b,k,H,C)` does not need a separately retained word even for internal
phase recovery: it reconstructs that word. Representing a word by one exact
integer can be algebraically convenient without reducing its information.
The finite-menu residue state in section 3 really does discard information,
and its guarantee is correspondingly limited to the declared observations.
This sharpens the inherited C3 wording that allowed retaining the word as an
additional phase sidecar; exact carry already encodes it.

All nonempty `C` and `Delta` are odd, and `gcd(Delta,b)=1`. Consequently the
reduced denominator `d=|Delta|/gcd(C,Delta)` is coprime to `2b`.

## 2. Exact all-suffix denominator law and local centres

Write `A=A_u,B=B_u,C=C_u`. For each suffix `v`, let

```
N_v=A_v C+B C_v,
D_v=|B B_v-A A_v|,
d_v=D_v/gcd(N_v,D_v).
```

This is an iff integrality test: the reclosed cycle is integral exactly when
`N_v=0 mod D_v`; positivity additionally requires `B B_v>A A_v`.
For prime `p` with `e=v_p(D_v)>0`, `A_v` is a unit modulo `p^e`. Define

```
xi_(v,p,e)=-B C_v A_v^(-1) mod p^e.
```

Then the exact prime-power law is

```
v_p(d_v)=e-min(e,v_p(C-xi_(v,p,e))).
```

The capped valuation is independent of the integer lift chosen for `xi`.
Proof: `N_v=A_v(C-xi)` modulo `p^e`, and multiplication by the unit `A_v`
preserves divisibility through precision `p^e`. This is a nested residue-ball
profile, not just a list of which primes occur.

For two prefixes with the same `(k,H)` and carries `C,C'`, their continuation
numerators differ by `A_v(C-C')`. Thus, if `t=v_p(C-C')`, their cancellation
depths agree when capped at `min(e,t)`. If one cancellation depth is strictly
below `min(e,t)`, both depths are equal. If `p` does not divide `C-C'`, at most
one continuation numerator is divisible by `p`. Differences beyond precision
`t` are exactly where phase-sensitive cancellation can split the outputs.
These are immediate consequences of the ultrametric divisibility law; they
are valid for every suffix, not a census extrapolation.

For the anchor, `D_(1)=175=5^2*7` and `N_(1)=3C+128`.
The cancellation centre modulo 175 is `xi=74`. At `C=19`,
`C-xi=-55` has one factor of 5 and none of 7, leaving `35`; at `C=31`,
`C-xi=-43` is coprime to 175, leaving `175`. For the rotation with `C=49`,
`C-xi=-25` cancels both factors of 5, leaving `7`. The original 101 test
has no way to predict these distinct residues modulo 25 and 7.

## 3. Finite menus: exact residue modulus, honest denominator quotient

Fix `(b,k,H)` and a finite nonempty suffix menu `V`. Set

```
L_V=lcm {D_v:v in V}.
```

Retaining `C mod L_V` suffices to compute every denominator in the menu.
More strongly, for **arbitrary integer** carries `C,C'`, the complete
numerator-residue responses `N_v mod D_v` agree for all `v` iff
`C=C' mod L_V`. Indeed the difference is `A_v(C-C')`, and `A_v` is a unit
modulo each `D_v`. This is an exact modulus characterization for full residue
responses. It does not assert that every integer carry is realized by a word.

**Do not overclaim minimality.** Denominators observe only capped valuations,
not complete residues. Their exact menu-equivalence is equality of the
vectors `min(e,v_p(C-xi_(v,p,e)))` over all suffix-prime pairs. This can merge
many different residue vectors. `C mod L_V` is a convenient sufficient state,
not proved minimal for denominator observations. An all-suffix minimal
quotient is not identified here. This distinction parallels the difference
between retaining a labelled attachment response and retaining only its rank.

## 4. Rotation is a symmetry of closure, not of marked attachment

Let `u=(a)r`. Then `M_u=M_r M_a`, while
`M_(r(a))=M_a M_r=M_a M_u M_a^(-1)` over the rational numbers.
The fixed point transforms by the actual branch `(b*x+1)/2^a`.
The integer identity

```
b C_u+Delta_u=2^a C_(r(a))
```

shows `gcd(C_u,Delta_u)=gcd(C_(r(a)),Delta_u)`, since `2b` is coprime to
`Delta_u`. Therefore the closure denominator is rotation-invariant.

But appending `v` means left multiplication by **the fixed marked matrix**
`M_v`. In general `M_v M_a M_u M_a^-1` is not conjugate to `M_v M_u` by
`M_a`: conjugating also changes the attachment matrix unless it commutes.
The explicit 175/175/7 witness proves failure without relying on a generic
noncommutativity slogan. The repair is a marked attachment site, or the full
response indexed over the orbit of sites. A multiset of site-responses is a
well-defined necklace invariant, though it loses correlations needed for
multiple successive, site-dependent operations.

**Connection contract (structural, not an arithmetic-to-tournament map):**
source: marked objects with extension operators; quotient: forgetting the
mark modulo symmetry; preserved predicate: the closed-object invariant;
lost data: where the new piece attaches; sidecar: orbit of marked sites and
its response function; decisive test: apply the same labelled extension to
two representatives of one symmetry orbit. The canonical Collatz test fails.
For tournaments this instructs retaining attachment subsets modulo the
parent automorphism action, not just the parent's unrooted invariant. No
specific tournament enumeration theorem or Collatz-to-tournament bijection
is claimed by this connection.

## 5. Three different operations, three denominator behaviours

1. **Actual forward odd branch on the cycle:** denominator remains `d`.
   Its numerator is `bp+q` for reduced `p/q`; `gcd(bp+q,q)=1` because
   `gcd(b,q)=1`. Dividing out the exact power of two does not alter odd `q`.
2. **Finite inverse prehistory:** `(2^a*p-q)/(bq)` has numerator coprime to
   odd `q`, so its reduced denominator remains a multiple of `q`. It can
   acquire factors of `b`; the inherited 101 obstruction cannot disappear.
3. **Append letters and reclose:** the new fixed point has new gap `D_v`
   and new numerator `N_v`; neither denominator divisibility nor monotonicity
   relative to the old `d` holds. The 101 to 35/175/7 examples concern only
   this operation. They do not describe temporal steps of one rational orbit.

## 6. Exact reproducibility and limits

Run from this directory:

```
python3 -B check_continuations.py > checks.json
python3 -B -O check_continuations.py > optimized.json
cmp checks.json optimized.json
```

Universe: multipliers 3,5,7; all 340 nonempty words of lengths 1--4 on
letters 1--4; suffixes empty or length 1--2 on letters 1--3 (13 choices).
There are 13,260 composition comparisons. Gates verify affine recurrence
against independent Fraction composition, fixed points, exact valuations,
rotation denominators, lossless marked-word decoding, every prime-power continuation formula, the complete
15-word `(3,7)` fibre, the finite-menu modulus law and cancellation precision.
No positivity filtering is hidden: signed rational cycles are included in
the broad box; all words in the 101 fibre are positive. Normal/optimized
outputs agree byte-for-byte. Factorization is trial division and all
arithmetic exact. The infinite statements are established by the proofs,
not by the finite checks.

Remaining frontier: minimal all-suffix denominator equivalence, efficient
finite-menu residue compression respecting successive changing clocks,
and integer-selective operations that preserve forcing 1. None is supplied
by mere denominator spectra or symmetry quotients.
