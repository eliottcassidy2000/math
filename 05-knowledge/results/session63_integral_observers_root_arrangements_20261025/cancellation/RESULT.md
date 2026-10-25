# Integral cancellation spectra: odd actions, every prime, and a 50-vertex tournament

2026-10-25. **PROVED, independently audited:** the scoped action refinements
and verification below. **RECOVERED CONCURRENT RESULT:** the complete regular
elementary-abelian spectrum (Theorem3) and the exact 50-vertex witness (§5)
already occur in incoming checkpoint `a971f9b4`; they are reproduced here with
independent proofs and checks, not claimed as new results. **FINITE-EXACT:**
the declared arithmetic/action/tournament universes in `check.py`. No historical novelty,
new black move, enumeration speedup, HN bound, or LRC conclusion is claimed.

## 0. Inheritance and the retained coordinates

The immediate proved input is
`05-knowledge/results/session62_agent_theme_bridges_20261025/clone/README.md`
at repository checkpoint `eaa81fd80089cc3b0a56d346f59403a0cfb730c2`.
It proves the rational minimum-cycle formula, exact invariant-toggle rank,
and canonical clone-fibre interpretation. Its underlying module theorem is
THM-4463, `01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md`,
§2, and its fibre source is
`05-knowledge/results/ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md`.
Those sources supply the identification with actual weighted inflation fibres.

### Recovered source and explicit identification

The initial search at `eaa81fd8` preceded the incoming checkpoint
`a971f9b4fbdf409c31ba61ad9e8b3b00648536cf`. Reading that checkpoint recovered:

- `05-knowledge/results/unprotected_ports_cancellation_20261025/valuation/README.md`:
  the complete regular elementary-abelian spectrum, its integral expansion,
  exceptional-prime boundaries, and cyclic regular-action comparison;
- `05-knowledge/results/unprotected_ports_cancellation_20261025/counting/POSITIVE_JET.md`,
  §4: the transcript's exact 50-vertex witness;
- `05-knowledge/results/unprotected_ports_cancellation_20261025/valuation/INDEPENDENT_REVIEW.md`:
  its prior independent algebraic and tournament checks.

The recovered connection set is **exactly the negative** of the one in §5.
The map `u -> -u` identifies the 25-vertex cores, and `(u,b) -> (-u,b)`
identifies the ordered-pair inflations. `verify_source_recovery.py` checks
every one of the 625 core and 2500 inflated adjacency entries. Thus the
transcript's example is **recovered up to this explicit isomorphism**, not
merely matched numerically. Its source convention `Q(u)=-P(-1-u)` differs
from our local shift `P(-1+x)` by `x=-u` and an overall minus sign; root
orders agree. Both conventions and the mod5 coefficients are checked.

Theorem3 and §5 retain their proofs as **independent rederivation and
verification of recovered concurrent results**. Theorem1's universal odd-action
parity/jump interface, Theorem2's exact denominator-precision criterion, and
Theorem4's diagonal all-rank spectrum are scoped extensions not stated in the
inspected incoming packet; this is a bounded comparison, not a historical
novelty claim. Independent order-join products retain the operation sidecar.
The incoming cyclic-action depth, sparse carry expansion and common ordinary
derivative are credited below and are not reintroduced as new theorems.

- **Anchor:** exact root order after reduction of the integral orbit polynomial.
- **Niche:** denominator precision and modular invariant spaces.
- **Wildcard:** a one-direction perturbation of a Cayley tournament produces
  a compact, independently checkable intrinsic realization.
- **Inherited hostile:** regular `C3²` already separates rational root order
  from one legal invariant toggle.
- **Corrected near miss:** reducing a Burnside numerator before dividing by a
  nonunit group order destroys the information needed by the quotient.
- **Least-used sidecar:** integral shifted coefficients, not just their
  rational leading degree or the rank of common fixed vectors.

| Live concept | Exact map / preserved information | Lost information / decisive control |
|---|---|---|
| Subset quotient | `C -> [C]`, retaining cardinality | Literal subset multiplicity; orbit-size parity restores mod2 comparison |
| Integral averaging | cycle-index numerator divided in `Z[t]` | Mod-p numerator alone loses p-adic precision; `C5²` needs division by25 |
| Root specialization | `P_Z -> P_Fell`, retaining all reduced coefficients | Positive real leading coefficient can vanish; keep its prime factors |
| Legal toggles | invariant subsets act by XOR | Higher root order supplies no additional XOR directions |
| Tournament realization | homogeneous clone pairs -> gap action | Abstract actions are not automatically full automorphism actions; certify stabilizer |

The methods instantiate the existing meta-patterns “Search the statement before
its method”, “Type every analogy and every implication”, and “Controlled
forgetting and unlabeled quotients require a sidecar”. No new repository-wide
meta-pattern is inferred from this single lane.

## 1. All odd actions: universal dyadic collapse and exact prime support

Let a finite odd-order group `A` act on a finite set `G`, with `n=|G|`, and put

```
P(t) = sum_[C] in (2^G/A) t^|C| in Z[t].
d_0 = ord_(t=-1) P over Q,
d_ell = ord_(t=-1) (P mod ell) over F_ell.
```

Empty-set conventions are `P=1` and all depths zero. The leading coefficient
of `P` is one, so no reduction is the zero polynomial. Let

```
a_j = [x^j] P(-1+x),       kappa = a_(d_0) > 0.
```

**Theorem 1 (all-prime interface).**

1. `d_0 = min_(g in A) c(g)`, with
   `kappa = |A|^(-1) sum_(c(g)=d_0) product_(D cycle of g) |D|`.
2. For every prime `ell`, `d_ell >= d_0`, and
   `d_ell > d_0` **iff** `ell` divides `kappa`.
3. For every such action,
   `P(t) == (1+t)^n (mod 2)`, and therefore **`d_2=n`**.
4. For odd `ell`, `d_ell == n (mod 2)`.

**Proof.** Part 1 is the inherited positive-leading-term theorem, recalled
here to specify the dependency rather than rebrand it. Since the coefficients
below `d_0` vanish integrally, reduction cannot introduce an earlier nonzero
coefficient. The coefficient at `d_0` survives exactly when `ell` does not
divide `kappa`, proving part 2. For part 3, each subset orbit has odd size,
since it divides `|A|`. In each cardinality the number of labelled subsets is
therefore congruent modulo2 to the number of orbits. Its enumerator is
`(1+t)^n`. For part 4, complement gives `P(t)=t^nP(1/t)` in every field.
Writing `P=(1+t)^d Q`, where `Q(-1) != 0`, yields
`Q(-1)=(-1)^(n-d)Q(-1)`. In odd characteristic this forces `n-d` even. ∎

In particular, any nontrivial effective odd action has `d_0<n` and even
`kappa`: characteristic2 always raises its depth. The identity action has
`P=(1+t)^n`, `kappa=1`, and no characteristic jump.

**Theorem 2 (integral averaging precision).** Write

```
N_j = [x^j] sum_(g in A) product_(D cycle of g) (1+(-1+x)^|D|),
e = v_ell(|A|),                  v_ell(0)=infinity.
```

Then `v_ell(N_j)>=e` for every `j`, and

```
d_ell = min { j : v_ell(N_j)=e }.
```

To determine the quotient coefficients modulo `ell`, it is sufficient to know
the numerator coefficients modulo `ell^(e+1)` together with the integer
`|A|`; numerator reduction modulo `ell` alone generally does not suffice.

**Proof.** Integral Burnside counting gives `N_j=|A| a_j`. Therefore
`v_ell(N_j)=e+v_ell(a_j)`. The first unit quotient coefficient is exactly the
first coefficient whose numerator has valuation `e`. Dividing after retaining
`e+1` digits yields `a_j mod ell`, after inversion of the prime-to-ell part of
`|A|`. ∎

This is a complete finite-coefficient criterion for every odd action and
prime, not a claim that minimum-cycle counts alone determine modular depths.
For the actual clone fibre the signed binomial moments satisfy

```
M_j = sum_[C] (-1)^(n-|C|) binom(|C|,j) = (-1)^(n+j) a_j.
```

Thus modulo `ell` they vanish exactly through degree `d_ell-1`, and the
degree-`d_ell` moment is nonzero. This is the target-bearing cancellation
statement; it uses binomial/Hasse coefficients rather than ordinary
derivatives that may vanish for factorial reasons.

The polynomial coefficients remain **integer orbit counts**; reducing those
integers is not the same as changing their numerical dimensions. The
specialization framework in the companion lane treats the invariant-lattice
base-change statement separately. No modular Reynolds projector is assumed.

## 2. Elementary abelian regular actions: complete characteristic spectrum

Let `p` be an odd prime, `k>=1`, `n=p^k`, `m=p^(k-1)`, and let
`E=(C_p)^k` act regularly on itself. Every nonidentity element consists of
exactly `m` cycles of length `p`. Hence

```
P_(p,k)(t) = ((1+t)^n + (n-1)(1+t^p)^m)/n.                 (1)
```

**Theorem 3 (recovered concurrent spectrum; independently rederived).** Its invariant-XOR rank is
one, and its depths are

```
d_0   = m,
d_ell = n          if ell divides n-1,
d_p   = n-p+1,
d_ell = m          otherwise.                             (2)
```

The cases in (2) are disjoint, since `p` does not divide `n-1`.
For `k=1`, `n-p+1=1=m`: characteristic `p` does **not** cause a jump.
For `k>=2` it does. In particular the rational leading coefficient is

```
kappa = (n-1) p^(m-k).                                    (3)
```

The exponent is nonnegative for odd `p` and `k>=1`; at `k=1` it is zero.
The prime-support criterion in Theorem1 agrees exactly with (2).

**Proof away from p.** Substitute `t=-1+x` and put
`f(x)=1+(x-1)^p`. In characteristic `ell != p`, `f` has simple zero at0.
If `ell` divides `n-1`, (1) reduces to `x^n`, since `n=1` in that field.
Otherwise the degree-`m` term has nonzero coefficient
`(n-1)p^m/n`, while the identity summand begins only at degree `n>m`.
This proves the two good-denominator cases and the rational statement.

**Proof at p: divide integrally before reducing.** Write

```
f(x)=x^p+p h_p(x),
h_p(x)=sum_(j=1)^(p-1) ((-1)^(p-j)/p) binom(p,j) x^j.
```

Thus `h_p` is an integral polynomial of exact order1 with linear coefficient1.
Expanding (1), the `j=0` terms combine into `x^n`, leaving

```
P(-1+x) = x^n + (n-1) sum_(j=1)^m
             (p^j binom(m,j)/p^k) x^(n-pj) h_p(x)^j.       (4)
```

Every displayed scalar is integral. Indeed, for `1<=j<=m=p^(k-1)`,

```
v_p(binom(m,j)) = k-1-v_p(j),
v_p(p^j binom(m,j)/p^k) = j-1-v_p(j).                      (5)
```

For clarity, the first identity follows from
`binom(m,j)=(m/j)binom(m-1,j-1)`. The latter binomial is a p-adic unit:
modulo `p`, `(1+z)^(p^(k-1)-1)=sum_(i=0)^(m-1) (-z)^i`, obtained by
multiplying by `1+z` and using `(1+z)^m=1+z^m`. For `j=1` the second
valuation in (5) is zero. For `j>=2` it is at least one, because `p` is odd.
Consequently

```
P(-1+x) == x^n - x^(n-p) h_p(x) (mod p).                  (6)
```

Its first nonzero term is `-x^(n-p+1)`. Its coefficient is `-1`, so there
is no hidden later cancellation. This proves the bad-prime case for every
`k>=1`, including `m=1`. Transitivity supplies the invariant-XOR rank1. ∎

### The 25-gap specialization

For `p=5,k=2`, `n=25,m=5`, the rational leading coefficient is `3000`.
Since `h_5=x-2x²+2x³-x⁴`, formula (6) is explicitly

```
P(-1+x) == -x^21+2x^22-2x^23+x^24+x^25 (mod 5).
```

Thus the complete spectrum is

| Coefficients | Root depth |
|---|---:|
| Q | 5 |
| F2 or F3 | 25 |
| F5 | 21 |
| every other prime field | 5 |

A leading-coefficient argument alone would identify exceptional primes
`2,3,5` but not their actual depths. Formula(6) supplies the missing precision.

## 3. All ranks: distinguish synchronized from independent orbits

**Theorem 4.** Let the same `E=(C_p)^k` act diagonally on `s>=1` disjoint
regular copies, with `n=p^k,m=n/p`. The point-orbit count and invariant-XOR
rank are `s`, and

```
P_diag(t) = ((1+t)^(sn)+(n-1)(1+t^p)^(sm))/n,
d_0 = sm,
d_ell = sn                    if ell divides n-1,
d_p = sn-p+1                  if p does not divide s,
d_p = sn                      if p divides s,
d_ell = sm                    otherwise.                  (7)
```

**Proof.** The cycle formula and all characteristics other than `p` are as
above. In the expansion of `(x^p+p h)^(sm)`, the normalized `j`-coefficient
has valuation at least

```
v_p(p^j binom(sm,j)/p^k) >= j-1+v_p(s)-v_p(j),
```

because `binom(sm,j)=(sm/j)binom(sm-1,j-1)`. This bound proves integrality
for all `j>=1` and divisibility by `p` for all `j>=2`; it does not require
`sm` itself to be a prime power. The exact `j=1` coefficient is `s`.
Therefore

```
P_diag(-1+x) == x^(sn)-s x^(sn-p)h_p(x) (mod p).            (8)
```

If `p` divides `s`, only `x^(sn)` remains; otherwise the leading degree is
`sn-p+1`. This proves all cases, with no assumption on `v_p(s)`. ∎

By contrast, if independent copies of `E` act on the separate copies of `G`,
the acting group is `E^s`, subset orbits form Cartesian products, and

```
P_ind(t)=P_(p,k)(t)^s,       d_ell(ind)=s*d_ell(p,k).        (9)
```

Both actions have `s` point orbits and the same number `sn` of points; the
shared versus independent group element is the missing coordinate. A sharp
small control is three 3-point orbits: diagonal `C3` gives characteristic3
depth9, whereas the independent `C3³` action gives depth3. This is an actual
operation difference, not a rank-renaming analogy. The diagonal all-rank
statement is about abstract actions; §5 realizes the independent product
construction intrinsically in tournaments.

## 4. Scope boundaries and hostile tests

1. **Even group order is forbidden in Theorem1's dyadic law.** A transposition
   on two points has `P=1+t+t²`, whose characteristic2 depth at `-1=1` is0,
   not2. Its rational minimum-cycle formula also fails.
2. **Ordinary derivatives are not root multiplicities in small characteristic.**
   Use coefficients of `P(-1+x)` (Hasse/binomial coefficients), or repeated
   monic division. For example `x^5` has zero ordinary derivative in F5
   but is not the zero polynomial.
3. **A larger modular depth does not manufacture toggles.** The 25-gap
   action has exactly one fixed symmetric-difference generator in every
   calculation above. Root depth25 is not a 25-dimensional legal toggle group.
4. **The dimension polynomial is reduced, not recounted.** Its coefficient
   integers count orbits. For odd permutation actions, oriented subset-orbit
   sums in exterior powers remain valid over every field; bad-prime failure
   of a division formula must not be promoted to failure of that basis.
5. **Regular elementary abelian is not regular cyclic.** The incoming
   `valuation/README.md` proves that regular `C_(p^k)` has rational depth1
   and characteristic-p depth `p^(k-1)`, with exactly `1+k(p-1)` nonzero
   terms after reduction. Its ordinary derivative agrees with the elementary
   polynomial's derivative modulo `p`, despite different depths. At order25
   the characteristic5 depths are5 versus21. Those are recovered results,
   not claims newly derived here; they strengthen the ordinary-derivative
   hostile while requiring the distinct cycle inventory to remain visible.
6. **An abstract subgroup is not the full tournament automorphism group.**
   The compact witness below certifies the entire group, not merely a regular
   subgroup. Replacing it by an unperturbed lexicographic `C5[C5]` would admit
   a larger wreath-product action and would not justify this regular-action
   orbit polynomial.

## 5. The recovered intrinsic 50-vertex witness, in negated coordinates

This is the witness in incoming `counting/POSITIVE_JET.md` §4, transported
by `u -> -u`. Take `X=F5²`, lexicographically ordered for our independent
certificate. Define

```
S = {(0,1),(0,2),(1,1),(1,2),(1,3),(1,4),
     (2,0),(2,1),(2,2),(2,3),(2,4),(4,0)}.
u -> v  iff  v-u in S.                                    (10)
```

Exactly one of every pair `{a,-a}` of nonzero vectors lies in `S`, so (10)
is a genuine tournament `R` on25 vertices. It is the lexicographic cyclic
5-tournament connection set with just `(1,0)` replaced by `(-1,0)`. Every
translation preserves arcs, and every outdegree is12. A homogeneous pair
would have outdegrees differing by one, so `R` is pair-free.

### Full automorphism certificate, not a subgroup guess

Fix the vertex `(0,0)`. Start with colors `c_0(0)=0` and `c_0(v)=1` otherwise.
Given `c_i`, let the key of `v` be its previous color together with its
outneighbor counts in each previous color class; rank distinct keys
lexicographically to obtain `c_(i+1)`. Every automorphism fixing0 preserves
every such color. The numbers of color classes are

```
2 -> 3 -> 15 -> 21 -> 25.
```

The final colors, in vertex order `(0,0),(0,1),...,(4,4)`, are

```
[0,12,11,23,24,20,9,10,7,6,1,3,2,4,5,13,15,14,16,17,8,19,18,21,22].
```

They are all distinct, forcing the root stabilizer to be trivial. Since the
25 translations act transitively, every automorphism is a unique translation:
`Aut(R)=C5²`. This is a finite exact certificate with an elementary proof of
why the certificate implies the group statement. The checker records every
intermediate coloring and also performs an independent exhaustive
arc-consistent bijection search fixing0; it finds exactly the identity.
Neither method assumes linearity or affineness of an unknown automorphism.

Now inflate each `u in X` to an ordered pair `(u,0)->(u,1)`, and put every
cross-pair arc in the direction prescribed by `R`. Call this tournament `T`.
It has50 vertices. Each intended pair is homogeneous. Distinct core vertices
have an external separator because `R` is pair-free; any lift of the separator
rules out a homogeneous pair crossing distinct blocks. Thus `T` has exactly
25 homogeneous pairs, all disjoint. No larger transitive module is possible:
a transitive module of size at least3 would contain overlapping adjacent
homogeneous pairs. Hence these are precisely its canonical transitive modules.

Every automorphism of `T` permutes the canonical pairs and induces an
automorphism of `R`; the ordered pairs have no internal automorphism.
Conversely each translation of `R` lifts uniquely. Therefore

```
|V(T)|=50,       Aut(T)=C5²,
G(T) is its regular 25-point action,
invariant-XOR rank=1,
(d_Q,d_F5,d_F2,d_F3)=(5,21,25,25).                         (11)
```

The inherited canonical clone-fibre theorem now applies to this **actual**
tournament, not just to a formal group polynomial. There are `1,342,208`
subset orbits in its gap action (sum of the integer coefficients), without
any claimed enumeration of all `2^25` subsets.

### Intrinsic arbitrary-toggle-rank family

For `s>=1`, take the order join of `s` copies of `T`. The tournament `R` is
strongly connected: in a regular tournament a proper source strong component
of size `h` would have average outdegree `25-h+(h-1)/2>12`. Inflating preserves
strong connectivity, so each copy of `T` is a strong component of the join.
Their total order fixes each component setwise. A pair crossing two components
is separated by an outneighbor of its earlier vertex within the earlier
component. Thus no homogeneous pair crosses components, and the full
automorphism group is `(C5²)^s`, with exactly `s` gap orbits and `25s` gaps. The invariant-XOR rank is exactly `s`,
the polynomial is `P_(5,2)^s`, and

```
|V|=50s,    d_Q=5s,    d_F5=21s,    d_F2=d_F3=25s,
d_Fell=5s for every other prime ell.                        (12)
```

This is an all-rank family of actual tournaments. It does not assert that
every abstract diagonal action in Theorem4 has been realized by this join.

## 6. Reproduction and declared finite universe

Python3.10+, standard library only, any working directory:

```bash
python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json
```

Every gate raises an explicit exception; there are no bare `assert` statements.
Both modes passed and their outputs are byte-identical (`normal.json`,
`optimized.json`): **80 arithmetic cases** plus the action and tournament
controls below. The checker uses two independent
root-depth calculations (Horner shift/binomial moments and repeated monic
division), literal subset orbits for its small controls, the generic cycle
index for the25-gap action, and two independent root-stabilizer certificates.

- Arithmetic universe: `p in {3,5,7,11,13,17}`, `1<=k<=5`, `p^k<=343`,
  `1<=s<=6`, `sp^k<=500`. In each case check every prime at most41, all prime
  factors of `p^k-1`, and `p`. These are finite controls of the all-parameter
  proofs, not the scope of those proofs.
- Full subset-orbit universe: empty/trivial actions; regular `C3,C5,C7,C9,C3²`;
  diagonal `C3` on2,3,4 regular copies and `C5` on2 regular copies. The largest
  labelled subset universe here has `2^12` elements.
- Hostiles: even transposition; three synchronized versus independent
  `C3` orbits; the `k=1` no-bad-prime-jump boundary.
- Tournament universe: the literal25-vertex core, all300 core unordered pairs,
  all1225 pairs of the50-vertex inflation, all25 translations, complete root
  stabilizer backtracking, full rooted color refinement, and its25-gap cycle
  index. **Not** a `2^25` full subset census and **not** a census of all
  25-vertex tournaments.

The source-recovery checker is separately reproducible:

```bash
python3 verify_source_recovery.py > source_recovery_normal.json
python3 -O verify_source_recovery.py > source_recovery_optimized.json
cmp source_recovery_normal.json source_recovery_optimized.json
```

Both modes passed with byte-identical output. This check is about provenance
and isomorphism; it neither changes nor duplicates the 80 arithmetic cases.
`SOURCE_PINS.json` retains the initial dependency checkpoint and adds the exact
incoming source files and their byte hashes.

The source of any future connection must preserve the integral action,
coefficient prime and the actual operation. Equal dimensions or equal depths
alone identify none of those objects.
