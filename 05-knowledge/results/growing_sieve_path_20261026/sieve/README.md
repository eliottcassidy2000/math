# Growing twin wheels: exact counts, self-hits, and the height boundary

**Status: PROVED elementary identities and boundary statements;
independently audited [PASS](../audit_sieve/REPORT.md).
Computation FINITE-EXACT.** This is a self-contained specialization of
classical Eratosthenes–Legendre inclusion-exclusion to the established
Pythagorean ray. No historical-priority, twin-prime infinitude, asymptotic,
or improved analytic sieve bound is claimed.

## 1. Inheritance and the precise carrier

The [previous tree packet](../../twin_pythagorean_124_20261026/tree/RESULT.md)
proves that, for odd x>=3,

    T(x)=(x(x+2),2(x+1),(x+1)^2+1)

is primitive, with inradius x and even-leg exradius x+2. Twin primes select
exactly its height-two, total-prime-multiplicity-two locus. Its Berggren word
is `D^((x-1)/2)`, and D sends x to x+2. The
[fixed-residue obstruction](../../twin_pythagorean_124_20261026/residue_obstruction/README.md)
already shows that every fixed modulus admits both-composite false positives
with the same finite observations as a twin pair.

The [prime-pair structural note, §§2 and 5](../../prime_pairs_local_symmetry_20261019/prime-pairs-structural-note.md)
already distinguishes roughness from primality, explicitly warns about
prime self-hits, and proves periodic four-form representation laws. Its
[power-depth continuation](../../prime_pairs_local_symmetry_20261019/power-depth-rough-pairs.md)
has a separate analytic four-form almost-prime theorem. The present carrier
has two forms x,x+2 and a finite-height exact count; it does not improve that
analytic theorem. [THM-2447 / local prime harmonics](../../../../01-canon/theorems/THM-2447-twin-gap-local-prime-harmonics-and-detrended-spectrum.md)
is further local-factor ancestry. The square-root threshold and Möbius
indicator are classical; see [source receipt](SOURCES.md).

Anchor: growing prime checks on T(x). Niche: an exact finite counting
recurrence with the self-hit correction. Wildcard: coherent wheel paths
that have no positive bounded integer lift. Canonical hostile: a prime is
divisible by itself. Corrected near miss: positive counts at every finite
wheel do not give a positive integer surviving all strict wheels. Least-used
sidecars: integer height, exact prime equality, and both endpoint factors.

| Concept | Representation | Preserved / missing coordinate |
|---|---|---|
| Integer triangle | x, T(x), D depth `(x-1)/2` | all arithmetic retained |
| Strict wheel | residue x modulo a primorial | small-prime avoidance; self-primes lost |
| Restored sieve | divisibility plus equality to tested prime | twin predicate at sufficient height |
| Finite count | residue floors or Möbius sum | exact endpoints and signed errors |
| Compatible path | nested residue classes | no bounded positive integer guaranteed |

Throughout B is an integer at least 2. Let

    P_B={primes ell<=B}, M_B=product_(ell in P_B) ell,
    W_B={a in {0,...,M_B-1}: gcd(a(a+2),M_B)=1},
    R_B(H)=#{x in Z: 3<=x<=H, gcd(x(x+2),M_B)=1},
    pi_2(H)=#{primes p: 3<=p<=H and p+2 is prime}.

Set `pi_2(H)=0` for H<3. The upper member of a counted twin pair may exceed
H by two. Including prime 2 makes all strict survivors odd; the explicit
lower bound x>=3 excludes x=1, which gives the non-twin root `(3,4,5)`.

## 2. Nested wheels, their lifts, and their exact counts

**Theorem W1.** Put `nu(2)=1` and `nu(ell)=2` for odd primes. Then

    |W_B|=product_(ell<=B)(ell-nu(ell)).                 (1)

When a new odd prime ell is adjoined to a wheel of modulus M, every allowed
residue a modulo M has exactly ell-2 allowed lifts modulo M*ell. Explicitly,
the lifts are `a+tM`, `0<=t<ell`, except the two distinct indices

    t=-a/M and t=-(a+2)/M modulo ell.                  (2)

*Proof.* A tested odd prime forbids precisely residues 0 and -2. At prime
2 these coincide. CRT multiplies the allowed local cardinalities. In a
lift, M is invertible modulo the new prime, giving (2). QED.

The reduction map from each finer wheel onto the preceding wheel is thus
surjective. In integer triangle coordinates the lift changes D depth by
`tM/2`. It is not a single D child or a ternary Berggren refinement. The
prime-five step has three children per residue; the prime-seven step has
five. The branch count ell-2 belongs to a residue-refinement tree.

For every H>=3, counting each residue class exactly gives

    R_B(H)=sum_(a in W_B)
      [floor((H-a)/M_B)-floor((2-a)/M_B)].             (3)

There is a second finite formula. For squarefree d define

    Z(d)={r in {0,...,d-1}: r(r+2)=0 modulo d},
    Z(1)={0}.

CRT gives `|Z(d)|=product_(ell|d)nu(ell)`. If mu is the Möbius function,

    R_B(H)=sum_(d|M_B) mu(d) sum_(r in Z(d))
      [floor((H-r)/d)-floor((2-r)/d)].                 (4)

*Proof.* For each integer x, the finite divisor sum
`sum_(d|gcd(x(x+2),M_B))mu(d)` is one if the gcd is one, and zero otherwise.
Sum over the interval and split each divisibility condition into its roots.
This is classical inclusion-exclusion with its actual interval endpoints.
QED.

Writing `delta_B=|W_B|/M_B`, two immediate uniform, but possibly large,
error bounds are

    |R_B(H)-(H-2)delta_B|
      <= min(|W_B|, product_(ell<=B)(1+nu(ell)))
      = min(|W_B|,2*3^(|P_B|-1)).                    (5)

Each residue-interval count differs from its interval length divided by
the modulus by less than one. Applying this once to (3) and once to (4)
proves the bounds. Formula (4) is exact, but its term count grows rapidly;
neither (4) nor (5) supplies a positive lower bound at the primality scale
merely by discarding the signed floor errors.

## 3. Restoring prime self-hits exactly

Define `E_B(x)` for integer x>=3 by the following local tests for every
ell<=B:

    (ell does not divide x OR x=ell), AND
    (ell does not divide x+2 OR x+2=ell).              (6)

This is trial division by proper prime factors, not strict roughness.

**Theorem W2.** For every B>=2, the set accepted by (6) is the disjoint union

    {strict wheel survivors x>=3}
       disjoint union {twin lower primes p<=B}.       (7)

*Proof.* A strict survivor passes (6). Every twin pair passes (6), and if
its lower prime is at most B it fails strict avoidance. Conversely suppose
(6) passes but strict avoidance fails. Some tested odd prime ell must equal
one of x,x+2. The other entry is at most ell+2. If composite, it has a
proper prime factor at most `sqrt(ell+2)<ell`, which is tested and violates
(6). Thus the other entry is prime too, and the lower prime is at most B.
Prime 2 cannot occur in this argument because an accepted x>=3 is odd. QED.

The sets in (6) are nested as B grows. In particular, the exact number of
self-hit survivors in a height-H box is `pi_2(min(H,B))`. This finite patch
is forced by the local equality conditions; it is not permission to admit
arbitrary small exceptional composites.

For example, strict sieving by prime 3 deletes the genuine pair (3,5).
At B=5 it also deletes (5,7). Test (6) keeps both. The same tests reject
(25,27): equality with a tested prime cannot rescue a composite multiple.

## 4. Height completeness and a recursive exact count

Let q be the least prime greater than B.

**Theorem W3.** If H>=3 and `H+2<q²`, then (6) accepts x in `[3,H]` if and
only if x,x+2 are prime. Consequently

    pi_2(H)=R_B(H)+pi_2(min(H,B)).                     (8)

In particular, the choice `B=floor(sqrt(H+2))` gives

    pi_2(H)=R_floor(sqrt(H+2))(H)
                 +pi_2(floor(sqrt(H+2))).           (9)

The latter argument is strictly smaller than H for H>=3. Equations (3) or
(4), (9), and the base values below 3 therefore give a finite exact counting
recurrence.

*Proof.* A composite integer at most H+2 has a prime factor at most its
square root, which is less than q. All such primes are at most B and are
tested. A proper divisor cannot use the equality exception. Twin pairs
always pass. Equation (7) gives (8). For the selected B, the next prime is
strictly greater than sqrt(H+2), and B<H. QED.

This theorem goes beyond a fixed observer by growing the tested primes
with the height. It does not turn the finite counting recurrence into an
infinitude proof or a guaranteed efficient algorithm.

**First square boundary.** No false accepted pair can have lower entry
less than `q²-2`. At exactly `x=q²-2`, a false pair survives if and only if
`q²-2` is prime: its upper entry is q², while a composite lower entry has
a prime factor less than q and is removed. Thus the universal threshold
is attained at B=2,3,5 by `(7,9),(23,25),(47,49)`. It need not be attained
at every stage: at B=7, q=11 and `119=7*17` is already removed.

## 5. Exact experiment: why checking only squares fails

The checker finds the first false strict survivor by inspecting **every**
smaller x, with complete least-prime-factor data and independent trial
division. Selected outcomes are:

| B | next prime q | first possible x=q²-2 | actual first false pair |
|---:|---:|---:|---|
|2|3|7|(7,9)|
|3|5|23|(23,25)|
|5|7|47|(47,49)|
|7|11|119|(167,169)|
|11|13|167|(167,169)|
|13|17|287|(359,361)|
|59|61|3719|(3719,3721)|
|61|67|4487|(4757,4759)|

At the final row, `4487=7*641` is removed and
`4757=67*71`, while 4759 is prime. Thus the first composite entry need not
be an upper prime square. The checker covers all prime cutoffs through 61
and proves this is the first stage in that finite sequence at which the
displayed square-upper pattern fails. Between prime cutoffs the wheel is
unchanged. No all-stage square pattern is proposed.

The recurrence is independently checked through lower height 10000. For
example, `pi_2(100)=6+2=8`, `pi_2(1000)=30+5=35`, and
`pi_2(10000)=197+8=205`; in each expression the first term is the strict
count at the square-root cutoff and the second is its exact self-hit patch.
These are finite counts, not fitted density estimates.

## 6. The quantifiers do not interchange

Every finite wheel has infinitely many positive strict survivors. In fact,
`x=-1 modulo M_B` always survives. The earlier CRT construction, applied
to any allowed residue and two new primes, also supplies infinitely many
both-composite strict survivors at that fixed stage.

Nevertheless no positive integer x>=3 survives **all** strict wheels:
a prime factor of x eventually enters the wheel and deletes it, even when
x itself is prime. Thus

    for every B there are arbitrarily large strict survivors,
    but there is no x>=3 that survives every B.

The compatible residue path `x=-1 modulo M_B` makes the missing coordinate
particularly explicit. It has the ordinary integer representative -1,
center `m=x+1=0`, and signed null triple `T(-1)=(-1,0,1)`, outside the
positive Berggren ray. Its nonnegative representatives `M_B-1` grow with
the modulus; choosing such a representative afresh at each stage does not
give one bounded positive integer. An inverse-limit existence statement
therefore supplies neither positive height nor the twin predicate.

In contrast, the intersection of the exception-restored sets (6) is
exactly the set of actual twin lower primes. Every fixed x stabilizes
correctly once the tested primes reach sqrt(x+2). This pointwise statement
does not establish that the intersection is infinite.

## 7. Connection contract and reproduction

Source: the positive integer D ray, with x>=3 and its full factorization.
Target: nested allowed residues plus exact self-hit conditions.
Map: reduction modulo M_B; branch lifts are (2), not Berggren child moves.
Preserved: small-prime divisibility, refined residue ancestry, and the twin
predicate only after the height gate in W3.
Lost: integer height and prime equality under a bare residue quotient.
Needed sidecars: height H, equality x=ell or x+2=ell, and both factors.
Decisive tests: deletion of (3,5), the boundary (23,25), the delayed boundary
119, the nonsquare failure 4757, and the coherent negative lift -1.

`check.py` is standard-library-only and uses explicit runtime gates. It
compares direct gcd, wheel lifts, residue-floor counts, Möbius root counts,
and independent trial division on declared finite universes. Run:

```bash
python3 -B check.py
python3 -B -O check.py
```

The two outputs must match `normal.json` and `optimized.json`. Source pins,
classical attribution, and the later independent audit are recorded beside
this note. No numerical census is promoted to a global twin-prime result.
