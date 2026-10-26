# Integer ancestry separates reset examples from the marked orbit

2026-10-25. **PROVED:** the elementary inverse, root, congruence and debt
statements below for the explicit integer map. **FINITE-EXACT:** the 33
reset parameters and bounded controls in the executable block. **OPEN:**
whether the guarded orbit from `(a,x)=(6,27)` ever halts. No Collatz
convergence or native-machine formalization claim is made.

The inherited inverse mechanism is
[Session64, section 4](../../session64_lca_reverse_valuation_20261026/session64-reverse-valuation/REPORT.md),
at incoming commit `9e79465c2c`. The independent derivation here connects
it to this packet's [reset family and guarded debt](README.md#4-a-marked-lca-barrier-and-an-explicit-precision-reset-family).
The hostile is that many integer resets with dangerous regenerated depth
belong to a different integer component from 27. The useful new sidecar is
the *integer root*, alongside the labelled valuation prefix and credit.
The root restores an orbit-membership predicate, not a survival certificate.

## 1. An inverse algorithm that always terminates

For odd integers `x>1`, write

    k=v_2(3x-1),    U(x)=3^(k-1)*(3x-1)/2^k.

If `y=U(x)`, then `3x-1` is prime to 3. Consequently

    r=v_3(y)=k-1,    m=y/3^r,
    k=r+1,          x=(2^(r+1)*m+1)/3.             (N1)

This candidate is an odd integer precisely when

    m = (-1)^r (mod 3).                            (N2)

Indeed `(N2)` is exactly divisibility of the candidate numerator by 3;
the numerator is odd. Conversely, for this candidate,
`3x-1=2^(r+1)m`, and `m` is odd, so the required valuation is exactly
`r+1`. This proves necessity, sufficiency and uniqueness. For `y>1` a
legal candidate also satisfies `x>1`: the only preimage of 1 is 1, and
the displayed positive candidate cannot be a nonpositive integer.

Moreover

    U(x)=(3/2)^k*(x-1/3) >= (3x-1)/2 > x.

Every legal inverse step therefore has
`1<x<y` and `x-1 <= (2/3)(y-1)`. Repeating `(N1)` until `(N2)` fails
terminates after finitely many strictly decreasing positive integers.
Call the terminal integer `root(y)`. Its forward orbit is exactly the
fiber of this root label. In particular, for odd `x,y>1`,

    root(x)=root(y)
      iff their forward U-orbits meet
      iff one lies on the forward U-orbit of the other.              (N3)

Proof: inverse chains are unique; reversing the finite chain from any
point reconstructs its forward segment from the root. Thus every fiber
is an unbranched ray, with no positive nontrivial cycle. Root computation
does not query an infinite continuation.

For `27=3^3*1`, condition `(N2)` fails: `1 != -1 (mod 3)`. Hence

    root(27)=27,
    y is on the unguarded orbit of 27 iff root(y)=27.                 (N4)

This is ordinary integer ancestry, not the 2-adic prefix LCA. For example
7 and 23 both have first valuation 2, but have roots 5 and 23. Conversely,
5 and 7 have different first valuations and the same root 5.

## 2. The reset family and an all-parameter modular exclusion

For `R=6t+5`, `t>=0`, the earlier reset is

    x_R=(2^(R+2)+7)/9=(128*64^t+7)/9,
    U(x_R)=1+2^R,    k(x_R)=2.                    (N5)

Since `2^R=5 (mod 9)`, the target has `v_3(1+2^R)=1` and its unit
part is `2 (mod 3)`. Formula `(N1)` therefore recovers **exactly**
`x_R` as its unique integer predecessor. Thus

    root(x_R)=root(1+2^R),
    this reset belongs to the unguarded 27 ray iff root(x_R)=27.      (N6)

No `x_R` equals 27: `x_5=15 < 27 < x_11=911`, and the family is strictly
increasing. Therefore any `x_R` which itself has no predecessor is
excluded from the 27 ray for every possible credit assignment.

Here is an exact modular description of those self-rooted resets. Put
`X(t)=(128*64^t+7)/9`. For distinct nonnegative integers `s,t`,

    v_3(X(t)-X(s)) = v_3(t-s).                     (N7)

To prove this without assuming a distribution law, write the difference
as `128*64^min(s,t)*(64^abs(t-s)-1)/9`, up to sign. The elementary
3-adic lifting identity is
`v_3(64^h-1)=2+v_3(h)`. For `3` prime to `h`, factor the geometric
sum and reduce it modulo 3; cubing an element `1+3^e*u`, `e>=1`,
increases the valuation of its difference from 1 by exactly one.
These two facts prove the identity and `(N7)`.

It follows that `t -> X(t)` permutes the residues modulo every `3^j`:
it is well-defined modulo `3^j`, and `(N7)` makes it injective, hence
bijective on the finite residue set.

For each `r>=0`, let `alpha_r` be the unique residue modulo `3^(r+1)`
satisfying

    X(alpha_r) = (-1)^(r+1)*3^r (mod 3^(r+1)).      (N8)

Then `(N1)` fails at `x_R` precisely on the disjoint union

    t = alpha_r (mod 3^(r+1)),    r>=0.             (N9)

Indeed `(N8)` says that the 3-adic valuation of `x_R` is exactly `r`
and that its unit has the unique residue forbidden by `(N2)`. Every
positive `x_R` has finite 3-adic valuation, so no integer parameter is
omitted from this failure/success classification.

The first six failure classes are exact, infinite exclusions:

| `r=v_3(x_R)` | `alpha_r` | Excluded reset parameters |
| --- | ---: | --- |
| 0 | 1 | `R = 11 (mod 18)` |
| 1 | 3 | `R = 23 (mod 54)` |
| 2 | 6 | `R = 41 (mod 162)` |
| 3 | 78 | `R = 473 (mod 486)` |
| 4 | 24 | `R = 149 (mod 1458)` |
| 5 | 429 | `R = 2579 (mod 4374)` |

All moduli are divisible by 6 and all residues are 5 modulo 6, as
required. The classes for different `r` are disjoint. Their union has
**relative natural density `1/2` among parameters `R=5 (mod 6)`**.
For a rigorous infinite-union limit, the classes with `r<M` have
density `sum_(r<M)3^(-r-1)=(1-3^-M)/2`. The uncounted classes are
contained in `v_3(X(t))>=M`, which by the permutation property is one
class modulo `3^M`, of density `3^-M`. The lower and upper densities
therefore both tend to `1/2`.

This density concerns reset *parameters*, not natural-number density of
the exponentially sparse values `x_R`, nor density on the marked 27
orbit. A successful first inverse is not evidence for root 27; its
remaining finite inverse chain must still be computed.

There is also no uniform bound on how many literal predecessor steps
are needed to reach a root in this family. The same lifting calculation
gives

    X(t)-1=2*(64^(t+1)-1)/9,
    v_3(X(t)-1)=v_3(t+1).                           (N9b)

When a target `y` is 1 modulo 3, its inverse has `r=0`, `k=1` and
`P(y)=(2y+1)/3`. Thus `P(y)-1=(2/3)(y-1)`: each such inverse step
reduces `v_3(y-1)` by exactly one. Starting at `X(t)`, **exactly**
`v_3(t+1)` consecutive initial inverse steps have `r=0`. At that
endpoint the new value is no longer 1 modulo 3; the next step either
uses `r>=1` or fails. Taking `t=3^j-1` produces at least `j` predecessor
steps, for arbitrary `j`. This obstructs a fixed budget of literal
backward steps, not every conceivable finite encoding or closed formula.

## 3. Attach exact chronological debt to the integer root

Let the finite inverse algorithm return the chronological chain

    z_0=root(y), z_1,...,z_J=y,
    U(z_i)=z_(i+1), k_i=v_2(3z_i-1), K_j=sum_(i<j)k_i.

Starting with formal credit `a_0`, every retained odd-section credit is

    a_j=a_0+2K_j-3j.                               (N10)

The least initial integer credit keeping **every** odd section, including
both endpoints, positive is the computable integer

    B(y)=1-min_(0<=j<=J)(2K_j-3j).                 (N11)

Thus `a_0>=B(y)` is necessary and sufficient for this entire segment to
be live also at its returned endpoint. Executing the last block only
requires the guards at `j<J`; a returned credit zero is already a halt.
Inside a block from credit `a_j>=1`, its minimum credit is `a_j-1`.
Zero is permitted there only when the intermediate scalar is even;
every subsequent even macro increases credit by 2.

For `root(y)=27`, the depth `J` and sum `K_J` are uniquely determined
finite integers, so define

    a_27(y)=6+2K_J-3J,    d(y)=v_2(y-1).

The exact fixed-source reformulation is

    the guarded orbit from (6,27) never halts
      iff every odd y>1 with root(y)=27 has a_27(y)>=1
      iff every such y has d(y)<=a_27(y).            (N12)

The first equivalence is `(N4)` with all chronological guards retained.
For the second, the consecutive valuation-1 run from `y` has length
`d(y)-1`. A positive credit `a<d(y)` therefore reaches zero after
`a` such blocks. Conversely `d(y)>=1`, so the last universal inequality
implies every credit is positive. If an orbit has already halted, later
points on its **unguarded** ray are not native states reached by that
run; the earlier zero remains a witness against both universal claims.

In fact the first formal point with `d(y)>a_27(y)` is reached while all
earlier guards, and its own initial credit, are positive. Otherwise the
preceding block would have had credit 1 and valuation 1, hence already
`d>=2>1`. Such a first failure provides a finite halt certificate: its
exact ancestry is executable, and its ensuing ones-run exhausts its
credit. The predicate at each individual candidate `y` is decidable by
the terminating inverse algorithm. The unbounded universal quantifier
in `(N12)` has **not** been eliminated.

The reset family samples only a narrow subset of possible bad states.
Excluding even every member of this family from the 27 ray would not
prove `(N12)`. Root membership also does not determine credit: both
the chronological word and the chosen initial credit remain necessary.

## 4. Finite classification and exact executable control

For all 33 parameters `R=5,11,...,197`, the checker computes the complete
inverse chain using exact integers. None has root 27. The inverse-chain
length distribution is `{0:17, 1:9, 2:6, 3:1}`. This is
**FINITE-EXACT**, not an extrapolation for the remaining parameters.
The emitted JSON contains every root, source and chronological word;
large integers are decimal strings to preserve them through JSON clients.

Two separate credit experiments are retained explicitly:

1. Start credit 6 at the computed root and transport it along the finite
   chain to `x_R`. All 33 chains have positive odd-section credits; their
   endpoint credits range from 3 to 12. These are generally other roots,
   not statements about the blank source 27.
2. Start credit 6 afresh at `x_R`. The reset gives `(7,1+2^R)`. For every
   sampled `R>=11` (and uniformly every `R=5 mod6` with `R>=11`), seven
   subsequent valuation-1 blocks reach credit zero. The halt occurs
   after eight odd blocks in total, at scalar `1+3^7*2^(R-7)`. For `R=5`,
   the reset followed by its four available ones ends at `(3,163)`;
   this local control makes no claim about its later continuation.

The block below is standalone standard-library Python. Run it directly,
or extract this sole Python fence and execute it with normal and `-O`
Python; their JSON must agree byte for byte. It imports no producer code,
uses explicit exceptions rather than removable assertions, and checks
the inverse formula against literal odd/even scalar macros. Its finite
universes, expected results and hostile controls are inside the block.
Both normal and optimized Python returned `PASS` with identical complete
JSON, SHA-256
`28df9d61d3ddaa4094c5b76fc60543a602ed1fa01a3f34104ef1bbb6578e69cb`.

```python
import itertools
import json
from collections import Counter


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuation(n, p):
    require(n > 0, "valuation domain")
    r = 0
    while n % p == 0:
        n //= p
        r += 1
    return r


def literal_return(x):
    require(x > 1 and x % 2 == 1, "odd scalar domain")
    y = (3*x - 1)//2
    k = 1
    while y % 2 == 0:
        y = 3*y//2
        k += 1
    return y, k


def predecessor(y):
    require(y > 1 and y % 2 == 1, "inverse domain")
    r = valuation(y, 3)
    m = y // (3**r)
    if m % 3 != pow(-1, r, 3):
        return None
    x = ((1 << (r+1))*m + 1)//3
    require(1 < x < y and x % 2 == 1, "inverse decrease")
    require(literal_return(x) == (y, r+1), "literal inverse replay")
    return x, r+1


def ancestry(y):
    reverse_values, reverse_letters = [y], []
    while True:
        edge = predecessor(reverse_values[-1])
        if edge is None:
            break
        x, k = edge
        reverse_values.append(x)
        reverse_letters.append(k)
    return list(reversed(reverse_values)), list(reversed(reverse_letters))


def credit_summary(word, a=6):
    changes = [0]
    for k in word:
        changes.append(changes[-1] + 2*k - 3)
    credits = [a+d for d in changes]
    return {
        "end": credits[-1], "minimum_odd": min(credits),
        "endpoint_inclusive_required_credit": 1-min(changes),
        "first_nonpositive_odd_index": next(
            (j for j, v in enumerate(credits) if v <= 0), None),
    }


def guarded_block(a, x):
    require(a > 0, "block attempted after halt")
    a -= 1
    y = (3*x - 1)//2
    k = 1
    minimum_macro_credit = a
    while y % 2 == 0:
        a += 2
        y = 3*y//2
        k += 1
        minimum_macro_credit = min(minimum_macro_credit, a)
    return a, y, k, minimum_macro_credit


def x_mod(t, modulus):
    # Division by 9 is retained in the power's modulus.
    return ((128*pow(64, t, 9*modulus) + 7)//9) % modulus


# Complete inverse census for odd targets 3..1999. A predecessor is
# smaller than its target, so the source census 3..1999 is complete here.
images = {}
for x in range(3, 2000, 2):
    y, k = literal_return(x)
    require(y > x and y not in images, "growth or injectivity")
    images[y] = (x, k)
for y in range(3, 2000, 2):
    require(predecessor(y) == images.get(y), "complete inverse census")
    values, word = ancestry(y)
    require(values[-1] == y and len(values) == len(word)+1,
            "chain dimensions")
    require(predecessor(values[0]) is None, "terminal root")
    for i, k in enumerate(word):
        require(literal_return(values[i]) == (values[i+1], k),
                "chronological reconstruction")

require(predecessor(27) is None, "27 must be an integer root")
require(literal_return(7)[1] == literal_return(23)[1] == 2,
        "shared-letter hostile")
require(ancestry(7)[0][0] == 5 and ancestry(23)[0][0] == 23,
        "same first letter does not identify a root")
require(ancestry(5)[0][0] == ancestry(7)[0][0] == 5 and
        literal_return(5)[1] != literal_return(7)[1],
        "same root does not identify a first letter")
require(credit_summary([1, 2])["endpoint_inclusive_required_credit"] == 2
        and credit_summary([2, 1])["endpoint_inclusive_required_credit"] == 1,
        "endpoint drift forgets internal debt")
require(guarded_block(1, 5)[:3] == (0, 7, 1)
        and guarded_block(6, 5)[:3] == (5, 7, 1),
        "root label does not specify credit")

# Residue permutations through 3^8, and six exact failure classes.
permutation_sizes = []
for exponent in range(1, 9):
    modulus = 3**exponent
    require(len({x_mod(t, modulus) for t in range(modulus)}) == modulus,
            "3-adic permutation")
    permutation_sizes.append(modulus)
expected_alpha = [1, 3, 6, 78, 24, 429]
failure_classes = []
for r, alpha in enumerate(expected_alpha):
    modulus = 3**(r+1)
    target = ((-1)**(r+1)*3**r) % modulus
    found = [t for t in range(modulus) if x_mod(t, modulus) == target]
    require(found == [alpha], "exact failure residue")
    for j in range(4):
        t = alpha + j*modulus
        x = (128*pow(64, t) + 7)//9
        require(valuation(x, 3) == r and predecessor(x) is None,
                "failure class integer lift")
    failure_classes.append({"r": r, "t_residue": alpha,
                            "t_modulus": modulus,
                            "R_residue": 6*alpha+5,
                            "R_modulus": 6*modulus})

# Unbounded inverse-depth construction: a bounded control of its formula,
# including the endpoint where an r=0 inverse is no longer available.
inverse_prefix_controls = 0
for t in list(range(100)) + [3**j-1 for j in range(1, 7)]:
    x = (128*pow(64, t) + 7)//9
    depth = valuation(t+1, 3)
    require(valuation(x-1, 3) == depth, "exact initial inverse depth")
    y = x
    for _ in range(depth):
        edge = predecessor(y)
        require(edge is not None and edge[1] == 1, "r=0 inverse prefix")
        y = edge[0]
    require(y % 3 != 1, "r=0 prefix has exact endpoint")
    edge = predecessor(y)
    require(edge is None or edge[1] >= 2, "next inverse branch or failure")
    inverse_prefix_controls += 1

# Positive marked-root controls, with every chronological guard retained.
x, a, prefix = 27, 6, []
for j in range(65):
    values, word = ancestry(x)
    require(values[0] == 27 and word == prefix, "marked ray reconstruction")
    require(credit_summary(word)["end"] == a, "marked credit transport")
    require(credit_summary(word)["minimum_odd"] > 0,
            "finite marked prefix guards")
    require(valuation(x-1, 2) <= a, "finite marked prefix fuel")
    if j < 64:
        a, x, k, _ = guarded_block(a, x)
        prefix.append(k)

rows = []
for R in range(5, 198, 6):
    x = ((1 << (R+2)) + 7)//9
    values, word = ancestry(x)
    root = values[0]
    require(root != 27, "finite reset unexpectedly meets marked ray")
    require(valuation(x-1, 2) == 1, "old reset depth")
    target = 1 + (1 << R)
    require(predecessor(target) == (x, 2), "unique reset predecessor")
    require(ancestry(target)[0][0] == root, "reset preserves integer root")
    debt = credit_summary(word)
    require(debt["minimum_odd"] > 0, "finite root-credit6 path guards")

    # Fresh credit6 at x_R is distinct from credit6 at its ancestry root.
    a, y, k, macro_minimum = guarded_block(6, x)
    require((a, y, k, macro_minimum) == (7, target, 2, 5),
            "reset credit and internal minimum")
    ones = min(R-1, 7)
    for _ in range(ones):
        a, y, k, _ = guarded_block(a, y)
        require(k == 1, "reset ones-run")
    if R == 5:
        require((a, y, ones) == (3, 163, 4), "short reset boundary")
        outcome = "LOCAL_WINDOW_SURVIVES_LATER_UNDECIDED"
    else:
        require(a == 0 and ones == 7 and
                y == 1 + 3**7*(1 << (R-7)), "exact reset halt")
        outcome = "HALT_AFTER_8_ODD_BLOCKS"
    rows.append({
        "R": R, "x": str(x), "root": str(root),
        "J": len(word), "K": sum(word), "root_to_x_word": word,
        "credit6_at_root": debt,
        "credit6_fresh_at_reset": {
            "outcome": outcome, "blocks_checked": 1+ones,
            "final_credit": a, "final_x": str(y)},
    })

depths = Counter(row["J"] for row in rows)
require(dict(depths) == {0: 17, 1: 9, 2: 6, 3: 1}, "finite root depths")
require(len(rows) == 33, "reset census size")
require(min(row["credit6_at_root"]["end"] for row in rows) == 3
        and max(row["credit6_at_root"]["end"] for row in rows) == 12,
        "root-credit endpoint range")

# All finite words of length<=4 over 1..4: exact required-credit gate.
word_controls = 0
for length in range(5):
    for word in itertools.product(range(1, 5), repeat=length):
        required = credit_summary(word)["endpoint_inclusive_required_credit"]
        for a in range(0, required+2):
            summary = credit_summary(word, a)
            require((summary["minimum_odd"] > 0) == (a >= required),
                    "endpoint-inclusive debt iff")
        word_controls += 1

print(json.dumps({
    "status": "PASS", "scope": "ELEMENTARY_PROOFS_PLUS_FINITE_CONTROLS",
    "inverse_odd_target_count": 999, "marked_root_prefix_returns": 64,
    "residue_permutation_sizes": permutation_sizes,
    "exact_failure_classes": failure_classes,
    "initial_inverse_prefix_controls": inverse_prefix_controls,
    "finite_word_controls": word_controls,
    "reset_count": len(rows), "reset_root27_count": 0,
    "reset_inverse_depth_histogram": dict(sorted(depths.items())),
    "reset_rows": rows,
    "open": "No universal marked-root debt bound or blank halting decision",
}, sort_keys=True, separators=(",", ":")))
```

The new map is `odd integer -> (root, exact finite forward word from root)`.
It preserves and decides integer-ray membership and reconstructs the input;
discarding the word would lose position and credit transport. The cheapest
decisive controls are a failed predecessor congruence, a literal forward
replay, and an internal-prefix debt check. The next unresolved target is a
uniform debt/fuel inequality on the fiber `root(y)=27`, not a density
statement across unrelated roots.
