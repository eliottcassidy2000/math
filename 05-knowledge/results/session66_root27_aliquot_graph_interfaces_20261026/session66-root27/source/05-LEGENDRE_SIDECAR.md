# Legendre-prime tournament: known extension theorem, separate prefix question

2026-10-26. **CITED known result + independently rederived elementary
argument:** the prime Legendre tournament is the countable random
tournament. **FINITE-EXACT:** eight prescribed extension patterns below.
The numerically ordered prime-prefix hypothesis is **not settled by this
audit**. No novelty, effective prime-size bound or enumeration-speedup
claim is made.

## Sources and inheritance

The incoming [QR_TOURNAMENT.md](../../analogy_reduction_atlas_20261025/QR_TOURNAMENT.md)
defines the tournament correctly and labels its prime-prefix frequency
claim HYP-9375. Its structural conclusion can be strengthened using
established literature:

- Peter J. Cameron, *The Random Graph Revisited* (2001),
  [primary paper PDF](https://www.math.uni-bielefeld.de/~rehmann/ECM/cdrom/3ecm/pdfs/pant3/camer.pdf),
  PDF pp. 1–3: the prime construction is proved using CRT and Dirichlet;
  p. 3 explicitly gives the analogous random tournament using primes
  congruent to minus one modulo four.
- Gareth A. Jones, *Infinite Paley graphs*, arXiv:1912.02709v1 (2019),
  [primary paper PDF](https://arxiv.org/pdf/1912.02709), section 4 and
  section 6, Remark 4: the same construction and tournament conclusion
  are explicit. These prime vertices should not be confused with the
  vertices of a single finite-field Paley tournament.

Both sources were inspected on 2026-10-26. The construction is known;
no claim of earliest priority is needed here. The supplied
[Session26 synthesis](../../session26_random_bridges_20261022/SYNTHESIS.md)
does not state this Legendre-tournament theorem. Its section 5 provides
a useful logical parallel: complete congruence realization need not
give the desired numerical spectrum or density. This is an analogy of
quantifiers, not a transported proof.

## Exact extension property, with the sign checked

Vertices are the primes congruent to 3 modulo 4. For distinct vertices,
orient `p -> q` iff `(p/q)=+1`, where parentheses denote the Legendre
symbol. Quadratic reciprocity gives `(p/q)=-(q/p)`, so exactly one arc
occurs.

Fix distinct existing primes `p_1,...,p_m`, prescribe signs
`epsilon_i=(p_i/q)` in `{+1,-1}`, and impose any lower bound on q.
Choose a nonzero residue `a_i mod p_i` with `(a_i/p_i)=-epsilon_i`.
CRT gives one class

    q=3 (mod4),    q=a_i (mod p_i),
    modulo M=4*product_i p_i.

Its representative is coprime to M. Dirichlet supplies infinitely many
prime q in that class, hence one above the chosen lower bound. Reciprocity
then gives every requested sign. This is the full finite extension
property, not merely the appearance of a few tournament classes.

Induction embeds any finite tournament. Back-and-forth between two
countable tournaments having this extension property extends finite
partial isomorphisms one vertex at a time; hence the isomorphism type is
unique. This is the countable random tournament. These deductions use
the stated classical results and are the known construction cited above.

The minus sign matters: to require `3 -> q`, one needs `q=2 (mod3)`;
q=11 works. The opposite residue gives q=7 and the opposite arc.
Allowing primes congruent to 1 modulo 4 destroys the orientation rule:
for example `(5/29)=(29/5)=+1` gives two proposed arcs, not a tournament.

## Why numerical prime prefixes remain a different obligation

The extension theorem makes every finite pattern available, with
infinitely many witnesses. It does not estimate the size of a witness,
the joint frequency of growing sets of prime vertices, or the error in
the increasing-prime exhaustion. Calling the abstract isomorphism type
random does not make its prime-labelled edges independent coin flips.
This audit therefore neither promotes nor refutes HYP-9375. It does not
claim that this separate analytic question is new or unknown in all
literature.

Here is a structural hostile to deriving prefix frequencies from the
extension property alone. Enumerate a countable random tournament in
stages. First insert the least unused vertex from a fixed exhaustive
list. Then append a very large transitive block, each new vertex losing
to every earlier vertex. The extension property permits every insertion.
No cyclic triangle uses a vertex in this appended block. Make the block
large enough that the cyclic-triangle density at the stage endpoint is
below `1/j` at stage j. The next mandatory vertex adds at most
`binom(N,2)` cyclic triangles to a prefix of size N; its density increment
is at most `3/(N+1)`. During the subsequent transitive block the density
only decreases. The resulting exhaustive enumeration therefore has
cyclic-triangle density tending to zero, despite the same extension
property and abstract random tournament.

This is a change of enumeration, **not** a refutation for increasing
prime order. The numerical order is precisely the additional sidecar
that abstract universality forgets.

Typed connection: a prescribed finite tournament maps, by successive
CRT extensions, to a prime-labelled induced tournament with all arcs
preserved. The construction loses any small numerical label bound and
any canonical isomorphism-class encoding. It supplies no smaller A000568
recurrence: all prescribed arc information remains in the growing CRT
constraints.

## Small reproducible control

This standalone block checks all eight extensions of the fixed vertices
`(3,7,11)`, using exact trial primality and Euler's criterion. It also
checks a cyclic and a transitive triple, the reciprocity sign hostile,
and the prime-1-mod4 failure. Its bounded search witnesses do not replace
Dirichlet's theorem or test asymptotic prefix frequencies.
Normal and optimized Python both returned `PASS` with identical output
SHA-256 `47009071eef4603b5426dd46da7da609d03540a1b94d57c07d05d37fda4297c3`.

```python
import json
from itertools import product
from math import gcd, isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def legendre(a, p):
    require(prime(p) and p % 2 == 1, "odd prime modulus")
    a %= p
    if a == 0:
        return 0
    value = pow(a, (p-1)//2, p)
    require(value in (1, p-1), "Euler criterion")
    return 1 if value == 1 else -1


anchors = (3, 7, 11)
receipts = []
for signs in product((-1, 1), repeat=len(anchors)):
    residue, modulus = 3, 4
    for p, sign in zip(anchors, signs):
        local = next(a for a in range(1, p) if legendre(a, p) == -sign)
        residue += modulus*((local-residue)*pow(modulus, -1, p) % p)
        modulus *= p
    require(modulus == 924 and gcd(residue, modulus) == 1,
            "reduced CRT class")
    candidates = (residue+k*modulus for k in range(1000))
    q = next((q for q in candidates if q > max(anchors) and prime(q)), None)
    require(q is not None and q % 4 == 3, "finite prime witness")
    require(tuple(legendre(p, q) for p in anchors) == signs,
            "prescribed outgoing signs")
    require(tuple(legendre(q, p) for p in anchors) == tuple(-s for s in signs),
            "reciprocity sign")
    receipts.append({"signs": signs, "residue": residue,
                     "modulus": modulus, "prime": q})

require([r["prime"] for r in receipts] ==
        [463, 211, 199, 2719, 2003, 827, 2663, 563], "closed finite universe")
require(all(legendre(a, b) == 1 for a, b in ((3,11),(11,7),(7,3))),
        "cyclic triple")
require(all(legendre(a, b) == 1 for a, b in ((7,19),(19,3),(7,3))),
        "transitive triple")
require(legendre(7,3) == 1 and legendre(3,7) == -1,
        "wrong CRT sign reverses the requested arc")
require(legendre(5,29) == legendre(29,5) == 1,
        "prime1mod4 pair cannot be tournament-oriented this way")
print(json.dumps({"status":"PASS", "extension_patterns":receipts,
                  "scope":"FINITE_SIGN_AND_CRT_CONTROLS_ONLY"},
                 sort_keys=True, separators=(",", ":")))
```

A targeted check of the October 5 tree/credit note found the ternary
grammar and its ballot interpretation, but no stated Doob transform or
Galton–Watson conditional-extinction law. This limited negative search
is not a claim of novelty or absence elsewhere in the repository.
