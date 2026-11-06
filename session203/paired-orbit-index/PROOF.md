# A single integer index for every paired O/Z carrier

Status: **PROVED, independently audited in the stated scope.** This is an exact all-order arithmetic description of the four physical paired carrier orbits from `dependencies/CARRIERS.md`, not universal tournament coverage. The all-target staggered-profile criterion remains open.

## 1. Exact parametrization

The inherited paired-carrier row is `(X,e,V,d)` and its physical moves are

    O: (X,e,V,d) -> (X,0,2V,d),
    Z: (X,e,V,d) -> (X+V-e,d,2V-d,d).

Fix an initial row `(X0,e0,V0,d)` with `0 <= e0 <= d`, `d>=1`, and `A=V0-d>0`. Let `g(1)=1`, `g(2h)=g(h)`, and `g(2h-1)=g(h)+h-(h mod 2)` for `h>=2`, as in `dependencies/CENTRAL.md`.

**Theorem.** The set of X coordinates obtained by all finite O/Z words is exactly

    X(u) = X0-V0 + A*u + d*g(u) + (d-e0)*(1-(u mod 2)),  u>=1.       (1)

For a physical word, read its letters chronologically as binary digits from least to most significant, with O=0 and Z=1; their integer value is `k=u-1`. Conversely the binary expansion of every `u-1` gives such a physical word. For `u=1` use the empty word. Leading zeroes in this binary integer are trailing O letters in the word and do not change X. Different integer indices are not asserted to yield different X values.

The index 1 is allowed here because it denotes the empty word on an already physical carrier. It does not add the forbidden singleton unit source to the two-index central family.

**Proof.** For binary `k`, define `q(0)=0` and

    q(2h)=q(h),
    q(2h+1)=q(h)+h+(h mod 2).

Induction on the word length gives its final first coordinate as

    X0+V0*k-d*q(k)-e0*(k mod 2).                              (2)

Indeed an initial O leaves X fixed, replaces V by 2V and e by zero, and leaves the suffix integer h, yielding the even recurrence. An initial Z adds V-e, replaces V by 2V-d and e by d, and yields the odd recurrence. The d coordinate never changes. An empty word gives k=0.

The displayed recurrence for g gives, for every k>=0,

    q(k) = k+1-(k mod 2)-g(k+1).                              (3)

At k=0 this is immediate. For k=2h, h>=1, use
`g(2h+1)=g(h+1)+h+(h mod 2)`; for k=2h+1 use `g(2h+2)=g(h+1)`. Substitution verifies both recurrences for q. Equation (1) now follows from (2)-(3), with `u=k+1`. Every finite binary integer specifies an actual word, proving equality of value sets rather than just containment. QED.

## 2. The four explicit arithmetic families and a complete inverse

The four initial rows are `(677,27,1061,49)`, `(759,29,927,29)`, `(1547,101,1691,109)`, `(821,33,1079,33)`. Thus their complete X spectra are respectively

    1012*u + 49*g(u) - 384 + 22*(1-(u mod 2)),
     898*u + 29*g(u) - 168,
    1582*u +109*g(u) - 144 +  8*(1-(u mod 2)),
    1046*u + 33*g(u) - 258,                         u>=1.    (4)

Every value of each formula supplies both staggered residue certificates, with the inherited fixed starts and j bounds. Therefore both `3X+2` and `3X+4` are actual tournament Hamilton-path counts. This conclusion inherits the literal physical carriers and the four-vertex compiler; the proof here introduces no new freely chosen profile coordinates.

Since g is a positive odd integer by its recurrence, (1) gives

    X(u) >= X0 + A*(u-1).                                    (5)

Thus membership of a requested integer X in any one orbit has the complete finite inverse

    1 <= u <= 1+floor((X-X0)/A),                             (6)

with an empty search if X<X0. Testing (4) is still necessary. This is stronger than a bound indexed by the length of a word: it tests each canonical integer encoding exactly once. The coefficient A, not V0, is the safe inverse denominator. The formula need not be monotone in u.

## 3. A sharper rigorous limitation of this carrier bank

Let C(N) be the number of distinct odd bases X<=N supplied by the union of the four orbits. Equation (5) proves

    C(N) <= sum_i max(0, 1+floor((N-X_i)/(V_i-d_i))).         (7)

Hence its upper density relative to the positive odd integers is at most

    2*(1/1012+1/898+1/1582+1/1046)
      = 693609731/93988567442 < 3/400.                       (8)

In particular this bank cannot certify more than 0.75 percent of all odd bases asymptotically; the previously recorded bound for this same bank was 1.5 percent. This is an upper bound, not an assertion that the density exists or that orbit values are distinct. It shows that these four orbits cannot by themselves prove the full staggered criterion. It does not exclude any target from the unrestricted tournament spectrum or even from the full bounded-start profile spectrum.

The arithmetic parameterization is genuinely all-order. Finite verification checks the formula and its boundaries but does not establish any unproved occupancy assertion. It also does not show that either prime atom track is covered infinitely often. Exponential subfamilies alone do not justify assertions about infinitely many primes.

## 4. Reproduction

`verify.py` compares (1) with direct row dynamics for every word of length at most14 from each physical carrier and two boundary-control rows. It separately checks the recurrence identity (3), oddness and positivity of g, the complete inverse bound for each tested word, the explicit four formulas, and the exact rational density bound. Run in normal and optimized Python and compare the JSON outputs. No input graph census or target coverage assumption is used by this arithmetic verifier.
