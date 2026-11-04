# An all-order failure of one-letter ternary closure at H=13

**PROVED, with exact finite certificates and a cited size cutoff.** In the
four-letter two-apex grammar, there is no word w with H(T(w))=13 for which
one-letter appending can produce both

    H(T(wc)) = 41 = 3*13+2,
    H(T(wd)) = 43 = 3*13+4.                            (1)

Moreover no word with H=13 has a one-letter appended child of count 39.
The statements permit arbitrarily long parent words, neutral padding,
non-strong parents, and arbitrary choices of the two apexes initially.
They retain that same apex pair and chain order in the children. They do not
exclude re-marking between branches, inserting at other positions, appending
multiple letters, or changing old incidences.

## 1. Every H=13 presentation has an exact finite-core normal form

By the SCC product formula and primeness of 13, any tournament with H=13 has
exactly one nontrivial strong component C with H(C)=13. Every other component
has H=1 and hence is a singleton: a nontrivial strong tournament has at least
three Hamilton paths, for example by rotating a Hamiltonian cycle, or by the
cited Busch lower bound.

Busch's inequality H(C)>=5^((|C|-1)/3) forces |C|<=5, since

    5^5 = 3125 > 13^3 = 2197.                          (2)

The component C must contain both named apexes. Otherwise it has a transitive
one-vertex deletion and order at most five. The complete one-apex presentations
at these orders have H-values {1,3,5,9,11}, not 13. `verify_append.py` checks
this directly by full vertex-subset path counting on all 31 binary words of
length at most four. This is consistent with, but needs less than, the
inherited S191 all-order one-apex exclusion of 13.

The chain names in C form an interval: if x_i,x_j belong to C with i<l<j,
then x_i->x_l->x_j brings x_l into C. Every earlier chain name belongs to an
earlier singleton SCC and therefore dominates both apexes, so its letter is
0. Every later chain name lies after C and is dominated by both apexes, so its
letter is 3. Thus every H=13 word is exactly

    w = 0^p u 3^q,          p,q>=0,                    (3)

where u presents C. No word of length at most two has count 13, and the
complete length-three list is

    u in {120, 212, 230, 302, 310, 321}.                (4)

All six words are strong. They are named presentations, not claimed to be
six different isomorphism classes. Enumerating all 4^3 words proves (4),
with independent full vertex-subset DP.

The initial 0^p is a source chain even after any future letters are appended,
so it can be erased without changing any appended child's H. Therefore the
entire all-order question reduces exactly to six cores and one nonnegative
tail clock q. In particular arbitrarily many neutral sinks have not been
silently bounded by the strong size cutoff.

## 2. Exact children at every tail clock

The four columns below refer to appended letters c=0,1,2,3. At q=0 the
children are respectively

| u | c=0 | c=1 | c=2 | c=3 |
|---|---:|---:|---:|---:|
| 120 | 33 | 23 | 29 | 13 |
| 212 | 37 | 29 | 25 | 13 |
| 230 | 37 | 29 | 29 | 13 |
| 302 | 33 | 29 | 25 | 13 |
| 310 | 37 | 23 | 37 | 13 |
| 321 | 43 | 25 | 37 | 13 |

For every q>=1 put x=2^q. The exact all-clock children are

| u | c=0 | c=1 | c=2 | c=3 |
|---|---|---|---|---|
| 120 | 7+32x | 13+14x | 13+18x | 13 |
| 212 | 9+33x | 13+18x | 13+15x | 13 |
| 230 | 5+36x | 13+18x | 13+18x | 13 |
| 302 | 5+33x | 13+18x | 13+15x | 13 |
| 310 | 3+40x | 13+14x | 13+26x | 13 |
| 321 | 7+41x | 13+15x | 13+26x | 13 |

These are not fitted extrapolations. On the invariant reachable state space,
the exact operator identity in `TWO_APEX.md` gives

    M_3^4 - 6M_3^3 + 11M_3^2 - 6M_3 = 0.

Hence every positive-clock child sequence has the form A+B*2^q+C*3^q.
The exact q=1,2,3 values uniquely determine these coefficients, yielding the
displayed rows with C=0. The verifier reconstructs an invariant reachable
span and checks the operator identity again; all 96 child values for
q=0,1,2,3 are independently recomputed by full subset DP. The separate q=0
row is essential: it generally does not equal the positive-clock formula
at x=1.

The c=1 and c=2 columns can also be read as the converse of the inherited
sole-entry wrapper: the new last vertex beats only b or only a among old
names, and the nonempty sink chain supplies the dyadic clock. The explicit
three-column certificate does not assume that those endpoint responses
coincide for different physical core presentations.

## 3. The obstruction and its strongest survivor

For q>=2, every child except c=3 exceeds 43, because every displayed positive
slope is positive and already exceeds 43 at x=4. Letter 3 always leaves H=13.
For q=0 the first table contains no 39 or 41. At q=1 the child rows are

    120: (71,41,49,13)     212: (75,49,43,13)
    230: (77,49,49,13)     302: (71,49,43,13)
    310: (83,41,65,13)     321: (89,43,65,13).

None contains 39, and none contains both 41 and 43. Together with (3), this
proves both all-order claims.

**Individual child existence survives.** Count 41 is an appended child of
1203 or 3103. Count 43 is an appended child of 321, 2123, 3023, or 3213.
The failure is simultaneous return from the same named parent, not global
attainability of either scalar target.

**Allowing both ends defeats this particular obstruction.** The actual
parent word 1203 has H=13. Appending 1 produces 12031 with H=41, whereas
prepending 3 produces 31203 with H=43. These three counts are independently
checked by subset DP. The lower child is strong, while the parent and upper
child have a final neutral sink and are nonstrong. Thus this example does not
supply two strong children on those literal names. It nevertheless shows why
the one-ended hypothesis cannot be dropped from the theorem. Repeated
bidirectional closure remains open.

## 4. Reproduction and dependencies

    python3 verify_append.py
    python3 -O verify_append.py

Both runs use explicit exception checks and write the deterministic
`append_receipt.json`. No bank or third-party Python package is used.
Dependencies are the proved run-assignment transfer, its exact invariant-space
certificate, SCC multiplicativity (proved above by ordered concatenation),
and the same cited Busch bound used by the main packet. No literature-priority
claim, new global spectrum gap, or impossibility of a richer grammar follows.
