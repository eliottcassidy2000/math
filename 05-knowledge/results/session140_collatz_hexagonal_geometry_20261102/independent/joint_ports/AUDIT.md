# Bounded independent review of joint native ports

**Accepted:** J1–J7 reconstruction identities for the supplied finite rooted
odd inverse tree, with the root as a **stopping boundary**. The final
leaf-Gram rank statement is accepted for **h≥1**, not h=0. This does not
promote the status of the full concurrent Collatz packet.

Source: [`JOINT_PORTS.md` at commit 2bf016747f984c148248e4e7cefa5cef2dfa9f22](https://github.com/eliottcassidy2000/math/blob/2bf016747f984c148248e4e7cefa5cef2dfa9f22/05-knowledge/results/collatz_hexagonal_growth_20261102/JOINT_PORTS.md).
Git blob: `676ef198ac8f4bfe754230a2fb41edda7d485645`.
SHA256: `012e51e6a918163cfbe0c5a66c259e3de3f9b24c411c20174aa858d48793d83d`.
The exact small source is retained as `SOURCE.md`; no sibling concurrent
packet is needed to replay this self-contained audit.

## Proof review

- **Shared tails (J1–J3):** the even states terminating at an odd parent q
  have precisely odd part q. Distinct parents' chains cannot overlap; sibling
  tails overlap exactly to their smaller exponent. Attaching each odd child
  at its own exponent produces the actual native path union, not independent
  macro-edge subdivisions. The vertex, edge and duplicate-count formulas follow.
- **Meetings (J4):** if neither source is their odd LCA q, the final two
  tails first coincide at `2^min(a,b) q`, exactly min(a,b) steps before q on
  both routes. Any earlier coincidence would force an earlier common odd
  ancestor. If one source is q, the stopped-tree meeting is q. The clocks
  are asynchronous first-hit times, not a synchronous collision assertion.
- **Generating recurrence (J5–J6):** the parent/even chain and the two child
  subtree unions are disjoint by odd-part classes. Their native first-hit
  times shift by a+1 and b+1. Thus the recurrence counts each integer state
  exactly once and retains the source label, rather than a bare residue phase.
- **Tail Gram (J7):** nested indicator vectors have inner product min(a,b).
  Splitting the stage coordinates into consecutive exponent blocks yields
  the lower-triangular factorization and determinant product. Nonnegative
  weights retain full rank exactly when every required block has positive
  total weight. Arithmetic states can remain distinct when zero weights
  erase a measured distinction.

**Boundary correction:** the subsequent leaf-Gram rank `2^h` uses a private
first outgoing edge for each leaf, which exists only for h≥1. At h=0 the
sole leaf is the stopping root 5; its edge-indicator vector is zero, so its
Gram matrix is `[0]` and has rank 0, not 1. No J1–J7 formula is affected.
The positive-depth identity-submatrix argument is valid.

## Exact hostiles and controls

The paths from 13 and 53 to 5 first meet at 40, with clocks (1,3), giving
native distance 4 rather than independently subdivided distance 10.
The paths from 853 and 53 first meet at 160, with clocks (5,1).
Yet 13 and 853 are legal 3-unit odd sources of parent 5 with exponents 3
and 9, equal modulo six. A fixed sibling therefore distinguishes responses
lost by the periodic exponent summary. At root 5 and depth one the exact
polynomial is `1+t+t²+t³+2t⁴+t⁵+t⁶`, with eight native vertices rather than
eleven independent-subdivision vertices.

The independent verifier constructs a valid 127-node odd tree by the two
safe exponents in 1 through 6, then checks every source pair against literal
native paths and native-tree BFS. It additionally compares source-dependent
polynomials at four roots/depths 0–4, exact min-kernel determinants, weighted
rank collapse, the two stated hostiles, the positive-depth private-edge
identity minor, and the depth-zero exception. This family is not claimed
to reproduce the concurrent packet's complete SF5 census.

Normal and Python `-O` outputs agree, as do relocated copies. All checks use
explicit exceptions. No convergence conclusion, new lattice realization,
whole-packet audit, Git mutation, or remote write follows from this review.
