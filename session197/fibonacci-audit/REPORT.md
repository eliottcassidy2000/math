# Independent review of the Fibonacci start-one enumeration

Verdict: **PASS**, after the inverse-map wording was repaired. The count formulas and their stated subclass scope are sound.

## Reviewed snapshot

- `session197-lowcycle/FIBONACCI_CLASS.md`: SHA256 `d95a133e360b96589009ebd8a1c8449525c4a2847a88b33a4b3f2b839739e63e`.
- `session197-lowcycle/verify_fibonacci_class.py`: SHA256 `01154015e10dd720a1ba217e6411b1bb5b8a9b211884b151145bf818621917b0`.
- `session197-lowcycle/fibonacci_class_certificate.json`: SHA256 `a121e50a85c73f074a2a22d2154fb3722d183009695ce1bc7c7bdf492dea34f4`.

The stated structural dependency is `ROOTED.md`, SHA256 `ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e`. Its canonical source-SCC decomposition and the bridge/start-one equivalence were read and checked as needed. The source checker imports only the named `verify.py` graph utilities. No unrelated repository search was used.

## Mathematical findings

1. Deleting the marked root makes the source strong component, its unique vertex beaten by the root, and the transitive tail intrinsically recoverable. Thus the recurrence counts rooted isomorphism classes, not construction histories. Tail words have exactly `2^(w-1)` possibilities when nonempty, and the sole two-vertex nonstrong output is correctly excluded.
2. The generating function is `B(z)=z+z^3/(1-3z+z^2)`. Its initial values and even-Fibonacci indexing are correct: `b_n=F_(2n-4)` for `n>=3`, with `b_1=1,b_2=0`.
3. Requiring the marked root to have outdegree one leaves one tail orientation for each rooted source piece. Hence `c_n=sum_(m=1)^(n-1)b_m`, for `n>=3`.
4. The distinction between marked and unmarked counts is essential and correctly handled. At order at least four there are one or two outdegree-one vertices, and every automorphism fixes them individually. The unique directed Hamiltonian cycle then forces rigidity. At order three the directed triangle is separately handled; its three possible marks form one rooted isomorphism class.
5. The two-low-degree-vertex construction is bijective on isomorphism classes. The direction of the edge between the two low-degree vertices identifies which one to delete, so no factor of two is missing.
6. The inverse proof originally claimed path avoidance for all targets, which failed for a path ending at the other low-degree vertex. It was repaired: first target a vertex outside `{u,v,w}`, then append its arc into `v`. The neighbor sentence now correctly excludes both `v` and `w`. This establishes strongness of the deleted graph and preserves `d_n=c_(n-1)`.
7. Consequently `a_n=c_n-c_(n-1)=b_(n-1)=F_(2n-6)` for `n>=4`, with `a_3=1`. These formulas count the C=1 **bridge subclass**, not all C=1 tournaments. No unproved bridge-exhaustion assertion enters the enumeration.

## Independent computation

- Replayed the original checker in an isolated directory, normally and with `python -O`. Outputs agreed exactly and reproduced the frozen certificate hash.
- Implemented a different canonicalizer based on the unique directed Hamiltonian cycle. For rooted objects it anchors the cycle at the root; for unrooted objects it minimizes over cyclic rotations. Its applicability is checked by a direct cycle DFS, rather than assuming a unique cycle.
- Independently generated and checked all 610 rooted class representatives through order nine. The rooted, outdegree-one-marked, unrooted, and two-outdegree-one counts all agree with the author's inventory.
- Checked automorphism rigidity directly by testing all rotations of each unique cycle, and checked both directions of the two-vertex bijection.
- Exhaustively examined all 33,867 labelled tournaments of orders one through six, independently counted root-start paths using a reverse-graph path dynamic program, and recovered exactly the constructed rooted classes. This is independent of the author's unlabelled extension generator and its color-refinement canonicalizer.

Detailed outcomes are in `INDEPENDENT_RESULTS.json`; the independent checker is `independent_check.py`. The author's replay additionally covers the complete unlabelled universe through order seven and the formal recurrence through order one hundred.

No unresolved mathematical issue remains. The author removed the harmless duplicated strongness sentence; that editorial-only final snapshot is bound to the proof hash above.
