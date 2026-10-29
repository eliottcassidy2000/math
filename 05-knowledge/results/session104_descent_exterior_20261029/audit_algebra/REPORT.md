# Independent algebra audit — Session 104

**Verdict: PASS for the stated algebraic claims, with provenance and abc boundaries retained.**

Read targets: `../exterior/K4_EXTERIOR_BOUNDARY.md` and
`../summands/SUMMAND_SQUARES.md`, followed by the newly recovered
`../summands/THREE_CHAINS_RECOVERED.md`. The independent checker imports only Python
standard-library modules, independently encodes the printed transition table,
and does not import or execute either producer checker. Explicit exception
guards remain active under optimization. No publication was performed.

## Exterior/automaton findings

* Moore partition refinement reproduces six minimal states for the
  safe/singular/illegal output and four for acceptance alone. The six-state
  letters have rank three, but `ee` has rank two. The quotient's letters have
  rank two and their exterior squares rank one: these are different models.
* The obstruction is universal, not experimental: for every linear map T over
  any field, `im(Λ²T)=Λ²(im T)`, so its rank is `r(r−1)/2`. A four-dimensional
  map gives only ranks 0, 1, 3, 6. Functoriality and conjugacy invariance then
  rule out a transition-compatible exterior-square representation of this
  six-state machine, because the word `ee` would require rank two.
* Independent oriented K4 incidence has rank three over Q and F2. Over Q,
  the computed Λ²H cycle basis has dimension three and complements the cut
  space; `BBᵗ=4I−J`. Over F2, `rank(BBᵗ)=1`, so the cut/cycle intersection
  has dimension `rank(Bᵗ)−rank(BBᵗ)=2`. Thus the characteristic-zero splitting
  must not be transferred unchanged to characteristic two.
* The Euclidean four-dimensional Hodge star is an involution with two
  three-dimensional eigenspaces and exchanges cuts with cycles. Those
  eigenspaces are not the cut/cycle summands. The familiar vector-valued
  cross product requires the metric/orientation identification in dimension
  three, unlike the canonical wedge. These checks do not manufacture a
  grammar-compatible map or a base space for a bivector field.
* Plücker cancellation is verified symbolically in four indeterminates, not
  merely sampled. Ordered potentials make the three opposite-edge products
  positive; `e01+e23` is a nondecomposable negative control. No coprimality or
  radical/height estimate follows. This supplies no abc theorem.
* For `C=[[1,0,1],[1,1,1],[1,1,0]]`, exact determinant expansion independently
  gives `det(C)=−1`, `char_C(z)=z³−2z²−z+1`, and
  `char_(Λ²C)(z)=z³+z²−2z−1`. The explicit volume-coordinate identity
  `Cᵗ J (Λ²C)=−J` checks the contragredient twist. The exterior square of
  `1⊕C` is exactly `C⊕Λ²C` in the chosen ordered basis. This is synthetic,
  not the actual legality or acceptance representation.

### Proof audit: the native observer's genuine wedge carrier

The final §3 addition is approved independently against the native definition
in Session103's README: entries count spanning covers by two ordered,
nonempty, disjoint directed paths, with singleton paths allowed. Simultaneously
swapping path labels swaps both ordered start and end ports and is an
involution of covers. Consequently `R M2 R=M2`, equivalently `M2 R=R M2`.

For n=4 the twelve ports comprise six transposed pairs. Each pair contributes
one rational +1 eigenvector and one −1 eigenvector; hence two invariant
six-dimensional rational sectors. The map
`ei∧ej ↦ e_(i,j)−e_(j,i)` identifies the antisymmetric sector with Λ²Q⁴.
It intertwines vertex re-marking with Λ²Pπ, including the correct sign when
the order of a pair reverses. On each integral two-port block, the sum/difference
basis matrix has determinant −2, so its full twelve-port lattice index is
`2^6`, exactly as stated. The integer antisymmetric lattice itself has the
printed difference basis; the index warning concerns combining both sectors.

The expanded bridge vector `(ea−es)∧er` agrees termwise with the native
ordered-pair certificate stated in Session103. At general n the dimension is
`binomial(n,2)`, and injections of labelled vertex sets induce injective wedge
maps. These maps alone supply no intertwining equation with the new observer
and therefore no automatic kernel persistence. The additional last-forbidden
bridge/sink proof is essential. This native observer sector is a real wedge
carrier, **not** the six-state Moore machine, and run-swap commutation does
not imply that its restricted operator is Λ² of a vertex operator. No new
computation was needed: this addition is certified by the explicit bijection,
two-by-two lattice calculation, and functorial definition, while the existing
101 coded checks and receipts remain unchanged.

### Independently audited integral dual / seventh-cyclotomic addition

The subsequently added §5 is correct. The checker independently verifies the
printed matrix S has determinant one and `W S=S(C−I)`, with `W=Λ²C`.
Both stated cyclic-basis matrices have determinant −1 and S maps the D=C−I
basis to the W basis. Thus this is integral similarity, stronger than merely
matching characteristic polynomials.

The exact polynomial identity

```
z³ q(z+z⁻¹) = 1+z+z²+z³+z⁴+z⁵+z⁶
```

is checked coefficientwise. A primitive seventh root therefore yields a q-root
`ζ+ζ⁻¹`; pairing the six primitive roots by inversion gives three distinct
real roots `2cos(2πk/7)`, k=1,2,3. Their shifts by one are the C-eigenvalues;
the positive largest is `1+2cos(2π/7)`, consistent with inherited growth.
No new counting theorem or numerical recognition is needed for this identity.

Exact matrix multiplication confirms `(D²−2I)(D+I)=−I`, and three successive
applications of `X↦X²−2I` return D, with all three intermediate matrices
distinct. At each q-root, `φ(x)=−1/(x+1)=x²−2` follows by rearranging q(x)=0.
The polynomial's values at 0 and −1 are nonzero, excluding the poles of
φ and φ². Its values at 2 and −1 exclude fixed points; thus the real roots
have exact period three. Equivalently, squaring ζ cycles their inversion
classes. This is a lawful cyclotomic / integral exterior-dual connection.
**Only the synthetic counting lift is integrally C⊕(C−I); the six-state
legality machine remains excluded by the earlier rank-two obstruction.**

A quick scope review of the current packet README found the exterior, square,
recovered-three-chain, and abc summaries consistent with these audited claims.
The descent/positivity portions remain the other independent audit's scope.

## Consecutive-block findings

The checker proves the all-k coefficient identities using exact polynomial
arithmetic for the sum formulas, independently of producer code. Since the
signed block has one more left term, its constant-in-power count is one.
Its first moment is `−k²`, and second moment is `−k²(2k+1)`. Hence

```
D1(a;k) = a−k²
D2(a;k) = a²−2k²a−k²(2k+1)
         = (a+k)(a−k(2k+1)).
```

For k≥1 the linear root k² is positive. The quadratic roots are −k<0 and
k(2k+1)>0, proving uniqueness among positive real bases (and therefore
positive integer bases). With step d>0, replace k² by dk² and the roots by
−dk and dk(2k+1). Direct enumeration supplies 2,196 independent arithmetic
progression controls, alongside positive/negative root and invalid-domain
controls.

Subtracting the displayed quadratic at a+1 and a yields `2(a−k²)+1`;
at a+1/2 and a−1/2 it yields `2(a−k²)`. Thus the forward and centered
finite-difference identities hold identically. Rational-value tests corroborate
this derivation. The balanced bases differ by k(k+1), so transporting a
pointwise square-balanced equality directly to the linear-balanced base
would be invalid. Only the defect-polynomial relation has been proved.

## Source and scope boundary

The second source-recovery pass now supersedes the earlier open status. Its
new note identifies the historical condensed framing and distinguishes closure
holes {1,4,6} from the square-sum graph's founders {1,2,4}. This audit read
that note and independently verified Q12's three component vertex sets,
three-path structure for 4≤n≤12, component counts 3/2/1 at n=12/13/14,
and exact merger neighbors. It did not independently fetch the remote
historical blobs; source attribution remains supported by the source-recovery
workstream's retrieval rather than by this checker.

The all-n attachment proof is valid: for n≥6, the next square above n is
at most n+2√n+1<2n, so a previous neighbor exists. Thus connectivity persists
from 14. The excluded-square-diagonal criterion is also exact: 2m square iff
m=2k². This equals the even-balance first term, but does not entail a recurring
leaf phenomenon: independently, 32 has earlier neighbors exactly {4,17}.
The ordered weak-fibre split remains a candidate construction, not the
historically recovered partition.

### Additional abc boundary strengthening

Every primitive positive additive triple is realized after normalization by
the difference-Plücker construction. Up to swapping A and B, take A<B,
d=B−A>0, and potentials (0,d,2d,2B). Strict ordering follows from A>0.
Their matching products are exactly (2dA,2dB,2d(A+B)); their gcd is 2d
because gcd(A,B)=1. For A=B, primitivity forces A=B=1; potentials
(0,2,3,6) give (6,6,12), normalized to (1,1,2). The checker verifies these
identities across all primitive pairs 1≤A≤B≤50. This supports the symbolic
proof, not a restriction to those pairs. Therefore this geometric family
already includes every primitive abc triple and cannot impose an additional
radical restriction merely by admitting the encoding. No abc consequence is
claimed.

## Reproduction and receipts

From this packet directory:

```
python3 audit_algebra/check_independent.py
python3 -O audit_algebra/check_independent.py
```

Both PASS with 101 explicit labeled checks. Receipts are
`receipt.normal.json` and `receipt.optimized.json`; each also records all check
labels. Finite controls support the stated universal proofs but do not replace
them. No blocker encountered for the algebra audit.
