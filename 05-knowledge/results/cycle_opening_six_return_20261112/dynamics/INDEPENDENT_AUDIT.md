# Independent audit of the six-state arithmetic and strong contextual return

**VERIFIED / PROOF AUDITED.** The full-count lane independently reviewed
[PROOF.md](PROOF.md), including the positive integral classification,
autonomous marked-cut controller, strong wrapper, contextual observer,
near-square boundary, and period statements. The elementary new theorems
pass with the domain guards recorded below. The review does not promote
Chan–Pak, Landau, a classification of all positive Collatz cycles, or an
indefinitely iterable ready-return grammar.

## Integral states and exact periods

For L(a,b)=(b,b/a), direct iteration gives L³(a,b)=(1/a,1/b) and L⁶=id.
The only rational fixed point is (1,1). L² has no further rational fixed
point because it would require a³=1. L³ fixes exactly the four sign pairs;
the three nonidentity pairs form one exact three-cycle. All remaining
rational nonzero states therefore have exact period six. These are statements
about the displayed two-dimensional map, not quadratic polynomials.

The integral classification correctly concerns **positive** h_i and C.
The relation h_i h_(i+3)=C² first forces rational C to be an integer.
For p^c dividing C, the six valuation inequalities are precisely

    c±x≥0, c±y≥0, c±(x−y)≥0.

They are equivalent to the stated radius-c hexagon. Primes outside C have
x=y=0, so independent local choices give every positive integral state
exactly once. The shell count 6j yields 1+3c(c+1), and the Chinese-independent
prime choices multiply. Apart from the constant positive state, all cycles
have length six; division by six is therefore justified. At p=2, requiring
all six values odd forces c=x=y=0. This verifies both directions of the
oddness statement, including C=1 and its single constant state.

The positivity guard was made explicit during review. If signs were allowed,
C=1 with normalized state (-1,-1) gives an integral three-cycle, so the
positive classification must not be advertised as a classification of all
signed integral states.

## Autonomous marked cuts and full count

After centering at (m,m), the state domain is
`|x|,|y|,|x+y|≤m` and the controller is Q(x,y)=(-y,x+y).
Its third power is -I; neither Q−I nor Q²−I has a kernel. Thus only the
center has a shorter orbit, and the lattice-point and orbit counts follow.
The marked cut and named resource copies are essential: they retain the
state that an unmarked count would forget. The displayed sorting rule is
an autonomous operation on those actual objects.

The ordinal-product Hamilton-path bijection proves the same full count
at every phase, independently of the prefix readout. The prefix recurrence
follows from Q+Q⁻¹=I. Its asserted exact period remains valid when a=b,
when a and b are multiplicatively dependent, or when exactly one is1:
the nonzero functional (log a,log b) cannot vanish on the orbit of a
nonzero vector under the real irreducible rotation Q. Only a=b=1 makes
all scalar prefixes constant. Repeated individual values in a six-cycle
do not contradict its least period.

A leading anchor retains its root indegree and multiplies full, start,
and relevant deletion counts by one fixed resource factor. This proves
phase neutrality but does not make a weak ordinal carrier strong. The
prime-target obstruction is correctly confined to this raw product:
any nonunit a or b supplies a nontrivial square divisor.

## Complete contextual observer and strong wrapper

The ordered-bin identity P_(A→B)=P_A P_B follows by splitting each bin at
the ordinal cut, including empty pieces. Hence all resource permutations
have exactly the same polynomial, not merely the same H.

The expansion `P_M(z)=sum_j K_j(M) binom(z,j)` uses ordered nonempty path
covers K_j. A uniform module context factors through these coefficients
because each global path has an ordered list of runs through M. Endpoints
or deletions outside M restrict the exterior skeleton without changing
that argument. The external-root and uniform-interface qualifications are
necessary and are explicitly retained.

The separating contexts C3[M,TT_t,TT_t] have at most t+1 module runs.
At the maximum, every exterior vertex is used as a singleton in the
alternating pattern M A B M ... A B M. Independent assignments of the A
and B vertices give coefficient (t!)². Together with the t=0 observation,
this is a triangular recovery of every K_j, validating the inherited iff
observer classification. It is not a classification of individually
addressable module vertices or roots inside M.

For W(M)=C3[1,1,M], three one-run patterns contribute3H(M), and the unique
two-run skeleton contributes K2(M)=P_M(2)−2H(M). Therefore

    H(W)=H(M)+P_M(2).

The physical root Y has the same sole predecessor X. Its Hamilton paths
must begin Y, traverse all of M, and end X; deleting Y also leaves M→X.
Both readouts equal H(M), proving defect zero. Strongness follows from
u→X→Y→v for module vertices u,v, even if M is weak.
This is an actual indefinitely repeatable strong return at a fixed full
target. It is contextually neutral, so it does not repair a deficient
external-root response. The equal-H9 modules with P2=144 and168 are an
appropriate hostile control for replacing the polynomial by H alone.

## Near-square and Fermat scope

For nonempty A,B, a,b are odd and their two-bin counts alpha,beta are even:
bin swapping has no fixed partition of a nonempty vertex set. Thus the
strong target is an odd square plus an even square. The common-divisor
obstruction is immediate. If m has an odd divisor d>1, the sum factors
as U^d+V^d; its factors are nontrivial because alpha beta≥4. The resulting
conditions on m and coprimality are necessary, never stated sufficient.

For the structural first-square-equals-one branch, ab=1 forces both
modules to be transitive. This uses the repaired all-order odd-cycle
collection theorem, as routed by the current correction banner of
[THM-209](../../../../01-canon/theorems/THM-209-independence-polynomial.md).
A transitive t-vertex module has P(2)=2^t by assigning each vertex to one
of two bins. Hence the displayed branch is exactly
`1+2^[2m(|A|+|B|)]`. The odd-factor argument forces its exponent to be a
power of two if the value is prime. Accidental alternative representations
as1+n² are explicitly left open. Neither general Landau coverage nor new
Fermat primes follow. The 5/17 example's factor29 is a decisive positive-
construction/negative-primality control.

## Signed storage, clocks, and literature boundary

The signed matrix claim uses A³=-I, which immediately excludes a whole
nonzero orbit from a pointed positive cone and excludes a constant positive
linear readout. For a nonnegative integral N with N⁶=I, both N and N⁻¹ are
nonnegative, forcing a monomial matrix; integrality of the inverse forces
a permutation. An A quotient requires a primitive sixth root among its
eigenvalues, hence a six-cycle and dimension at least six. The finite-order,
linear, integral hypotheses are necessary; the statement does not bound
nonlinear or affine encodings or prove compiler minimality.

The raw/Terras/odd-only Collatz lengths are mutually consistent. The
17-cycle has seven odd steps and valuation sum11, giving raw length18
and Terras length11. A phase-equivariant quotient of cycles requires and
suffices that the target length divide the source length. These are clock
maps, not value-preserving identifications or a classification of all
positive cycles.

The circulant assertion now explicitly assumes q≥3. At q=1 its zero parity
word has fixed point0, so the former unqualified positive-numerator phrase
had a degenerate endpoint. For q≥3 the numerator is positive and
2^(2k+1)−3^k>0, so reflection produces a negative3x−1 point. The marked
word and base-point guard remain indispensable; q7's two H values show
that the arithmetic shape alone is not a full-count state.

For f(x)=x²−29/16, odd-prime negative valuations double, while at2 only
valuation-2 can remain bounded. The real escape bound then leaves exactly
the eight odd-quarter candidates listed in the proof, all preperiodic.
This verifies the specific parameter's complete affine rational set and
its absence of exact six-cycles. The normalized interpolation matrix has
linear order six but projective order three, so it creates no extra scalar
quadratic orbit.

The primary [Flynn–Poonen–Schaefer paper](https://people.maths.ox.ac.uk/flynn/arts/art11.pdf)
was rechecked for its exact-period-five theorem and its account of Morton's
period-four exclusion. [Stoll's introduction](https://www.mathe2.uni-bayreuth.de/stoll/papers/Xdyn06-revised.pdf)
was rechecked for the analytic continuation, functional equation, and
specific genus-four Jacobian BSD rank assumptions. The packet accurately
preserves that conditional boundary. This is an audit of the cited theorem
statements and assumptions, not a new proof of their arithmetic geometry.

The positive five-period exchange identity also passes. Its only positive
fixed value is irrational, and a direct odd-count input produces an even
numerator over an odd denominator, preventing its immediate interpretation
as successive full tournament counts. Deleting the additive one changes
the recurrence, rather than conjugating its five-cycle to the six-cycle.

## Verification scope and outcome

The existing exact receipts were already replayed by the coordinating
lane; this review did not repeat the expensive20-/16-vertex counts.
It inspected the literal strong-wrapper implementation, independently
rederived its full and rooted formulas, and supplied the cheap signed-state
and q=1 endpoint hostiles above. No unresolved mathematical defect remains
in the reviewed stated scope.

The final stopping boundary now distinguishes the separately proved
[one-successor ready return](../compiler/MODULE_RETURN.md) from the still
open indefinitely iterable construction. The context-neutral six-phase
controller is valuable but cannot by itself change the rooted readout
needed by that larger problem.
