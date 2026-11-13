# A guarded core-changing descent for all seventy inherited charts

**PROVED within the stated infinite grammar. Chan–Pak and the473 repair obligations remain OPEN.**

## 1. What is inherited and what changes

The [birth packet](../../birth_grammar_repairs_20261112/PROOFS.md), B1–B2,
constructs24 even-order synchronized charts and46 odd-order charts, all
strong unit seeds. Their existence, birth status, and quotient separation
are inherited. Here the formerly forward-only family grammar is supplied
with exact scalar inversion, guarded core-changing deletion, a strict
Hamilton-count rank, and a finite terminal bank. None of those algorithms
is claimed to recognize arbitrary tournaments or arbitrary attainable counts.

Let E_j(m), j=0,...,23, replace the two coupled singleton ports of the
inherited10-vertex source by TT_m; its free stable port is fixed at length1.
Let A_i(m), i=0,...,45, replace the two specified members of the inherited
13-vertex source by TT_m; those members originally belong to two TT_2 blocks.
Thus |E_j(m)|=8+2m and |A_i(m)|=11+2m. Every m is a positive integer.
O is the inherited selected set and neither growing port is selected.
The unit defect is delta=H−W(O)=1. Write U for the selected-start slope,
and C_t(S) for the inherited physical tail of length t≥1. Then

    H(C_t(S))=1+2^t U(S).

Adjacency rows, selected labels, ports, chart identity, parameter m and t
are the state. We retain those actual data; we do **not** assert numerical
invariance of the full endpoint table. Indeed a full endpoint table includes
H, so keeping it numerically fixed would contradict strict H descent.

## 2. Exact physical core-change map

Every A_i has a distinguished unselected vertex v, the last vertex12 in
its13-vertex presentation. Its deletion leaves E_j(m+1), where
j=`parent_synchrony_family`. This is an actual marked graph identity up to
the order-preserving relabelling within the twin blocks, not a conclusion
from equal scalar slopes.

`compile_descent.py` checks for every i that deleting vertex12 at m=1
literally yields the even source with its two ports expanded to length2.
The two growing representatives map to distinct inherited ports under the
explicit expansion labels. All selected labels also agree. At general m,
each of these two blocks contains m new copies of one member and the other
old member: a chain of length m+1. Their external orientations and selected
status are unchanged. Hence the same graph identity holds at every m≥1.
This deletion changes the marked twin quotient in the inherited mixed-birth
construction. It is not an eligible stable-twin contraction of the odd birth.

The full independent endpoint/skeleton compiler gives finite exponential
profiles u_j(m)=U(E_j(m)), v_i(m)=U(A_i(m)) and

    Delta_i(m)=v_i(m)−u_j(m+1)>0.                    (1)

The positive assertion is certified for all m, not fitted to a finite grid.
Each profile is computed from all compressed Hamilton/two-cover skeletons
on the acyclic two-port set. For an exponential polynomial
P(m)=sum_b e_b b^m, the exact identity

    P(m)=sum_(r=1)^R [sum_(b>=r) e_b binom(b,r)] r! S(m,r), m≥1,

gives a finite positivity certificate. Every coefficient on the right is
nonnegative and the r=1 coefficient is strictly positive for every Delta_i.
The same test certifies strict growth u_j(m+1)>u_j(m), v_i(m+1)>v_i(m),
and u_j(m),v_i(m)≥2^m. All coefficients and source graphs are in profiles.json.gz.

For example i=0 has j=0 and

    u_0(m)=−100+200·2^m+212·4^m+200·6^m,
    v_0(m)=−260+1040·2^m+2176·4^m+3120·6^m,
    Delta_0(m)=−160+640·2^m+1328·4^m+1920·6^m.

The last profile has positive Stirling coefficients
17952,37408,43712,30128,11520,1920. Hence at m=1 the source slopes
are v_0(1)=29244 and u_0(2)=11292, a strict drop17952. The tail targets
at t=1 are58489 and22585. These initial values are included in literal checks.

## 3. Strict target-decreasing recursion

The allowed inverse rules are:

1. (odd i,m,t) -> (even j,m+1,t), using the actual vertex deletion above.
2. (even j,m,t) -> (even j,m−1,t), for m>1, simultaneously removing one
   vertex from each coupled transitive block.
3. (even j,1,t) -> (even j,1,t−1), for t>1, shortening the external tail.

Every state remains a genuine unit seed followed by a positive-length tail.
Rule2 is a joint synchronized contraction, not an inference from individual
stable-port contraction; the unit identity on the diagonal is proved on the
whole chart. Rule3 is justified only while t−1≥1. Each rule strictly lowers
the actual Hamilton count N. Their terminal states are the24 explicit
(even j,1,1) presentations, which may have repeated scalar values.
The24 terminal chart presentations share12 actual marked source graphs,
each appearing twice with different coupled-port choices. In chart order
the pairs are(0,1),(2,3),...,(22,23); their adjacency rows and selected masks
are literally identical at m=1. Thus24→12 is a precise forgetful-map fiber,
not a consequence of a shared number in another subject or of converse
identification. The inherited birth proof separates each pair at m≥2.
These12 terminal graph targets have12 distinct counts:
2009,2297,3475,4025,4697,6241,7087,9375,11041,12741,16599,21525.
The rule counts are m+t−2 from an even state and m+t from an odd state.
Each rule is also an induced deletion of the complete tail tournament:
source deletion commutes with tail attachment because all tail arcs depend
only on the retained selected marking. Reversing the recorded rule sequence
is a literal graph construction, and
its selected labels/endpoints are always available from the retained rows.

This proves well-foundedness and physical validity **inside the70-chart
language**. It does not prove that each allowed odd N admits one of these
states. The compiled language deliberately holds each free port at length1;
adding that further clock is a different language and is not silently counted.

## 4. Exact arbitrary-target recognition within that language

For an odd N>1 enumerate t=1,...,v2(N−1), and set U=(N−1)/2^t.
In each of the70 charts, strict profile growth permits exact binary search
for m in1,...,floor(log2 U). The bound follows from U(m)≥2^m.
Equality supplies an actual graph and the strict decreasing trace above;
no equality means only that N lies outside this declared language.
All calculations use integers. No approximate logarithm, recurrence fitting,
or prime heuristic supplies a success condition.

This recognizes every scalar value of the declared infinite language,
including all multiple representations, and gives witness order O(log N).
There are at most70 floor(log2 X)^2 parameter states below X, hence at most
that many distinct values. The target language is sparse despite having
complete inverse rules. This is a positive, fully costed model for the
needed wider recursion, not an eventual-tail theorem.

In particular, a source with U=c2^j pays s=t+j≥j+1. The unpaid heads s≤j
have no positive tail-length state in this grammar. Changing to normalized
odd c while omitting j would erase exactly the remaining endpoint debt.
None of the473 infinite repair obligations is declared closed by this lane.

## 5. Independent checks and a general monotonicity lemma

For completeness, U admits an injection under induced marked extensions:
if T is induced in T' and O remains selected (extra selected vertices allowed),
then U(T,O)≤U(T',O'). It suffices to add one vertex x. An object counted by U
is an ordered cover (P,Q), with Q nonempty starting in O and P possibly empty.
When P is nonempty, insert x at its first valid Hamilton insertion position.
When P is empty, insert x into Q at its first valid non-prefix position if
one exists; otherwise take P=(x) and leave Q unchanged. Removing x identifies
the old object and its case. Selected starts are retained. Thus this is an
injection, not merely a comparison of equal-dimensional numerical observers.
Equality can occur for neutral additions, so the strict coefficients in(1)
are indispensable. The general injection does not prove unit-defect retention;
that guard comes from the inherited/computed full signed chart kernels.

`literal.cpp` independently counts actual subset endpoints, ordered two-path
covers and Hamilton paths. `verify.py` tests280 source/tail graphs at m=1,2,
including every70 chart, plus210 exact large parameter targets. Source and
optimized Python runs use explicit gates and have identical deterministic
compressed receipts. The finite graph checks are controls; all-length claims
come from complete kernels and finite exact coefficient identities.

Reproduction from repository root:

    python3 -B 05-knowledge/results/session290_core_descent_20261113/descent/compile_descent.py
    python3 -B -O 05-knowledge/results/session290_core_descent_20261113/descent/compile_descent.py
    python3 -B 05-knowledge/results/session290_core_descent_20261113/descent/verify.py
    python3 -B -O 05-knowledge/results/session290_core_descent_20261113/descent/verify.py
    python3 -B 05-knowledge/results/session290_core_descent_20261113/descent/decode.py 58489 --rows

The compiler imports the inherited finite skeleton algorithm directly; the
literal C++ verifier is an independent counting path. No claim of a new
complete census or novel inheritance of the original70 charts is made.
