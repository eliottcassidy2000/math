# Independent audit: weak-allowed native preparation reduction

**Analytic and inherited-operation audit: PASS after two explicit candidate
repairs.** This accepts [PREPARATION_REDUCTION.md](PREPARATION_REDUCTION.md)
as a necessary-and-sufficient reduction under its stated operation contract.
It does not establish membership of an infinite ternary orbit in either
ready spectrum. Status promotion still requires the separately audited
all-tournament joint classification used as its load-bearing input.

## 1. Repairs and exact dependencies

The primary originally attributed the entire weak prefix decomposition to
JOINT_CLASSIFICATION, whose corresponding corollary only gave the strong
SCC zero table. The producer authorized a repair: I inserted the direct
terminal-SCC argument into the necessity proof. It is now self-contained
after the positive joint iff and does not depend on an unwritten corollary.

I also checked both native-compiler section headings. The METAGRAPH link
fragment was incorrect; I changed it to
`#5-one-temporary-vertex-compiles-any-old-arc-reversal`. The NATIVE_SIX
section4 link is correct as written. No producer computation or frozen
joint-classification file was changed.

I read the full inherited METAGRAPH §5 theorem and NATIVE_SIX §4 guard
construction. The former is explicitly for **any named tournament**, with
no strongness or Hamilton-count hypothesis. It reverses one prescribed
unordered pair, keeps every other old pair and every old name, and uses one
temporary helper, two legal complements, then helper-only deletion. Traced
append and prepend cases, inverse gadgets for a guard legal in the reversed
orientation, and the order-two boundary are paid. NATIVE_SIX supplies the
direct pair-block construction of such a guard. The general compiler does
not promise strong intermediate states; the reduction does not assume it.

## 2. Necessity and all weak-source cases

The all-tournament joint guard gives `N^+(a)={r}` and `e_r(V−a)=s` with
s=1 or2. Positive endpoint count puts r in the terminal SCC C of U=V−a.
Hamilton paths pass through its SCCs in their unique order, hence

    e_r(U)=product_(earlier SCC K) H(K) * e_r(C).

Every nontrivial strong SCC has at least three Hamilton paths. Therefore
the earlier SCCs are singletons, with total product one. They form TT_l
and all beat a, since a's only exit is r. Inside C+a, the root a reaches
C through r; if C is nontrivial, C reaches a through any vertex other
than r. Thus C+a is strong. If C is a singleton, the whole V is transitive,
with a penultimate and r last; H=1. The hypothesis h>1 excludes exactly
this boundary. No parity theorem is needed for this argument.

This proves V=TT_l→Q with a genuine strong ready Q and H(Q)=h. Removing
the prefix is only a proof of spectral membership. It is not an allowed
deletion of old physical names. In particular the necessity proof does
not quietly authorize discarding a port outside Q.

## 3. Sufficiency, physical names, and exact cost

For an actual witness Q of order m and an actual input T of order n, take
N=max(n,m) and W=TT_(N−m)→Q. There are N positions and n old names. One
can inject all old names into those positions while fixing the designated
root at Q's root; the remaining positions receive distinct fresh names.
This preserves identities of every additional marked port, but does not
promise to preserve its old response. That is precisely the contract.

This naming assignment specifies a target matrix on the actual names. It
does not act as a relabelling operation on T. Birth N−n permanent vertices
with declared cuts, then compile each differing pair of the current and
target matrices. The inherited one-pair theorem applies after every prior
flip, because it has no special-class hypotheses. If d pairs differ, the
cost is exactly N−n permanent births, d helper births, 2d guarded native
complements, and d helper-only deletions. The preparation peak is at most
N+1, or N if d=0. Subsequent C2 adds two permanent vertices; its completed
order is N+2. The primary correctly labels its displayed smaller peak as
the preparation budget, not the whole branch budget.

The transitive prefix leaves H and sigma unchanged. For the retained-root
two-bin cover F, each prefix vertex can independently belong to either
ordered bin, and the vertices within each bin occur in their transitive
order before its core part. The first bin has a nonempty core part because
it ends at r. This gives the exact bijection

    F_r(W)=2^(N−m) F_r(Q),

including an empty core part in the second bin. The sole-exit condition
then pays tau=0. The C2 identity gives the primary's exact child H, F and
changing root defect. There is no unproved claim that H is constant at
intermediate native states, or that a marginal observer determines the
next guards: the full named relation is retained.

## 4. Both directions of the unbounded criterion

If an integer set contains the actual seed, lies in both ready spectra,
and is closed under h→3h+2 and h→3h+4, apply the proved preparation
equivalence at each actual output. This is valid even as order and root
defect change: the equivalence quantifies over every actual input with the
given h and a retained root. Padding pays old-name retention at the next
stage. Repeated selection of witnesses and compilation therefore produces
every finite requested word, coherently under a fixed deterministic
witness-ordering convention if desired.

Conversely, the existence of both contract-compliant preparations at every
reachable word endpoint puts the full scalar orbit in both spectra by
necessity. The orbit is closed under the two scalar maps. This proves the
claimed equivalence, not merely one sufficient direction.

Weak targets are essential to this exact reduction: if a transitive prefix
is needed, both it and the subsequent C2 child remain weak. An all-strong
version cannot use padding to absorb arbitrary retained names and is not
proved here. The finite induced obstruction and the example surgeries in
the primary are contextual controls, not hidden premises of this logical
equivalence.

## 5. Finite complete witness selection

The primary-source Busch bound was independently checked during the
[order-budget audit](../frame/ORDER_AUDIT.md): a strong n-vertex witness
with H=h satisfies `5^(n−1)<=h³`. Thus only finitely many n>=3 can occur.
Enumerating all labelled tournaments at those orders, all actual marks,
testing strongness and adjacency-derived ending counts, is finite and
complete for membership in R_s. The inequality uses integers only and
includes the possible equality boundary. It does not assume the 34-core
selector is complete, or promote a bounded core bank to the whole spectrum.

The exact rooted unique/two-start grammar in TWO_ENTRY is an alternative
enumerator, in converse and followed by the sole-exit birth. The
two-exit scalar-ready exception is not incorrectly included as a joint
guard. No efficient general spectrum algorithm is claimed: the complete
labelled witness enumeration may be extremely large.

## 6. Final audited scope

The repaired primary is 7,803 bytes, SHA256
`1fda22d20eb9260805a0468f452d8a9114866480d5929df1735fc3f7ec6fe14b`.
A later status/link-only promotion does not change this proof audit.
No computation was rerun because the new theorem is an analytic reduction
from the audited joint guard and inherited all-order native compiler.

The reduction identifies an exact remaining integer-membership question.
It does not assert that any specific seed's entire orbit stays in both
spectra, that every odd integer is covered, that preparation can preserve
F, or that all intermediate/output graphs can remain strong. Those
distinctions are explicit in the primary and are part of this acceptance.
