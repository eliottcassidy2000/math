# A quantitative run-budget boundary

**PROVED deduction from the inherited fixed-apex deletion kernel.**
This is a bound on a specified presentation class, not an exclusion from the
unrestricted tournament spectrum.

Fix a labelled apex tournament A on k>=1 vertices and an alphabet of q>=2
allowed chain-vertex neighbourhoods against A. Let C_(A,q,r)(N) be the number
of DISTINCT Hamilton-path counts h<=N supplied by tape words having at most
r nonempty constant runs. Use the inherited Session221
`automata/PROOF.md`, Theorem A:

    L(N)=max {2a+3b: a,b>=0, 3^a5^b<=N} <= floor(3 log_5 N).

For every finished word of count at most N, deleting its singleton chain
strong components preserves A and H and leaves a subsequence of length at
most L(N). Deleting letters cannot increase the number of constant runs;
it also cannot introduce a new letter. Thus a short representative remains
inside the same alphabet and run budget. This is a statement about finished
words, NOT a permissible state identification during continuation.

Consequently

    C_(A,q,r)(N) <= 1 + sum_(j=1)^min(r,L(N))
                            q(q-1)^(j-1) binom(L(N),j).       (1)

Indeed a nonempty word with j runs first chooses its j successive, distinct
run letters, in q(q-1)^(j-1) ways. Its positive run lengths with total at
most L are counted by binom(L,j), by summing binom(ell-1,j-1) over j<=ell<=L.
The empty word contributes at most one scalar. Different words may have the
same H, so this is only an upper bound, which is the required direction.

In particular:

- For fixed k,q,r, this class supplies O((log N)^r) values through N.
- More generally, if r(N)=o(log N), then it supplies N^{o(1)} values through N.

For the latter, replace L(N) by L0=floor(3 log_5 N). If r stays bounded the
first statement suffices. Otherwise r=o(L0), so eventually 1<=r<L0/2 and

    sum_(j<=r) q(q-1)^(j-1) binom(L0,j)
       <= q r (q-1)^r (e L0/r)^r.

Its logarithm is at most O(log r)+r log(q-1)+r log(e L0/r)=o(L0)=o(log N),
because x log(1/x)->0 as x->0. Cases where r is zero are immediate. This
proves the claim without assuming that the family has any positive density.

Thus fixed apexes and merely longer constant runs cannot furnish a
positive-density covering under a sublogarithmic run budget. The unrestricted
fixed-apex language may use Theta(log N) or more runs and is NOT excluded.
For the two-initial-star carrier q=3, (1) applies verbatim. The spectral
calculus makes each constant run simple; coverage must still contend with
an increasing number of ordered run boundaries.

## Connection and provenance

The source is the inherited count-preserving SCC deletion kernel, not the
letter spectrum alone. The target is an upper bound on the number of attained
scalars with limited run complexity. The map is induced deletion of singleton
chain SCCs; it preserves the physical core, H, alphabet and run budget, but
not arbitrary continuation states. The restored sidecar is the reduced word.
No new Busch bound, fixed-apex inverse, or general scalar gap theorem is claimed.
The deduction combines that prior kernel with elementary word enumeration.

For a finite set of fixed apex cores and alphabets, sum their bounds; the same
asymptotic conclusions hold. Allowing the number of apexes or alphabet size to
grow with N requires retaining those factors and is outside this statement.
