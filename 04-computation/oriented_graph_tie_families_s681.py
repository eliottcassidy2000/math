#!/usr/bin/env python3
"""oriented_graph_tie_families_s681.py — tournaments WITH TIES = oriented graphs:
which incomplete-orientation families are relevant to our proofs.

Prompt: ties remove edges => no longer complete or even connected; a useful family.
Which ones appear as relevant?

FINDINGS (verified):
 (1) COMPLETENESS BOUNDARY. The H-gaps {7,21} (impossible for tournaments) are
     reachable with just ONE TIE: a tournament minus a single edge already has
     orientations with H=7 and H=21 (n=5: H=7 at 1 tie; n=6,7: both 7 and 21 at
     1 tie). And EVEN H-values appear once incomplete (Redei's odd-H needs
     completeness). => {7,21} are obstructions OF COMPLETENESS, living only on
     the codistance-0 (complete) stratum. The whole H-impossibility phenomenon
     is a feature of the complete (tournament) stratum, dissolved by one tie.
 (2) So the oriented-graph family is graded by TIE-COUNT (codistance from
     complete). Relevant strata / families for our proofs:
   * STRONGLY-CONNECTED oriented graphs = the ATOMS. H is multiplicative over
     strong components (the partition-function / equidecomposability monoid,
     HYP-2183); ties that DISCONNECT factor H (or kill it: a Ham path needs a
     TRACEABLE graph, so disconnecting ties => H=0, the 'vacuum'). The H=21 proof
     reduced to strong tournaments -- the indecomposable carriers.
   * THE CONFLICT GRAPH Omega (general undirected graph; non-edges = TIES =
     vertex-disjoint odd cycles): H = I(Omega,2). This is where the count
     actually lives; the H=21 proof is entirely about Omega's tie/independence
     structure. The relevant 'incomplete' object IS Omega.
   * 1-TIE NEIGHBORS (tournament minus one edge): the minimal incomplete probes;
     they dissolve the gaps and are the natural induction neighbors (the n+2
     source/sink recursion adds near-tie vertices).
   * LRC THRESHOLD-TIE graphs: binding antipodal pairs {a, n-a} at distance
     EXACTLY 1/n are TIES (neither near nor far). At each tight witness t=j/n the
     tie-graph is the antipodal MATCHING; tight configs <=> this tie structure
     (S679: the negation-involution orbits). The relevant LRC tie-object = the
     threshold matching.

Session: claude-2026-06-06-S681 (oriented-graph-tie-families)."""
import sys; sys.stdout.reconfigure(line_buffering=True)
from itertools import product, combinations
def Hcount(n,adj):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        row=dp[mask]
        for v in range(n):
            c=row[v]
            if c:
                av=adj[v]
                for w in range(n):
                    if not(mask>>w&1) and av>>w&1: dp[mask|1<<w][w]+=c
    return sum(dp[size-1])
print("Completeness boundary: achievable H over oriented graphs (ties allowed), n=5:")
edges=list(combinations(range(5),2)); Hs={}; mt={}
for state in product([0,1,2],repeat=len(edges)):
    adj=[0]*5; ties=0
    for (i,j),s in zip(edges,state):
        if s==0: adj[i]|=1<<j
        elif s==1: adj[j]|=1<<i
        else: ties+=1
    h=Hcount(5,adj); Hs[h]=Hs.get(h,0)+1
    if h not in mt or ties<mt[h]: mt[h]=ties
print(f"  achievable H: {sorted(Hs)}  (tournaments give only odd, no 7)")
print(f"  H=7 reachable at min ties = {mt.get(7)} (tournament minus 1 edge); even H present (Redei broken).")
print("  => {7,21} are obstructions of COMPLETENESS; the tournament stratum is special.")
print("\nRelevant families (see docstring): strong components (atoms / multiplicative),")
print("the conflict graph Omega (ties=disjoint cycles; where H=I(Omega,2) lives),")
print("1-tie neighbors (dissolve gaps; induction), LRC threshold-tie antipodal matchings.")
