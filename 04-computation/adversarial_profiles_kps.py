"""
ADVERSARIAL check of cluster-profile claims (D) and (E):
 - H=7 has EXACTLY ONE graph-realizable independent-set profile = [1,3,0] = K3.
 - H=21 has FOUR graph-realizable profiles [1,4,3],[1,6,2],[1,8,1],[1,10].

A 'graph-realizable profile' here = an alpha-vector (alpha_0=1, alpha_1, alpha_2, ...)
that is the independent-set-count sequence of SOME simple graph G.
We enumerate ALL simple graphs on alpha_1 vertices (up to iso, but for safety we
enumerate ALL labeled graphs and dedup by alpha-vector), compute their
independent-set count sequences, evaluate I(G,2)=sum alpha_k 2^k, and bucket.

We bound alpha_1: since I(G,2)>=1+2*alpha_1 (the k=0,1 terms), H=7 => alpha_1<=3,
H=21 => alpha_1<=10. We enumerate graphs on up to 10 vertices? That's 2^45 — too
many. Instead: I(G,2)=sum alpha_k 2^k, and the EMPTY graph on alpha_1 vertices gives
3^alpha_1 which is the MAX. So for a target T we need 3^alpha_1 >= T, i.e.
alpha_1 >= log_3(T). And 1+2*alpha_1 <= T gives alpha_1 <= (T-1)/2.
For H=7: alpha_1 in [2,3]. For H=21: alpha_1 in [3,10].
We enumerate labeled graphs on v vertices for v up to 10 — but 2^45 infeasible at v=10.

KEY REDUCTION: the independent-set sequence depends only on the graph up to iso.
For a target value T with I(G,2)=T and v vertices, note alpha_1=v (all single
vertices are independent sets of size 1). The number of edges only REMOVES
independent sets. We can enumerate by number of vertices v and brute-force over
labeled graphs for v<=7 (2^21 ~ 2M, feasible), and handle v=8,9,10 by the fact
that with v vertices the minimum value of I(G,2) is achieved by complete graph K_v
giving 1+2v (since alpha_k=0 for k>=2). So I(G,2) ranges in [1+2v, 3^v].
For H=21: v can be 3..10 but we must check each.
 v=10: only K_10 gives 1+20=21 (the [1,10] profile = EMPTY graph? wait).

CAREFUL: [1,10] profile means alpha_0=1, alpha_1=10, alpha_2=0,...
 alpha_2=0 means NO independent pair => every pair is an edge => COMPLETE graph K_10.
 I(K_10,2)=1+20=21. YES.
So [1,10] = K_10 (complete). Good.

So we enumerate v from 3..10. For each v we want all distinct alpha-vectors with
I(G,2)=21 (and =7). For v<=7 brute force labeled graphs. For v=8,9,10:
 v=10: max independent set size limited; to get I=21 with v=10 we need alpha_2..=0
   beyond — 1+20+4*alpha_2+...=21 => alpha_2=0 and higher 0 => K_10. one profile.
 v=9: 1+18=19, need +2 from higher => 4*alpha_2+...=2 impossible (multiple of 4
   plus 8.. ) => alpha_2 must satisfy 4 alpha_2 + 8 alpha_3+...=2 -> impossible. So
   NO v=9 graph gives 21? Actually 1+18=19, remaining 2 not reachable. So none.
 v=8: 1+16=17, remaining 4 => 4*alpha_2=4 => alpha_2=1, rest 0 => profile [1,8,1].
   Realizable? need a graph on 8 vertices with exactly ONE independent pair and no
   independent triple. That's K_8 minus one edge. Check.
 v=7: 1+14=15, remaining 6 => 4 a2 + 8 a3 +... =6 -> no integer solution (6 not
   reachable by 4,8,16..). So NO v=7 profile for 21? Wait 4*a2=6 impossible. Hmm.
   But worker lists [1,6,2] (v=6) and [1,4,3] (v=4). Let's just brute force v<=7
   and special-case v=8,9,10.
"""
import itertools

def indep_profile_from_adj(v, edges):
    """edges = set of frozenset({i,j}); compute alpha-vector via enumerating
    independent sets (v<=10 ok via subset enumeration over 2^v)."""
    adj=[set() for _ in range(v)]
    for e in edges:
        a,b=tuple(e)
        adj[a].add(b); adj[b].add(a)
    alpha=[0]*(v+1)
    for mask in range(1<<v):
        verts=[i for i in range(v) if (mask>>i)&1]
        # independent?
        ind=True
        for i in range(len(verts)):
            for j in range(i+1,len(verts)):
                if verts[j] in adj[verts[i]]:
                    ind=False;break
            if not ind:break
        if ind:
            alpha[len(verts)]+=1
    return tuple(alpha)

def Iz(alpha,z=2):
    return sum(a*(z**k) for k,a in enumerate(alpha))

def enumerate_target(T, vmax_brute=7):
    """Find all distinct alpha-profiles (as tuples, trailing zeros stripped) that are
    graph-realizable with I(G,2)=T. Brute force labeled graphs on v<=vmax_brute;
    for larger v use complete-graph-minus-edges reasoning is NOT general, so we
    also brute force v=8 by enumerating graphs that are dense (complement sparse)."""
    found={}  # stripped_profile -> example (v, edges)
    # brute force v up to vmax_brute over ALL labeled graphs
    for v in range(2, vmax_brute+1):
        pairs=[frozenset((i,j)) for i in range(v) for j in range(i+1,v)]
        m=len(pairs)
        if m>21:
            continue
        for bits in range(1<<m):
            edges=set(p for idx,p in enumerate(pairs) if (bits>>idx)&1)
            alpha=indep_profile_from_adj(v,edges)
            if Iz(alpha)==T:
                # strip trailing zeros
                a=list(alpha)
                while len(a)>1 and a[-1]==0: a.pop()
                key=tuple(a)
                if key not in found:
                    found[key]=(v, sorted(tuple(sorted(e)) for e in edges))
    # handle v=8,9,10 via COMPLEMENT enumeration (complete minus a sparse graph)
    # We enumerate the COMPLEMENT graph (non-edges) since for I close to 1+2v the
    # graph is near-complete (few independent sets). Enumerate complement with up to
    # C(v,2) edges but realistically the complement is sparse for large v.
    for v in range(vmax_brute+1, 11):
        pairs=[frozenset((i,j)) for i in range(v) for j in range(i+1,v)]
        m=len(pairs)
        # complement (the non-edges of G) = the edges present in Gbar. G near complete
        # means few non-edges. We enumerate non-edge sets of size up to some bound.
        # I(G,2) - we just need profiles equal to T. Since 1+2v could exceed T we skip.
        if 1+2*v>T:
            continue
        # enumerate complement edge sets of size 0..min(m, some bound)
        # The independent sets of G of size>=2 correspond to cliques in complement.
        # alpha_2(G) = #edges of complement. Since I(G,2)=1+2v+4*alpha_2+8*alpha_3+...
        # = T, we need 4*alpha_2 <= T-1-2v, so #complement-edges <= (T-1-2v)/4.
        bound=(T-1-2*v)//4
        if bound<0:
            continue
        maxne=min(m, bound)
        for ne in range(0, maxne+1):
            for nonedges in itertools.combinations(pairs, ne):
                # G edges = all pairs except nonedges
                ne_set=set(nonedges)
                edges=set(p for p in pairs if p not in ne_set)
                alpha=indep_profile_from_adj(v,edges)
                if Iz(alpha)==T:
                    a=list(alpha)
                    while len(a)>1 and a[-1]==0: a.pop()
                    key=tuple(a)
                    if key not in found:
                        found[key]=(v, "complement_size_%d"%ne)
    return found

if __name__=="__main__":
    for T in (7,21):
        print(f"=== Graph-realizable independent-set profiles with I(G,2)={T} ===")
        res=enumerate_target(T, vmax_brute=7)
        for prof in sorted(res):
            print(f"  profile {prof}  (example v={res[prof][0]})")
        print(f"  TOTAL distinct profiles: {len(res)}")
        print()
