"""Exact Johnson x K2 neutral fibers in the physical A/B macro language."""
from itertools import combinations
from collections import defaultdict,deque
from math import comb
import json
import constructor_language as c
MACROS=['TTTLL','TTLTL']

def diameter(adj,metric=None):
    maximum=0
    for start in adj:
        dist={start:0};queue=deque([start])
        while queue:
            v=queue.popleft()
            for u in adj[v]:
                if u not in dist:dist[u]=dist[v]+1;queue.append(u)
        c.need(len(dist)==len(adj),'fiber graph connectivity')
        if metric is not None:
            c.need(all(distance==metric(start,end) for end,distance in dist.items()),'every exact pairwise graph distance')
        maximum=max(maximum,max(dist.values()))
    return maximum

def peel_births(rows,source_order):
    ids=list(range(len(rows)));current=tuple(rows);removed=[]
    while len(current)>source_order:
        n=len(current);full=(1<<n)-1
        incoming=[full^row^(1<<i) for i,row in enumerate(current)]
        low=[i for i in range(n) if incoming[i].bit_count()==1]
        c.need(len(low)==2,'exactly two indegree-one vertices after each birth')
        latest=[i for i in low if incoming[i].bit_length()-1 in low]
        c.need(len(latest)==1,'canonical youngest vertex')
        v=latest[0];removed.append(ids[v])
        current=c.a.induced(current,[i for i in range(n) if i!=v]);ids.pop(v)
    return list(reversed(removed))

def main():
    receipt=[];seed=(2,4,1);root=0
    weak,_=c.construct((2,0),1,'L')
    incoming=[7^row^(1<<i) for i,row in enumerate(weak)]
    c.need(all(mask.bit_count()==1 for mask in incoming),'minimal weak-source peeling hostile control')
    for m in range(1,8):
        groups=defaultdict(list);data={}
        for bits in range(1<<m):
            word=''.join(MACROS[(bits>>i)&1] for i in range(m))
            graph,r=c.construct(seed,root,word);matrix=c.word_matrix(word)
            Hrow=tuple(matrix[0][i]+matrix[2][i] for i in range(3))
            recovered=sum(((graph[3+8*i+4]>>(3+8*i+5))&1)<<i for i in range(m))
            c.need(recovered==bits,'private arc recovers every macro choice')
            c.need(peel_births(graph,3)==list(range(3,len(graph))),'canonical peeling recovers chronological labels')
            reverse=list(reversed(range(len(graph))))
            relabeled=c.a.induced(graph,reverse)
            peeled=peel_births(relabeled,3)
            c.need(peeled==[len(graph)-1-i for i in range(3,len(graph))],'intrinsic peeling after nontrivial relabeling')
            b=(bits&((1<<(m-1))-1)).bit_count();a=m-1-b
            groups[a].append(bits);data[bits]=(Hrow,graph)
        for a,vertices in groups.items():
            b=m-1-a;adj={v:[] for v in vertices};near={v:[] for v in vertices}
            c.need(len(vertices)==2*comb(a+b,a),'labeled fiber vertex census')
            for u,v in combinations(vertices,2):
                c.need(data[u][0]==data[v][0],'full-H matrix identity inside fiber')
                delta=u^v;positions=[i for i in range(m) if delta>>i&1]
                edges=4*len(positions)
                if m<=5:
                    measured=sum((x^y).bit_count() for x,y in zip(data[u][1],data[v][1]))//2
                    c.need(measured==edges,'exact four changed arcs per private macro bit')
                if len(positions)<=2:
                    c.need(positions==[m-1] or (len(positions)==2 and positions[-1]<m-1),'all neutral at-most-two-bit moves classified')
                    adj[u].append(v);adj[v].append(u)
                    if positions==[m-1] or positions[1]==positions[0]+1:
                        near[u].append(v);near[v].append(u)
            c.need(all(len(adj[v])==a*b+1 for v in vertices),'Johnson x K2 degree')
            prefix_mask=(1<<(m-1))-1
            johnson_metric=lambda u,v: ((u^v)&prefix_mask).bit_count()//2 + ((u^v)>>(m-1)&1)
            def adjacent_metric(u,v):
                up=[i for i in range(m-1) if not(u>>i&1)]
                vp=[i for i in range(m-1) if not(v>>i&1)]
                return sum(abs(i-j) for i,j in zip(up,vp)) + ((u^v)>>(m-1)&1)
            c.need(diameter(adj,johnson_metric)==min(a,b)+1,'Johnson x K2 diameter')
            c.need(diameter(near,adjacent_metric)==a*b+1,'adjacent-swap refinement diameter')
            expected=[0]*(min(a,b)+2)
            for j in range(min(a,b)+1):
                value=comb(a,j)*comb(b,j);expected[j]+=value;expected[j+1]+=value
            observed=[sum(johnson_metric(vertices[0],v)==j for v in vertices) for j in range(len(expected))]
            c.need(observed==expected,'Johnson distance generating polynomial')
            adjacent_poly=[sum(adjacent_metric(vertices[0],v)==j for v in vertices) for j in range(a*b+2)]
            receipt.append(dict(macros=m,prefix_A=a,prefix_B=b,labeled_vertices=len(vertices),Johnson_edges=sum(map(len,adj.values()))//2,degree=a*b+1,Johnson_diameter=min(a,b)+1,adjacent_diameter=a*b+1,Johnson_distance_coefficients=expected,adjacent_extreme_distance_coefficients=adjacent_poly))
        # Distinct exponent strata do not accidentally share the full-count row.
        representatives=[data[v[0]][0] for v in groups.values()]
        c.need(len(set(representatives))==len(groups),'distinct a,b strata have distinct H rows')
    (c.HERE/'fibre_graph.json').write_text(json.dumps(dict(status='VERIFIED exact matrix and labeled-graph finite controls for the all-order combinatorial proof',max_macros=7,weak_peeling_hostile=dict(source=[2,0],root=1,word='L',output=weak,eligible_newborns=3),records=receipt),indent=2)+'\n',encoding='utf-8')
    print('PASS every macro word through depth7, private-bit recovery and canonical peeling, Johnson x K2 degree/diameter and adjacent refinement')
if __name__=='__main__':main()
