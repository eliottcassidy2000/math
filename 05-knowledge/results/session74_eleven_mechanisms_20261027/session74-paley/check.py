#!/usr/bin/env python3
"""Independent exact certificate checker; no removable assertions or dependencies."""
from itertools import combinations, permutations
from pathlib import Path
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def score(rows): return [r.bit_count() for r in rows]
def L(rows):
    s=score(rows); a=sum(abs(v-i) for i,v in enumerate(sorted(s)))
    b=sum(abs(sum(v<=k for v in s)-(k+1)) for k in range(len(s)-1))
    need(a==b and a%2==0,'independent score transports disagree')
    return a//2

def tour(n,edges):
    # Ascending baseline scores 0,...,n-1: higher index beats lower.
    rows=[(1<<i)-1 for i in range(n)]
    for a,b in edges:
        need(0<=a<b<n,'bad edge');rows[a]^=1<<b;rows[b]^=1<<a
    return rows

def components(n,edges):
    adj=[set() for _ in range(n)]
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    rem=set(range(n));out=[]
    while rem:
        reached={min(rem)};todo=list(reached)
        while todo:
            v=todo.pop()
            for w in adj[v]-reached:reached.add(w);todo.append(w)
        out.append(sorted(reached));rem-=reached
    return out

def frame(rows,order):
    n=len(rows);need(sorted(order)==list(range(n)),'not permutation')
    ed=[(order[i],order[j]) for i in range(n) for j in range(i+1,n) if (rows[order[j]]>>order[i])&1]
    cs=components(n,ed)
    return {'edges':ed,'components':cs,'m':len(ed),'c':len(cs),'beta':len(ed)-n+len(cs)}

def legal(rows,path):return sorted(path)==list(range(len(rows))) and all(rows[a]>>b&1 for a,b in zip(path,path[1:]))

def tau(rows,path):
    need(legal(rows,path),'illegal named path')
    keep={frozenset((a,b)) for a,b in zip(path,path[1:])};n=len(rows);out=[0]*n
    for a,b in combinations(range(n),2):
        ab=bool(rows[a]>>b&1)
        if frozenset((a,b)) not in keep:ab=not ab
        if ab:out[a]|=1<<b
        else:out[b]|=1<<a
    return out

def fas(rows):
    # Independent direction from search: select FIRST vertex, price incoming arcs.
    n=len(rows);N=1<<n;inc=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
    dp=[0]*N
    for S in range(1,N):
        vals=[];todo=S
        while todo:
            low=todo&-todo;todo-=low;v=low.bit_length()-1;R=S^low
            vals.append(dp[R]+(inc[v]&R).bit_count())
        dp[S]=min(vals)
    return dp[-1]

def small_checks():
    graph_rows=[];move_rows=[]
    for n in range(2,7):
        es=list(combinations(range(n),2));maximum=-1;forests=0
        for mask in range(1<<len(es)):
            ed=[e for k,e in enumerate(es) if mask>>k&1];rows=tour(n,ed);cs=components(n,ed);b=len(ed)-n+len(cs);l=L(rows)
            need(l<=n-2+b,f'cycle bound failed n={n} mask={mask}')
            if len(cs)==1:need(l<=len(ed)-1,'connected strictness failed')
            if b==0:forests+=1;maximum=max(maximum,l)
            if n<=4:
                best=min(sum(abs(v-t) for v,t in zip(score(rows),pi)) for pi in permutations(range(n)))
                need(best==2*l,'assignment comparison failed')
        graph_rows.append({'n':n,'graphs':1<<len(es),'forests':forests,'max_forest_L':maximum})
    for n in range(2,6):
        es=list(combinations(range(n),2));moves=0
        for mask in range(1<<len(es)):
            rows=tour(n,[e for k,e in enumerate(es) if mask>>k&1])
            for path in permutations(range(n)):
                if legal(rows,path):
                    after=tau(rows,path);moves+=1
                    need(abs(L(after)-L(rows))<=1,'inherited move bound failed')
                    need(tau(after,path)==rows,'same-path involution failed')
        move_rows.append({'n':n,'legal_moves':moves})
    return graph_rows,move_rows

def main():
    p=Path(__file__).resolve().parent;cert=json.loads((p/'certificate.json').read_text());n=cert['n'];need(n==11,'wrong order')
    qr={x*x%n for x in range(1,n)};rows=[sum(1<<y for y in range(n) if (y-x)%n in qr) for x in range(n)]
    need(score(rows)==[5]*n and L(rows)==15 and fas(rows)==20,'initial Paley controls failed')
    outputs=[]
    for k,step in enumerate(cert['steps'],1):
        need(step['depth']==k,'out-of-order certificate');before=rows;rows=tau(before,step['path']);need(rows==step['rows'],'literal rows mismatch')
        need(tau(rows,step['path'])==before,'witness involution failure')
        predicted=[n-1-d for d in score(before)];predicted[step['path'][0]]+=1;predicted[step['path'][-1]]-=1
        need(score(rows)==predicted,'endpoint score update mismatch')
        fr=frame(rows,step['forward_order']);l=L(rows);f=fas(rows)
        need(l==15-k==step['L'],'L mismatch');need(f==step['f'],'f mismatch');need(fr['beta']==6-k==step['beta'],'beta mismatch');need(fr['c']==step['c'],'component mismatch')
        need(fr['m']==f,'listed frame not feedback optimal')
        outputs.append({'depth':k,'L':l,'f':f,'frame_beta':fr['beta'],'frame_edges':fr['edges'],'components':fr['components']})
    need(len(outputs)==6 and outputs[-1]['frame_beta']==0 and len(outputs[-1]['components'])==1,'final tree missing')
    # Independent sharp forest-score witness at n=11, not taken from move search.
    path=[4,7,0,8,1,9,2,10,3,5,6];ed=[tuple(sorted(e)) for e in zip(path,path[1:])];sharp=tour(11,ed)
    need(L(sharp)==9 and len(ed)-11+len(components(11,ed))==0,'sharp forest failed')
    staircase=[]
    for r in range(1,16):
        nn=2*r+1;ed=[(i,j) for i in range(r) for j in range(i+r+1,2*r+1)];t=tour(nn,ed);b=len(ed)-nn+len(components(nn,ed))
        need(score(t)==[r]*nn,'staircase not regular');need(b==(r-1)*(r-2)//2==max(0,L(t)-nn+2),'staircase equality failed')
        staircase.append({'n':nn,'L':L(t),'frame_beta':b})
    # Hostile: edge count <= n-1 does NOT imply a forest or L <= n-2.
    ed=[(i,j) for i in range(3) for j in range(i+4,7)];hostile=tour(7,ed)
    need(len(ed)==6 and L(hostile)==6 and len(ed)-7+len(components(7,ed))==1,'edge-count hostile missing')
    graphs,moves=small_checks()
    print(json.dumps({'status':'PASS','universe':'All simple ordered backward graphs n=2..6; all labelled tournament/legal Hamilton-path pairs n=2..5; six named P11 transitions; staircase r=1..15','initial':{'L':15,'f':20},'witness_steps':outputs,'all_ordered_graph_checks':graphs,'all_legal_move_checks':moves,'staircase_controls':staircase,'sharp_forest':{'scores':score(sharp),'L':L(sharp)},'hostile_edgecount_not_forest':{'n':7,'edges':6,'L':6,'beta':1},'theorem_dependent_conclusions':{'exact_P11_forest_distance':6,'exact_min_beta_at_depth_k_1_through_6':[5,4,3,2,1,0]}},indent=2,sort_keys=True))
if __name__=='__main__': main()
