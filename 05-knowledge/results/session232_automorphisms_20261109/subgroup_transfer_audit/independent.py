#!/usr/bin/env python3
"""Independent finite controls; no producer or repository imports."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations, combinations, product
from pathlib import Path
import hashlib, json, sys

def need(ok, message):
    if not ok: raise RuntimeError(message)

def mat(n, code):
    t=[[False]*n for _ in range(n)]
    for i,(u,v) in enumerate(combinations(range(n),2)):
        t[u][v]=bool((code>>i)&1);t[v][u]=not t[u][v]
    return t

def cayley(n,steps): return [[(v-u)%n in steps for v in range(n)] for u in range(n)]
def join(*ts):
    t=[]; before=0; n=sum(map(len,ts))
    for b in ts:
        for row in b:t.append([False]*before+row+[True]*(n-before-len(b)))
        before+=len(b)
    return t

def cyclic_join(a,b,c):
    ts=[a,b,c]; n=sum(map(len,ts));offsets=[0,len(a),len(a)+len(b)];t=[[False]*n for _ in range(n)]
    for i,x in enumerate(ts):
        for u in range(len(x)):
            for v in range(len(x)):t[offsets[i]+u][offsets[i]+v]=x[u][v]
            for v in range(len(ts[(i+1)%3])):t[offsets[i]+u][offsets[(i+1)%3]+v]=True
    return t

def aut(t):
    n=len(t); scores=list(map(sum,t)); cells=[[v for v in range(n) if scores[v]==d] for d in sorted(set(scores))]
    out=[]
    for choices in product(*(permutations(c) for c in cells)):
        p=[0]*n
        for cell,im in zip(cells,choices):
            for x,y in zip(cell,im):p[x]=y
        if all(t[u][v]==t[p[u]][p[v]] for u,v in combinations(range(n),2)):out.append(tuple(p))
    return out

def subs(group):
    identity=tuple(range(len(group[0])))
    def mul(a,b):return tuple(a[b[v]] for v in range(len(a)))
    result={frozenset((identity,))}; work=list(result)
    while work:
        h=work.pop()
        for g in group:
            if g in h:continue
            k=set(h); k.add(g)
            while True:
                z={mul(a,b) for a in k for b in k}
                if z<=k:break
                k|=z
            fk=frozenset(k)
            if fk not in result:result.add(fk);work.append(fk)
    return sorted(result,key=lambda h:(len(h),tuple(sorted(h))))

def hcount(t):
    n=len(t);dp={}
    for size in range(1,n+1):
        for vs in combinations(range(n),size):
            s=frozenset(vs)
            for v in vs:
                dp[s,v]=1 if size==1 else sum(dp[s-{v},u] for u in s-{v} if t[u][v])
    return sum(dp[frozenset(range(n)),v] for v in range(n))

def sccs(t):
    n=len(t); reach=[list(row) for row in t]
    for u in range(n):reach[u][u]=True
    for k in range(n):
        for u in range(n):
            for v in range(n):reach[u][v]=reach[u][v] or (reach[u][k] and reach[k][v])
    unseen=set(range(n));out=[]
    while unseen:
        u=min(unseen);b={v for v in unseen if reach[u][v] and reach[v][u]};out.append(sorted(b));unseen-=b
    return out

def induced(t,b):return [[t[u][v] for v in b] for u in b]
def val(x,p):
    n=0
    while x%p==0:n+=1;x//=p
    return n

def allowed_stack(t):
    factors=Counter()
    for b in sccs(t):
        u=induced(t,b); key=(len(b),hcount(u),len(aut(u)))
        if key not in ((1,1,1),(3,3,3),(5,15,5),(7,189,21)):return False,dict(factors)
        factors[len(b)]+=1
    return True,dict(factors)

def subgroup_check(t,group=None):
    group=aut(t) if group is None else group; h=hcount(t); target,factors=allowed_stack(t); count=0; equalities=0
    for a in subs(group):
        need(h%len(a)==0,'subgroup divisibility');q=h//len(a);count+=1
        need(q not in (7,21),'subgroup gaps')
        need(q==(len(group)//len(a))*(h//len(group)),'index law')
        if q<23:
            z=len(a)
            for p in (3,5,7):
                while z%p==0:z//=p
            need(z==1,'low q prime support')
            cost=3**(val(len(a),5)+2*val(len(a),7))
            need(q>=cost,'prime cost')
            need(q==cost or q>=3*cost,'strengthened equality gap')
            if len(a)%7==0:need(q==9,'sharp seven stratum')
            if q<9 and len(a)%5==0:need(q==3,'sharp five stratum')
            expected=target and len(a)==len(group)
            need((q==cost)==expected,'equality iff stack and full group')
            equalities+=q==cost
    return count,equalities

def prefix_transfer(t,g=None):
    n=len(t);g=aut(t) if g is None else g
    state_id={}; states=[]; allstates={(frozenset(vs),u) for r in range(1,n+1) for vs in combinations(range(n),r) for u in vs}
    while allstates:
        state=min(allstates,key=lambda st:(len(st[0]),tuple(sorted(st[0])),st[1]));s,u=state
        orb={(frozenset(p[v] for v in s),p[u]) for p in g}
        i=len(states);states.append(state)
        for a in orb:state_id[a]=i
        allstates-=orb
    k=[];state_stab=[]
    for s,u in states:
        k.append(sum(all(p[v]==v for v in s) for p in g))
        state_stab.append(sum(p[u]==u and {p[v] for v in s}==s for p in g))
        if len(s)>=n-2:need(k[-1]==1,'terminal pointwise stabilization')
    c=[0]*len(states);mass=[Fraction(0)]*len(states);naive=[0]*len(states);bad=[Fraction(0)]*len(states)
    correction=0;type_difference=0;edges=[]
    for b,(s,u) in enumerate(states):
        if len(s)==1:c[b]=naive[b]=1;mass[b]=Fraction(1,k[b]);bad[b]=1
        count=Counter(state_id[s|{v},v] for v in range(n) if v not in s and t[u][v])
        for a,N in count.items():
            w=Fraction(N*k[a],k[b]);need(w.denominator==1,'integral orbit transfer')
            c[a]+=c[b]*int(w);mass[a]+=mass[b]*N;naive[a]+=naive[b]*N
            bad[a]+=bad[b]*Fraction(N*state_stab[a],state_stab[b])
            correction+=w!=N;edges.append((b,a,N,k[b],k[a],int(w)))
    paths=defaultdict(set)
    def visit(p,s):
        u=p[-1];b=state_id[frozenset(s),u]
        canonical=min(tuple(x[v] for v in p) for x in g)
        paths[b].add(canonical)
        need(sum(all(x[v]==v for v in p) for x in g)==k[b],'path stabilizer equals pointwise')
        for v in range(n):
            if v not in s and t[u][v]:visit(p+(v,),s|{v})
    for u in range(n):visit((u,),{u})
    for b,st in enumerate(states):
        need(c[b]==len(paths[b]),'literal prefix orbit count')
        need(mass[b]==Fraction(c[b],k[b]),'groupoid measure')
        type_difference+=state_stab[b]!=k[b]
    full=[i for i,(s,u) in enumerate(states) if len(s)==n]
    final=sum(c[i] for i in full);need(final*len(g)==hcount(t),'final free quotient')
    return dict(n=n,group=len(g),H=hcount(t),orbits=final,states=len(states),corrections=correction,
                naive=sum(naive[i] for i in full),wrong_state_stabilizer_result=str(sum(bad[i] for i in full)),
                state_vs_path_stabilizer_differences=type_difference)

def two_path_covers(t):
    total=0
    for p in permutations(range(len(t))):
        for split in range(1,len(t)):
            total+=all(t[p[i]][p[i+1]] for i in range(len(t)-1) if i!=split-1)
    by_subsets=sum(hcount(induced(t,vs))*hcount(induced(t,[v for v in range(len(t)) if v not in vs])) for r in range(1,len(t)) for vs in combinations(range(len(t)),r))
    need(total==by_subsets,'two path independent subset and literal counts')
    return total

def partition_lemma(t):
    n=len(t);ans=0
    def parts(v,blocks):
        if v==n:yield tuple(tuple(b) for b in blocks);return
        for j in range(len(blocks)):
            blocks[j].append(v);yield from parts(v+1,blocks);blocks[j].pop()
        blocks.append([v]);yield from parts(v+1,blocks);blocks.pop()
    for bs in parts(0,[]):
        if any(len(sccs(induced(t,b)))!=1 for b in bs):continue
        block={v:i for i,b in enumerate(bs) for v in b}
        cross=any((t[u][v] and t[v][w] and t[w][u] or t[v][u] and t[u][w] and t[w][v]) and len({block[u],block[v],block[w]})>1 for u,v,w in combinations(range(n),3))
        if not cross:
            need({frozenset(b) for b in bs}=={frozenset(b) for b in sccs(t)},'no cross triangle iff SCC blocks')
            ans+=1
    return ans

def main():
    if '--sentinel' in sys.argv:need(False,'deliberate independent sentinel')
    counts=Counter(); subgroup_pairs=eq=state_count=corrected=partition_cases=0
    for n in range(1,6):
        for mask in range(1<<(n*(n-1)//2)):
            t=mat(n,mask);g=aut(t);a,b=subgroup_check(t,g);subgroup_pairs+=a;eq+=b
            out=prefix_transfer(t,g);state_count+=out['states'];corrected+=out['corrections'];counts[n]+=1
            partition_cases+=partition_lemma(t)
    c3=cayley(3,{1});r5=cayley(5,{1,2});p7=cayley(7,{1,2,4});q7=cayley(7,{1,2,3});one=[[False]]
    selected={}
    for name,t in [('source_C3',join(one,c3)),('C3_sink',join(c3,one)),('two_C3',join(c3,c3)),('R5',r5),('P7',p7),('carousel7',q7)]:
        a,b=subgroup_check(t);selected[name]=dict(prefix_transfer(t),subgroups=a,equality_pairs=b)
    need(selected['source_C3']['naive']==3 and selected['source_C3']['orbits']==1,'naive hostile')
    need(selected['C3_sink']['wrong_state_stabilizer_result']=='3' and selected['C3_sink']['orbits']==1,'setwise stabilizer hostile')
    diag=[tuple((v//3)*3+(v+j)%3 for v in range(6)) for j in range(3)]
    selected['diagonal_C3']=prefix_transfer(join(c3,c3),diag)
    need(selected['diagonal_C3']['orbits']==3,'diagonal subgroup')
    for name,t,expected in [('two_R5',join(r5,r5),(225,25,9)),('P7_C3',join(p7,c3),(567,63,9))]:
        g=aut(t);a,b=subgroup_check(t,g);actual=(hcount(t),len(g),hcount(t)//len(g));need(actual==expected,'large equality control')
        selected[name]=dict(H=actual[0],group=actual[1],orbits=actual[2],subgroups=a,equality_pairs=b)
    for name,t,f2 in [('P7',p7,1638),('R5',r5,90)]:
        result=two_path_covers(t);need(result==f2,'nonempty two-path covers');selected[name]['F2_nonempty_ordered']=result
    for name,t,expected in [('cyclic_P7_singletons',cyclic_join(p7,one,one),(2205,21,105)),('cyclic_R5_singletons',cyclic_join(r5,one,one),(135,5,27))]:
        g=aut(t);actual=(hcount(t),len(g),hcount(t)//len(g));need(actual==expected,'nonordinal positive triangle gap')
        selected[name]=dict(H=actual[0],group=actual[1],orbits=actual[2])
    # Exhaustive arithmetic control for legal exceptional-block product orders.
    arith=0
    for c,r,p in product(range(7),range(3),range(2)):
        P=3**(c+p)*5**r*7**p
        for a3,a5,a7 in product(range(c+p+1),range(r+1),range(p+1)):
            A=3**a3*5**a5*7**a7;lower=(P//A)*3**(r+2*p)
            if lower<23:
                cost=3**(a5+2*a7);need(lower>=cost,'arithmetic cost');need((lower==cost)==(A==P),'arithmetic equality');arith+=1
    print(json.dumps(dict(status='PASS',universe=dict(counts),subgroup_pairs=subgroup_pairs,equality_pairs=eq,
        endpoint_state_orbits=state_count,corrected_edges=corrected,no_cross_triangle_partitions=partition_cases,
        arithmetic_cases=arith,selected=selected,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
if __name__=='__main__':main()
