#!/usr/bin/env python3
"""Independent named-graph / literal-path audit of realization.md.
Imports no producer code and uses no subset DP for small/contextual checks.
"""
from itertools import combinations,permutations,product
from math import comb,isqrt
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def require(ok,msg):
    if not ok:raise RuntimeError(msg)

def paths(vertices,arcs):
    if not vertices:return [()]
    out=[]
    def walk(word,remaining):
        if not remaining:out.append(tuple(word));return
        for v in remaining:
            if not word or (word[-1],v) in arcs:walk(word+[v],remaining-{v})
    walk([],set(vertices));return out

def build(blocks,internal,cut):
    vertices=[v for b in blocks for v in b]
    arcs=set(internal)
    for i,b in enumerate(blocks):
        for c in blocks[i+1:]:arcs|={(v,w) for v in b for w in c}
    for v in vertices:arcs.add(('x',v) if v in cut else (v,'x'))
    return tuple(vertices)+('x',),arcs

def check_transport(u,v,block):
    vertices,a=u;other,b=v
    require(set(vertices)==set(other),'lost vertex name')
    expected={(z,w) for w in block for z in ('z',)}|{(w,z) for w in block for z in ('z',)}
    require(a^b==expected,'edge changes beyond the named singleton/block cut')
    require({e for e in a if 'x' in e}=={e for e in b if 'x' in e},'changed apex ports')

def profile(vertices,arcs,cut):
    n=len(vertices);byset={}
    for k in range(n+1):
        for subset in combinations(vertices,k):byset[frozenset(subset)]=paths(subset,arcs)
    whole=frozenset(vertices);hs=byset[whole]
    H=len(hs);s=sum(p[0] in cut for p in hs);t=sum(p[-1] not in cut for p in hs)
    E=P=Q=0
    for subset,left in byset.items():
        for p in left:
            for q in byset[whole-subset]:
                E+=1;P+=bool(p) and p[-1] not in cut;Q+=bool(q) and q[0] in cut
    child_arcs=set(arcs)|{('x',v) if v in cut else (v,'x') for v in vertices}
    f=len(paths(tuple(vertices)+('x',),child_arcs))
    return H,E,P,Q,s,t,f

def mul(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def subtract(a,b):return [[x-y for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]

def check_matrices(prof):
    H,E,P,Q,s,t,f=prof;D=E-H-P
    M=[[H,P,H-s,f-H],[0,E,0,Q],[0,0,0,H-t],[0,0,0,H]]
    S0=[[1,1,1,0],[0,2,0,0],[0,0,0,0],[0,0,0,1]]
    S1=[[1,0,0,0],[0,2,0,1],[0,0,0,1],[0,0,0,1]]
    want0=[[0,D,-s,H+Q-t],[0,0,0,Q],[0,0,0,t-H],[0,0,0,0]]
    want1=[[0,P,s-H,H+P-s],[0,0,0,E-H-Q],[0,0,0,-t],[0,0,0,0]]
    require(subtract(mul(S0,M),mul(M,S0))==want0,'zero commutator')
    require(subtract(mul(M,S1),mul(S1,M))==want1,'one commutator')

def cover_profile(vertices,arcs):
    # A path cover is one permutation plus cuts containing every failed arc.
    n=len(vertices);result=[0]*(n+1)
    for p in permutations(vertices):
        bad=sum((p[i],p[i+1]) not in arcs for i in range(n-1));good=n-1-bad
        for k in range(bad+1,n+1):result[k]+=comb(good,k-1-bad)
    return result[1:]

def strong(vertices,arcs):
    for start in vertices:
        seen=set();todo=[start]
        while todo:
            v=todo.pop()
            if v in seen:continue
            seen.add(v);todo.extend(w for u,w in arcs if u==v and w not in seen)
        if seen!=set(vertices):return False
    return True

def family(r,a):
    R=[f'r{i}' for i in range(r)];C=['c0','c1','c2'];z=['z']
    internal={('c0','c1'),('c1','c2'),('c2','c0')}
    cut=set(R)|set(C[:a])
    pre=[[v] for v in R]
    return build(pre+[z,C],internal,cut),build(pre+[C,z],internal,cut)

def main():
    profiles=0;contexts=0;boundaries=0
    for n in range(1,5):
        B=[f'b{i}' for i in range(n)];pairs=list(combinations(B,2))
        for mask in range(1<<len(pairs)):
            internal={(v,w) if mask>>i&1 else (w,v) for i,(v,w) in enumerate(pairs)}
            for bits in range(1<<n):
                cut={v for i,v in enumerate(B) if bits>>i&1};prof=profile(B,internal,cut)
                check_matrices(prof);H,E,P,Q,s,t,f=prof
                if n>=2:require(Q-s>=len(cut) and P-t>=n-len(cut),'proper-cover injection')
                profiles+=1
                if n not in (2,3) and not(n==4 and mask==4):continue
                if n==4:require(strong(B,internal),'named asymmetric strong4 control')
                # None=empty context; 0/1=one named singleton with fixed cut.
                for pe,se,e in product((None,0,1),(None,0,1),(0,1)):
                    pre=[] if pe is None else [['pre']];suf=[] if se is None else [['post']]
                    totalcut=cut|({'pre'} if pe else set())|({'post'} if se else set())|({'z'} if e else set())
                    u=build(pre+[['z'],B]+suf,internal,totalcut)
                    v=build(pre+[B,['z']]+suf,internal,totalcut)
                    check_transport(u,v,B)
                    hu,hv=len(paths(*u)),len(paths(*v))
                    if e==0:
                        require((hv>hu) if cut else (hv==hu),'zero sorting sign')
                    else:
                        require((hu>hv) if len(cut)<n else (hu==hv),'one sorting sign')
                    if not cut or len(cut)==n:boundaries+=1
                    contexts+=1
    families=[]
    for r in range(1,6):
        for a in (1,2):
            u,v=family(r,a);check_transport(u,v,['c0','c1','c2'])
            # DFS enumerates only valid literal vertex words, not subset DP.
            hu,hv=len(paths(*u)),len(paths(*v))
            want=(15*2**r-7,18*2**r-7) if a==1 else (9*2**r-1,15*2**r-1)
            require((hu,hv)==want,'family exact counts')
            require(strong(*u) and strong(*v),'family not strong')
            row={'r':r,'a':a,'vertices':r+5,'H':[hu,hv],'changed_named_edges':3}
            if r<=3:
                # A second literal route checks every permutation, including rejected words.
                literal=[]
                for vertices,arcs in (u,v):
                    literal.append(sum(all((p[i],p[i+1]) in arcs for i in range(len(p)-1)) for p in permutations(vertices)))
                require(literal==[hu,hv],'full permutations disagree with adjacency DFS')
                parents=[]
                for vertices,arcs in (u,v):
                    old=tuple(w for w in vertices if w!='x');old_arcs={e for e in arcs if 'x' not in e}
                    parents.append(cover_profile(old,old_arcs))
                require(parents[0]==parents[1] and parents[0][0]==3,'parent cover profile changes')
                row['parent_ordered_nonempty_cover_profile']=parents[0]
            families.append(row)
    require(all(all(p%d for d in range(2,isqrt(p)+1)) for p in (5,17,23,29,41,59,71,479)),
            'displayed prime lane endpoint not prime')
    require(35==7*5 and 119==7*17 and 287==7*41,'lane factors wrong')
    # The r>=1 clock guard is necessary, not cosmetic.
    zero=[]
    for a in (1,2):
        u,v=family(0,a);zero.append({'a':a,'actual_H':[len(paths(*u)),len(paths(*v))]})
    require(zero==[{'a':1,'actual_H':[5,9]},{'a':2,'actual_H':[5,13]}],'r0 hostile drift')
    result={'status':'PASS','profile_universe':'all labelled tournaments n1..4 with every apex cut',
            'literal_profiles':profiles,'named_contextual_swaps':contexts,'empty_full_context_controls':boundaries,
            'families':families,'r0_outside_formula_scope':zero,
            'independence':'No producer imports; adjacency DFS, full permutations, and cut-of-permutation cover counts.'}
    out=json.dumps(result,indent=2,sort_keys=True)+'\n';(HERE/'audit_realization.json').write_text(out);print(out,end='')

if __name__=='__main__':main()
