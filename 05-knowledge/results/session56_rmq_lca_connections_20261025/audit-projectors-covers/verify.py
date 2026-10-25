#!/usr/bin/env python3
"""Independent exact audit. Uses block intersections and integer triangular solve.
No producer imports, removable assertions, floats, or repository mutation.
"""
import hashlib,json,sys
from pathlib import Path
from itertools import combinations,product
from fractions import Fraction as Q
from functools import lru_cache
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def partitions(n):
    if not n:return [frozenset()]
    ans=[]
    for p in partitions(n-1):
        ans.append(p|{frozenset({n-1})})
        for b in p:ans.append((p-{b})|{b|{n-1}})
    return ans

def refine(a,b):return all(any(x<=y for y in b) for x in a)
def meet(a,b):return frozenset(x&y for x in a for y in b if x&y)
def algebra(states):
    ix={p:i for i,p in enumerate(states)};N=len(states)
    z=lambda v:[sum(v[j] for j in range(N) if refine(states[i],states[j])) for i in range(N)]
    def conv(v,w):
        out=[0]*N
        for i,a in enumerate(v):
            for j,b in enumerate(w):
                if a and b:out[ix[meet(states[i],states[j])]]+=a*b
        return out
    # Solve Z e_i = standard basis from coarse to fine, not Mobius recurrence.
    es=[]
    for k in range(N):
        v=[0]*N
        for i in reversed(range(N)):
            v[i]=int(i==k)-sum(v[j] for j in range(i+1,N) if refine(states[i],states[j]))
        need(z(v)==[int(i==k) for i in range(N)],'triangular inverse')
        es.append(v)
    for i,a in enumerate(es):
        for j,b in enumerate(es):need(conv(a,b)==(a if i==j else [0]*N),'projector products')
    return z,conv,es

def compositions(n,k):
    if k==1:yield (n,);return
    for x in range(n+1):
        for tail in compositions(n-x,k-1):yield (x,)+tail

def cross(a,b):return bool(a&b and a-b and b-a)
def edges(F):return {(i,j) for i,j in combinations(range(len(F)),2) if cross(F[i],F[j])}
def coloring(n,E):
    # Exhaustive canonical color assignments, fixed original vertex order.
    if not n:return 0,[]
    for cap in range(1,n+1):
        def solve(colors):
            i=len(colors)
            if i==n:return colors
            for c in range(min(cap,1+max(colors,default=-1)+1)):
                if all(colors[j]!=c for j in range(i) if (j,i) in E):
                    answer=solve(colors+[c])
                    if answer is not None:return answer
            return None
        answer=solve([])
        if answer is not None:return cap,answer

def cover(F):
    n=len(F)
    good=[]
    for m in range(1,1<<n):
        chosen=[F[i] for i in range(n) if m>>i&1]
        if all(not (a&b) or a<=b or b<=a for a,b in combinations(chosen,2)):good.append(m)
    @lru_cache(None)
    def rec(mask):
        if not mask:return 0
        return 1+min(rec(mask^m) for m in good if m&mask==m and m&(mask&-mask))
    return rec((1<<n)-1)
def tree(U,F):
    need(bool(U),'empty universe excluded')
    need(all(a and a<=U for a in F),'empty or out-of-universe cluster excluded')
    H=set(F)|{U}|{frozenset({x}) for x in U}
    parents={}
    for a in H-{U}:
        supers=[b for b in H if a<b]
        mins=[b for b in supers if not any(a<c<b for c in supers)]
        need(len(mins)==1,'unique inclusion parent');parents[a]=mins[0]
    def descend(a):
        children=[b for b in H if parents.get(b)==a]
        if not children:
            need(len(a)==1,'singleton leaves');return a
        need(all(not b&c for b,c in combinations(children,2)),'disjoint children')
        labels=frozenset().union(*(descend(b) for b in children))
        need(labels==a,'actual descendant cluster');return labels
    need(descend(U)==U,'fixed labels survive')

def check(F,U):
    k,c=coloring(len(F),edges(F));need(k==cover(F),'coloring versus laminar cover')
    for j in range(k):tree(U,[a for i,a in enumerate(F) if c[i]==j])
    return {'vertices':len(F),'edges':len(edges(F)),'trees':k}
def realize(n,E):return [frozenset([('v',v)]+[('e',a,b) for a,b in E if v in (a,b)]) for v in range(n)]

def main():
    if '--self-test-failure' in sys.argv:need(False,'intentional exception guard')
    pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
    for p in pins:need(hashlib.sha256((ROOT/p['snapshot']).read_bytes()).hexdigest()==p['sha256'],'source pin '+p['snapshot'])
    out={'source_pins':len(pins),'projector_state_counts':[]}
    for n in range(1,6):
        ps=sorted(partitions(n),key=lambda p:-len(p));z,conv,es=algebra(ps)
        out['projector_state_counts'].append(len(ps))
        for i in range(len(ps)):
            for j in range(len(ps)):
                a=[int(x==i) for x in range(len(ps))];b=[int(x==j) for x in range(len(ps))]
                need(z(conv(a,b))==[u*v for u,v in zip(z(a),z(b))],'basis zeta multiplication')
    ps=sorted(partitions(3),key=lambda p:-len(p));z,conv,es=algebra(ps);N=len(ps)
    positive=0
    for bits in product((0,1),repeat=N):
        v=[sum(bits[j]*es[j][i] for j in range(N)) for i in range(N)]
        need(conv(v,v)==v,'all signed idempotents')
        positive+=all(x>=0 for x in v)
    need(positive==6,'positive signed idempotent count')
    count=stable=0
    for d in range(1,9):
        for nums in compositions(d,N):
            p=[Q(x,d) for x in nums];same=conv(p,p)==p
            need(same==(sum(x>0 for x in p)==1),'probability iff Dirac')
            count+=1;stable+=same
    p=[Q(int(i in (1,2)),2) for i in range(N)]
    need(conv(p,p)!=p and conv(p,p)[0]==Q(1,2),'independent mixture')
    coupled=[sum(p[j] for j in range(N) if meet(ps[j],ps[j])==ps[i]) for i in range(N)]
    need(coupled==p,'coupled diagonal')
    need([sum(Q(1,4) for a,b in product(range(2),repeat=2) if (a+b)%2==i) for i in range(2)]==[Q(1,2)]*2,'C2 group countercontrol')
    # Non-top meet-semilattice: bottom plus two incomparable atoms.
    z0,c0,e0=algebra(ps[:3]);identity=[sum(e[i] for e in e0) for i in range(3)]
    need(identity==[-1,1,1],'identity need not be a Dirac without top')
    out['probability_grid']={'count':count,'stable':stable,'nonnegative_signed_idempotents':positive}
    out['no_top_identity']=identity
    U=frozenset(range(3));subsets=[frozenset(i for i in U if m>>i&1) for m in range(1,8)]
    for mask in range(128):check([s for i,s in enumerate(subsets) if mask>>i&1],U)
    out['cluster_families_n3']=128
    out['graph_realizations']={}
    for n in range(6):
        poss=list(combinations(range(n),2))
        for mask in range(1<<len(poss)):
            E={e for i,e in enumerate(poss) if mask>>i&1};F=realize(n,E)
            need(edges(F)==E,'private graph realization')
            check(F,frozenset().union(*F) if F else frozenset({'dummy'}))
        out['graph_realizations'][str(n)]=1<<len(poss)
    U=frozenset(product(range(3),repeat=2))
    F=[frozenset((x,y) for x,y in U if x==b) for b in range(3)]+[frozenset((x,y) for x,y in U if (y-m*x)%3==b) for m in range(3) for b in range(3)]
    out['affine']=check(F,U);need(out['affine']=={'vertices':12,'edges':54,'trees':4},'affine K3333')
    U=frozenset(range(6));F=[frozenset(i for i in U if i%p==r) for p in (2,3) for r in range(p)]
    out['six_label_CRT']=check(F,U);need(out['six_label_CRT']['trees']==2,'six label CRT')
    D={a&b for a in F[:2] for b in F[2:] if a&b}|{U}
    out['diagonal']=check(list(D),U);need(not any(a in D for a in F),'independent blocks absent')
    A=[frozenset(s) for s in ({0,4,8},{2,3,7},{1,5,6})];B=[frozenset(s) for s in ({2,4,6},{0,5,7},{1,3,8})]
    need(all(len(a&b)==1 for a in A for b in B),'Session55 orthogonality')
    out['Session55_nine_cell_AB']=check(A+B,frozenset(range(9)))
    need(out['Session55_nine_cell_AB']=={'vertices':6,'edges':9,'trees':2},'actual K33')
    E={(0,1),(1,2),(2,3),(3,4),(0,4)};F=realize(5,E)
    out['C5']=check(F,frozenset().union(*F))
    clique=max(len(s) for r in range(6) for s in combinations(range(5),r) if all(tuple(sorted(e)) in E for e in combinations(s,2)))
    need(clique==2 and out['C5']['trees']==3,'clique not chromatic');out['C5']['clique']=clique
    need(check([],frozenset({0}))['trees']==0,'empty family')
    need(check([frozenset({0})],frozenset({0}))['trees']==1,'singleton universe')
    try:tree(frozenset({0}),[frozenset()])
    except RuntimeError:pass
    else:raise RuntimeError('empty cluster should reject')
    out['status']='PASS';print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
