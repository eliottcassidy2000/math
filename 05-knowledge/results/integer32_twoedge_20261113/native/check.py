#!/usr/bin/env python3
"""Independent small-order replay of the marked native cut presentation.

No imports from prior producers. Integer codes use increasing-arc bit1;
the proof's reversal mask is its bitwise complement. beta ignores that
global complement. The only prior data input is the pinned matrix census.
"""
from collections import Counter
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
DATA=ROOT/'05-knowledge/results/incidence_closure_20261113/closure/data.json'
EXPECTED='c0484a99712edff88e8158b5f6000097f32f88874150d059da5aefa20f3de8f7'
CHECKS=0

def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

class Universe:
    def __init__(self,n):
        self.n=n
        self.E=tuple(combinations(range(n),2))
        self.pos={e:k for k,e in enumerate(self.E)}
        self.top=(1<<len(self.E))-1
        self.perms=tuple(permutations(range(n)))
        self.W=tuple((v,u,w) for v in range(n) for u,w in combinations([z for z in range(n) if z!=v],2))
        self.links=tuple((self.bit(v,u),self.bit(v,w)) for v,u,w in self.W)
        self.base=self.beta(self.top)
        # Independent orbit discovery: fill the entire labelled universe by
        # actual permutation orbits, taking the least unused code each time.
        self.owner=[-1]*(self.top+1)
        self.reps=[]
        self.aut={}
        for t in range(self.top+1):
            if self.owner[t]>=0:continue
            orbit={self.relabel(t,p) for p in self.perms}
            need(min(orbit)==t,'least-unseen canonical representative')
            need(all(self.owner[c]<0 for c in orbit),'disjoint permutation orbits')
            for c in orbit:self.owner[c]=t
            self.reps.append(t)
            self.aut[t]=factorial(n)//len(orbit)
        self.dual={t:self.owner[t^self.top] for t in self.reps}
        self.merged=sorted({min(t,self.dual[t]) for t in self.reps})

    def bit(self,u,v):return 1<<self.pos[min(u,v),max(u,v)]
    def arc(self,t,u,v):return bool(t&self.bit(u,v))==(u<v)
    def relabel(self,t,p):return sum(self.arc(t,p[u],p[v])<<k for k,(u,v) in enumerate(self.E))
    def beta(self,t):return sum((self.arc(t,v,u)!=self.arc(t,v,w))<<k for k,(v,u,w) in enumerate(self.W))
    def cut(self,s):return sum((bool(s&a)!=bool(s&b))<<k for k,(a,b) in enumerate(self.links))
    def hp(self,t,p):return all(self.arc(t,u,v) for u,v in zip(p,p[1:]))
    def support(self,p):return sum(self.bit(u,v) for u,v in zip(p,p[1:]))
    def rho(self,t):return self.relabel(t^self.top,tuple(reversed(range(self.n))))
    def merged_code(self,t):return min(self.owner[t],self.owner[t^self.top])

    def cut_inverse(self,b):
        # Integrate relative to edge0; verify all links, hence reject noncuts.
        values={1:False};todo=[1]
        adj={1<<k:[] for k in range(len(self.E))}
        for k,(a,c) in enumerate(self.links):
            d=bool(b>>k&1);adj[a].append((c,d));adj[c].append((a,d))
        while todo:
            a=todo.pop()
            for c,d in adj[a]:
                value=values[a]^d
                if c in values:
                    if values[c]!=value:return None
                else:values[c]=value;todo.append(c)
        return sum(k for k,v in values.items() if v) if len(values)==len(self.E) else None

    def directed_path_support(self,t,s):
        # No permutation enumeration: degrees, connected traversal, orientation.
        adj=[[] for _ in range(self.n)]
        for k,(a,b) in enumerate(self.E):
            if s>>k&1:adj[a].append(b);adj[b].append(a)
        ends=[v for v in range(self.n) if len(adj[v])==1]
        if len(ends)!=2 or any(len(a) not in (1,2) for a in adj):return False
        p=[ends[0]];previous=-1
        while True:
            nxt=[v for v in adj[p[-1]] if v!=previous]
            if not nxt:break
            previous=p[-1];p.append(nxt[0])
            if len(p)>self.n:return False
        return len(p)==self.n and (self.hp(t,p) or self.hp(t,p[::-1]))

    def guard_supports(self,t,jump):
        s=self.cut_inverse(jump)
        if s is None:return set()
        return {q for q in (s,s^self.top) if self.directed_path_support(t,q)}

def sparse(counts):return [list(k)+[v] for k,v in sorted(counts.items())]

def run(n,prior):
    U=Universe(n);M=len(U.E)
    need(U.reps==[r['code'] for r in prior['ordinary_nodes']],'independent ordinary census')
    need(U.merged==[r['code'] for r in prior['merged_nodes']],'independent merged census')
    need(len({U.beta(t) for t in range(U.top+1)})==1<<(M-1),'labelled converse fields')
    for t in range(U.top+1):
        need(U.beta(t)==U.base^U.cut(t),'coordinate equation')
    O=Counter();OB=Counter();A=Counter();B=Counter()
    spine=U.support(tuple(range(n)));free=[1<<k for k in range(M) if not spine>>k&1]
    blue_count=0;rho_orbits=set();distance=2*(n-2)**2
    for w in range(1<<len(free)):
        t=spine|sum(b for k,b in enumerate(free) if w>>k&1)
        target=t^U.top^spine
        source_code=U.owner[t];target_code=U.owner[target]
        blue=U.rho(t)==t
        need(U.rho(U.rho(t))==t,'reflection involution')
        need(U.rho(target)==U.rho(t)^U.top^spine,'native/reflection commute')
        need(U.beta(target)==U.beta(t)^U.cut(spine),'native wedge update')
        need((U.beta(target)^U.beta(t)).bit_count()==distance,'native jump weight')
        need(U.beta(target)!=U.beta(t),'no labelled native loops')
        O[source_code,target_code]+=1
        A[U.merged_code(t),U.merged_code(target)]+=1
        if blue:
            blue_count+=1;OB[source_code,target_code]+=1
            B[U.merged_code(t),U.merged_code(target)]+=1
        rho_orbits.add(min(t,U.rho(t)))
    for name,counts in [('ordinary_A',O),('ordinary_blue',OB),('merged_A',A),('merged_blue',B)]:
        need(sparse(counts)==prior[name],('entrywise prior matrix',n,name))
    need(blue_count==1<<((n-1)**2//4),'blue fixed-word formula')
    need(len(rho_orbits)==((1<<len(free))+blue_count)//2,'premature phase quotient')
    total_paths=0;max_supports=0;double_states=0
    for t in U.reps:
        jumps={};counts=Counter();bcounts=Counter()
        hp=[p for p in U.perms if U.hp(t,p)];total_paths+=len(hp)
        for p in hp:
            s=U.support(p);j=U.cut(s);jumps.setdefault(j,set()).add(s)
            target=t^U.top^s
            need(U.beta(target)==U.beta(t)^j,'every representative path update')
            need(j.bit_count()==distance,'every representative path weight')
            counts[U.owner[target]]+=1
            # Relabelling the chosen path to the fixed spine reads reflection.
            in_spine=U.relabel(t,p)
            if U.rho(in_spine)==in_spine:bcounts[U.owner[target]]+=1
        for j,supports in jumps.items():
            need(U.guard_supports(t,j)==supports,'exact candidate guard on every move')
            need(len(supports)<=(2 if n==4 else 1),'sharp jump-support multiplicity')
            max_supports=max(max_supports,len(supports));double_states+=len(supports)==2
        if n<=5:
            for target in range(0,U.top+1,2):
                j=U.beta(target)^U.beta(t)
                need(U.guard_supports(t,j)==jumps.get(j,set()),'exhaustive target realization guard')
        need(all(v%U.aut[t]==0 for v in counts.values()),'path orbit divisibility')
        need({k[1]:v for k,v in O.items() if k[0]==t}=={k:v//U.aut[t] for k,v in counts.items()},'independent physical-path source weights')
        need({k[1]:v for k,v in OB.items() if k[0]==t}=={k:v//U.aut[t] for k,v in bcounts.items()},'independent physical-path blue weights')
        merged=U.merged_code(t);m=1+(U.dual[t]!=t)
        need(sum(v for (a,b),v in A.items() if a==merged)*U.aut[t]==m*len(hp),'merged source mass')
    need(max_supports==(2 if n==4 else 1),'support ambiguity boundary attained')
    for counts in (A,B):
        need(all(counts[b,a]==v for (a,b),v in counts.items()),'symmetric integer transport')
        need(all(v%2==0 for (a,b),v in counts.items() if a==b),'doubled native loops')
    need(all((v-B[k])%2==0 for k,v in A.items()),'black twin multiplicity')
    return {'n':n,'labelled_tournaments':1<<M,'ordinary_classes':len(U.reps),'merged_classes':len(U.merged),
            'fixed_spine_darts':1<<len(free),'blue_darts':blue_count,'premature_reflection_quotient':len(rho_orbits),
            'native_jump_weight':distance,'physical_paths':total_paths,'maximum_supports_per_jump':max_supports,
            'double_support_jumps_in_representatives':double_states,'exhaustive_target_guards':n<=5}

def hostile_controls():
    U=Universe(4);p=(0,1,2,3);q=(2,0,3,1)
    t=0
    for order in (p,q):
        for a,b in zip(order,order[1:]):
            if a<b:t|=U.bit(a,b)
    need(U.hp(t,p) and U.hp(t,q),'complementary path witness orientation')
    need(U.support(p)^U.support(q)==U.top,'complementary path witness covers K4')
    need(U.cut(U.support(p))==U.cut(U.support(q)),'order4 identical jumps')
    need(len(U.guard_supports(t,U.cut(U.support(p))))==2,'order4 exact double guard')
    mixed=U.support((0,2,1,3))
    need(not U.directed_path_support(U.top,mixed),'orientation guard hostile')
    V=Universe(5);disconnected=V.bit(0,1)|V.bit(1,2)|V.bit(0,2)|V.bit(3,4)
    degrees=[sum(bool(disconnected&V.bit(a,b)) for b in range(5) if a!=b) for a in range(5)]
    need(sorted(degrees)==[1,1,2,2,2],'degree-only false positive')
    need(not V.directed_path_support(V.top,disconnected),'connectivity guard hostile')
    return {'order4_tournament_code':t,'complementary_paths':[p,q],
            'degree_hostile':'triangle plus disjoint edge','orientation_hostile':[0,2,1,3]}

def main():
    need(digest(DATA)==EXPECTED,'immutable matrix input')
    old=json.loads(DATA.read_text());rows=[]
    for n in range(3,7):rows.append(run(n,next(r for r in old['rows'] if r['n']==n)))
    hostile=hostile_controls()
    print(json.dumps({'status':'FINITE-EXACT','checks':CHECKS,'producer_sha256':digest(Path(__file__)),
                      'input_sha256':EXPECTED,'rows':rows,'hostiles':hostile},sort_keys=True,indent=2))

if __name__=='__main__':main()
