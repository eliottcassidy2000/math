"""Independent exact checks for the nonconsecutive-port theorems.

Standard-library only. All gates are explicit and survive python -O.
"""
from itertools import combinations
from collections import Counter
import json
import random


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def pair(a,b):
    return (min(a,b),max(a,b))


def pathset(q):
    return {pair(a,b) for a,b in zip(q,q[1:])}


class Equations:
    def __init__(self,n):
        self.n=n
        self.edges=list(combinations(range(n),2))
        self.index={e:i for i,e in enumerate(self.edges)}
        self.const=len(self.edges)
        self.parent=list(range(self.const+1))
        self.depth=[0]*(self.const+1)
        self.delta=[0]*(self.const+1)
        self.ok=True

    def clone(self):
        out=Equations(self.n)
        out.parent=self.parent[:];out.depth=self.depth[:]
        out.delta=self.delta[:];out.ok=self.ok
        return out

    def root(self,i):
        parity=0
        while self.parent[i]!=i:
            parity ^= self.delta[i]
            i=self.parent[i]
        return i,parity

    def equal(self,i,j,value):
        a,u=self.root(i);b,v=self.root(j)
        if a==b:
            self.ok &= ((u^v)==value)
            return
        if self.depth[a]<self.depth[b]:
            a,b=b,a
        self.parent[b]=a
        self.delta[b]=u^v^value
        if self.depth[a]==self.depth[b]:
            self.depth[a]+=1

    def force(self,a,b,value=1):
        self.equal(self.index[pair(a,b)],self.const,value^(a>b))

    def forced(self,a,b,value):
        r,u=self.root(self.index[pair(a,b)])
        c,v=self.root(self.const)
        return r==c and (u^v)==(value^(a>b))

    def anti(self,q,flips=frozenset()):
        r=dict(zip(q,reversed(q)))
        for a,b in self.edges:
            x,y=r[b],r[a]
            self.equal(self.index[(a,b)],self.index[pair(x,y)],
                       (x>y)^((a,b) in flips)^(pair(x,y) in flips))

    def model(self,seed):
        need(self.ok,'inconsistent source has no model')
        rng=random.Random(seed)
        c,v=self.root(self.const);values={c:v}
        arcs=set()
        for (a,b),i in self.index.items():
            r,u=self.root(i)
            if r not in values:
                values[r]=rng.randrange(2)
            arcs.add((a,b) if values[r]^u else (b,a))
        return arcs


def source(n,interior,guards):
    q=tuple(range(n));s=Equations(n)
    s.anti(q)
    for a,b in zip(q,q[1:]):
        s.force(a,b)
    incoming=set(interior)|{n-1}
    for v in range(1,n):
        s.force(v,0,int(v in incoming))
    for a,b in guards:
        s.force(a,b)
    return s


def old(n,ell):
    d=tuple(range(3,ell+3));u=tuple(range(ell+3,n-1))
    return [(0,2)+d+u+(n-1,1),d+(0,2)+u+(n-1,1),d+(0,1,2)+u+(n-1,)]


def one_hole(n,branch):
    h=n-5;end=tuple(range(h,n))
    if branch=='first':
        return old(n,1),[(4,2)]
    if branch=='second':
        return old(n,3),[(2,4),(6,2)]
    if branch=='single':
        v=tuple(range(5,h))
        return [(0,1,2,3)+v+end+(4,), (1,2,3)+v+(0,)+end+(4,),
                (1,2,3,4)+v+(0,)+end],[(2,4),(2,6),(3,5)]
    u=tuple(range(6,h))
    return [(0,1,2)+u+end+(4,5,3), (1,2)+u+(0,)+end+(4,5,3),
            (1,2,3,4,5)+u+(0,)+end],[(2,4),(2,6),(5,3)]


def general(n,h,forward):
    if not forward:
        return old(n,1),[(4,2)]
    u=tuple(range(4,h));end=tuple(range(h,n))
    return [(0,1,2)+u+end+(3,), (1,2)+u+(0,)+end+(3,),
            (1,2,3)+u+(0,)+end],[(2,4)]


def symbolic(n,interior,words,guards,strict=True):
    s=source(n,interior,guards)
    need(s.ok,('empty claimed source cell',n,interior,guards))
    flips=set();allpairs=set(s.edges);illegal=[];blue=[]
    for step,q in enumerate(words):
        need(sorted(q)==list(range(n)),('not a Hamilton word',n,q))
        for a,b in zip(q,q[1:]):
            if not s.forced(a,b,1^(pair(a,b) in flips)):
                illegal.append([step,a,b])
        trial=s.clone();trial.anti(q,flips);blue.append(trial.ok)
        flips ^= allpairs^pathset(q)
    target=(flips == allpairs^pathset(tuple(range(n))))
    if strict:
        need(not illegal,('unforced join',n,interior,guards,illegal))
        need(not any(blue),('possible blue move',n,interior,guards,blue))
        need(target,('wrong literal target',n,interior,guards))
    return s,{'illegal':illegal,'blue_feasible':blue,'target':target}


def literal(arcs,words):
    n=len(words[0]);original=set(arcs)
    for q in words:
        joins=set(zip(q,q[1:]))
        need(joins<=arcs,('literal illegal join',joins-arcs))
        r=dict(zip(q,reversed(q)))
        need(any((r[b],r[a]) not in arcs for a,b in arcs),'literal blue')
        arcs=joins|{(b,a) for a,b in arcs-joins}
    p=set(zip(range(n-1),range(1,n)))
    need(arcs==p|{(b,a) for a,b in original-p},'literal target mismatch')


def permutation_controls():
    count=0
    for n in range(8,121):
        for h in range(4,n-3):
            q=tuple(range(1,h))+(0,)+tuple(range(h,n))
            rev=dict(zip(q,reversed(q)));g={i:rev[n-1-i] for i in range(n)}
            seen=set();cycles=[]
            for a in range(n):
                if a in seen:
                    continue
                c=[];v=a
                while v not in seen:
                    c.append(v);seen.add(v);v=g[v]
                cycles.append(c)
            orbit=next(set(c) for c in cycles if 0 in c)
            if 2*h<=n:
                need(orbit==set(range(h)),'disjoint insertion orbit')
            elif n%2:
                need(len(cycles)==1,'odd overlapping insertion orbit')
            else:
                offset=n-h
                expected=set(range(offset+1))|set(range(offset+2,h-1,2))
                need(orbit==expected and sorted(map(len,cycles))==[n//2,n//2],
                     'even overlapping insertion orbit')
            count+=1
    return count


def main():
    counts=Counter()
    for n in range(12,61):
        interior={n-5,n-3,n-2}
        for branch in ('first','second','single','triple'):
            words,guards=one_hole(n,branch)
            s,_=symbolic(n,interior,words,guards)
            counts['first_hole_affine_cells']+=1
            for seed in (17,30030):
                literal(s.model(seed),words);counts['first_hole_literal_sources']+=1
    for k in range(4,7):
        for n in range(4*k-1,4*k+2):
            required={n-4,n-2}
            possible=set(range(5,n-1))-required
            for extra in combinations(sorted(possible),k-3):
                interior=required|set(extra);h=min(interior)
                for forward in (False,True):
                    words,guards=general(n,h,forward)
                    s,_=symbolic(n,interior,words,guards)
                    counts['protected_port_affine_cells']+=1
                    if sum(extra)%23==0:
                        literal(s.model(n+k),words);counts['protected_literal_sources']+=1
    for k in range(5,17):
        for n in range(2*k+5,2*k+15):
            interior={n-k-1}|set(range(n-k+1,n-1))
            for forward in (False,True):
                words,guards=general(n,min(interior),forward)
                s,_=symbolic(n,interior,words,guards)
                literal(s.model(k*100+n),words)
                counts['higher_one_hole_affine_and_literal_cells']+=1
    counts['permutation_cases']=permutation_controls()
    # Two real source controls for the removed insertion-predecessor mistake.
    words,guards=one_hole(11,'triple')
    s,head=symbolic(11,{6,8,9},words,guards,strict=False)
    need([2,5,0] in head['illegal'] and not head['target'],'n11 overlap hostile missing')
    rows=[sum(1<<b for a,b in s.model(0) if a==v) for v in range(11)]
    words,guards=general(15,4,True)
    s,h4=symbolic(15,{4,11,13},words,guards,strict=False)
    need([2,3,0] in h4['illegal'] and not h4['target'],'h4 overlap hostile missing')
    print(json.dumps({'status':'PASS','counts':dict(counts),
        'universes':{'first_hole':'n12..60, all four branch cells and every completion',
            'protected':'k4..6,n4k-1..4k+1,every admissible port set with min>=5 and both chord cells',
            'higher_one_hole':'k5..16,n2k+5..2k+14,both chord cells',
            'permutations':'n8..120,all4<=h<=n-4'},
        'hostiles':{'n11':head,'n11_rows':rows,'protected_h4':h4},
        'scope':'affine cells quantify all source completions; literal models are deterministic corroboration'},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
