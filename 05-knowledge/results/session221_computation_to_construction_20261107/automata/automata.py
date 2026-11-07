#!/usr/bin/env python3
"""Exact fixed-apex tape inverse with a Busch-certified finite horizon.
No external packages. A mask bit i is 1 iff apex i beats the chain vertex.
"""
from itertools import product, permutations

def tournament(apex, word):
    k=len(apex); n=k+len(word); out=list(apex)+[0]*len(word)
    for j,mask in enumerate(word):
        v=k+j
        for i in range(k):
            if mask>>i&1: out[i]|=1<<v
            else: out[v]|=1<<i
        for u in range(k,v): out[u]|=1<<v
    return tuple(out)

def subset_dp(out):
    n=len(out)
    if not n:return 1,(),()
    ends=[[0]*n for _ in range(1<<n)]
    for v in range(n):ends[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(ends[mask]):
            if not c:continue
            rest=out[v]&~mask
            while rest:
                bit=rest&-rest;rest-=bit;u=bit.bit_length()-1
                ends[mask|bit][u]+=c
    rev=tuple(sum(1<<u for u in range(n) if out[u]>>v&1) for v in range(n))
    starts=[[0]*n for _ in range(1<<n)]
    for v in range(n):starts[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(starts[mask]):
            if not c:continue
            rest=rev[v]&~mask
            while rest:
                bit=rest&-rest;rest-=bit;u=bit.bit_length()-1
                starts[mask|bit][u]+=c
    return sum(ends[-1]),tuple(starts[-1]),tuple(ends[-1])

def strong_components(out):
    n=len(out);reach=[out[v]|1<<v for v in range(n)]
    for u in range(n):
        for v in range(n):
            if reach[v]>>u&1:reach[v]|=reach[u]
    left=(1<<n)-1;components=[]
    while left:
        v=(left&-left).bit_length()-1
        comp=tuple(u for u in range(n) if reach[v]>>u&1 and reach[u]>>v&1)
        components.append(comp)
        for u in comp:left&=~(1<<u)
    return tuple(components)

def normalize(apex,word):
    """Delete precisely singleton chain SCCs. FINAL words only, not prefixes."""
    k=len(apex);full=(1<<k)-1;n=len(word)
    suffix=[full]*(n+1)
    for j in range(n-1,-1,-1):suffix[j]=suffix[j+1]&word[j]
    union=0;kept=[];deleted=[]
    for j,mask in enumerate(word):
        apex_cut=all(not(apex[a]>>b&1) for a in range(k) if not(mask>>a&1)
                     for b in range(k) if mask>>b&1)
        singleton=not(union&~mask) and not(mask&~suffix[j+1]) and apex_cut
        (deleted if singleton else kept).append(j)
        union|=mask
    return tuple(word[j] for j in kept),tuple(deleted)

def normalize_endpoints(apex,word):
    """Also preserve every start/end label that belongs to the apex core."""
    reduced,deleted=normalize(apex,word)
    guards=set()
    if word and word[0]==0:guards.add(0)
    if word and word[-1]==(1<<len(apex))-1:guards.add(len(word)-1)
    deleted=tuple(i for i in deleted if i not in guards)
    remove=set(deleted)
    return tuple(x for i,x in enumerate(word) if i not in remove),deleted

def length_cap(m):
    """max{2a+3b : 3**a*5**b <= m}, computed without floating point."""
    if m<1:raise ValueError('The target/census upper bound must be positive')
    best=0;b=0;p5=1
    while p5<=m:
        a=0;p=p5
        while p*3<=m:a+=1;p*=3
        best=max(best,2*a+3*b);b+=1;p5*=5
    return best

class Transfer:
    """Inherited positive apex-permutation/gap-state transfer, reimplemented."""
    def __init__(self,apex):
        self.apex=tuple(apex);self.k=k=len(apex)
        states=list(product(range(3),repeat=k))
        states=[s+(last,) for s in states for last in (0,2)]
        perms=list(permutations(range(k)))
        ix={s:i for i,s in enumerate(states)};d=len(states)
        self.dimension=len(perms)*d
        init=[0]*self.dimension
        self.edges=[[] for _ in range(1<<k)]
        self.accept=[];self.start=[[] for _ in range(k)];self.end=[[] for _ in range(k)]
        self.joint=[[[] for _ in range(k)] for _ in range(k)]
        for pi,p in enumerate(perms):
            init[pi*d+ix[(0,)*(k+1)]]=1
            for si,s in enumerate(states):
                source=pi*d+si
                ok=all(x==2 or (x==0 and (g==0 or g==k or apex[p[g-1]]>>p[g]&1)) for g,x in enumerate(s))
                if ok:
                    self.accept.append(source)
                    if k and s[0]==0:self.start[p[0]].append(source)
                    if k and s[k]==0:self.end[p[-1]].append(source)
                    if k and s[0]==0 and s[k]==0:self.joint[p[0]][p[-1]].append(source)
                for mask in range(1<<k):
                    targets=[]
                    for g in range(k+1):
                        if s[g]==0 and g>0 and not(mask>>p[g-1]&1):continue
                        t=list(s);t[g]=2 if g==k or not(mask>>p[g]&1) else 1
                        targets.append(pi*d+ix[tuple(t)])
                    if targets:self.edges[mask].append((source,tuple(targets)))
        self.initial=tuple(init)
    def step(self,state,mask):
        nxt=[0]*self.dimension
        for source,targets in self.edges[mask]:
            c=state[source]
            if c:
                for dest in targets:nxt[dest]+=c
        return tuple(nxt)
    def output(self,state):return sum(state[i] for i in self.accept)
    def endpoints(self,state):
        return (tuple(sum(state[i] for i in inds) for inds in self.start),
                tuple(sum(state[i] for i in inds) for inds in self.end))
    def joint_endpoints(self,state):
        return tuple(tuple(sum(state[i] for i in inds) for inds in row) for row in self.joint)
    def evaluate(self,word):
        state=self.initial
        for mask in word:state=self.step(state,mask)
        return self.output(state),self.endpoints(state)

def census(apex,upper):
    """Complete all-length census <= upper, with literal shortest witnesses.
    All prefixes remain unnormalized. Depth cap is the proved length_cap.
    """
    tf=Transfer(apex);limit=length_cap(upper);witness={};visited=pruned=0
    def dfs(state,word):
        nonlocal visited,pruned
        visited+=1;h=tf.output(state)
        if h>upper:pruned+=1;return
        old=witness.get(h)
        if old is None or (len(word),word)<(len(old),old):witness[h]=word
        if len(word)==limit:return
        for mask in range(1<<len(apex)):dfs(tf.step(state,mask),word+(mask,))
    dfs(tf.initial,())
    return {'upper':upper,'tape_cap':limit,'visited':visited,'pruned':pruned,
            'witnesses':{str(h):list(w) for h,w in sorted(witness.items())}}

def compile_target(apex,target,wanted=None):
    """Decide full fixed-core membership and return a bounded physical witness.
    Search never rewrites a prefix. DFS uses O_k(log^2 target) working bits.
    """
    k=len(apex);limit=length_cap(target)+(2 if wanted else 0)
    if wanted:
        parsed={}
        for key,value in wanted.items():
            if key not in ('start','end','joint'):raise ValueError('Unknown endpoint coordinate: '+str(key))
            if key=='joint':
                value=tuple(tuple(row) for row in value)
                if len(value)!=k or any(len(row)!=k for row in value):raise ValueError('joint must be k by k')
                flat=tuple(x for row in value for x in row)
            else:
                value=tuple(value)
                if len(value)!=k:raise ValueError(key+' must have k entries')
                flat=value
            if any(not isinstance(x,int) for x in flat):raise ValueError('Endpoint counts must be integers')
            if any(x<0 or x>target for x in flat):
                return {'target':target,'tape_cap':limit,'visited':0,'pruned':0,
                        'status':'ABSENT','reason':'endpoint count outside [0,target]'}
            parsed[key]=value
        wanted=parsed
    tf=Transfer(apex);visited=pruned=0
    def accepted(state):
        if not wanted:return True
        start,end=tf.endpoints(state)
        actual={'start':start,'end':end,'joint':tf.joint_endpoints(state)}
        return all(actual[key]==value for key,value in wanted.items())
    def dfs(state,word):
        nonlocal visited,pruned
        visited+=1;h=tf.output(state)
        if h==target and accepted(state):return word
        if h>target:pruned+=1;return None
        if len(word)==limit:return None
        for mask in range(1<<len(apex)):
            answer=dfs(tf.step(state,mask),word+(mask,))
            if answer is not None:return answer
        return None
    word=dfs(tf.initial,())
    result={'target':target,'tape_cap':limit,'visited':visited,'pruned':pruned,
            'status':'FOUND' if word is not None else 'ABSENT'}
    if word is not None:
        result.update(word=list(word),outgoing=list(tournament(apex,word)))
        if wanted:result['endpoint_constraints']={key:value for key,value in wanted.items()}
    return result

if __name__=='__main__':
    import sys,json
    for target in map(int,sys.argv[1:] or [7,21,23,63]):
        print(json.dumps(compile_target((2,0),target),sort_keys=True))
