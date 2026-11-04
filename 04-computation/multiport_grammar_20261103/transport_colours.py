#!/usr/bin/env python3
"""Exact coloured endpoint-forest operators and independent small graph checks.
No producer imports, external data, native-complement claims or assertion gates.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
import hashlib
import json

gates=0
counts=Counter()
def need(name,value):
    global gates
    if not value: raise RuntimeError(name)
    gates+=1
def same(name,a,b):
    need(name,a==b)
def cleaned(c):return Counter({k:v for k,v in c.items() if v})
def vertex_token(v,assignment):return v if assignment[v]<0 else -assignment[v]-1

def canonical(exponents,anchors,q):
    ex=list(exponents); an=[]
    for a,b in anchors:
        if a<0 and b<0:ex[(-a-1)*q+(-b-1)]+=1
        else:an.append((a,b))
    return (tuple(ex),tuple(sorted(an)))

@lru_cache(None)
def covers(rows):
    """Independent actual forests: component of least remaining vertex first.
    Each path is counted by literal permutations, not by transport operators.
    Returned pairs are physical endpoints; distinct internal paths retain weight.
    """
    n=len(rows);paths={}
    for mask in range(1,1<<n):
        vs=[v for v in range(n) if mask>>v&1];p=Counter()
        for word in permutations(vs):
            if all(rows[a]>>b&1 for a,b in zip(word,word[1:])):
                p[(word[0],word[-1])]+=1
        paths[mask]=p
    @lru_cache(None)
    def rec(mask):
        if not mask:return (((),1),)
        first=mask&-mask;ans=Counter();sub=mask
        while sub:
            if sub&first:
                for edge,k in paths[sub].items():
                    for rest,m in rec(mask^sub):
                        ans[tuple(sorted((edge,)+rest))]+=k*m
            sub=(sub-1)&mask
        return tuple(ans.items())
    return rec((1<<n)-1)

def direct(rows,assignment,q):
    ans=Counter()
    for pairs,mult in covers(tuple(rows)):
        translated=[(vertex_token(a,assignment),vertex_token(b,assignment)) for a,b in pairs]
        ans[canonical((0,)*(q*q),translated,q)]+=mult
    return ans

def birth(state,q,new,in_active,in_colours):
    """Actual positive operator. Anonymous repeated selections use falling counts."""
    def enters(tok):return tok in in_active if tok>=0 else -tok-1 in in_colours
    ans=Counter()
    for (ex,anch),weight in state.items():
        options=[('a',i,pair) for i,pair in enumerate(anch)]
        options += [('z',i,(-(i//q)-1,-(i%q)-1)) for i,k in enumerate(ex) if k]
        left=[None]+[o for o in options if enters(o[2][1])]
        right=[None]+[o for o in options if not enters(o[2][0])]
        for aa in left:
            for bb in right:
                if aa is not None and bb is not None and aa[0]=='a' and bb[0]=='a' and aa[1]==bb[1]:continue
                nex=list(ex);removed=set();mult=weight
                for opt in (aa,bb):
                    if opt is None:continue
                    kind,i,_=opt
                    if kind=='a':removed.add(i)
                    else:
                        mult*=nex[i]
                        if not mult:break
                        nex[i]-=1
                if not mult:continue
                start=aa[2][0] if aa is not None else new
                end=bb[2][1] if bb is not None else new
                pairs=[pair for i,pair in enumerate(anch) if i not in removed]+[(start,end)]
                ans[canonical(nex,pairs,q)]+=mult
    return cleaned(ans)

def join(state,q,a,b):
    ans=Counter()
    for (ex,anch),weight in state.items():
        left=next((i for i,pair in enumerate(anch) if pair[1]==a),None)
        right=next((i for i,pair in enumerate(anch) if pair[0]==b),None)
        if left is None or right is None or left==right:continue
        pairs=[pair for i,pair in enumerate(anch) if i not in (left,right)]
        pairs.append((anch[left][0],anch[right][1]))
        ans[canonical(ex,pairs,q)]+=weight
    return ans

def flip(state,q,a,b):
    ans=Counter(state)
    ans.update(join(state,q,b,a))
    ans.subtract(join(state,q,a,b))
    return cleaned(ans)

def forget(state,q,name,colour):
    ans=Counter()
    for (ex,anch),weight in state.items():
        pairs=[tuple(-colour-1 if v==name else v for v in pair) for pair in anch]
        ans[canonical(ex,pairs,q)]+=weight
    return ans

def colour_map(state,q,newq,phi):
    ans=Counter()
    for (ex,anch),weight in state.items():
        nex=[0]*(newq*newq)
        for i,k in enumerate(ex):nex[phi[i//q]*newq+phi[i%q]]+=k
        pairs=[tuple(v if v>=0 else -phi[-v-1]-1 for v in pair) for pair in anch]
        ans[canonical(nex,pairs,newq)]+=weight
    return ans

def converse(state,q):
    ans=Counter()
    for (ex,anch),weight in state.items():
        nex=[0]*(q*q)
        for a in range(q):
            for b in range(q):nex[b*q+a]=ex[a*q+b]
        ans[canonical(nex,[(b,a) for a,b in anch],q)]+=weight
    return ans

def graph_birth(rows,incoming):
    n=len(rows);out=list(rows)+[0]
    for v in range(n):
        if v in incoming:out[v]|=1<<n
        else:out[n]|=1<<v
    return tuple(out)
def graph_flip(rows,a,b):
    out=list(rows);out[a]^=1<<b;out[b]^=1<<a;return tuple(out)
def graph_converse(rows):
    full=(1<<len(rows))-1
    return tuple(full^row^(1<<i) for i,row in enumerate(rows))

@lru_cache(None)
def reach(rows):
    n=len(rows);r=[rows[i]|1<<i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if r[i]>>k&1:r[i]|=r[k]
    return tuple(r)
def strong(rows):return all(v==(1<<len(rows))-1 for v in reach(rows))
def birth_guard(rows,incoming):
    if not rows:return True
    rr=reach(rows);full=(1<<len(rows))-1
    first={i for i,v in enumerate(rr) if v==full}
    last={i for i in range(len(rows)) if all(v>>i&1 for v in rr)}
    return bool(last.intersection(incoming) and first.difference(incoming))

@lru_cache(None)
def endpoint_dp(rows):
    n=len(rows);size=1<<n
    def endings(a):
        dp=[[0]*n for _ in range(size)]
        for v in range(n):dp[1<<v][v]=1
        for mask in range(1,size):
            for v in range(n):
                if mask>>v&1:
                    old=mask^(1<<v)
                    if old:dp[mask][v]=sum(dp[old][w] for w in range(n) if old>>w&1 and a[w]>>v&1)
        return dp
    end=endings(rows);start=endings(graph_converse(rows));h=[sum(a) for a in end];h[0]=1
    full=size-1;f=[0]*n;joint=[[0]*n for _ in range(n)]
    for mask in range(1,size):
        rest=full^mask
        for a in range(n):
            f[a]+=end[mask][a]*h[rest]
            if rest:
                for b in range(n):joint[a][b]+=end[mask][a]*start[rest][b]
    p2=sum(h[mask]*h[full^mask] for mask in range(size))
    return h[full],tuple(f),tuple(tuple(a) for a in joint),p2

def observe(state,n,active):
    h=1 if n==0 else 0;f={v:0 for v in active};j={(a,b):0 for a in active for b in active};p2=1 if n==0 else 0
    for (ex,anch),weight in state.items():
        degree=sum(ex)+len(anch)
        if degree==1:h+=weight
        if degree not in (1,2):continue
        p2+=2*weight
        ends={b:i for i,(_,b) in enumerate(anch) if b>=0}
        starts={a:i for i,(a,_) in enumerate(anch) if a>=0}
        for a in ends:f[a]+=weight
        if degree==2:
            for a,i in ends.items():
                for b,k in starts.items():
                    if i!=k:j[a,b]+=weight
    return h,f,j,p2

def tournaments(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for word in range(1<<len(pairs)):
        rows=[0]*n
        for k,(i,j) in enumerate(pairs):
            if word>>k&1:rows[i]|=1<<j
            else:rows[j]|=1<<i
        yield tuple(rows)

# Entire named finite universe, no isomorphism filter or source imports.
for n in range(5):
    for rows in tournaments(n):
        counts['graphs']+=1;expected_h,expected_f,expected_j,expected_p2=endpoint_dp(rows)
        for assignment in product((-1,0,1),repeat=n):
            active=tuple(v for v in range(n) if assignment[v]<0)
            state=direct(rows,assignment,2);counts['states']+=1
            h,f,j,p2=observe(state,n,active)
            same('H observer',h,expected_h);same('P2 observer',p2,expected_p2)
            for a in active:
                same('F observer',f[a],expected_f[a])
                same('positive defect',f[a]-h,sum(expected_j[a][v] for v in range(n) if rows[v]>>a&1))
                for b in active:same('named joint observer',j[a,b],expected_j[a][b])
            for abit in range(1<<len(active)):
                ia={v for k,v in enumerate(active) if abit>>k&1}
                for cbit in range(4):
                    ic={c for c in range(2) if cbit>>c&1}
                    incoming=ia|{v for v in range(n) if assignment[v]>=0 and assignment[v] in ic}
                    target=graph_birth(rows,incoming)
                    same('coloured birth',birth(state,2,n,ia,ic),direct(target,assignment+(-1,),2))
                    same('actual SCC guard',birth_guard(rows,incoming),strong(target))
                    counts['births']+=1
            for a in active:
                for colour in range(2):
                    other=list(assignment);other[a]=colour
                    same('forget',forget(state,2,a,colour),direct(rows,tuple(other),2));counts['forgets']+=1
                for b in active:
                    if rows[a]>>b&1:
                        result=flip(state,2,a,b)
                        same('named flip',result,direct(graph_flip(rows,a,b),assignment,2))
                        need('actual flip positivity',all(v>=0 for v in result.values()))
                        same('flip inverse',flip(result,2,b,a),state);counts['flips']+=1
            same('converse',converse(state,2),direct(graph_converse(rows),assignment,2));counts['converses']+=1
            merged=tuple(-1 if c<0 else 0 for c in assignment)
            same('merge',colour_map(state,2,1,(0,0)),direct(rows,merged,1));counts['merges']+=1
            swapped=tuple(-1 if c<0 else 1-c for c in assignment)
            same('colour permutation',colour_map(state,2,2,(1,0)),direct(rows,swapped,2))
            if n<=2:same('add unused colour',colour_map(state,2,3,(0,1)),direct(rows,assignment,3))
            # Distinct named joins commute, including same-component rejection.
            if len(active)>=3 and n<=3:
                a,b,c=active[:3]
                same('join commutation',join(join(state,2,a,b),2,b,c),join(join(state,2,b,c),2,a,b))

# Explicit two-colour distinction and the invalid split of its merged state.
redblue=(0,1);trans=(2,0);reverse=(0,1)
s1=direct(trans,redblue,2);s2=direct(reverse,redblue,2)
need('coloured hostile genuinely differs',s1!=s2)
same('merged hostile aliases',colour_map(s1,2,1,(0,0)),colour_map(s2,2,1,(0,0)))
b1=birth(s1,2,2,set(),{0});b2=birth(s2,2,2,set(),{0})
same('two colour H1',observe(b1,3,(2,))[0],1);same('two colour H3',observe(b2,3,(2,))[0],3)
# A single (blue,red) component cannot be selected twice; two can.
one=Counter({((0,0,1,0),()):1});two=Counter({((0,0,2,0),()):1})
one_after=birth(one,2,0,set(),{0});two_after=birth(two,2,0,set(),{0})
need('same anonymous component not selected twice',not any(not an for ex,an in one_after))
same('two distinct identical types falling multiplicity',two_after[((0,0,1,0),())],2)
# Paired-colour hostile realized as covers of transitive order a,c,b,d.
paired1=Counter({((0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0),()):1})
ex=[0]*16;ex[0*4+3]=1;ex[2*4+1]=1;paired2=Counter({(tuple(ex),()):1})
pa=birth(paired1,4,0,set(),{0,1});pb=birth(paired2,4,0,set(),{0,1})
same('paired first anonymous join count',sum(v for (ex,an),v in pa.items() if not an),1)
same('paired second cycle rejected',sum(v for (ex,an),v in pb.items() if not an),0)

# An actual indefinitely iterable strong two-colour construction; no digit claim.
rows=(2,4,1);assignment=(-1,0,1);state=direct(rows,assignment,2);family=[]
for step in range(17):
    h,f,j,p2=observe(state,len(rows),(0,))
    same('family H all-length formula',h,2**(step+1)+1)
    same('family F all-length formula',f[0],2**(step+2)-1)
    same('family changing defect',f[0]-h,2**(step+1)-2)
    need('family literal strongness',strong(rows))
    if len(rows)<=6:same('family full-profile independent',state,direct(rows,assignment,2))
    if len(rows)<=12:
        hh,ff,jj,pp=endpoint_dp(rows)
        same('family independent H',h,hh);same('family independent F',f[0],ff[0]);same('family independent P2',p2,pp)
    family.append([step,h,f[0],f[0]-h,len(state)])
    n=len(rows);incoming={v for v in range(n) if assignment[v]==0}
    state=birth(state,2,n,set(),{0});rows=graph_birth(rows,incoming)
    colour=step%2;state=forget(state,2,n,colour);assignment=assignment+(colour,)

# Recurrence is a formal coordinate count, not a reachable-rank claim.
for q in (1,2,3):
    vals=[1,2+2*q]
    for p in range(2,7):vals.append((2+2*q)*vals[-1]+2*(p-1)*vals[-2])
    if q==1:same('inherited signature counts',vals[:7],[1,4,18,88,460,2544,14776])
record={'universe':dict(sorted(counts.items())),'family':family,'gates':gates,
        'hostiles':{'one_bulk_outputs':[1,3],'same_type_falling_factor':2,'paired_join':[1,0]}}
encoded=json.dumps(record,sort_keys=True,separators=(',',':')).encode()
print('coloured endpoint-forest transport: PASS')
print('exact gates:',gates)
print('finite universe:',json.dumps(dict(sorted(counts.items())),sort_keys=True))
print('long word last [births,H,F0,b0,profile terms]:',family[-1])
print('two-colour / split / same-component / paired-colour hostiles: PASS')
print('native guard and two-branch closure claims: NONE')
print('semantic sha256:',hashlib.sha256(encoded).hexdigest())
