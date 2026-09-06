#!/usr/bin/env python3
"""Exact inherited-profile graph compiler; sharp actual balanced e<=6 control."""
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import combinations
from math import gcd,prod
from pathlib import Path
import json

Q=91**6
GATES=0
TRACE=sha256()
ROOT=Path(__file__).resolve().parents[1]
PROFILE_HASH='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'

def gate(ok,label,data=None):
    global GATES
    if not ok:raise RuntimeError(f'{label}: {data!r}')
    GATES+=1
    TRACE.update((label+':'+repr(data)+'\n').encode())

def load_profiles():
    raw=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    gate(sha256(raw).hexdigest()==PROFILE_HASH,'pinned-full-profile-supplier')
    data=json.loads(raw)
    return {int(k):{(c,tuple(w)) for c,w in L['profiles']} for k,L in data['levels'].items()}

def finite_graph(P,threshold):
    # A partial complement word must occur as a submultiset of at least one
    # complete inherited word. This retains exactly a necessary extension test.
    subwords={}
    for k,rows in P.items():
        for c,w in sorted(rows):
            for j in range(k+1):
                subwords.setdefault((k,c,j),set()).update(tuple(w[i] for i in S) for S in combinations(range(k),j))
    @cache
    def valid(S):
        n=len(S)
        for k in range(1,min(6,n)+1):
            for I in combinations(range(n),k):
                c=gcd(*(S[i] for i in I))
                w=tuple(sorted(gcd(c,S[j]) for j in range(n) if j not in I))
                if w not in subwords.get((7-k,c,n-k),()):return False
        return n<7 or gcd(*S)==1
    values=tuple(sorted({c for c,w in P[6] if c>=threshold}))
    neighbors={v:{w for w in values if gcd(v,w)>=threshold} for v in values}
    current={(v,) for v in values if valid((v,))};counts=[len(current)]
    for n in range(2,8):
        candidates=set()
        for S in sorted(current):
            for x in sorted(set.union(*(neighbors[v] for v in S))):
                candidate=tuple(sorted(S+(x,)))
                if valid(candidate):candidates.add(candidate)
        current=candidates;counts.append(len(current))
    expected={7:[36,71,153,249,256,111,0],6:[37,111,381,825,1052,575,26]}
    gate(counts==expected[threshold],'exhaustive-connected-prefix-counts',(threshold,counts))
    for S in sorted(current):
        for k in range(1,7):
            for I in combinations(range(7),k):
                c=gcd(*(S[i] for i in I));w=tuple(sorted(gcd(c,S[j]) for j in range(7) if j not in I))
                gate((c,w) in P[7-k],'direct-complete-survivor-profile',(S,I))
    print('EXHAUSTIVE_PROFILE_GRAPH threshold='+str(threshold)+' layer_counts='+str(counts)+' survivors='+json.dumps(sorted(current),separators=(',',':')))
    return current

def hostile(P):
    D=(14,48,50,63,75,80,84);caps=(0,90,30,9,4,2,1,1)
    for k in range(1,8):
        for S in combinations(D,k):gate(gcd(*S)<=caps[k],'scalar-cap-hostile-all-subsets',(k,S))
    expected=((0,6,14),(1,5,16),(1,6,12),(2,4,25),(2,5,10),(3,6,21))
    edges=tuple((i,j,gcd(D[i],D[j])) for i,j in combinations(range(7),2) if gcd(D[i],D[j])>=10)
    gate(edges==expected,'scalar-cap-hostile-connected-tree')
    word=tuple(sorted(gcd(D[0],D[j]) for j in range(1,7)))
    gate(word==(1,2,2,2,7,14) and (14,word) not in P[6],'full-word-destroys-scalar-hostile')
    print('SCALAR_CAP_HOSTILE states='+str(D)+' rejected_full_word='+str((14,word)))
    D=(8,8,9,9,10,60,72)
    for k in range(1,7):
        cores={c for c,w in P[7-k]}
        for S in combinations(D,k):gate(gcd(*S) in cores,'exact-scalar-projection-hostile',(k,S))
    gate(gcd(*D)==1,'exact-scalar-projection-hostile-primitive')
    edges=[(i,j) for i,j in combinations(range(7),2) if gcd(D[i],D[j])>=7]
    seen={0}
    while True:
        new={j for i,j in edges if i in seen}|{i for i,j in edges if j in seen}
        if new<=seen:break
        seen|=new
    gate(len(seen)==7 and min(gcd(D[i],D[j]) for i,j in edges)==8,'exact-scalar-projection-hostile-connected')
    word=tuple(sorted(gcd(60,D[j]) for j in range(7) if j!=5))
    gate(word==(3,3,4,4,10,12) and (60,word) not in P[6],'exact-scalar-projection-loses-complement-word')
    print('EXACT_SCALAR_PROJECTION_HOSTILE states='+str(D)+' rejected_full_word='+str((60,word)))

def inert(n):
    if not 2<=n<=356:return False
    p=2
    while p*p<=n:
        power=0
        while n%p==0:
            n//=p;power+=1
        if power and (p%3!=2 or power>2):return False
        p+=1
    return n==1 or n%3==2

def rank(A):
    A=[list(map(F,row)) for row in A];r=0
    for j in range(len(A[0])):
        i=next((i for i in range(r,len(A)) if A[i][j]),None)
        if i is None:continue
        A[i],A[r]=A[r],A[i];v=A[r][j];A[r]=[x/v for x in A[r]]
        for i in range(r+1,len(A)):
            v=A[i][j]
            if v:A[i]=[x-v*y for x,y in zip(A[i],A[r])]
        r+=1
        if r==len(A):break
    return r

def no_crossing(A,B,Y):
    d=gcd(A,B);p,q=sorted((A//d,B//d));delta=gcd(d,Y)
    gate(q<=Q,'internal-pair-height-for-signed-box')
    if d//delta>Q:return True
    x=Y//delta;r=pow(p,-1,q)*x%q;s=(x-p*r)//q
    ceil=lambda a,b:-((-a)//b)
    return max(ceil(-Q-r,q),ceil(s-Q,p))>min((Q-r)//q,(s+Q)//p)

def actual_control(P,survivors):
    V=(1,2,3,4,5,6)
    U=(12584,14872,117,9999,98890,132990,10296)
    t=360*(1000*Q+1);g=1
    D=tuple(gcd(t,u) for u in U)
    row=tuple(t*v for v in V)+U
    gate(D==(8,8,9,9,10,30,72) and tuple(sorted(D)) in survivors,'actual-sharp-sheet-states')
    gate(gcd(*V)==gcd(*U)==gcd(*row)==1 and len(set(row))==13,'actual-primitive-distinct-labels')
    gate(t==204432930734760360 and sum(row)==4293091545430247308<Q**2,'actual-physical-box')
    gate(gcd(1000*Q+1,prod(U))==1,'scale-does-not-change-sheet-states')
    gate(not inert(125) and inert(25),'strict-atlas-exponent-hostile')
    atlas=[(p,q) for p in range(1,178) for q in range(p+1,357-p) if gcd(p,q)==1 and inert(p+q)]
    gate(len(atlas)==5855,'literal-strict-atlas-size')
    edges=[];relations=[]
    for i,j in combinations(range(13),2):
        d=gcd(row[i],row[j]);p,q=row[i]//d,row[j]//d
        if inert(p+q):
            edges.append((i,j));relation=[0]*13;relation[i],relation[j]=q,-p
            gate(max(map(abs,relation))<=355 and sum(a*b for a,b in zip(relation,row))==0,'literal-decoder-relation')
            relations.append(relation)
    unseen=set(range(13));components=[]
    while unseen:
        reached={min(unseen)}
        while True:
            extra={j for i,j in edges if i in reached}|{i for i,j in edges if j in reached}
            if extra<=reached:break
            reached|=extra
        components.append(tuple(sorted(reached)));unseen-=reached
    gate(components==[tuple(range(6)),tuple(range(6,13))] and rank(relations)==11,'actual-two-components-rank11')
    Uedges=[]
    for i,j in combinations(range(7),2):
        d=gcd(U[i],U[j]);p,q=U[i]//d,U[j]//d
        if inert(p+q):Uedges.append((i,j,p,q,gcd(t,d)))
    expected=[(0,6,11,9,8),(1,6,13,9,8),(2,6,1,88,9),(3,6,101,104,9),(4,5,29,39,10),(5,6,155,12,6)]
    gate(Uedges==expected and min(e[-1] for e in Uedges)==6,'actual-forced-edge-bound-sharp')
    gate(t>Q*sum(sorted(U)[-2:]),'both-mixed-orientations-separated-by-amplitude')
    mixed=0
    for inside,outside in ((row[:6],row[6:]),(row[6:],row[:6])):
        for A,B in combinations(inside,2):
            for Y in outside:
                gate(no_crossing(A,B,Y),'full-literal-mixed-signed-box')
                mixed+=1
    gate(mixed==231,'all-mixed-supports')
    profile_count=0;maxima=[]
    for size in range(7,13):
        maximum=0
        for I in combinations(range(13),size):
            c=gcd(*(row[i] for i in I));w=tuple(sorted(gcd(c,row[j]) for j in range(13) if j not in I))
            gate((c,w) in P[13-size],'all-actual-inherited-profile-words')
            maximum=max(maximum,c);profile_count+=1
        maxima.append(maximum)
    gate(profile_count==4095 and maxima==[72,10,9,3,2,1],'full-profile-universe-and-maxima')
    phase=F(1,7);clearance=min(min((v*phase)%1,(-v*phase)%1) for v in row)
    gate(clearance==F(1,7)>F(1,14),'literal-full-row-safe-phase')
    print('ACTUAL_SHARP_CONTROL '+json.dumps({'V':V,'U':U,'t':t,'g':g,'sheet_states':D,'sum':sum(row),'U_edges':Uedges,'decoder_rank':11,'mixed_supports':mixed,'profiles':profile_count,'max_subset_gcds_7_to_12':maxima,'phase':str(phase),'clearance':str(clearance)},sort_keys=True,separators=(',',':')))

if __name__=='__main__':
    P=load_profiles();hostile(P);finite_graph(P,7);survivors=finite_graph(P,6);actual_control(P,survivors)
    print('CLAIM: every full-profile surviving actual balanced entry has some connected seven-component atlas edge with sheet gcd e<=6; sharp in that actual domain')
    print('SCOPE: finite exact graph quotient plus all-height actual-entry reduction; no claimed strict failure or LRC14 solution')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())
