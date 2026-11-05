import hashlib
import itertools
import json
from pathlib import Path

def need(x,message):
    if not x: raise RuntimeError(message)

def labelled(n):
    pairs=list(itertools.combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        rows=[0]*n
        for i,(a,b) in enumerate(pairs):
            if bits>>i&1: rows[a]|=1<<b
            else: rows[b]|=1<<a
        yield rows

def profile(rows,r):
    size=1<<len(rows); full=size-1
    hs=[0]*size; ss=[0]*size; es=[0]*size; cycles=0
    def dfs(mask,last,first):
        nonlocal cycles
        hs[mask]+=1; ss[mask]+=first==r; es[mask]+=last==r
        if mask==full and first==r and rows[last]>>r&1: cycles+=1
        more=rows[last]&(full^mask)
        while more:
            bit=more&-more; more-=bit
            dfs(mask|bit,bit.bit_length()-1,first)
    for v in range(len(rows)): dfs(1<<v,v,v)
    hs[0]=1; deletion=full^(1<<r)
    return dict(h=hs[full],e=es[full],s=ss[full],d=hs[deletion],
        R=sum(ss[S]*hs[full^S] for S in range(size)),
        L=sum(es[S]*hs[full^S] for S in range(size)),
        b=sum(hs[S]*hs[deletion^S] for S in range(size) if not S>>r&1),C=cycles)

def join(A,a,B,b):
    # Keep A labels; append all nonroot B labels in their original order.
    mapping={b:a}
    for v in range(len(B)):
        if v!=b: mapping[v]=len(A)+len(mapping)-1
    rows=A+[0]*(len(B)-1)
    for v in range(len(B)):
        for w in range(len(B)):
            if B[v]>>w&1: rows[mapping[v]]|=1<<mapping[w]
    other=sum(1<<mapping[v] for v in range(len(B)) if v!=b)
    for v in range(len(A)):
        if v!=a: rows[v]|=other
    return rows,a

def converse(A):
    full=(1<<len(A))-1
    return [full^row^(1<<v) for v,row in enumerate(A)]

M=[[[0,0,0,3],[0,0,0,-3],[0,0,12,-3],[3,0,-12,3]],
   [[0,0,0,3],[0,0,0,-2],[0,0,12,-6],[3,-1,-9,5]],
   [[0,0,0,3],[0,0,0,-1],[0,0,12,-9],[3,-2,-6,5]],
   [[0,0,0,3],[0,0,0,0],[0,0,12,-12],[3,-3,-3,3]]]

def main():
    bank=[]; digest=hashlib.sha256(); general=family=0
    for n in range(1,4):
        for A in labelled(n):
            for a in range(n): bank.append((A,a,profile(A,a)))
    for A,a,p in bank:
        for B,b,q in bank:
            T,r=join(A,a,B,b); got=profile(T,r)
            expected=dict(h=(p['R']-p['d'])*(q['L']-q['d'])+(p['h']-p['e'])*q['d']+(q['h']-q['s'])*p['d']+p['e']*q['s'],
                R=(p['R']-p['d'])*q['b']+p['d']*q['R'],
                L=(q['L']-q['d'])*p['b']+q['d']*p['L'],
                d=p['d']*q['d'],e=(p['e'] if len(B)==1 else p['d']*q['e']),s=(q['s'] if len(A)==1 else p['s']*q['d']))
            need(all(got[k]==v for k,v in expected.items()),'general root amalgamation')
            general+=1
    triangle=[2,4,1]; constants=[]
    for j in range(4):
        K=triangle+[sum(1<<v for v in range(j))]
        for v in range(j,3): K[v]|=8
        q=profile(K,3)
        constants.append([q[k] for k in ('h','s','e','R','L','d','b')])
        need(M[j]==[list(x) for x in zip(*M[3-j])],'converse matrix')
        for A,a,p in bank:
            if p['s']!=1: continue
            for B,b,z in bank:
                if z['s']!=1: continue
                if len(A)==1 and (len(B)==1 or j!=3): continue
                if len(B)==1 and j!=0: continue
                T,r=join(A,a,K,3); T,r=join(T,r,converse(B),b)
                got=profile(T,r)
                v=[p[k] for k in ('h','e','R','d')]
                w=[z[k] for k in ('h','e','R','d')]
                scalar=sum(v[x]*M[j][x][y]*w[y] for x in range(4) for y in range(4))
                need((got['h'],got['C'],got['s'],got['e'],got['d'])==(scalar,3,3*z['d'],3*p['d'],3*p['d']*z['d']),'physical three-cycle family')
                Kopp=triangle+[sum(1<<v for v in range(3-j))]
                for v in range(3-j,3): Kopp[v]|=8
                S,sr=join(B,b,Kopp,3); S,sr=join(S,sr,converse(A),a)
                perm={a:b}
                for k,vv in enumerate(vv for vv in range(len(A)) if vv!=a): perm[vv]=len(B)+3+k
                for k in range(3): perm[len(A)+k]=len(B)+2-k
                for k,vv in enumerate(vv for vv in range(len(B)) if vv!=b): perm[len(A)+3+k]=vv
                need(len(set(perm.values()))==len(T),'converse bijection')
                need(all(((T[x]>>y)&1)==((S[perm[y]]>>perm[x])&1) for x in range(len(T)) for y in range(len(T)) if x!=y),'literal converse arcs')
                digest.update(f'{j},{A},{a},{B},{b},{scalar}\n'.encode())
                family+=1
    need(constants==[[3,0,3,3,12,3,12],[5,1,2,6,9,3,12],[5,2,1,9,6,3,12],[3,3,0,12,3,3,12]],'triangle table')
    result=dict(status='PASS',general_pairs=general,three_cycle_pairs=family,triangle_constants=constants,digest=digest.hexdigest())
    print(json.dumps(result,indent=2,sort_keys=True))
    (Path(__file__).parent/'root-amalgamation-certificate.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')

if __name__=='__main__': main()
