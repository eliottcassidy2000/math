"""Independent literal deleted-card and family-completion controls.

All labelled quotient tournaments through order five. For quotients through
order three, exhaust every deletion and every incident replacement cut.
"""
from itertools import combinations
from hashlib import sha256
import json

def require(test,message):
    if not test:
        raise RuntimeError(message)

def blowup(k,code):
    q=[0]*k
    for i,j in combinations(range(k),2):
        if code&1:q[i]|=1<<j
        else:q[j]|=1<<i
        code>>=1
    a=[0]*(3*k)
    for b in range(k):
        for i in range(3):
            a[3*b+i] |= 1<<(3*b+(i+1)%3)
            for c in range(k):
                if q[b]>>c&1:a[3*b+i]|=7<<(3*c)
    return a

def homogeneous(a,vertices):
    full=(1<<len(a))-1
    mask=sum(1<<v for v in vertices)
    outside=full^mask
    return len({a[v]&outside for v in vertices})==1

def triangles(a):
    answer=[]
    for vs in combinations(range(len(a)),3):
        mask=sum(1<<v for v in vs)
        if all((a[v]&mask).bit_count()==1 for v in vs) and homogeneous(a,vs):
            answer.append(mask)
    return answer

def in_family(a):
    modules=triangles(a)
    def cover(left):
        if not left:return True
        v=left&-left
        return any(cover(left^m) for m in modules if m&v and m&left==m)
    return cover((1<<len(a))-1)

def run():
    cases=cards=cuts=0
    digest=sha256()
    for k in range(1,6):
        n=3*k
        for code in range(1<<(k*(k-1)//2)):
            a=blowup(k,code)
            require(triangles(a)==[7<<(3*i) for i in range(k)],'canonical triangle modules')
            cases+=1
            for deleted in range(n):
                vs=[i for i in range(n) if i!=deleted]
                b=[sum(((a[i]>>j)&1)<<t for t,j in enumerate(vs)) for i in vs]
                pairs=[pair for pair in combinations(range(n-1),2) if homogeneous(b,pair)]
                require(len(pairs)==1,'unique homogeneous pair')
                u,w=pairs[0]
                if not b[u]>>w&1:u,w=w,u
                recovered=(1<<u)|(b[u]&~(1<<w))
                original=sum(((a[deleted]>>j)&1)<<t for t,j in enumerate(vs))
                require(recovered==original,'literal canonical cut reconstruction')
                cards+=1
                digest.update(str((k,code,deleted,recovered)).encode())
                if k<=3:
                    family_cuts=[]
                    for cut in range(1<<(n-1)):
                        child=[b[i]|((not(cut>>i&1))<<(n-1)) for i in range(n-1)]+[cut]
                        if in_family(child):family_cuts.append(cut)
                        cuts+=1
                    require(family_cuts==[original],'unique completion into entire substitution family')
    print(json.dumps({'scope':'all labelled quotients k1..5, including nonstrong;all cuts for each deletion at k1..3',
                      'quotients':cases,'deleted_cards':cards,'literal_completion_cuts':cuts,
                      'semantic_sha256':digest.hexdigest(),
                      'verdict':'PASS unique homogeneous pair, intrinsic original cut, canonical triangle modules, family-wide completion uniqueness'},sort_keys=True,indent=2))

if __name__=='__main__':run()
