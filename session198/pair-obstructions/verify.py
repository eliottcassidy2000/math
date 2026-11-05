"""Independent hostile controls: literal gadget counts, exact source bank, word counts."""
from pathlib import Path
from functools import lru_cache
from itertools import permutations
from math import gcd
import hashlib,json
import source_oracle
import c3
ROOT=Path(__file__).resolve().parent

def require(test, message):
    if not test: raise RuntimeError(message)

@lru_cache(None)
def hwords(T, vertices):
    return tuple(p for p in permutations(vertices)
                 if all((T[a]>>b)&1 for a,b in zip(p,p[1:])))

def count(T):
    T=tuple(T); ws=hwords(T,tuple(range(len(T))))
    C=sum(p[0]==0 and bool(T[p[-1]]&1) for p in ws) if len(T)>=3 else 0
    return len(ws),C

def profile(T,c):
    T=tuple(T);n=len(T);full=(1<<n)-1
    ws=hwords(T,tuple(range(n)))
    H=len(ws); E=sum(p[-1]==c for p in ws);S=sum(p[0]==c for p in ws)
    d=len(hwords(T,tuple(v for v in range(n) if v!=c)))
    R=L=0
    for mask in range(1<<n):
        if not(mask>>c)&1:continue
        q=hwords(T,tuple(v for v in range(n)if mask>>v&1))
        p=hwords(T,tuple(v for v in range(n)if not(mask>>v&1)))
        R+=len(p)*sum(word[0]==c for word in q)
        L+=len(p)*sum(word[-1]==c for word in q)
    return (H,E,R,d),S,L

@lru_cache(None)
def gadget(L):
    require(L>=1,'positive gadget index')
    m=(L-1).bit_length();T=[0]*(m+1)
    for a in range(1,m+1):
        for b in range(a+1,m+1):T[a]|=1<<b
        if not((L-1)>>(a-1))&1:T[0]|=1<<a
        else:T[a]|=1
    return tuple(T)

@lru_cache(None)
def constants(L):
    p,s,end_cover=profile(gadget(L),0)
    require(end_cover==L and p[1]==1 and p[3]==1,'literal gadget compiler')
    require(s==L%2,'literal gadget starting parity')
    return p[0],p[2],1<<(len(gadget(L))-1)

def independent_bank(cap):
    o=(1,1,1,1);out={o:()};queue=[o]
    for D,E,U,d in queue:
        choices=[(('O',None),(D,0,2*U,d))]
        for L in range(1,cap//U+1):
            g,Q,B=constants(L)
            choices.append((('M',L),(g*D-(L%2)*E+L*U,D-(E if L==1 else 0),Q*D+B*U,D)))
        for move,nxt in choices:
            if nxt[2]<=cap and nxt not in out:
                out[nxt]=out[(D,E,U,d)]+(move,);queue.append(nxt)
    return out

def glue(V,c,G,entries=None):
    n=len(V);out=list(V)+[r<<n for r in G]
    entries={c}if entries is None else set(entries)
    for a in range(n):
        for b in range(len(G)):
            if b==0 and a in entries:out[n]|=1<<a
            else:out[a]|=1<<(n+b)
    return out,n

def build(recipe):
    T=[0];c=0
    for kind,L in recipe:
        if kind=='O':
            n=len(T);T=[row|(1<<n) for row in T]+[0]
        else:T,c=glue(T,c,gadget(L))
    return T,c

def form(kind,v,w):
    D,E,U,d=v;g,s,L,f=w
    if kind=='N':h,Ea,R,K=D-E+U,D-E+d,D+3*U-d,U+E
    else:h,Ea,R,K=D+2*U,U+d,2*D-E+6*U-d,U+E
    return h*g-Ea*s+R*L-K*f

def content(v):
    g=0
    for x in v:g=gcd(g,x)
    return g

def cone(v):
    D,E,U,d=v
    return 1<=d<=D<=U and 0<=E<=d and D%2==1 and d%2==1 and (U-E)%2==0 and D-E>=d and U>=D+d-E and U>=2*d

def main():
    b=independent_bank(22)
    require(set(b)==set(source_oracle.source_bank(22)),'independent exact source-bank comparison')
    for v,recipe in b.items():
        T,c=build(recipe);actual,s,L=profile(T,c)
        require(s==1 and actual==v,'literal source profile')
        require(v[0]%2 and v[3]%2 and (v[2]-v[1])%2==0,'source parity')
    small=independent_bank(11);nontrivial=[v for v in small if v!=(1,1,1,1)]
    c2_values={15};tests=0
    for kind in ('N','M'):
        for i,v in enumerate(nontrivial):
            for w in nontrivial[i:]:
                z=form(kind,v,w);tests+=1
                require(z==source_oracle.value(kind,v,w),'independent C2 formula')
                if z<=67:c2_values.add(z)
    require(67 not in c2_values and len(small)==15,'complete C2 exclusion67')
    v=(3,1,3,1);fake=(1,0,6,1)
    require(v in b and fake not in b and cone(v)and cone(fake),'parity-respecting cone hostile')
    require(form('N',v,fake)==67,'formal67')
    # Independently assembled B template. Both sources are literal transitive graphs.
    V=[6,4,0];W=[2,0]
    T,a=glue(V,0,[0]);T,r=glue(T,0,[0],{a,0})
    # Converse W is [0,1], with the root still vertex0.
    T,root=glue(T,0,[0,1],{r,0})
    require(count(T)==(67,3),'actual C3 witness by permutation count')
    expected=c3.construct('B',V,0,W,0)[0]
    require(T==expected,'independent versus inherited B construction')
    require(c3.value('B',(1,0,4,1),(1,0,2,1))==67,'C3 scalar')
    primitive=(3,1,3,1);trap=(9,3,9,3);w=(1,0,2,1)
    require(trap in b and content(primitive)==1 and content(trap)==3,'gcd may increase')
    require(form('N',trap,w)==69 and 209%3==2 and 211%3==1,'ternary-child content trap')
    divisibility=0
    for v in b:
        D,E,U,d=v;g0=content(v)
        successors=[(D,0,2*U,d)]
        for L in range(1,33):
            g,Q,B=constants(L)
            successors.append((g*D-(L%2)*E+L*U,D-(E if L==1 else 0),Q*D+B*U,D))
        for u in successors:
            require(all(x%g0==0 for x in u),'persisting source divisor');divisibility+=1
    a=(9,0,18,1);b2=(9,0,18,3)
    require(a in b and b2 in b and form('N',a,w)==133 and form('N',b2,w)==129,'deletion-sidecar collision')
    require(content(a)==1 and content(b2)==3 and gcd(gcd(a[0],a[1]),a[2])==9,'false three-coordinate content')
    # Congruence hostiles for many combined moduli, each directly between powers of2.
    for m in range(1,101):
        k=(2*m).bit_length();actual=1<<k;fakeU=actual+2*m
        require(actual<fakeU<2*actual and fakeU%m==actual%m and cone((1,0,fakeU,1)),'finite-modulus hostile')
    out=dict(status='PASS',source_cap=22,source_states=len(b),complete_c2_cutoff=11,c2_source_states=len(small),c2_pair_tests=tests,c2_support_through67=sorted(c2_values),formal67=dict(left=[3,1,3,1],right=[1,0,6,1],value=67),actual67=dict(rows=T,H=67,C=3),content_trap=dict(left=list(trap),right=list(w),H=69,excluded_children=[209,211]),deletion_collision=dict(a=list(a),b=list(b2),outputs=[133,129]),source_divisibility_controls=divisibility,source_files={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in ROOT.glob('*.py')})
    (ROOT/'VALIDATION.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
