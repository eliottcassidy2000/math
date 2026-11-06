#!/usr/bin/env python3
"""Exact controls for bounded-offset physical normal forms; no tail census."""
from collections import deque, Counter
from functools import lru_cache
import hashlib, json
from pathlib import Path

T=(1,0,2,1)

def require(ok, msg):
    if not ok: raise RuntimeError(msg)

def step(s,ch):
    D,E,U,d=s
    if ch=='O':return D,0,2*U,d
    if ch=='Z':return D-E+U,d,2*U-d,d
    if ch=='J':return D-E+U,D-E,D+U,D
    raise ValueError(ch)

@lru_cache(None)
def profile(w):
    return T if not w else step(profile(w[:-1]),w[-1])

def cone(s):
    D,E,U,d=s
    return D-E-d,U-D+E-d,E,d-E

def cone_step(y,ch):
    a,f,e,p=y
    if ch=='O':return a+e,a+2*f+e+2*p,0,e+p
    if ch=='Z':return 2*a+f+e+p,f,e+p,0
    if ch=='J':return f+p,0,a+e+p,e
    raise ValueError(ch)

def F(s):D,E,U,d=s;return D-E+U

def sides(s):
    D,E,U,d=s
    return D-2*E,U-D+E-2*d,U-D+E

def strip(w,ch):
    i=0
    while w.endswith(ch):w=w[:-1];i+=1
    return w,i

def d_normal(w,C):
    require(profile(w)[0]<=C,'D normalization domain')
    b,i=strip(w,'O')
    require(profile(b)[2]<=4*C,'D core cutoff')
    return b,i

def c_normal(w,C):
    require(sides(profile(w))[2]<=C,'c normalization domain')
    v,k=strip(w,'Z')
    if not v or v[-1]=='O':
        require(profile(v)[2]<=4*C,'c direct core cutoff')
        return v,'Z',k
    require(v[-1]=='J','c terminal class')
    b,i=d_normal(v[:-1],C)
    return b,'OJ',i+k

def normalize(w,C,which):
    s=profile(w);ell,h,c=sides(s)
    require((ell<=C if which=='ell' else abs(h)<=C),'normalization domain')
    if which=='ell':
        if not w:return '', 'OZ',0,0
        if w[-1]=='O':
            b,i=d_normal(w[:-1],C);return b,'OZ',i+1,0
        if w[-1]=='Z':return w,'OZ',0,0
        b,typ,i=c_normal(w[:-1],C)
        return b,'ZJ' if typ=='Z' else 'OJJ',i,0
    v,k=strip(w,'Z')
    if h<0:
        if not v:return '', 'OZ',0,k
        require(v[-1]=='J','negative h terminal class')
        b,i=d_normal(v[:-1],C)
        return b,'OJ',i+k,0
    require(h>0 and v[-1]=='O','positive h terminal class')
    q=v[:-1]
    if profile(q)[1]==0:
        b,i=d_normal(q,C);return b,'OZ',i+1,k
    if q[-1]=='Z':return q,'OZ',1,k
    require(q[-1]=='J','positive h precursor class')
    b,typ,i=c_normal(q[:-1],C)
    return b,'ZJOZ' if typ=='Z' else 'OJJOZ',i,k

def suffix(typ,i,j):
    return {'OZ':'O'*i+'Z'*j,
            'ZJOZ':'Z'*i+'JO'+'Z'*j,
            'OJJOZ':'O'*i+'JJO'+'Z'*j,
            'OJ':'O'*i+'J','ZJ':'Z'*i+'J','OJJ':'O'*i+'JJ'}[typ]

def inventory(Umax):
    out={T:''}; todo=deque([T])
    while todo:
        s=todo.popleft(); w=out[s]
        for ch in 'OZJ':
            z=step(s,ch)
            if z[2]<=Umax and z not in out:
                out[z]=w+ch;todo.append(z)
    return out

def graph(w):
    out=[2,0];root=0
    for ch in w:
        n=len(out);out.append(0)
        for i in range(n):
            if ch!='O' and i==root:out[n]|=1<<i
            else:out[i]|=1<<n
        if ch=='J':root=n
    return out

def count_graph(out):
    n=len(out);end=(1<<n)-1
    dp=[{} for _ in range(1<<n)]
    cyc=[{} for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    cyc[1][0]=1
    for mask in range(1,1<<n):
        for table in (dp,cyc):
            for i,v in table[mask].items():
                avail=out[i]&~mask
                while avail:
                    bit=avail&-avail;avail-=bit;j=bit.bit_length()-1
                    nxt=table[mask|bit];nxt[j]=nxt.get(j,0)+v
    return sum(dp[end].values()),sum(v for i,v in cyc[end].items() if out[i]&1)

def length_bound(N):
    L=0
    while 3**(L+2)<=N*2**(L+1):L+=1
    return L

def fib(n):
    a,b=0,1
    for _ in range(n):a,b=b,a+b
    return a

def main():
    Umax=512
    inv=inventory(Umax)
    hist=Counter(); controls=0; samples={}
    for s,w in inv.items():
        a,f,e,p=cone(s);ell,h,c=sides(s)
        require(min(a,f,e,p)>=0 and e+p>=1,'physical cone')
        require(ell>0 and ell%2==1 and h%2==1 and c>=1,'sidecar positivity/parity')
        for ch in 'OZJ':
            z=step(s,ch)
            require(cone(z)==cone_step((a,f,e,p),ch),'independent coordinate transfer')
            require(2*F(z)>=3*F(s),'readout growth')
            require(z[2]>s[2],'strict cover growth')
        require(step(step(s,'J'),'Z')==step(step(s,'O'),'J'),'JZ=OJ')
        require(F(step(s,'O'))==F(step(s,'J'))==2*F(s)-ell,'O/J offset')
        require(F(step(s,'Z'))==2*F(s)+h,'Z offset')
        require(sides(step(s,'Z'))[1:]==(h,c),'Z invariants')
        require(sides(step(s,'J'))[1:]==(-s[0],s[0]),'J sidecars')
        require(sides(step(s,'O'))[1:]==(a+2*f+p,2*s[2]-s[0]),'O sidecars')
        for C in range(1,13):
            if s[0]<=C:
                b,i=d_normal(w,C)
                require(profile(b+'O'*i)==s,'D core reconstruction')
            if c<=C:
                b,typ,i=c_normal(w,C)
                require(profile(b+('Z'*i if typ=='Z' else 'O'*i+'J'))==s,'c core reconstruction')
            for which in ('ell','h'):
                if not (ell<=C if which=='ell' else abs(h)<=C):continue
                b,typ,i,j=normalize(w,C,which)
                require(profile(b)[2]<=4*C,'normal core size')
                ww=b+suffix(typ,i,j)
                require(profile(ww)==s,'normal form profile equality')
                require(3**(len(ww)+1)<=F(s)*2**len(ww),'normal form length bound')
                hist[typ]+=1;controls+=1
                if typ not in samples or len(ww)<len(samples[typ]['normal_word']):
                    samples[typ]={'bound':C,'sidecar':which,'source_word':w,'normal_word':ww,'profile':s}
    # Separate growing exponents, including every six normal forms.
    exp_checks=0
    for b in ('','Z','OZ','JJO'):
        for i in (0,1,2,5,12,40):
            for j in (0,1,3,17,40):
                require(profile(b+'O'*i+'J'+'Z'*j)==profile(b+'O'*(i+j)+'J'),'all-height rewrite test')
                exp_checks+=1
    # Literal graph controls include a nontrivial rewrite and hostile parent 25.
    gw=sorted(set(['','O','Z','JO','OJOO','OJJ','ZJO','OJJJO']+[x['source_word'] for x in samples.values()]+[x['normal_word'] for x in samples.values()]))
    graphs=[]
    for w in gw:
        if len(w)+3>15:continue
        hp,hc=count_graph(graph(w+'J'))
        require(hp==F(profile(w)) and hc==1,'literal final-J count')
        graphs.append({'word':w+'J','vertices':len(w)+3,'paths':hp,'cycles':hc})
    require(F(profile('OJOO'))==25,'hostile parent25')
    hostile=[]
    for n in (3,5,9,17,25,33,65,129,257):
        by_n=[s for s in inv if F(s)==n]
        # U<F, so inventory is complete for these prescribed parent fibers.
        require(n<=Umax,'complete parent fiber scope')
        offs=sorted(set([-sides(s)[0] for s in by_n]+[sides(s)[1] for s in by_n]))
        hostile.append({'parent':n,'profile_count':len(by_n),'minimum_absolute_one_letter_offset':min(map(abs,offs)) if offs else None})
    require(next(x for x in hostile if x['parent']==25)['minimum_absolute_one_letter_offset']>1,'binary hostile25')
    # Count rewrite-normal words by an independent two-state automaton.
    nonJ,J=1,0
    fib_controls=[]; cumulative=1
    for l in range(31):
        if l:
            nonJ,J=2*nonJ+J,nonJ+J
            cumulative+=nonJ+J
        require(nonJ+J==fib(2*l+2),'forbidden-bigram count')
        require(cumulative==fib(2*l+3)-1,'cumulative Fibonacci count')
        fib_controls.append({'length':l,'words':nonJ+J,'cumulative':cumulative})
    # Independently count sidecar support in the complete physical inventory.
    finite=[]
    for C in (1,2,3,5,8,12):
        vals=sorted({F(s) for s in inv if F(s)<=Umax and (sides(s)[0]<=C or abs(sides(s)[1])<=C)})
        bcount=sum(s[2]<=4*C for s in inv)
        L=length_bound(Umax)
        fibbound=fib(2*length_bound(8*C)+3)-1
        require(bcount<=fibbound,'polynomial core count')
        bound=bcount*(3*(L+1)**2+3*(L+1))
        require(len(vals)<=bound,'quantitative finite control')
        finite.append({'C':C,'bounded_core_profiles':bcount,'Fibonacci_core_bound':fibbound,'qualifying_parent_counts_through_512':len(vals),'proven_profile_bound':bound})
    print(json.dumps({'status':'PASS; unreviewed author controls','scope':'physical normal forms and bounded-offset obstruction, not target occupancy','cover_inventory_bound':Umax,'source_profiles':len(inv),'normal_form_reconstructions':controls,'normal_form_usage':dict(hist),'large_exponent_rewrites':exp_checks,'normal_form_examples':samples,'literal_graph_controls':graphs,'parent_fiber_controls':hostile,'finite_bound_controls':finite,'Fibonacci_controls':fib_controls,'proof_sha256':hashlib.sha256(Path(__file__).with_name('PROOF.md').read_bytes()).hexdigest()},indent=2,sort_keys=True))

if __name__=='__main__':main()
