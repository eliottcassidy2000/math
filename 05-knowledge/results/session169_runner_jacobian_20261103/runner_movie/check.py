"""Exact LRC wall movie, inverse-edge transports, and tournament lift.

Standard library only. All checks survive python -O. No repository reads.
"""
from fractions import Fraction as Q
from itertools import permutations, combinations
from math import factorial
from random import Random
import hashlib, json

def require(c, message):
    if not c: raise RuntimeError(message)

def frac(x): return x - x.numerator // x.denominator
def distance(x):
    y=frac(x)
    return min(y,1-y)

def status(v,p,d,t):
    a=distance(p+v*t)
    return 0 if a<d else 1 if a==d else 2

CODE=((0,0),(0,1),(1,1))

def walls(v,p,d):
    require(0<d<Q(1,2), 'threshold range')
    require(all(a and isinstance(a,int) for a in v),'nonzero integral speeds')
    require(all(distance(b)!=d for b in p),'cut must not be a wall')
    out={}
    for i,(a,b) in enumerate(zip(v,p)):
        b=frac(b)
        for k in range(-abs(a)-2,abs(a)+3):
            for sign in (-1,1):
                t=(Q(k)+sign*d-b)/a
                if 0<t<1: out.setdefault(t,[]).append(i)
    result=sorted(out.items())
    require(sum(len(inds) for t,inds in result)==2*sum(abs(a) for a in v),'wall incidence count')
    return result

def threshold_path(v,p,d,reflected=False):
    ev=walls(v,p,d); times=[Q(0)]+[t for t,_ in ev]+[Q(1)]
    initial=tuple(b for a,q in zip(v,p) for b in CODE[status(a,q,d,Q(0))])
    state=list(initial); states=[initial]; steps=[]; cells=[]; atoms=[]
    for l,r in zip(times,times[1:]):
        st=tuple(status(a,q,d,(l+r)/2) for a,q in zip(v,p))
        cells.append((l,r,st))
    for h,(t,inds) in enumerate(ev):
        atom=tuple(status(a,q,d,t) for a,q in zip(v,p)); atoms.append((t,atom))
        left=cells[h][2]; right=cells[h+1][2]
        order=sorted(inds,reverse=bool(reflected and t>Q(1,2)))
        for phase,target in [('atom',atom),('cell',right)]:
            for i in order:
                for b in (0,1):
                    k=2*i+b; desired=CODE[target[i]][b]
                    if state[k]!=desired:
                        old=state[k]; state[k]=desired
                        steps.append({'t':str(t),'phase':phase,'runner':i,'bit':b,'direction':desired-old})
                        states.append(tuple(state))
            require(tuple(state)==tuple(b for s in target for b in CODE[s]),'macrostate decoding')
    require(states[-1]==initial,'periodic threshold path')
    return states,steps,cells,atoms

def omega(n,a,b):
    if not n:return 0
    bit=1<<(n-1); h=bool(a&bit); k=bool(b&bit); a%=bit; b%=bit
    if not h and not k:return omega(n-1,a,b)
    if not h and k:return omega(n-1,b,a)
    if h and not k:return omega(n-1,a,b)^bool(b)
    return omega(n-1,b,a)^1^bool(b)

def low_digit_omega(n,a,b):
    w=v=A=B=0
    for j in range(n):
        h=(a>>j)&1;k=(b>>j)&1
        w,v=((v if k else w)^(h&k)^(h&B),
             (w if h else v)^(h&k)^(k&A))
        A|=h;B|=k
    return w

def transport(n,a):
    return [(s^omega(n,a,x))*(1<<n)+(a^x) for s in (0,1) for x in range(1<<n)]

def inverse(p):
    q=[0]*len(p)
    for x,y in enumerate(p):q[y]=x
    require(len(set(p))==len(p),'basis transport permutation')
    return q

def holonomy(n,assign,steps):
    pos=[transport(n,a) for a in assign]; neg=[inverse(p) for p in pos]
    g=list(range(2<<n)); word=[]
    for step in steps:
        if step['bit']:continue
        i=step['runner']; direction=step['direction']
        p=(pos if direction==1 else neg)[i]
        g=[p[x] for x in g]
        word.append((i,direction))
    return g,word

def tournament(word):
    N=5*len(word); A=[[0]*N for _ in range(N)]
    for i,j in combinations(range(N),2):
        bi,bj=i//5,j//5
        if bi!=bj:A[i][j]=1
        else:
            pairs=list(combinations(range(5),2)); bit=pairs.index((i%5,j%5))
            if (8+2*word[bi])>>bit&1:A[i][j]=1
            else:A[j][i]=1
    return A

def encoded_movie(n,assign,states,steps):
    threshold=list(states[0]); unit=0; unitbits=[0]*(n+1)
    words=[tuple(threshold+unitbits)]; flips=[]
    pos=[transport(n,a) for a in assign]; neg=[inverse(p) for p in pos]
    for st,step in zip(states[1:],steps):
        k=2*step['runner']+step['bit']; threshold=list(st)
        words.append(tuple(threshold+unitbits)); flips.append(k)
        if step['bit']==0:
            p=(pos if step['direction']==1 else neg)[step['runner']]
            new=p[unit]; newbits=[(new>>j)&1 for j in range(n+1)]
            for j in range(n+1):
                if unitbits[j]!=newbits[j]:
                    unitbits[j]=newbits[j]
                    words.append(tuple(threshold+unitbits)); flips.append(len(threshold)+j)
            unit=new
    for old,new,k in zip(words,words[1:],flips):
        A,B=tournament(old),tournament(new)
        diff=[(i,j) for i in range(len(A)) for j in range(i+1,len(A)) if A[i][j]!=B[i][j]]
        require(diff==[(5*k,5*k+2)],'literal one-arc flip')
    digest=hashlib.sha256(json.dumps(words).encode()).hexdigest()
    return {'order':5*len(words[0]),'arc_flips':len(flips),'final_unit':unit,'state_hash':digest}

def profile(mask):
    A=tournament((0 if mask==8 else 1,)); n=5; H={0:1}
    for S in range(1,1<<n):
        vs=[i for i in range(n) if S>>i&1]
        H[S]=sum(all(A[a][b] for a,b in zip(p,p[1:])) for p in permutations(vs))
    cache={0:[1]+[0]*n}
    for S in range(1,1<<n):
        out=[0]*(n+1); first=S&-S; U=S
        while U:
            if U&first:
                for j,c in enumerate(cache[S^U][:-1]):out[j+1]+=H[U]*c
            U=(U-1)&S
        cache[S]=out
    # Independent mandatory-cut enumeration of all ordered path lists.
    raw=[0]*(n+1)
    for p in permutations(range(n)):
        for cuts in range(1<<(n-1)):
            if all((cuts>>j&1) or A[p[j]][p[j+1]] for j in range(n-1)):
                raw[1+cuts.bit_count()]+=1
    independent=[raw[j]//factorial(j) for j in range(n+1)]
    require(cache[(1<<n)-1]==independent,'independent profile')
    return independent

def safe_intersection(v,p,d):
    # Independent closed-interval construction; no wall/state function used.
    current=[(Q(0),Q(1))]
    for a,b in zip(v,p):
        b=b-b.numerator//b.denominator; own=[]
        for k in range(-abs(a)-2,abs(a)+3):
            ends=sorted(((Q(k)+d-b)/a,(Q(k+1)-d-b)/a))
            lo=max(Q(0),ends[0]);hi=min(Q(1),ends[1])
            if lo<=hi:own.append((lo,hi))
        new=[]
        for x,y in current:
            for z,w in own:
                lo=max(x,z);hi=min(y,w)
                if lo<=hi:new.append((lo,hi))
        current=sorted(set(new))
    return current

def case(v,p,d,n,assign,reflect):
    states,steps,cells,atoms=threshold_path(v,p,d,reflect)
    g,word=holonomy(n,assign,steps)
    size=1<<n
    require(all(g[x]%size==x%size for x in range(2*size)),'zero XOR return')
    require(all(g[g[x]]==x for x in range(2*size)),'two-lap ceiling')
    if reflect:
        require(states==list(reversed(states)),'literal reflected threshold path')
        require(g==list(range(2*size)),'reflection transport cancellation')
    safe=[(l,r) for l,r,s in cells if min(s)==2]
    closed=safe_intersection(v,p,d)
    cut_trits=[status(a,q,d,Q(0)) for a,q in zip(v,p)]
    cut={'time':'0','trits':cut_trits,'strict_safe':min(cut_trits)==2,'closed_safe':min(cut_trits)>=1}
    require(cut['strict_safe']==cut['closed_safe'],'nonwall cut strict/closed membership')
    require(cut['closed_safe']==any(l<=0<=r for l,r in closed),'independent cut membership at zero')
    require(cut['closed_safe']==any(l<=1<=r for l,r in closed),'periodic cut membership at one')
    require(cut['strict_safe']==(min(cells[0][2])==2)==(min(cells[-1][2])==2),'cut agrees with adjacent open cells')
    mass=sum((r-l for l,r in safe),Q(0))
    require(mass==sum((r-l for l,r in closed),Q(0)),'independent safe mass')
    for t,s in atoms:
        require((min(s)>=1)==any(l<=t<=r for l,r in closed),'closed endpoint membership')
    return {'speeds':v,'phases':list(map(str,p)),'threshold':str(d),'wall_times':[str(t) for t,_ in walls(v,p,d)],'strict_intervals':[[str(l),str(r)] for l,r in safe],'cut_atom':cut,'closed_intervals':[[str(l),str(r)] for l,r in closed],'strict_mass':str(mass),'transport_word':word,'holonomy':g,'threshold_flips':len(steps)},states,steps

def main():
    pair_checks=0
    for n in range(1,8):
        for a in range(1<<n):
            for b in range(1<<n):
                require(omega(n,a,b)==low_digit_omega(n,a,b),'independent doubling sign')
                pair_checks+=1
    profiles=[profile(x) for x in (8,10)]
    require(profiles==[[0,9,33,31,10,1]]*2,'inherited binary seed profile')
    examples={}
    for label,v,p,reflection in [('unshifted',[1,2],[Q(0),Q(0)],True),('shifted',[1,2],[Q(0),Q(1,4)],False)]:
        out,states,steps=case(v,p,Q(1,3),2,[1,2],reflection)
        out['physical_lift']=encoded_movie(2,[1,2],states,steps)
        examples[label]=out
    require(examples['unshifted']['strict_mass']=='0','unshifted boundary atom control')
    require(examples['shifted']['strict_mass']=='1/8','shifted strict window')
    require(examples['shifted']['holonomy']==list(range(4,8))+list(range(4)),'quaternion -identity')
    cut_regression,_,_=case([1],[Q(1,2)],Q(1,3),1,[1],False)
    require(cut_regression['cut_atom']['strict_safe'],'strict-safe cut atom retained')
    require(cut_regression['strict_intervals']==[['0','1/6'],['5/6','1']],'cut-crossing open intervals')
    require(cut_regression['strict_mass']=='1/3','cut-crossing mass')
    examples['strict_safe_cut']=cut_regression
    phase_chamber_checks=0
    for q in range(4,25):
        for r in range(1,q):
            theta=Q(r,q)
            if theta>=Q(1,3):continue
            out,_,_=case([1,2],[Q(0),theta],Q(1,3),2,[1,2],False)
            require(Q(out['strict_mass'])==theta/2,'phase chamber linear mass')
            require(out['holonomy']==list(range(4,8))+list(range(4)),'phase chamber constant return')
            phase_chamber_checks+=1
    # Simultaneous entries/exits: all exact atom states are retained, but tie
    # ordering still affects noncommuting transport. The asymmetry is supplied.
    reflected,_,_=case([1,4],[Q(0),Q(0)],Q(1,3),2,[1,2],True)
    unreflected,_,_=case([1,4],[Q(0),Q(0)],Q(1,3),2,[1,2],False)
    require(reflected['holonomy']!=unreflected['holonomy'],'serialization changes holonomy')
    examples['tie_reflected']=reflected; examples['tie_same_order']=unreflected
    count=0
    for v in combinations(range(1,9),3):
        for d in (Q(1,4),Q(1,3),Q(2,5)):
            case(list(v),[Q(0)]*3,d,3,[1,2,4],True);count+=1
    rng=Random(169);shifted=0
    for _ in range(400):
        v=rng.sample(list(range(-7,0))+list(range(1,8)),3)
        p=[Q(rng.randrange(11),11) for i in range(3)]
        d=rng.choice((Q(1,4),Q(1,3),Q(2,5)))
        if any(distance(q)==d for q in p):continue
        n=rng.choice((2,3,4,5));assign=[rng.randrange(1,1<<n) for i in range(3)]
        case(v,p,d,n,assign,False);shifted+=1
    return {'seed_profiles':profiles,'doubling_sign_pair_checks':pair_checks,'phase_chamber_checks':phase_chamber_checks,'unshifted_exact_cases':count,'shifted_exact_cases':shifted,'examples':examples,'claims':'Finite controls only; universal statements use NOTE.md proofs.'}

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
