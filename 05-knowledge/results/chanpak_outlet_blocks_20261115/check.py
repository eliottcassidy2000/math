from fractions import Fraction as F
from math import prod
from pathlib import Path

def rev(x, width):
    return int(f'{x:0{width}b}'[::-1], 2)

def J(e):
    bins = [i for i in range(e.bit_length()) if e >> i & 1]
    return sum(2 if j-i <= 2 else 1 << (j-i-1) for i in bins for j in bins if i < j)

def graph(e, B):
    m = e.bit_length()+1
    u = m+len(B)
    return tuple([sum(1<<j for j in range(i+1,u)) + ((1<<u) if i == m-1 or not ((e >> (m-2-i))&1) else 0) for i in range(m)] +
        [(row << m) | (1<<u) for row in B] +
        [sum(1<<i for i in range(m-1) if (e >> (m-2-i))&1)])

def subset_counts(rows):
    n = len(rows)
    hp = [[0]*n for _ in range(1<<n)]
    cnt = [1]*(1<<n)
    for s in range(1,1<<n):
        for v in range(n):
            if not (s>>v)&1: continue
            sub = s^(1<<v)
            hp[s][v] = 1 if not sub else sum(hp[sub][u] for u in range(n) if rows[u]>>v&1)
        cnt[s] = sum(hp[s])
    return cnt

def profile(rows):
    d = subset_counts(rows)
    full=len(d)-1
    return d[-1],sum(d[s]*d[full^s] for s in range(full+1))


def mixed_literal(cores, marks):
    root=sum(map(len,cores))
    assert root>0
    rows=[]
    contacts=0
    offset=0
    for B,w in zip(cores,marks):
        for i,row in enumerate(B):
            idx=offset+i
            if w: contacts|=1<<idx
            rows.append((row<<offset)|sum(1<<j for j in range(offset+len(B),root))|((1<<root) if not w else 0))
        offset+=len(B)
    rows.append(contacts)
    new=root+1
    rows_out=tuple([r|(1<<new) for r in rows[:-1]]+[rows[-1],1<<root])
    return tuple(rows),rows_out

def mixed_states(cores,marks):
    D,C,N,z=1,0,1,0
    for B,w in zip(cores,marks):
        h,P=profile(B)
        D,C,N,z=h*D,P*C+w*(P-h)*D,h*N+(2-w)*(P-h)*C+w*(P-h*z)*D,w
    return D,C,N,z

def main():
    data=iter(Path(__file__).with_name('profiles.txt').read_text().splitlines())
    profiles=[]
    counts=[]
    for line in data:
        n,k=map(int,line.split())
        counts.append(k)
        profiles.extend(tuple(map(int,next(data).split())) for _ in range(k))
    assert counts == [0,1,1,5,16,131]
    assert len(profiles)==len(set(profiles))==154
    a=prod(F(p,p-1) for h,p in profiles)
    b=prod(F(p*p,p*p-h) for h,p in profiles)
    assert all(4*h<=p for h,p in profiles)
    assert a<F(1398,1000) and b<F(1046,1000)
    # log(4/3) < 2/7 + 2/(3*7^3)/(1-1/49) < 288/1000
    assert F(2,7)+F(2,3*7**3)/(1-F(1,49)) < F(288,1000)
    c=prod(F(p**3,p**3-h**2) for h,p in profiles)
    log_upper=sum((F(2,(2*j+1)*7**(2*j+1)) for j in range(3)),F(0))+F(2,7**8)/(1-F(1,49))
    bound=log_upper+(a-1)/2+(b-1)/8+(c-1)/24+(a-1)/3584
    assert bound<F(493,1000)
    print(f'profiles {counts}; certified mass bound {float(bound):.9f}')
    B=(8,1,3,6)
    assert profile(B)==(5,32)
    checked=0
    for core in [(0,),(2,4,1),B,(24,17,19,22,0)]:
        h,p=profile(core)
        for e in range(1,30):
            expect=h+2*p*e-h*J(e)
            assert subset_counts(graph(e,core))[-1]==expect
            checked+=1
    print(f'{checked} independent literal enlarged tournament counts')
    ss=[]
    for x,y in [(106,17773),(103,18613)]:
        coef=y*(64*128-5*rev(x,7)); base=5+64*x-5*(J(x)+J(y))
        assert x.bit_count()+y.bit_count()==12 and y.bit_length()==15
        ss.append((coef,base))
        for k in range(8,24):
            e=x+(y<<k)
            assert 5+64*e-5*J(e)==coef*2**(k-7)+base
    assert ss==[(141775221,-56761),(141775221,-56763)]
    print('same-graph adjacent bands checked')
    C3=(2,4,1)
    checked=0
    for cores in [((0,),), (B,C3), (C3,B), ((0,),B,(0,),C3), ((0,),C3,(0,),(0,))]:
        for mask in range(1<<len(cores)):
            marks=[mask>>j&1 for j in range(len(cores))]
            s,t=mixed_literal(cores,marks)
            assert mixed_states(cores,marks)[2]==subset_counts(t)[-1]
            if C3 in cores:
                i=cores.index(C3)
                flat_cores=cores[:i]+((0,),(0,))+cores[i+1:]
                flat_marks=marks[:i]+[marks[i],marks[i]]+marks[i+1:]
                sf,tf=mixed_literal(flat_cores,flat_marks)
                assert 3*subset_counts(sf)[-1]==subset_counts(s)[-1]
                assert 3*subset_counts(tf)[-1]==subset_counts(t)[-1]
            checked+=1
    for cores,marks,W in [(((0,),)*5+(B,),[1,1,0,1,1,0],1623),((B,)+((0,),)*5,[1,0,0,0,0,0],1711)]:
        _,t=mixed_literal(cores,marks)
        assert subset_counts(t)[-1] == W
        assert mixed_states(cores,marks)==(5,864,W,0)
        checked+=1
    print(f'{checked} literal mixed-radix controls, continuation pair (1623,1711) checked')
    unit=[]
    for core,p,loc,mask,want in [
        ((34,8,3,5,15,30),272,10,30557,9971529),
        ((48,1,3,7,14,30),296,14,19413,9971531),
    ]:
        assert profile(core)==(23,p)
        cores=((0,),)*loc+(core,)+((0,),)*(15-loc)
        marks=[mask>>j&1 for j in range(16)]
        assert mixed_states(cores,marks)==(23,6510178,want,0)
        marked=sum(1<<i for i in range(21) if marks[i if i<loc else loc if i<loc+6 else i-5])
        full=(1<<21)-1
        H,E=0,0
        ds=subset_counts(core)
        # Independent literal split around apex, and enumeration of even paths from it.
        for right in range(full+1):
            left=full^right
            if right and not (marked & (right & -right)): continue
            t=ds[(right>>loc)&63]*ds[(left>>loc)&63]
            if not left or not ((marked>>(left.bit_length()-1))&1): H+=t
            if right.bit_count()&1: E+=t
        assert H+2*E==want
        unit.append((H,E))
    assert unit==[(3461351,3255089),(3461353,3255089)]
    print(f'physical continuation-locked unit pair, independent exhaustive apex splits: {unit}')
if __name__=='__main__':
    main()
