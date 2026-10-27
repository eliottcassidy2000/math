#!/usr/bin/env python3
"""Exact native-source, clock, paid-forcing, and exception-budget controls."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product, permutations
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
VPATH = HERE/'frontier.txt'
DEN = 10**12


def check(ok, why):
    if not ok:
        raise RuntimeError(why)


def edge(s):
    a = 2 if s % 8 == 1 else 1
    y, t = (3*s+1) >> a, a
    while y % 8 == 5:
        y = (y-1)//4
        t += 2
    return y, t, a


def norm(z):
    a,b = z
    return a*a-a*b+b*b


def charge(s, N):
    row = [0,0]
    lower = ((1,2),(-2,-1),(1,-1))
    upper = ((0,1),(-1,-1),(1,0))
    m = s
    real_rows = [0,0,0]
    while m < 4*N:
        v = (lower if m < 2*N else upper)[m % 3]
        row[0] += v[0]; row[1] += v[1]
        w = (0,-3,3) if m < 2*N else (-1,-1,2)
        for r in range(3):
            real_rows[r] += w[(m+r) % 3]
        m = 4*m+1
    a,b = row
    check(real_rows == [2*a-b,-a-b,-a+2*b], 'independent three-row source charge')
    check(norm(row) in (0,1,3), 'native charge palette')
    return tuple(row)


def rational_interval(numerator, denominator):
    q,r = divmod(numerator,denominator)
    return q,q+bool(r)


def prefix_debt(records, S, max_modulus, position):
    moduli=[];q=1
    while q<=max_modulus:moduli.append(q);q*=2
    bins=[defaultdict(lambda:[0,0]) for _ in moduli]
    energies=[0]*len(moduli)
    batches=defaultdict(list)
    for row in records:batches[row[position]].append(row)
    peak=(0,0);last=[]
    for level in sorted(batches):
        for s,u,t,d,p,a,b in batches[level]:
            if (a,b)==(0,0):continue
            sign=(-1)**t;a*=sign;b*=sign
            for j,q in enumerate(moduli):
                v=bins[j][u,((s-1)//2)%q];old=norm(v)
                v[0]+=a;v[1]+=b;energies[j]+=norm(v)-old
        debt=sum(max(energies[j]-energies[j+1],0) for j in range(len(energies)-1))
        if debt>peak[0]:peak=(debt,level)
        check(energies[0] <= energies[-1]+debt, 'prefix dyadic sufficient bound')
        check(energies[-1] <= S, 'pending sources retained in denominator')
    return dict(max_positive_debt=peak[0],at_cutoff=peak[1],source_norm=S,
                cutoffs_checked=len(batches),final_energies=energies)


def native_checks():
    V = {int(l.split()[0]) for l in VPATH.read_text().splitlines()
         if l.strip() and not l.startswith('#')}
    ranks = {int(l.split()[0]):int(l.split()[1]) for l in VPATH.read_text().splitlines()
             if l.strip() and not l.startswith('#')}
    check(len(V)==627 and max(V)==63665, 'frontier source pin')
    check(all(edge(s)[0] in V for s in V), 'frontier forward closed')
    check(ranks[1]==0 and all(s==1 or ranks[edge(s)[0]]<ranks[s] for s in V),
          'frontier grounded rank retained')
    check(all(s in V for s in range(1,1009,2) if s % 8 != 5), 'all low seeds stop')
    checked=0
    for s in range(1,100002,2):
        if s % 8 == 5:
            continue
        y,t,a=edge(s); X=3*s+1; Y=3*y+1
        check(2**t*Y == 3*X + 2**a, 'corrected native affine carry')
        b=F(3*X,3*X+2**a); d=1-b*b
        lam=F(2**a,6*X+2**a)
        check(lam==(1-b)/(1+b), 'phase-gauged forcing price')
        if s not in V:
            check(lam <= F(t,9139), 'sharp per-clock forcing bound')
        # Both native signs; vector arithmetic uses the Eisenstein norm.
        for A in ((1,2),(-2,-1),(1,-1),(0,1),(-1,-1),(1,0)):
            for sign in (-1,1):
                f=tuple((sign-b)*x for x in A)
                r=tuple(d*x-b*z for x,z in zip(A,f))
                check(norm(tuple(sign*x for x in A))+norm(r)/d == norm(A)+norm(f)/d,
                      'native two-port energy')
                price=lam if sign==1 else 1/lam
                check(norm(f)/d == price*norm(A), 'paid sign cost')
        checked+=1
    check(F(2,18*1015+8)==F(edge(1015)[1],9139), 'sharp clock witness')
    check(1015 not in V, 'sharp witness is exterior')
    word_controls=0
    for initial in (1009,1011,1015,4097,10007,65535,1000001):
        if initial%8==5:continue
        x=initial;total=0;carry=0;cost=F(0);gain=F(1)
        for dep in range(1,17):
            y,t,a=edge(x);X=3*x+1
            carry=3*carry+2**(total+a);total+=t
            lam=F(2**a,6*X+2**a);cost+=lam
            gain*=F(3*X,3*X+2**a)
            check(2**total*(3*y+1)==3**dep*(3*initial+1)+carry,'word carry')
            check(gain==F(3**dep*(3*initial+1),2**total*(3*y+1)), 'telescoped gain')
            check(cost<=F(carry,2*3**dep*(3*initial+1)), 'rational carry payment')
            word_controls+=1
            x=y
            if x in V:break
    # Route record: terminal, integer clock, native depth, entry port,
    # lower/upper integer certificates for gauged and ungauged forcing.
    memo={s:(s,0,0,0,0,0,0,0) for s in V}
    def route(s):
        x=s; path=[]; seen=set()
        while x not in memo:
            if x in seen or len(path)>=10000:
                raise RuntimeError(('unresolved source retained; stop finite test',s,x))
            seen.add(x); y,t,a=edge(x); path.append((x,t,a)); x=y
        u,t,d,p,lo,hi,ulo,uhi=memo[x]
        for x,dt,a in reversed(path):
            if d==0:p=x
            num,den=2**a,18*x+6+2**a
            dl,dh=rational_interval(DEN*num,den)
            if a==1:
                ul,uh=rational_interval(DEN*den,num)
            else:
                ul,uh=dl,dh
            lo+=dl;hi+=dh;ulo+=ul;uhi+=uh;t+=dt;d+=1
            memo[x]=(u,t,d,p,lo,hi,ulo,uhi)
        return memo[s]
    rows=[]
    for N in (4096,5997,11002,16384,65536):
        records=[]; S=lo=hi=ulo=uhi=weighted_clock=0
        q=defaultdict(lambda:[0,0]); ports=defaultdict(lambda:[0,0])
        digest=hashlib.sha256()
        maxclock=maxdepth=0
        for s in range(1,4*N,2):
            if s % 8==5:continue
            A=charge(s,N); w=norm(A)
            u,t,d,p,l,h,ul,uh=route(s)
            S+=w;lo+=w*l;hi+=w*h;ulo+=w*ul;uhi+=w*uh;weighted_clock+=w*t
            signed=tuple((-1 if t%2 else 1)*x for x in A)
            for j in range(2):q[u][j]+=signed[j];ports[u,p,t][j]+=signed[j]
            records.append((s,u,t,d,p,*A))
            digest.update((','.join(map(str,records[-1]))+'\n').encode())
            maxclock=max(maxclock,t);maxdepth=max(maxdepth,d)
        check(9139*lo <= DEN*weighted_clock, 'source-weighted clock payment')
        energies=[]
        modulus=1
        while modulus < 2*N:
            bins=defaultdict(lambda:[0,0])
            for s,u,t,d,p,a,b in records:
                sign=-1 if t%2 else 1
                v=bins[u,((s-1)//2)%modulus]
                v[0]+=sign*a;v[1]+=sign*b
            energies.append([modulus,sum(norm(v) for v in bins.values())])
            modulus*=2
        bins=defaultdict(lambda:[0,0])
        for s,u,t,d,p,a,b in records:
            sign=-1 if t%2 else 1
            v=bins[u,((s-1)//2)%modulus];v[0]+=sign*a;v[1]+=sign*b
        energies.append([modulus,sum(norm(v) for v in bins.values())])
        E=sum(norm(v) for v in q.values())
        check(energies[0][1]==E and energies[-1][1]==S,'source refinement endpoints')
        increments=[energies[j][1]-energies[j+1][1] for j in range(len(energies)-1)]
        check(sum(increments)==E-S,'coherence debt telescoping')
        # Separate pair-polarization formulation at the first six cuts.
        for k in range(min(6,len(energies)-1)):
            m=energies[k][0];fine=defaultdict(lambda:[0,0])
            for s,u,t,d,p,a,b in records:
                sign=(-1)**t;v=fine[u,((s-1)//2)%(2*m)]
                v[0]+=sign*a;v[1]+=sign*b
            cross=0
            for u in q:
                for r in range(m):
                    x=fine[u,r];y=fine[u,r+m]
                    cross+=norm((x[0]+y[0],x[1]+y[1]))-norm(x)-norm(y)
            check(cross==increments[k],'independent residue-cut polarization')
        prefix={}
        if N<=16384:
            for name,position in (('depth',3),('clock',2)):
                prefix[name]=prefix_debt(records,S,modulus,position)
                check(prefix[name].pop('final_energies')==[e for _,e in energies],
                      'streamed prefix versus separately rebuilt bins')
        rows.append(dict(N=N,sources=len(records),source_norm=S,nyquist_norm=E,
                         max_clock=maxclock,max_depth=maxdepth,route_sha256=digest.hexdigest(),
                         retained_port_clock_norm=sum(norm(v) for v in ports.values()),
                         forcing_ratio_interval=[str(F(lo,DEN*S)),str(F(hi,DEN*S))],
                         ungauged_ratio_interval=[str(F(ulo,DEN*S)),str(F(uhi,DEN*S))],
                         weighted_clock=weighted_clock,refinement_energies=energies,
                         sum_positive_refinement_debt=sum(max(x,0) for x in increments),
                         prefix_debt=prefix))
    # Target chosen before the complete scan: positive dyadic source debt <= S/2.
    # The complete integer-scale scan uses three threshold events, not B-chain sums.
    maxN=65536;moduli=[];q0=1
    while q0<2*maxN:moduli.append(q0);q0*=2
    moduli.append(q0)
    bins=[defaultdict(lambda:[0,0]) for _ in moduli];energies=[0]*len(moduli)
    sources=defaultdict(lambda:[0,0]);S=Flo=Fhi=0;violations=0
    debt_peak=(0,1,0);force_peak=(0,1,0);first_bad=None
    row_lookup={r['N']:r for r in rows}
    upper=((0,1),(-1,-1),(1,0));change=((1,1),(-1,0),(0,-1))
    for N in range(1,maxN+1):
        for m,delta in ((4*N-3,upper[(4*N-3)%3]),(4*N-1,upper[(4*N-1)%3]),
                        (2*N-1,change[(2*N-1)%3])):
            s=m
            while s%8==5:s=(s-1)//4
            u,t,d,p,l,h,ul,uh=route(s)
            A=sources[s];old=norm(A);A[0]+=delta[0];A[1]+=delta[1];dw=norm(A)-old
            S+=dw;Flo+=dw*l;Fhi+=dw*h
            sign=(-1)**t;a,b=(sign*x for x in delta)
            for j,q0 in enumerate(moduli):
                v=bins[j][u,((s-1)//2)%q0];old=norm(v)
                v[0]+=a;v[1]+=b;energies[j]+=norm(v)-old
        check(energies[-1]==S,'finest dyadic energy equals actual source norm')
        debt=sum(max(energies[j]-energies[j+1],0) for j in range(len(energies)-1))
        if N in row_lookup:
            r=row_lookup[N]
            check((S,energies[0],debt)==(r['source_norm'],r['nyquist_norm'],r['sum_positive_refinement_debt']),
                  'independent all-N threshold-event source and debt')
            check(str(F(Fhi,DEN*S))==r['forcing_ratio_interval'][1],'independent paid forcing')
        if N>=4096:
            if 2*debt>S:
                violations+=1
                if first_bad is None:first_bad=N
            if debt*debt_peak[1]>debt_peak[0]*S:debt_peak=(debt,S,N)
            if Fhi*force_peak[1]>force_peak[0]*S:force_peak=(Fhi,S,N)
    scan=dict(first_N=4096,last_N=maxN,integer_scales=maxN-4096+1,
              sufficient_half_source_debt_violations=violations,first_violation=first_bad,
              maximum_debt_ratio=str(F(debt_peak[0],debt_peak[1])),maximum_debt_N=debt_peak[2],
              maximum_forcing_upper_ratio=str(F(force_peak[0],DEN*force_peak[1])),
              maximum_forcing_N=force_peak[2])
    return dict(primitive_edges_checked=checked,denominator=DEN,rows=rows,all_N_scan=scan,
                word_carry_controls=word_controls,
                sharp_exterior_cost_per_clock='1/9139',frontier_sha256=hashlib.sha256(VPATH.read_bytes()).hexdigest())


def graph_beta(n,edges):
    p=list(range(n));components=n
    def root(x):
        while p[x]!=x:x=p[x]
        return x
    for x,y in edges:
        a,b=root(x),root(y)
        if a!=b:p[a]=b;components-=1
    return len(edges)-n+components


def colour_orders(sizes):
    word=[];remaining=list(sizes)
    def visit():
        if not any(remaining):
            yield tuple(word);return
        for c in range(3):
            if remaining[c]:
                word.append(c);remaining[c]-=1
                yield from visit()
                remaining[c]+=1;word.pop()
    yield from visit()


def cycle_checks():
    cases=[]
    for a in range(1,5):
        for b in range(a,5):
            for c in range(b,5):
                n=a+b+c;bestbeta=n*n;bestedges=n*n;bestbeta_at_edges=n*n;count=0
                for word in colour_orders((a,b,c)):
                    ds=[(i,j) for i in range(n) for j in range(i+1,n)
                        if word[i]!=word[j] and (word[j]-word[i])%3!=1]
                    beta=graph_beta(n,ds)
                    bestbeta=min(bestbeta,beta)
                    if len(ds)<bestedges:bestedges=len(ds);bestbeta_at_edges=beta
                    elif len(ds)==bestedges:bestbeta_at_edges=min(bestbeta_at_edges,beta)
                    count+=1
                check(bestbeta==max(a-2,0)*(b-1),'three-module exact cycle minimum')
                check(bestedges==a*b,'three-module minimum feedback edges')
                cases.append(dict(sizes=[a,b,c],orders=count,minimum_beta=bestbeta,minimum_edges=bestedges,
                                  minimum_beta_at_minimum_edges=bestbeta_at_edges))
    # All labelled orders for the 2+2+2 control, including reversed clone orders.
    n=6;edgesmin=99;betamin=99
    for p in permutations(range(n)):
        ds=[]
        for i,x in enumerate(p):
            for y in p[i+1:]:
                forward=(x<y) if x//2==y//2 else (y//2-x//2)%3==1
                if not forward:ds.append((x,y))
        edgesmin=min(edgesmin,len(ds));betamin=min(betamin,graph_beta(n,ds))
    check((edgesmin,betamin)==(4,0),'unrestricted clone-order control')
    modules=0
    for mask in range(1,128):
        inside=[x for x in range(7) if mask>>x&1]
        if not 2<=len(inside)<=6:continue
        modules+=1
        outside=[x for x in range(7) if not mask>>x&1]
        check(not all(len({(y-x)%7 in (1,2,4) for y in inside})==1 for x in outside),
              'Paley7 has no nontrivial homogeneous module')
    return dict(cases=cases,unrestricted_six_vertex_orders=720,paley_module_subsets_checked=modules)


def aliquot_checks():
    bound=400;s=[0]*(bound+1)
    for d in range(1,bound+1):
        for n in range(2*d,bound+1,d):s[n]+=d
    fibres={m:[n for n in range(1,(m-1)**2+1) if s[n]==m] for m in (2,5,7,21)}
    check(fibres[2]==[] and fibres[5]==[] and fibres[7]==[8], 'certified small fibres')
    check(18 in fibres[21], '21 is touchable')
    check(sum(d for d in (1,2,4))==7 and sum(d*d for d in (1,2,4))==21,
          'two weighted views of one divisor chain')
    return dict(complete_inverse_fibres=fibres,source_bound=400)


def norm_seven_checks():
    def mul(x,y):
        a,b=x;c,d=y
        return a*c-b*d,a*d+b*c-b*d
    alpha=(3,1);seen=0
    check(mul((1,-1),alpha)==(4,-1), 'norm21 factorization')
    for a,b in product(range(-7,8),repeat=2):
        x=(a,b);y=mul(alpha,x)
        check(y==(3*a-b,a+2*b),'twist matrix')
        check(norm(y)==7*norm(x),'norm multiplier')
        check((5*y[0]-y[1])==7*(2*a-b),'marked ROOT readout')
        check((2*y[0]+y[1])%7==0,'index7 image guard')
        z=x
        for k in range(1,7):
            z=mul(alpha,z)
            check(norm(z)==7**k*norm(x),'iterated normalization')
        seen+=1
    inverse_lifts=0
    for u in range(1,102,2):
        if u%8==5:continue
        for a in (1,2):
            j0=(3-a-u)%3
            for r in range(4):
                j=j0+3*r
                numerator=2**a*(3*u+1)*4**j-(2**a+3)
                check(numerator%9==0,'native inverse divisibility guard')
                s=numerator//9
                snext=(2**a*(3*u+1)*4**(j+3)-(2**a+3))//9
                check(snext==64*s+7*(2**a+3),'norm-seven inverse lift offset')
                check(edge(s)==(u,a+2*j,a),'inverse lift preserves full native clock')
                check(edge(snext)==(u,a+2*j+6,a),'six additional ticks')
                inverse_lifts+=1
    return dict(integer_vectors=seen,inverse_lifts=inverse_lifts,twist_matrix=[[3,-1],[1,2]],
                marked_readout_numerator=[5,-1],denominator=7,
                source_norm_shells_after_twist=[0,7,21])


if __name__=='__main__':
    print(json.dumps(dict(status='PASS',native=native_checks(),cycles=cycle_checks(),
                         aliquot=aliquot_checks(),norm_seven=norm_seven_checks()),indent=2,sort_keys=True))
