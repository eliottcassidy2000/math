"""Golden ramified lifts: exact ports/modules and interval spectral receipts."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'pinned_discrepancy_reset_tails_20261029'))
from check_pins import sine_centres, Q, require
from fractions import Fraction as F
from itertools import combinations, product
import numpy as np
import json


def mul(u,v,m,eta=1):
    a,b=u;c,d=v
    return ((a*c+b*d)%m,(a*d+b*c+eta*b*d)%m)


def bank(m,eta=1,p=5):
    points=list(product(range(m),repeat=2))
    left=set(points)-{(0,0)};orbits=[];owner={(0,0):(-1,0)}
    while left:
        first=min(left);o=[];v=first
        while v not in o:
            o.append(v);v=mul((eta*eta+1,eta*(eta*eta+2)),v,m,eta)
        neg=[((-x)%m,(-y)%m) for x,y in o]
        require(not set(o)&set(neg),'antipodal guard')
        i=len(orbits)
        for v in o: owner[v]=(i,1)
        for v in neg: owner[v]=(i,-1)
        orbits.append(o);left-=set(o)|set(neg)
    s=np.array(sine_centres(m),dtype=np.int64)
    reps=np.array([o[0] for o in orbits],dtype=np.int64)
    C=np.zeros((len(orbits),len(orbits)),dtype=np.int64)
    for i,o in enumerate(orbits):
        v=np.array(o,dtype=np.int64)
        phases=(reps@v.T)%m
        C[:,i]=2*s[phases].sum(axis=1)
    return {'m':m,'eta':eta,'p':p,'points':points,'orbits':orbits,'owner':owner,'C':C}


def hierarchy(up,down):
    groups=[[] for _ in down['orbits']]; kernel=[]
    for i,o in enumerate(up['orbits']):
        x,y=o[0];v=(x%down['m'],y%down['m'])
        if v==(0,0): kernel.append(i);continue
        j,orientation=down['owner'][v]
        groups[j].append((i,orientation))
        require(len(o)==5*len(down['orbits'][j]),'orbit length lifts by five')
    require(all(len(g)==5 for g in groups) and len(kernel)==4,'five children and four kernel bits')
    return groups,kernel


def assemble(groups,kernel,parent,choices,kerbits,b):
    eps=np.zeros(b,dtype=np.int64)
    for j,group in enumerate(groups):
        minority=set(choices[j])
        for t,(i,orientation) in enumerate(group):
            eps[i]=orientation*parent[j]*(-1 if t in minority else 1)
    for i,s in zip(kernel,kerbits): eps[i]=s
    return eps


def spectral_interval(data,eps):
    peak=int(np.max(np.abs(data['C']@eps)))
    err=2*(data['m']**2-1)
    return [str(F(max(0,peak-err),Q)),str(F(peak+err,Q))]


def literal_symbol(data,eps):
    return {v:0 if i<0 else int(eps[i])*s for v,(i,s) in data['owner'].items()}


def ideal_constraints(data):
    m=data['m'];eta=data['eta'];delta=((-eta)%m,2);v=(1,0);constraints=[]
    while True:
        v=mul(v,delta,m,eta)
        I={mul(v,x,m,eta) for x in data['points']}
        if len(I)==1: break
        seen=set(I);pairs=set()
        for x in data['points']:
            if x in seen:continue
            coset={((x[0]+h[0])%m,(x[1]+h[1])%m) for h in I}
            i,s=data['owner'][x]
            for y in coset:
                j,t=data['owner'][y]
                if i!=j:pairs.add((min(i,j),max(i,j),s*t))
                else:require(s==t,'no contradictory same-port module constraint')
            seen |= coset
        constraints.append((len(I),sorted(pairs)))
    return constraints


def exact_variance(data):
    m=data['m'];p=data['p'];eta=data['eta'];q=m//p;values={};checked=0
    for orbit in data['orbits']:
        x,y=orbit[0]
        if x%p==0 and y%p==0:continue
        v=1 if (2*x+eta*y)%p==0 else 0
        poly=np.zeros(m,dtype=np.int64)
        for O in data['orbits']:
            coeff=np.zeros(m,dtype=np.int64)
            for s,t in O:
                phase=(x*s+y*t)%m
                coeff[phase]+=1;coeff[-phase % m]-=1
            conv=-np.convolve(coeff,coeff)
            poly+=conv[:m]
            poly[:m-1]+=conv[m:]
        # Reduce exactly by Phi_(p^k)=sum_{j=0}^{p-1} X^(j p^(k-1)).
        for r in range(q):
            c=int(poly[r+(p-1)*q])
            for j in range(p-1):poly[r+j*q]-=c
            poly[r+(p-1)*q]=0
        expected=F((p**(2*v))*(m*m+p),p+1)
        require(expected.denominator==1 and int(poly[0])==expected.numerator and not np.any(poly[1:]),
                'exact primitive Fourier variance, cyclotomic reduction')
        values[str(v)]=int(poly[0]);checked+=1
    return {'modulus':m,'eta':eta,'p':p,'rows':checked,'variance_by_delta_valuation':values}


def exact_orbit_gauss(data):
    m=data['m'];p=data['p'];eta=data['eta'];q=m//p;k=0;power=1
    while power<m:power*=p;k+=1
    def vp(x):
        answer=0
        while answer<k and x%p==0:x//=p;answer+=1
        return answer
    def valuation(x,y):
        return min(2*vp(2*x+eta*y),1+2*vp(y))
    checked=zero=0
    for row in data['orbits']:
        x,y=row[0];v=valuation(x,y)
        if v>=2:continue
        for O in data['orbits']:
            ell=valuation(*O[0]);r=v+ell
            coeff=np.zeros(m,dtype=np.int64)
            for s,t in O:coeff[(x*s+y*t)%m]+=1
            neg=coeff[np.array([(-i)%m for i in range(m)])]
            raw=np.convolve(coeff,neg);poly=raw[:m].copy();poly[:m-1]+=raw[m:]
            for a in range(q):
                c=int(poly[a+(p-1)*q])
                for j in range(p-1):poly[a+j*q]-=c
                poly[a+(p-1)*q]=0
            if r>=2*k:expected=len(O)**2
            elif r%2==0:expected=0;zero+=1
            else:expected=p**(k+(r-1)//2+1-2*((ell+1)//2))
            require(int(poly[0])==expected and not np.any(poly[1:]),
                    'exact orbit Gauss norm via cyclotomic reduction')
            checked+=1
    return {'p':p,'eta':eta,'modulus':m,'exact_orbit_norms':checked,'zero_sums':zero}


def has_ideal_module(eps,constraints):
    return [size for size,pairs in constraints if all(int(eps[i])*int(eps[j])==s for i,j,s in pairs)]


def native_prime(data,eps):
    # Translation permits testing all pair closures containing the original zero.
    points=data['points'];m=data['m'];n=len(points);sym=literal_symbol(data,eps)
    out=[]
    for x in points:
        mask=0
        for j,y in enumerate(points):
            if sym[((y[0]-x[0])%m,(y[1]-x[1])%m)]==1:mask |= 1 << j
        out.append(mask)
    full=(1 << n)-1
    for v in range(1,n):
        closure=1 | (1 << v);todo=1 << v
        while todo and closure!=full:
            low=todo & -todo;todo-=low;i=low.bit_length()-1
            new=(out[i]^out[0]) & (full^closure)
            closure |= new;todo |= new
        if closure!=full:return False,closure.bit_count()
    return True,n


def verify_coarse(up,down,parent,eps,groups):
    sym=literal_symbol(up,eps);m=up['m']
    # Integer character polynomials, not merely numerical eigenvalues.
    for o in down['orbits']:
        u,v=o[0];coeff=[0]*down['m'];lower=[0]*down['m']
        for (x,y),s in sym.items():coeff[(u*x+v*y)%down['m']]+=s
        for (x,y),(i,s) in down['owner'].items():
            if i>=0:lower[(u*x+v*y)%down['m']]+=int(parent[i])*s
        require(coeff==[5*x for x in lower],'all inherited character polynomials scale by five')
    for j,group in enumerate(groups):
        require(sum(int(eps[i])*s for i,s in group)==int(parent[j]),'named-parent balance')


def exhaustive_25(up,down,parent):
    groups,kernel=hierarchy(up,down);b=len(up['orbits']);C=up['C']
    choices=list(combinations(range(5),2));best=None;worst=0;count=0;composite=0
    constraints=ideal_constraints(up)
    for choice in product(choices,repeat=len(groups)):
        signs=np.array([assemble(groups,kernel,parent,choice,ker,b) for ker in product((-1,1),repeat=4)])
        values=np.max(np.abs(signs@C.T),axis=1)
        j=int(np.argmin(values))
        if best is None or values[j]<best[0]:best=(int(values[j]),signs[j].copy())
        worst=max(worst,int(np.max(values)))
        alive=np.zeros(16,dtype=bool)
        for _,pairs in constraints:
            ok=np.ones(16,dtype=bool)
            for i,j,s in pairs:ok &= signs[:,i]*signs[:,j]==s
            alive |= ok
        composite+=int(alive.sum());count+=16
    eps=best[1];verify_coarse(up,down,parent,eps,groups)
    native=native_prime(up,eps)
    require(native[0] and composite==0,'balanced lifts prime by ideal and native checks')
    new_rows=[i for i,o in enumerate(up['orbits']) if any(x%5 for x in o[0])]
    new_peak=int(np.max(np.abs(up['C'][new_rows]@eps)))
    new_upper=F(new_peak+2*(up['m']**2-1),Q)
    inherited_lower=5*F(spectral_interval(down,parent)[0])
    require(new_upper<inherited_lower,'candidate attains the unavoidable inherited norm')
    return eps,{'labelled_balanced_completions':count,'composite_completions':composite,
                'best_norm_interval':spectral_interval(up,eps),'max_peak_numerator':worst,
                'new_mode_upper':str(new_upper),'forced_old_mode_lower':str(inherited_lower),
                'exact_optimum':'25*sqrt(5+2*sqrt(5))',
                'native_pair_closures':len(up['points'])-1,'native_prime':native[0],
                'signs':eps.tolist(),'parent_signs':parent.tolist(),
                'ideal_constraint_counts':[[s,len(p)] for s,p in constraints]}


def search_balanced(up,down,parent):
    groups,kernel=hierarchy(up,down);b=len(up['orbits']);C=up['C'];rng=np.random.default_rng(20261029)
    choices=list(combinations(range(5),2));best=None
    for run in range(12):
        choice=[choices[int(rng.integers(10))] for _ in groups]
        eps=assemble(groups,kernel,parent,choice,rng.choice([-1,1],4),b)
        vals=C@eps;peak=int(np.max(np.abs(vals)))
        for _ in range(150):
            moves=[];deltas=[]
            for group in groups:
                for (i,s),(j,t) in combinations(group,2):
                    if eps[i]*s==eps[j]*t:continue
                    moves.append((i,j));deltas.append(-2*int(eps[i])*C[:,i]-2*int(eps[j])*C[:,j])
            for i in kernel:
                moves.append((i,));deltas.append(-2*int(eps[i])*C[:,i])
            candidate=vals[:,None]+np.array(deltas).T
            costs=np.max(np.abs(candidate),axis=0);j=int(np.argmin(costs))
            if costs[j]>=peak:break
            for i in moves[j]:eps[i]*=-1
            vals=candidate[:,j];peak=int(costs[j])
        if best is None or peak<best[0]:best=(peak,eps.copy())
    eps=best[1];verify_coarse(up,down,parent,eps,groups)
    require(not has_ideal_module(eps,ideal_constraints(up)),'no proper ideal module')
    return eps,{'norm_interval':spectral_interval(up,eps),'signs':eps.tolist(),
                'parent_signs':parent.tolist(),'search_starts':12,'prime_via_ideal_criterion':True}


if __name__=='__main__':
    base=bank(5);up=bank(25);higher=bank(125)
    parent=np.ones(len(base['orbits']),dtype=np.int64)
    eps,small=exhaustive_25(up,base,parent)
    big,bigreceipt=search_balanced(higher,up,eps)
    print(json.dumps({'base_norm_interval':spectral_interval(base,parent),'level25':small,
                      'level125':bigreceipt,'exact_variances':[exact_variance(d) for d in
                          [base,up,higher,bank(13,eta=3,p=13)]],
                      'exact_orbit_gauss':[exact_orbit_gauss(d) for d in
                          [base,up,higher,bank(13,eta=3,p=13),bank(25,eta=4,p=5),
                           bank(169,eta=3,p=13)]]},indent=2,sort_keys=True))
