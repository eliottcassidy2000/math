"""Exact joint-response recurrence and finite-field clock, with literal controls."""
from bulk import *
from return_program import bin_values,golden,field
from sympy import Matrix,diag,eye,symbols,Poly,rem,factor


def certify(B,a,b,coefficients):
    powers=[b]
    for _ in coefficients[1:]:powers.append(B*powers[-1])
    residual=Matrix.zeros(B.rows,1)
    for c,p in zip(coefficients,powers):residual+=c*p
    row=a
    for _ in range(B.rows):
        need((row*residual)[0]==0,'reachable-space polynomial annihilation')
        row=row*B
    return list(residual)


def recurrence(seq,coefficient,count):
    out=list(seq)
    while len(out)<count:out.append(sum(c*out[-len(coefficient)+i] for i,c in enumerate(coefficient)))
    return out


def tick(v):
    h0,h1,h2,s0,s1,s2,e=v
    return (h1,h2,(h0+h1+h2)%7,s1,s2,(2*s0+s1+2*s2)%7,e)


def mul(a,b):
    x,y=a;u,v=b
    return ((x*u+2*y*v)%7,(x*v+y*u+4*y*v)%7)


def power(a,n):
    out=(1,0)
    while n:
        if n&1:out=mul(out,a)
        a=mul(a,a);n//=2
    return out


def spectral(v):
    h0,h1,h2=v[:3];r0=(h1-3*h0)%7;r1=(h2-3*h1)%7
    return ((r1+4*r0)%7,2*r0%7)


def mpower(A,n,modulus):
    d=len(A);out=[[int(i==j) for j in range(d)] for i in range(d)]
    def times(A,B):return [[sum(A[i][k]*B[k][j] for k in range(d))%modulus for j in range(d)] for i in range(d)]
    while n:
        if n&1:out=times(out,A)
        A=times(A,A);n//=2
    return out


def on(A,v,modulus):return [sum(a*b for a,b in zip(row,v))%modulus for row in A]


def run():
    records=[];x=symbols('x')
    D=[[0,1,0,0,0,0,0],[0,0,1,0,0,0,0],[1,1,1,0,0,0,0],
       [0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,2,1,2,2],[0,0,0,0,0,0,1]]
    need(int(Matrix(D).det())==2,'one binary inverse guard')
    need(mpower(D,2736,7)==[[int(i==j) for j in range(7)] for i in range(7)],'entire ambient clock closes modulo7')
    need(Poly(x**3-2*x**2-x-2,x,modulus=7).is_irreducible,'actual cubic is irreducible')
    f=Poly(x**3-2*x**2-x-2,x,modulus=7)
    def fpow(n):
        a=Poly(x,x,modulus=7);v=Poly(1,x,modulus=7)
        while n:
            if n&1:v=rem(v*a,f)
            a=rem(a*a,f);n//=2
        return v
    need(fpow(171)==Poly(1,x,modulus=7),'F343 root order divides171')
    need(all(fpow(171//p)!=Poly(1,x,modulus=7) for p in (3,19)),'F343 root order exactly171')
    beta=(2,2)
    need(power(beta,48)==(1,0) and power(beta,24)!=(1,0) and power(beta,16)!=(1,0),'F49 spectral root order48')
    need(power(beta,7)==(3,5),'explicit Frobenius conjugate')
    for kind in ('straight','revisit'):
        ds=[blueprint(3,kind,p) for p in (None,'q2','q5')]
        h,a,b=ds
        need(a['states']==b['states'] and a['transfer']==b['transfer'] and a['initial']==b['initial'],'common path frame for Q4 probes')
        B=Matrix(a['transfer']);alpha=Matrix([a['initial']]);terminal=(Matrix(b['terminal'])-Matrix(a['terminal']))/2
        certS=certify(B,alpha,terminal,[-2,-1,-2,1])
        H=Matrix(h['transfer']);AA=Matrix([h['initial']+a['initial']]);BB=diag(H,B)
        # Use exact rational values, including when a terminal coefficient is odd.
        from sympy import Rational
        EE=Matrix([5*z for z in h['terminal']]+[Rational(7*y-9*z,2) for y,z in zip(b['terminal'],a['terminal'])])
        certE=certify(BB,AA,EE,[-1,1])
        vv=[values(d,12) for d in ds];hh=vv[0];ss=[];ee=[];observations=[]
        for n,(hh0,qa,qb) in enumerate(zip(*vv),4):
            need((qb-qa)%2==0 and (7*qb-9*qa+10*hh0)%2==0,'integral probe guards')
            s=(qb-qa)//2;e=(7*qb-9*qa+10*hh0)//2
            need((s-e)%2==0,'native golden parity')
            v=(hh0,(s-e)//2,(s+e)//2);ss.append(s);ee.append(e)
            observations.append(dict(n=n,observer=v,phase=field(v),all_Q4_phases={q:field(golden.step(v,golden.BANK[q])) for q in ('q0','q1','q2','q4','q5')}))
            if n in (4,5):
                T,ad=target(3,n,kind);P,E,S=bin_values(T,ad.index((0,2)),2)
                need(v==(P[1],(S[2]-E[1])//2,(S[2]+E[1])//2),'independent native joint response');tables.cache_clear()
        need(ss==recurrence(ss[:3],[2,1,2],len(ss)) and len(set(ee))==1,'S cubic and constant E')
        detS=int(Matrix([[ss[i+j] for j in range(3)] for i in range(3)]).det());need(detS!=0,'minimal S cubic')
        # Any common scalar annihilator of the three coordinates contains
        # both coprime cubics and x-1, hence order exactly seven.
        state=tuple(z%7 for z in hh[:3]+ss[:3]+ee[:1]);first=state;seen=set()
        while state not in seen:
            seen.add(state)
            need(spectral(tick(state))==mul(beta,spectral(state)),'physical bulk growth transports spectral F49')
            state=tick(state)
        need(state==first and len(seen)==2736,'exact actual joint clock period')
        intstate=hh[:3]+[(s-ee[0])//2 for s in ss[:3]]+[ee[0]]
        moved=on(mpower(D,2736,49),intstate,49)
        lift=[((a-b)%49)//7 for a,b in zip(moved,intstate)]
        need(any(lift) and all((a-b)%7==0 for a,b in zip(moved,intstate)),'precision lift has exact valuation1')
        for a in (1,2,3):
            modulus=7**a;period=2736*7**(a-1)
            need(on(mpower(D,period,modulus),intstate,modulus)==[v%modulus for v in intstate],'exact lifted period returns')
            if a>1:need(on(mpower(D,period//7,modulus),intstate,modulus)!=[v%modulus for v in intstate],'no shorter precision period')
        longh=recurrence(hh[:3],[1,1,1],997);longs=recurrence(ss[:3],[2,1,2],997)
        long=[]
        for n in (20,100,1000):
            j=n-4;v=(longh[j],(longs[j]-ee[0])//2,(longs[j]+ee[0])//2)
            long.append(dict(n=n,vertices=n+8,H=str(v[0]),phase=field(v)))
        records.append(dict(kind=kind,probes=ds,S_certificate=[str(z) for z in certS],E_certificate=[str(z) for z in certE],
                            observations=observations,H_initial=hh[:3],S_initial=ss[:3],E=ee[0],S_hankel=detS,
                            clock_period=len(seen),clock_initial=first,integral_initial=intstate,precision_lift_mod7=lift,long=long))
        print(kind,'joint order7; modulo7 period',len(seen),'H/S matrices4/17',flush=True)
    out=dict(records=records,field_orders=dict(F49=48,F343=171),frobenius_beta=power(beta,7),integral_transfer=D)
    (HERE/'phase_bulk.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':run()
