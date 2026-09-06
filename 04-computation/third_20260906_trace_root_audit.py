#!/usr/bin/env python3
"""Independent formal quotient traces and resultants for the exact boundary jets."""
from math import factorial
from sympy import symbols, Poly, Rational as F, Matrix, expand, cancel, resultant, diff

Y,T=symbols('y t')
GATES=0

def need(ok,why):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(why)

def falling(a,n):
    r=1
    for j in range(n):
        r*=a-j
    return expand(r)

def main():
    traces=[]
    for h in range(1,6):
        p=sum(F(factorial(2*h+1),factorial(3*h-3*j)*factorial(1+2*j))*falling(Y,h-j)*T**j for j in range(h+1))
        q={e:falling(2*Y,2*h-e)/F(factorial(6*h-3*e)*factorial(2+2*e)) for e in range(-1,2*h+1)}
        ppoly=Poly(p,T)
        carry=cancel(q[-1]/ppoly.nth(0))
        qplus=sum(q[e]*T**e for e in range(2*h+1))
        remainder=Poly(qplus-carry*sum(ppoly.nth(j)*T**(j-1) for j in range(1,h+1)),T).rem(ppoly)
        cols=[Poly(remainder.as_expr()*T**j,T).rem(ppoly) for j in range(h)]
        trace=expand(-sum(cols[j].nth(j) for j in range(h)))
        traces.append(trace)
        dh=F(factorial(2*h-1),factorial(6*h))
        bh=F(2*factorial(2*h),factorial(6*h+3))
        ah=F((-1)**(h-1)*factorial(h-1)*factorial(2*h+1),factorial(3*h))
        if h==1:
            need(trace==F(1,90720)*(884*Y**2+123*Y+1),'height-one exact exception')
        else:
            vh=F(3*h*(3*h-1)*(3*h-2),6)
            loss=F(4*vh,(h-1)*(6*h+1)*(6*h+2)*(6*h+3))
            need(trace.subs(Y,0)==0,'trace vanishing')
            need(diff(trace,Y).subs(Y,0)==h*dh*(1-loss)>0,'exact positive trace slope')
            need(0<loss<F(1,12*(h-1)),'strict one-twelfth loss bound')
        if h<=4:
            # Product of t*q at the roots divided by p(0) is the characteristic constant.
            tq=Poly(T*qplus+q[-1],T).rem(ppoly).as_expr()
            norm=Poly(cancel(resultant(p,tq,T)/ppoly.nth(0)),Y)
            plusnorm=Poly((-1)**h*resultant(p,Poly(qplus,T).rem(ppoly).as_expr(),T),Y)
            need(all(norm.nth(j)==0 for j in range(h-1)),'carried lower jet orders')
            need(norm.nth(h-1)==bh**h/ah,'exact carried norm leading jet')
            need(all(plusnorm.nth(j)==0 for j in range(h)),'deleted-carry lower jet orders')
            need(plusnorm.nth(h)==dh**h,'exact deleted-carry leading jet')
            divisor=1
            for r in range(2,h+1):
                divisor*=(Y+r-h)**(r-1)
            deflated,rem=norm.div(Poly(divisor,Y))
            need(rem.is_zero,'complete forced norm divisor')
            need((1 if deflated.nth(0)>0 else -1)==(-1)**(h*(h-1)//2),'four-periodic deflated boundary sign')
            if h==2:
                hostile=F(1,100000)
                need(trace.subs(Y,hostile)>0 and norm.eval(hostile)<0,'positive trace negative norm hostile')
        need(all(v>=0 for v in Poly(trace,Y).all_coeffs()),'finite trace coefficient positive control')
        print(f'height={h} trace_degree={Poly(trace,Y).degree()} boundary_jet_controls=PASS')
    monics=[Poly(v,Y).monic().as_expr() for v in traces[:3]]
    v1,v2,v3=monics
    columns=[Y*v2,v2,Y**2*v1,Y*v1,v1,v3-Y**2*v2]
    C=Matrix([[Poly(v,Y).nth(j) for v in columns] for j in range(6)])
    mod=lambda v:int(v.p)*pow(int(v.q),-1,101)%101
    modular=C.applyfunc(mod)
    need(int(modular.det())%101==29,'independent rational recurrence rank hostile')
    need(C.det()!=0,'recurrence impossible over R, not only modulo101')
    print(f'PASS gates={GATES}')

if __name__=='__main__':
    main()
