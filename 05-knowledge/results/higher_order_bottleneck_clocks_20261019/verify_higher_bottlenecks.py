#!/usr/bin/env python3
"""Exact algebra guards plus directed-decimal escape controls; no dependencies."""
import json
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from pathlib import Path

def need(test, message):
    if not test:
        raise RuntimeError(message)

a = F(1,8)
checks = 0
for r in range(2,102,2):
    need(r*F(5,32)**(r-1)<=F(5,16), 'first derivative')
    need(r*(r-1)*F(5,32)**(r-2)<=2, 'second derivative')
    need(r*r*a**(r-2)<=4, 'square derivative')
    need(2*a**r<=F(1,32) and 2*a**r<=a/4, 'step bound')
    checks += 4
clock_constant = F(1,3)+(F(8,3)+F(1,80))*F(16,11)+F(1,2)+(2+F(1,64))*F(16,11)
need(clock_constant == F(5057,660) and clock_constant<8,'clock constant')
need(4*F(5,32)**3+6*F(5,32)**5<F(5,16),'quartic derivative')
need(12*F(5,32)**2+30*F(5,32)**4<2,'quartic second derivative')
need((4+6*a*a)**2*a*a<4,'quartic square derivative')
need(2*a**4+a**6<=a/4,'quartic step')
checks += 5

def even_power(y, r):
    q=y*y
    result=Decimal(1)
    for _ in range(r//2):
        result*=q
    return result

def count(r, power, lam):
    """Every iterate enclosed using monotonicity of the full step on J."""
    with localcontext() as ctx:
        ctx.prec=65
        eps=Decimal(1)/Decimal(1<<power)
        l=u=-Decimal(1)/8
        target=Decimal(1)/8
        lam=Decimal(lam)
        n=0
        while True:
            need(-target<=l<=u<target,'pre-exit enclosure left monotonicity domain; raise precision')
            ctx.rounding=ROUND_FLOOR
            l=l+eps+even_power(l,r)+lam*even_power(l,6)
            ctx.rounding=ROUND_CEILING
            u=u+eps+even_power(u,r)+lam*even_power(u,6)
            n+=1
            need(-target<=l<=u<=Decimal(5)/32,'enclosure left monotonicity domain; raise precision')
            if l>=target:
                return n
            need(u<target,'ambiguous exit enclosure; raise precision')
            need(n<1000000,'control budget exceeded')

def atan_series(x):
    t=x; total=t; k=1
    while True:
        t *= -x*x
        add=t/(2*k+1)
        total+=add
        if abs(add)<Decimal('1e-65'):
            return total
        k+=1

rows=[]
with localcontext() as ctx:
    ctx.prec=60
    pi=16*atan_series(Decimal(1)/5)-4*atan_series(Decimal(1)/239)
    for r,power,lam in [(2,10,0),(2,14,0),(4,12,0),(4,16,0),(4,20,0),(4,16,'0.5'),(4,20,'0.5'),(4,16,1),(4,20,1),(6,18,0)]:
        n=count(r,power,lam)
        eps=Decimal(1)/Decimal(1<<power)
        if r==2: leading=pi/eps.sqrt()
        elif r==4: leading=pi/Decimal(2).sqrt()/(eps.sqrt().sqrt()**3)
        else: leading=2*pi/3/(eps**(Decimal(5)/6))
        correction=Decimal(0)
        if r==4:
            correction=3*pi*Decimal(2).sqrt()/8*Decimal(lam)/eps.sqrt().sqrt()
        remainder=Decimal(n)-leading+correction
        if r==4 and Decimal(lam):
            bound=1+3*(Decimal(5)/4).ln()+18*Decimal(1)/8+Decimal(2)/(3*(Decimal(1)/8)**3)+16+Decimal(1)/4
        else:
            bound=2*(Decimal(1)/8)**(1-r)/(r-1)+Decimal(r)/2*(Decimal(5)/4).ln()+18*Decimal(1)/8
        need(abs(remainder)<bound,'explicit theorem bound on certified count')
        checks+=1
        rows.append({'r':r,'epsilon':'2^(-%d)'%power,'lambda':str(lam),'certified_escape_count':n,'two_term_remainder':format(remainder,'.12f'),'explicit_bound':format(bound,'.12f')})
result={'status':'PASS','algebra_checks':checks,'clock_constant':str(clock_constant),'count_method':'65-digit directed Decimal products/divisions; powers use nonnegative repeated products; ambiguous threshold comparisons are rejected','real_constant_comparisons':'60-digit numerical controls, not substitutes for the written analytic proof','rows':rows}
print(json.dumps(result,indent=2,sort_keys=True))
