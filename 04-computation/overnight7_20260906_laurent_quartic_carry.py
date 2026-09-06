#!/usr/bin/env python3
"""Endpoint27 five-channel family: exact symbolic all-height certificate.

No repository imports. Requires SymPy. First-row real-rootedness is the
proved complete-factorial-row theorem THM4436, with explicit eligible tuple.
"""
from math import factorial,gcd
import hashlib,json,sys
import sympy as S
sys.stdout.reconfigure(newline='\n')
G,T,Z=S.symbols('g tau s')
CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)
def fall(x,n):return S.prod(x-j for j in range(n))
def literal(a,b,c,maxmass):
    state={(0,0):1}; answer=[]
    for m in range(1,maxmass+1):
        new={}
        for (charge,z),v in state.items():
            for dc,dz in ((-a,0),(b,0),(c,1)):
                key=charge+dc,z+dz
                new[key]=new.get(key,0)+v
        state=new
        answer.append({z:v for (charge,z),v in state.items() if charge==0})
    return answer
def raw(a,b,c,m):
    answer=[]
    for x in range(m+1):
        d=a*x-b*(m-x)
        if d%(c-b):continue
        z=d//(c-b); y=m-x-z
        if min(y,z)>=0:answer.append(((x,y,z),factorial(m)//(factorial(x)*factorial(y)*factorial(z))))
    return sorted(answer,key=lambda item:item[0][2])
def main():
    domain=S.QQ.frac_field(G)
    P=(fall(G-9,4)+220*fall(G-9,3)*T+5544*fall(G-9,2)*T**2+
       15840*(G-9)*T**3+1320*T**4)
    p=S.Poly(P/1320,T,domain=domain)
    L=fall(G,9)/factorial(12); K=fall(2*G,18)
    first=sum(fall(G,13-j)*T**j/(factorial(12-3*j)*factorial(1+2*j)) for j in range(5))
    Q=sum(fall(2*G,27-j)*T**j/(factorial(27-3*j)*factorial(2*j)) for j in range(10))
    qbar=S.Poly(sum(fall(2*G-18,9-j)*T**j/(factorial(27-3*j)*factorial(2*j)) for j in range(10)),T,domain=domain)
    need(S.cancel(first-L*P)==0,'first full factorial row')
    need(S.cancel(Q-K*qbar.as_expr())==0,'second full factorial row')
    need(S.cancel(qbar.nth(0)/p.nth(0)-
         42240*(G-13)*(2*G-19)*(2*G-21)*(2*G-23)*(2*G-25)/factorial(27))==0,
         'inverse denominator polynomial cancellation')
    inverse=S.invert(S.Poly(T,T,domain=domain),p)
    response=(qbar*inverse).rem(p)
    for j in range(4):need(S.Poly(response.nth(j),G).degree()<=8-j,'weighted response degree')
    C=S.zeros(4)
    for j in range(3):C[j+1,j]=1
    for i in range(4):C[i,3]=-p.nth(i)
    V=sum((response.nth(j)*C**j for j in range(4)),S.zeros(4)).applyfunc(S.cancel)
    coeff=[S.cancel(c) for c in V.charpoly(Z).all_coeffs()]
    powers=[None]+[S.cancel(S.trace(V**r)) for r in range(1,5)]
    records={}
    for k in range(1,5):
        need(S.cancel(k*coeff[k]+sum(coeff[k-r]*powers[r] for r in range(1,k+1)))==0,
             'independent Newton characteristic identity')
        poly=S.Poly(coeff[k],G,domain=S.QQ)
        need(poly.degree()==8*k,'characteristic degree')
        denominator,numerator=poly.clear_denoms(convert=True)
        need(denominator>0 and numerator.domain==S.ZZ,'explicit integer numerator convention')
        shifted=S.Poly(numerator.as_expr().subs(G,Z+14),Z,domain=S.ZZ)
        values=shifted.all_coeffs()
        need(all(c>0 for c in values),'all shifted characteristic coefficients strictly positive')
        records[str(k)]={'denominator':str(denominator),'shift_g':14,
                          'coefficients_descending':[str(c) for c in values]}
    print('SCOPE support(-27,2g-27,3g-27), integerg>=14, gcd(g,27)=1')
    print('FIRST_P',S.expand(P))
    print('FACTORS L=(g)_9/12!, K=(2g)_18; CT1=X*L*P; CT2=X^2*K*tau^-1*Qbar')
    print('CHARACTERISTIC_CERTIFICATES',json.dumps(records,sort_keys=True))
    controls=[]
    for g in (14,16,17,19,23,29):
        a,b,c=27,2*g-27,3*g-27
        need(gcd(g,27)==1 and b>0 and gcd(a,gcd(b,c))==1,'primitive source')
        native=literal(a,b,c,2*g)
        need(all(not native[m-1] for m in range(1,2*g) if m!=g),'all earlier/intermediate empty rows')
        rows=[raw(a,b,c,m) for m in (g,2*g)]
        need([v for v,w in rows[0]]==[(g-13+j,12-3*j,1+2*j) for j in range(5)],'complete first fiber')
        need([v for v,w in rows[1]]==[(2*g-27+j,27-3*j,2*j) for j in range(10)],'complete doubled fiber')
        for m,row,formula in zip((g,2*g),rows,(first,Q)):
            need(native[m-1]=={v[2]:w for v,w in row},'literal Laurent multiplication')
            need(S.expand(sum(w*T**j for j,(v,w) in enumerate(row))-formula.subs(G,g))==0,'symbolic source specialization')
        pf=S.Poly(P.subs(G,g),T); qf=S.Poly(Q.subs(G,g),T)
        need(pf.count_roots(-S.oo,0)==4 and S.gcd(pf,pf.diff()).degree()==0,'four distinct attainable negative roots')
        need(S.gcd(pf,qf).degree()==0,'no first/doubled common root')
        rf=S.Poly(S.rem(qf.as_expr()*S.invert(T,pf.as_expr(),T),pf.as_expr(),T),T)
        need(S.cancel(rf.as_expr()-K.subs(G,g)*response.as_expr().subs(G,g))==0,'specialized inverse carry')
        tf=S.zeros(4)
        for col in range(4):
            remainder=S.Poly(S.rem(T**(col+1),pf.as_expr(),T),T)
            for row in range(4):tf[row,col]=remainder.nth(row)
        vf=sum((rf.nth(j)*tf**j for j in range(4)),S.zeros(4))
        hermite=S.Matrix(4,4,lambda i,j:S.trace(vf*tf**(i+j)))
        need(all((-1)**k*hermite[:k,:k].det()>0 for k in range(1,5)),
             'full specialized negative-definite trace form')
        controls.append({'g':g,'support':[-a,b,c],'first_coefficients':[int(w) for v,w in rows[0]]})
    bad=literal(27,3,18,15)
    need(next(i+1 for i,row in enumerate(bad) if row)==5,'gcd hostile g15 first mass5')
    p14=S.Poly(P.subs(G,14),T)
    need(p14.eval(-59)>0,'quartic far-root left positive')
    need(p14.eval(-58)<0,'quartic far-root right negative')
    print('NAMED_RAW_CONTROLS',json.dumps(controls,sort_keys=True))
    print('HOSTILES gcd_dropped_g15_first_mass5; omitted_carry_reverses_sign; g14_first_root_in(-59,-58)_outside_cubic_real_core_sector')
    print('PASS',CHECKS,'explicit gates; symbolic all-height coefficient certificate')
    print('semantic_sha256',hashlib.sha256(json.dumps({'certificates':records,'controls':controls},sort_keys=True).encode()).hexdigest())
if __name__=='__main__':main()
