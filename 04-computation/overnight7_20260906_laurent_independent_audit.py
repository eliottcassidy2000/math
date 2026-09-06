"""Independent raw-source/resultant and grouped-Hermite audits.

No producer code imports. Complete polynomial degree bounds in the audit
report justify interpolation; finite grouped signs are not generalized.
"""
from pathlib import Path
from math import factorial,comb,gcd,prod
from fractions import Fraction
import json,hashlib,sys
import sympy as S
T,X,G=S.symbols('t x g')
GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)
def raw(a,g,m):
    ans=[]
    for z in range(m+1):
        rem=a*m-3*g*z
        if rem<0 or rem%(2*g):continue
        y=rem//(2*g);x=m-y-z
        if x>=0:ans.append(((x,y,z),factorial(m)//(factorial(x)*factorial(y)*factorial(z))))
    return ans
def falling(x,k):return S.prod(x-j for j in range(k))
def monic_first(h):
    return S.Poly(sum(falling(G-2*h-1,h-j)*S.Rational(factorial(2*h+1),
         factorial(3*h-3*j)*factorial(1+2*j))*T**j for j in range(h+1)),T)

def family():
    base=Path(__file__).parent;file=base/'overnight7_20260906_laurent_quartic_carry.out'
    if not file.exists():file=base.parent/'05-knowledge/results'/file.name
    data=file.read_text(encoding='utf-8')
    rec=json.loads(next(s for s in data.splitlines() if s.startswith('CHARACTERISTIC_CERTIFICATES ')).split(' ',1)[1])
    polys={};samples={k:[] for k in range(1,5)}
    for k in range(1,5):
        row=rec[str(k)];co=list(map(S.Rational,row['coefficients_descending']));den=S.Integer(row['denominator'])
        need(row['shift_g']==14 and len(co)==8*k+1 and den>0 and all(c>0 for c in co),'full positive coefficient record')
        polys[k]=S.Poly.from_list(co,G).as_expr().subs(G,G-14)/den
    psym=monic_first(4)
    ratio=S.cancel(falling(2*G-18,9)/factorial(27)/psym.nth(0))
    need(not S.denom(ratio).has(G) and S.Poly(ratio,G).degree()==5,'inverse denominator cancellation')
    for g in range(14,47):
        one,two=raw(27,g,g),raw(27,g,2*g)
        need([v for v,c in one]==[(g-13+j,12-3*j,1+2*j) for j in range(5)],'complete first raw fiber')
        need([v for v,c in two]==[(2*g-27+j,27-3*j,2*j) for j in range(10)],'complete doubled raw fiber')
        p=S.Poly(sum(c*T**j for j,(v,c) in enumerate(one)),T).monic()
        need(S.expand(p.as_expr()-psym.as_expr().subs(G,g))==0,'factorial-normalized first row')
        K=prod(range(2*g-17,2*g+1))
        q=sum(S.Rational(c,K)*T**j for j,(v,c) in enumerate(two))
        characteristic=S.Poly(S.resultant(p.as_expr(),T*X-q,T),X).monic()
        need(characteristic.degree()==4,'four complete response values')
        for k in range(1,5):
            value=characteristic.nth(4-k);samples[k].append((g,value))
            need(value==polys[k].subs(G,g),'exact characteristic specialization')
        need(p.count_roots(-S.oo,0)==4 and S.gcd(p,p.diff()).degree()==0,'simple negative first roots')
        need(characteristic.count_roots(0,S.oo)==0 and characteristic.nth(0)>0,'negative response spectrum')
        if gcd(g,27)==1:need(all(not raw(27,g,m) for m in range(1,2*g) if m!=g),'first-return interpretation')
    for k in range(1,5):
        interpolated=S.interpolate(samples[k],G)
        need(S.Poly(interpolated,G).degree()==8*k and S.expand(interpolated-polys[k])==0,'degree-bounded identity')
    need(next(m for m in range(1,16) if raw(27,15,m))==5,'gcd-dropped hostile')
    print('QUARTIC complete raw fibers/resultant interpolation at33 parameters: all84 coefficients PASS')
    print('QUARTIC record_sha256',hashlib.sha256(json.dumps(rec,sort_keys=True,separators=(',',':')).encode()).hexdigest())

def conv(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,0)+x*y
    return {i:x for i,x in out.items() if x}
def add(*ds):
    out={}
    for d in ds:
        for i,x in d.items():out[i]=out.get(i,0)+x
    return {i:x for i,x in out.items() if x}
def shift(a,k,c=1):return {i+k:c*x for i,x in a.items()}
def star(a,b):return {i:x*b[i] for i,x in a.items() if i in b and x*b[i]}
def beta(x,y):return {j:comb(x+y-2*j,x+j) for j in range(-x,y//3+1)}
def expr(d):return sum(c*T**j for j,c in d.items())
def eval_at(d,t):return sum(S.Rational(c)*t**j for j,c in d.items())

def groups():
    controls=[(1,5),(1,7),(1,11),(2,8),(2,11),(2,16),
              (3,11),(3,13),(3,17),(4,14),(4,16),(4,19)]
    for h,g in controls:
        x=g-3*h-1;y=3*h
        O={j:comb(g,2*j+1) for j in range((g-1)//2+1)}
        E={j:comb(g,2*j) for j in range(g//2+1)}
        B,C,D=beta(x,y),beta(x,y-1),beta(x,y-2)
        oo,ee,bb,cd=conv(O,O),shift(conv(E,E),-1),conv(B,B),shift(conv(C,D),1)
        A2={j:comb(2*g,2+2*j) for j in range(-1,g)}
        B2=beta(2*x,2*y)
        need(add(oo,ee)==A2,'alpha exact parity split')
        need(add(bb,shift(cd,0,2))==B2,'beta complete hit and skipped-level split')
        p=S.Poly(expr(star(O,B)),T)
        blocks=[star(ee,bb),shift(star(oo,cd),0,2),shift(star(ee,cd),0,2)]
        actual=star(A2,B2);virtual=star(oo,bb)
        need(add(virtual,*blocks)==actual,'all three grouped responses')
        rawrow={j:factorial(2*g)//(factorial(2*x+j)*factorial(2*y-3*j)*factorial(2+2*j)) for j in range(-1,2*h+1)}
        need(actual==rawrow,'full anchored multinomial row')
        # Direct companion from polynomial division; no published matrix.
        mat=S.zeros(h)
        for col in range(h):
            r=S.Poly(S.rem(T**(col+1),p.as_expr(),T),T)
            for row in range(h):mat[row,col]=r.nth(row)
        need(p.count_roots(-S.oo,0)==h,'complete finite first-root bank')
        inv=S.invert(T,p.as_expr(),T)
        for block in blocks:
            # Every actual group has minimum exponent at least minus one.
            numerator=S.Poly(S.expand(T*expr(block)),T)
            rem=S.Poly(S.rem(numerator.as_expr()*inv,p.as_expr(),T),T)
            op=sum((rem.nth(j)*mat**j for j in range(h)),S.zeros(h))
            gram=S.Matrix(h,h,lambda i,j:S.trace(op*mat**(i+j)))
            for k in range(1,h+1):need((-1)**k*gram[:k,:k].det()>0,'independent Hermite group negativity')
        N=x+y
        need(S.expand(g*expr(E)-(1+(g-1)*T)*expr(O)-2*T*(1-T)*S.diff(expr(O),T))==0,'alpha contiguous identity')
        need(S.expand(N*expr(C)-2*T*S.diff(expr(C),T)-y*expr(B)+3*T*S.diff(expr(B),T))==0,'first beta Euler identity')
        need(S.expand((N-1)*expr(D)-2*T*S.diff(expr(D),T)-(y-1)*expr(C)+3*T*S.diff(expr(C),T))==0,'second beta Euler identity')
        if (h,g)==(1,5):
            need([eval_at(star(U,V),S.Integer(-2)) for U,V in ((O,B),(O,C),(O,D),(E,B))]==[0,15,10,-16],'contiguous change loses zero constraint')
    for h in range(1,5):
        p0=monic_first(h).nth(0)
        ratio=S.cancel(falling(2*G-4*h-2,2*h+1)/factorial(6*h+3)/p0)
        proposed=S.Rational(2**(h+1)*factorial(3*h),factorial(2*h+1)*factorial(6*h+3))*(G-3*h-1)*S.prod(2*G-4*h-3-2*r for r in range(h))
        need(S.expand(ratio-proposed)==0,'general cancellation formula controls')
    print('TRANSPORT complete Laurent identities at12 named sources; all90 grouped root signs via Hermite forms PASS')
    print('TRANSPORT all-h identities and degree bounds use analytic proofs; uniform grouped negativity stays OPEN')

if __name__=='__main__':
    sys.stdout.reconfigure(newline='\n')
    family();groups();print('PASS_OPTIMIZATION_LIVE_GATES',GATES)
