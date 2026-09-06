"""Independent even-ideal referee: residual modular layers and rational minors."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial
import sys
sys.stdout.reconfigure(newline='\n')
G=0
def need(ok,why):
    global G
    G+=1
    if not ok:raise ArithmeticError(why)
def det(a):
    b=[[F(x) for x in row] for row in a];ans=F(1)
    for i in range(len(b)):
        k=next((j for j in range(i,len(b)) if b[j][i]),None)
        if k is None:return 0
        if k!=i:b[i],b[k]=b[k],b[i];ans=-ans
        v=b[i][i];ans*=v
        for j in range(i+1,len(b)):
            t=b[j][i]/v
            for k in range(i+1,len(b)):b[j][k]-=t*b[i][k]
    need(ans.denominator==1,'integer determinant')
    return ans.numerator
def scalar(m,s):
    h=m-s;num=den=1
    for c in range(m,m+2*s):
        for t in range(h):num*=c-t
    for j in range(h,m):
        for t in range(h):den*=(j-t)**2
    need(num%den==0,'integral falling-factorial scalar')
    return num//den
def vp(n):
    if not n:return 10**9
    k=0
    while n%31==0:n//=31;k+=1
    return k
def layers(e,a):
    precision=48*e+5;mod=31**precision
    b=[[comb(c,j)*pow(31,e*(c-j),mod)*pow(x,c-j,mod)%mod for c in range(16,48)] for x in (1,a) for j in range(16)]
    vals=[];level=0
    while b:
        hit=next(((i,j) for i in range(len(b)) for j in range(len(b)) if b[i][j]%31),None)
        if hit is None:
            v=min(vp(x) for row in b for x in row)
            need(v<precision,'visible residual matrix')
            div=31**v;precision-=v;mod//=div;level+=v
            b=[[x//div for x in row] for row in b]
            continue
        i,j=hit;b[0],b[i]=b[i],b[0]
        for row in b:row[0],row[j]=row[j],row[0]
        inv=pow(b[0][0],-1,mod)
        b=[[(row[k]-row[0]*inv*b[0][k])%mod for k in range(1,len(b))] for row in b[1:]]
        vals.append(level)
    return vals
def main():
    literal=0
    for s in range(1,6):
        for m in (s,s+2,s+5):
            for a in (2,-2,4):
                h=m-s
                bank=[[comb(c,j)*x**(c-j) for c in range(m,m+2*s)] for x in (1,a) for j in range(h,m)]
                need(det(bank)==scalar(m,s)*(a*(a-1))**(s*s),'signed even determinant')
                literal+=1
    c28=scalar(16,14);c30=scalar(16,15)
    need(c28==factorial(43)*factorial(42)//(factorial(15)**3*factorial(14)**3),'C28 factorial form')
    need(c30==factorial(45)//factorial(15)**3,'C30 factorial form')
    need((vp(c28),vp(c30))==(2,1),'exact scalar valuations')
    roworders=list(range(16))*2
    for r in (28,30):
        h=16-r//2
        targetrows=2*sum(range(h,16));targetcols=sum(range(16,16+r))
        rb=[[],[]];cb=[[],[]]
        for omitted in combinations(range(32),32-r):
            loss=sum(roworders[i] for i in omitted)-(240-targetrows)
            need(loss>=0,'nonnegative row deficit')
            if loss<2:rb[loss].append(omitted)
            excess=1008-sum(16+i for i in omitted)-targetcols
            need(excess>=0,'nonnegative column excess')
            if excess<2:cb[excess].append(omitted)
        need(len(rb[0])==len(cb[0])==1,'unique minimum even minor')
        if r==28:
            need(tuple(map(len,rb))==(1,4) and tuple(map(len,cb))==(1,1),'all five next-band minors')
            for dr,dc in ((0,1),(1,0)):
                for ro in rb[dr]:
                    for co in cb[dc]:
                        keptj=[roworders[i] for i in range(32) if i not in ro]
                        need(15 not in co,'column31 retained')
                        need(all(j>=1 and comb(31,j)%31==0 for j in keptj),'polynomial zero column31')
            need(targetcols-targetrows==588,'rank28 minimum weight')
        else:need(targetcols-targetrows==675,'rank30 minimum weight')
    controls=[(e,a+31*e) for e in (0,1,2,3) for a in range(2,31)]
    controls += [(5,879),(4,-20),(7,3+31**4),(6,4-31**3)]
    for e,a in controls:
        v=layers(e,a)
        need(sum(v)==768*e,'whole determinant valuation')
        if e==0:need(v==[0]*32,'unimodular depth zero');continue
        k=int(a%31 in (3,11,15,17,21,29))
        need(sum(v[:28])==588*e+2,'E28 full residual matrix')
        need(sum(v[:29])==631*e+1+k,'inherited E29 full residual matrix')
        need(sum(v[:30])==675*e+1,'E30 full residual matrix')
        need(v[28:30]==[43*e-1+k,44*e-k],'two actual factors')
    for e in range(1,51):
        lo=(43*e-1,44*e);hi=(43*e,44*e-1)
        for b in range(44*e+3):
            delta=sum(min(b,x) for x in hi)-sum(min(b,x) for x in lo)
            need(delta==int(43*e<=b<=44*e-1),'exact precision window for pair contribution')
        need(sum(lo)==sum(hi)==87*e-1,'fixed determinant contribution')
        need((hi[0]==hi[1])==(e==1),'equal pair boundary')
    print('rational signed even minors',literal)
    print('contents',c28,c30,'valuations 2 1')
    print('complete rank28 bands 1 minimum 5 next; rank30 unique minimum')
    print('independent residual-layer matrices',len(controls),'all residues depths0..3 and four high/negative controls')
    print('all-depth formulas E28=588e+2 E30=675e+1; factors=(43e-1+kappa,44e-kappa)')
    print('pair-kernel window 43e<=b<=44e-1; no whole-kernel extrapolation')
    print('PASS',G,'always-active gates')
if __name__=='__main__':main()
