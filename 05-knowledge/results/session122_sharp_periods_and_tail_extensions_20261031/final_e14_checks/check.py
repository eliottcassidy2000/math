from fractions import Fraction as F
from math import gcd
import json

def need(x,msg):
    if not x:raise RuntimeError(msg)
lam=F(762939453,1220703125);alpha=F(134277301,1220703125)
count=nonone=0
# Direct integer affine coefficients, no potential/corridor/source import.
def visit(P,Q,B,A,S,changed,depth):
    global count,nonone
    if depth==8:return
    for a in range(1,5):
        p=3*P;q=Q*2**a;b=3*B+Q;clock=A+a
        g=gcd(p-q,b);h=0
        while g and g%5==0:h+=1;g//=5
        score=S+min(h,14);nn=changed or a!=1
        need(F(score)<=lam*clock-(alpha if nn else 0),'literal finite-word inequality')
        count+=1;nonone+=nn
        visit(p,q,b,clock,score,nn,depth+1)
visit(1,1,0,0,0,False,0)
print(json.dumps({'status':'PASS','direct_affine_words':count,'nonone_words':nonone,'universe':'every nonempty word of length<=8 over letters1..4','potential_imports':0,'scope':'finite hostile/implementation control, not replacement for all-word proof'},indent=2,sort_keys=True))
