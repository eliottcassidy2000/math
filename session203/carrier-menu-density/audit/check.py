from fractions import Fraction
from pathlib import Path
import json
SEEDS=[(677,27,1061,49),(759,29,927,29),(1547,101,1691,109),(821,33,1079,33)]
def need(b,s):
    if not b:raise RuntimeError(s)
count=nodes=0
for X0,e0,V0,d in SEEDS:
    A=V0-d;need(A>0 and 0<=e0<=d,'seed inequalities')
    def visit(n,X,e,V,subtraction):
        global count,nodes
        nodes+=1
        need(V==(1<<n)*V0-d*subtraction,'closed-form V from complete word')
        need(V>=(1<<n)*A+d and X>=X0 and 0<=e<=d,'all-prefix invariants')
        if n==16:return
        visit(n+1,X,0,2*V,2*subtraction)
        newX=X+V-e
        need(newX>=X0+(1<<n)*A,'last-Z lower bound')
        count+=1
        visit(n+1,newX,d,2*V-d,2*subtraction+1)
    visit(0,X0,e0,V0,0)
need(count==4*((1<<16)-1),'complete normalized word universe')
s=sum((Fraction(1,V-d) for X,e,V,d in SEEDS),Fraction())
need(s==Fraction(693609731,187977134884),'reciprocal sum')
need(2*s==Fraction(693609731,93988567442) and 4*s==Fraction(693609731,46994283721),'density arithmetic')
need(4*s<Fraction(3,200),'strict odd-relative bound')
for A in [V-d for X,e,V,d in SEEDS]:
    for M in [0,1,A-1,A,A+1,2*A-1,2*A,3*A,1_000_000]:
        length=(M//A).bit_length()
        words=1+((1<<length)-1 if length else 0)
        need(Fraction(words)<=1+Fraction(2*M,A),'geometric word bound, including endpoints')
print(json.dumps({'status':'PASS','last_Z_words':count,'all_prefix_nodes':nodes,'maximum_length':16,'upper_integer_density':str(2*s),'upper_odd_relative_density':str(4*s),'scope':'Only four fixed paired O/Z orbits; no claim for all seeds, J extensions or F.'},indent=2,sort_keys=True))
