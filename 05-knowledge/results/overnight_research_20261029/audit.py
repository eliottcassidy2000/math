"""Alternative paths, with no import of the producer under audit."""
from itertools import product,combinations
from collections import Counter
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
GATES=0
def gate(x,label):
    global GATES
    GATES+=1
    if not x:raise RuntimeError(label)

def matchings(ids):
    if not ids:return [(1,())]
    out=[]
    for j in range(1,len(ids)):
        rem=ids[1:j]+ids[j+1:]
        out.extend(((-1)**(j+1)*s,((ids[0],ids[j]),)+pairs) for s,pairs in matchings(rem))
    return out
MATCH={n:matchings(tuple(range(n))) for n in (2,4,6,8)}
def pf(a):
    total=0
    for sign,pairs in MATCH[len(a)]:
        for i,j in pairs:sign*=a[i][j]
        total+=sign
    return total

def det(a):
    a=[[F(x) for x in r] for r in a];n=len(a);d=1
    for j in range(n):
        p=next((p for p in range(j,n) if a[p][j]),None)
        if p is None:return 0
        if p!=j:a[p],a[j]=a[j],a[p];d=-d
        x=a[j][j];d*=x
        for i in range(j+1,n):
            q=a[i][j]/x
            for k in range(j+1,n):a[i][k]-=q*a[j][k]
    return int(d)
def sub(a,rows,cols):return [[a[i][j] for j in cols] for i in rows]

def cut_audit():
    m=[[0,1,1,1],[-1,0,-1,1],[-1,1,0,-1],[-1,-1,1,0]]
    gate(pf(m)==-3,"oriented base sign")
    directions=[(1,)+x for x in product((-1,1),repeat=3)]
    counts=Counter();clock=Counter();seen=set();perm=(0,2,3,1)
    for u,v in combinations(directions,2):
        for choices in product(range(4),repeat=4):
            if len({x//2 for x in choices})<2:continue
            b=[[(1 if k%2 else -1)*x for x in (u if k<2 else v)] for k in choices]
            mask=sum((x==1)<<(4*i+j) for i,row in enumerate(b) for j,x in enumerate(row))
            gate(mask not in seen,"unique alternative cut");seen.add(mask)
            a=[m[i]+b[i] for i in range(4)]+[[-b[j][i] for j in range(4)]+m[i] for i in range(4)]
            p=pf(a);g=0
            for ids in combinations(range(8),6):
                g=gcd(g,pf(sub(a,ids,ids)))
                if g==1:break
            gate(g==1,"six-Pfaffian ideal gives one nontrivial hyperbolic factor")
            counts[abs(p)]+=1
            if all(b[perm[i]][perm[j]]==b[i][j] for i in range(4) for j in range(4)):
                clock[abs(p)]+=1
    gate(len(seen)==6272,"complete rank-two generation")
    expected={1:576,3:768,5:1152,7:576,9:64,11:576,13:1152,15:384,17:576,21:384,27:64}
    gate(counts==expected,"full census agrees")
    gate(clock=={9:4,27:4},"clock census agrees")
    return {"rank_two_cuts":len(seen),"group_histogram":dict(sorted(counts.items())),"clock_histogram":dict(clock)}

def smith_audit():
    cases=0
    for d,e in product(range(1,5),repeat=2):
        for a,b,c,f in product(range(-2,3),repeat=4):
            t=d*e-a*f+b*c
            if t==0:continue
            k=[[0,d,a,b],[-d,0,c,f],[-a,-c,0,e],[-b,-f,-e,0]]
            s=gcd(d,e,a,b,c,f)
            # Ordinary (not Pfaffian) determinantal ideals.
            g2=0
            for rows in combinations(range(4),2):
                for cols in combinations(range(4),2):
                    g2=gcd(g2,k[rows[0]][cols[0]]*k[rows[1]][cols[1]]-k[rows[0]][cols[1]]*k[rows[1]][cols[0]])
            g3=0
            for rows in combinations(range(4),3):
                for cols in combinations(range(4),3):g3=gcd(g3,det(sub(k,rows,cols)))
            gate((g2,g3,det(k))==(s*s,s*abs(t),t*t),"ordinary Smith ideals")
            cases+=1
    return cases

def controller_audit():
    # Check the universal inequalities at the worst valuation, independently.
    for n in range(9,10000,2):
        gate(F(5*(3*n+1),16*n)<=F(35,36),"worst leaf inequality")
        gate(F(8*(3*n+1),25*n)<=F(224,225),"worst branch inequality")
    cores={1:[],3:[1,4],5:[4],7:[1,1,2,3,4]}
    for n,w in cores.items():
        x=n
        for a in w:
            gate((3*x+1)%(1<<a)==0 and ((3*x+1)>>a)%2==1,"core exact valuation")
            x=(3*x+1)>>a
        gate(x==1,"core hits one")
    tests=[]
    for source,qstart in product((9,11,25,27,57,343,1001,2**127-1,2**521-1), (1,2,11,50)):
        x=source;q=qstart;initial=F(source)*F(8,5)**q;v=initial;j=0;first_core=None
        while q:
            if x<9 and first_core is None:first_core=(j,q)
            y=3*x+1;a=0
            while y%2==0:y//=2;a+=1
            qnew=q+(2 if a%3==0 else -1);vnew=F(y)*F(8,5)**qnew
            if x>=9:gate(225*vnew<=224*v,"actual rational potential")
            x,q,v=y,qnew,vnew;j+=1
            gate(j<100000,"diagnostic finite guard")
        tests.append({"bits":source.bit_length(),"queue":qstart,"steps":j})
    gate(F(3,2)**3<F(2)**2,"27 less than 32 potential window")
    gate(not F(3,2)**3<F(2),"critical rule fails this window")
    return {"large_source_queue_cases":len(tests),"max_steps":max(r["steps"] for r in tests)}

def main():
    out={"cut_audit":cut_audit(),"ordinary_smith_cases":smith_audit(),"controller":controller_audit(),"gates":GATES}
    (ROOT/"audit.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))
if __name__=="__main__":main()
