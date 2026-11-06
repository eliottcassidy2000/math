"""Exact algebra, divisor inverse, and small hostile controls, not an occupancy census."""
from functools import lru_cache
from pathlib import Path
from math import isqrt
import hashlib
import json

ROOT=Path(__file__).resolve().parent
MENU_PATH=ROOT.parent/'session201-compact-central'/'MENUS.json'
MENU=json.loads(MENU_PATH.read_text())['5']
MATRICES={entry['word']:entry['matrix'] for entry in MENU}

def need(condition,message):
    if not condition: raise RuntimeError(message)

@lru_cache(None)
def g(u):
    need(u>=1,'positive source index')
    if u==1:return 1
    if u%2==0:return g(u//2)
    h=(u+1)//2
    return g(h)+h-h%2

def value(word,u,v):
    t=(g(u),u%2,u,1);s=(g(v),v%2,v,1)
    return sum(t[i]*MATRICES[word][i][j]*s[j] for i in range(4) for j in range(4))

@lru_cache(None)
def x_fiber(x):
    """Every actual u>=2 with u+g(u)=x; proved exact finite bounds."""
    if x<3:return ()
    if x%2:
        lo=(2*x+2)//3;hi=x-1;parity=0
    else:
        lo=(x+1)//2;hi=(2*x-3)//3;parity=1
    lo=max(lo,2)
    lo+=(parity-lo)%2
    return tuple(u for u in range(lo,hi+1,2) if u+g(u)==x)

def factor_pairs(n):
    """All ordered positive factor pairs, exact trial division."""
    out=[]
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            e=n//d;out.append((d,e))
            if d!=e:out.append((e,d))
    return sorted(out)

def inverse_factor(n,c):
    """All actual witnesses for c=3/JJZJ or c=5/JJZZJ, at all heights."""
    need(c in (3,5) and n>0,'inverse domain')
    out=set()
    # Both indices even: n=c*x*y, x,y odd source sums.
    if n%c==0:
        for x,y in factor_pairs(n//c):
            if x%2 and y%2:
                out.update((u,v) for u in x_fiber(x) for v in x_fiber(y))
    # One index even: n=x_even*(c*x_odd-1).
    for x,d in factor_pairs(n):
        if x%2==1 and d%(2*c)==2*c-1:
            y=(d+1)//c
            for u in x_fiber(x):
                for v in x_fiber(y):
                    out.add((u,v));out.add((v,u))
    # Both odd: c*n-(c-1)=(c*x-1)*(c*y-1).
    for d,e in factor_pairs(c*n-(c-1)):
        if d%(2*c)==e%(2*c)==2*c-1:
            x=(d+1)//c;y=(e+1)//c
            out.update((u,v) for u in x_fiber(x) for v in x_fiber(y))
    word={3:'JJZJ',5:'JJZZJ'}[c]
    need(all(value(word,u,v)==n for u,v in out),'every divisor witness is actual')
    return sorted(out)

def prime_trial(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def main():
    # Formal coefficient equality, no sampling of free source coordinates.
    for c,word in ((3,'JJZJ'),(5,'JJZZJ')):
        vec=(c,-1,c,0)
        product=[[vec[i]*vec[j] for j in range(4)] for i in range(4)]
        product[1][1]+=c-1
        need(product==[[c*z for z in row] for row in MATRICES[word]],'formal factorization')
    # The complete polynomial rank-one classification.
    dets={};rank_one=[]
    for word,M in MATRICES.items():
        if word=='':continue
        a,c=M[0][0],M[0][2];d,b=M[2][0],M[2][2]
        dets[word]=a*b-c*d
        if a*b==c*d and a*M[2][1]==d*M[0][1] and a*M[1][2]==c*M[1][0]:
            rank_one.append(word)
    need(set(rank_one)=={'JJZJ','JJZZJ'},'only two full polynomial products')
    # Independent table, as opposed to the recursive implementation above.
    N=4096;table=[0,1]
    for u in range(2,N+1):
        h=(u+1)//2
        table.append(table[u//2] if u%2==0 else table[h]+h-h%2)
    for u in range(1,N+1):
        need(g(u)==table[u] and g(u)%2==1 and 1<=g(u)<=u,'g invariants')
        need((g(u)==1)==(u&(u-1)==0),'unit fiber')
        if u%2 and u>=3:
            need(2*g(u)>=u+3,'odd core lower bound')
            need((2*g(u)==u+3)==((u+1)&u==0),'lower equality Mersenne')
            need((g(u)==u)==((u-1)&(u-2)==0),'upper equality binary-plus-one')
    brute={}
    for u in range(2,1025):brute.setdefault(u+table[u],[]).append(u)
    for x in range(3,1025):need(x_fiber(x)==tuple(brute.get(x,[])),'exact x fiber inverse')
    # Cheap exact checks for the complete divisor inverse. Since each c-form
    # dominates uv, all witnesses at n<=1000 lie in 2<=u,v<=500.
    inverse_checks=0
    for c,word in ((3,'JJZJ'),(5,'JJZZJ')):
        direct={}
        for u in range(2,501):
            for v in range(2,1000//u+1):
                n=value(word,u,v)
                if n<=1000:direct.setdefault(n,[]).append((u,v))
        targets=sorted(set(direct)|{1,3,5,7,21,97,169,997,999})
        for n in targets:
            need(inverse_factor(n,c)==sorted(direct.get(n,[])),'all-height divisor inverse control')
            inverse_checks+=1
    # Prime obstruction, deliberately paired with a positive witness elsewhere.
    hostile=2135743
    prime_triple=[hostile,3*hostile-2,5*hostile-4]
    need(all(prime_trial(n) for n in prime_triple),'exact trial division prime triple')
    need(inverse_factor(hostile,3)==inverse_factor(hostile,5)==[],'all-height two-form absence')
    need(value('J',7,191317)==hostile,'hostile is not a twelve-template gap')
    # Odd-core collapse checked on actual source indices, with both-core-one edge.
    for u in range(2,129):
        for v in range(2,129):
            a,b=u,v;r=s=0
            while a%2==0:a//=2;r+=1
            while b%2==0:b//=2;s+=1
            want=g(a)*g(b)+a*b-1 if r+s==0 else g(a)*g(b)+(1<<(r+s))*a*b
            need(value('J',u,v)==want,'J odd-core collapse')
            if a==b==1:need(r+s>=2,'singleton source exclusion')
    return {'status':'PASS','menu_sha256':hashlib.sha256(MENU_PATH.read_bytes()).hexdigest(),
            'formal_core_determinants':dets,'full_factorization_templates':rank_one,
            'g_and_equality_controls_through':N,'x_fibers_through':1024,
            'divisor_inverse_target_checks':inverse_checks,
            'prime_hostile_triple':prime_triple,'hostile_J_witness':[7,191317],
            'uniform_twelve_template_occupancy':'OPEN'}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
