"""Local exploratory identity scout. No target census or new occupancy claim."""
from pathlib import Path
from functools import lru_cache
import json

@lru_cache(None)
def g(n):
    if n==1:return 1
    if n%2==0:return g(n//2)
    h=(n+1)//2;return g(h)+h-h%2

def rev(k,A):
    R=0
    for _ in range(k):R=2*R+(A&1);A//=2
    return R

def Q(A):
    b=[int(x) for x in bin(A)[2:][::-1]]
    return sum(b[i]*b[j]*((1<<(j-i-1))+int(j==i+1)) for i in range(len(b)) for j in range(i+1,len(b)))

def need(c,m):
    if not c:raise RuntimeError(m)

menu=json.loads((Path(__file__).parent/'session201-compact-central/MENUS.json').read_text())['5']
words=0;row_checks=0
for k in range(11):
    B=1<<k
    for A in range(B):
        R=rev(k,A);u=2*B-R;v=2*B-A
        need(g(u)==2*A+1-Q(A),'exact source digit polynomial')
        need(Q(A)==Q(R),'reversal invariance')
        need(g(u)-2*u==g(v)-2*v,'actual reversal difference')
        if k:
            C=B-1-A
            need(Q(C)==Q(A)+B-2-A-R+(A&1)+((A>>(k-1))&1),'complement identity')
        if k<=6:
            p=0 if k==0 else (A>>(k-1))&1
            for left in (2,3,7,11):
                s=(g(left),left%2,left,1)
                for entry in menu:
                    M=entry['matrix'];row=[sum(s[i]*M[i][j] for i in range(4)) for j in range(4)]
                    alpha,gamma,beta,delta=row
                    want=alpha*(2*A+1-Q(A))+beta*(2*B-R)+gamma*p+delta
                    actual=sum(s[i]*M[i][j]*t for i in range(4) for j,t in enumerate((g(u),u%2,u,1)))
                    need(want==actual,'exact weighted binary-quadratic row equation')
                    row_checks+=1
        words+=1
for m in range(21):
    B=4**m
    need((4*B+2)%3==0 and (4*B+6*m+5)%9==0 and (28*B+6*m+17)%9==0,
         'periodic formulas are exact integers')
    v=(4*B+2)//3
    need(g(v)==(4*B+6*m+5)//9,'periodic backup g identity')
    need(2*v+g(v)==(28*B+6*m+17)//9,'periodic backup J identity')
nonempty_AB={e['word']:e['matrix'][0][0] for e in menu if e['word']}
need(all(a>0 and a%2 for a in nonempty_AB.values()),'uncancelled positive odd Q squared coefficient')
print(json.dumps({'status':'PASS','source_words':words,'row_identity_controls':row_checks,
                  'nonempty_Q_squared_coefficients':nonempty_AB,
                  'stopping_point':'No interval or congruence-class constructor; periodic family not advanced to a packet',
                  'universal_occupancy':'OPEN'},indent=2,sort_keys=True))
