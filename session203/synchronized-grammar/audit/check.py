"""Independent basis evaluation and exact rational observability certificate."""
from fractions import Fraction
from itertools import product
import json

def need(q,why):
    if not q: raise RuntimeError(why)
def step(v,c):
    h,s,e,L,R,d,C,b,T=v
    if c=='O':return [h,s,0,L,2*R,d,2*C,2*b,3*T]
    if c=='Z':return [h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b]
    if c=='J':return [h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C]
    raise RuntimeError('letter')
def continuation(v,w):
    for c in w:v=step(v,c)
    return v
unit=[[int(i==j) for i in range(9)] for j in range(9)]
def observation(w):
    ans=[]
    for v in unit:
        z=continuation(v,w)
        ans.append(z[0]+z[3]-v[0]-v[3])
    return ans

def echelon(rows):
    a=[[Fraction(x) for x in row] for row in rows];p=0;det=Fraction(1)
    for col in range(len(a[0])):
        pivot=next((i for i in range(p,len(a)) if a[i][col]),None)
        if pivot is None:continue
        if pivot!=p:a[p],a[pivot]=a[pivot],a[p];det=-det
        val=a[p][col];det*=val
        a[p]=[x/val for x in a[p]]
        for i in range(len(a)):
            if i!=p:
                factor=a[i][col]
                a[i]=[x-factor*y for x,y in zip(a[i],a[p])]
        p+=1
        if p==len(a):break
    return p,det,a
words=['Z','J','OZ','OJ','ZZ','ZJ','JJ','JOZ']
expected=[[0,0,-1,0,1,0,0,1,0],[0,0,-1,-2,1,0,1,0,0],[0,0,0,0,2,0,0,2,0],[0,0,0,-2,2,0,2,0,0],[0,0,-1,0,3,-2,0,3,0],[0,0,-1,-4,3,-2,2,-1,2],[0,0,0,-2,2,0,1,0,2],[2,0,-1,-2,3,0,3,0,0]]
rows=[observation(w) for w in words]
need(rows==expected,'eight literal rows')
rank,det,_=echelon([[x for i,x in enumerate(row) if i!=1] for row in rows])
need(rank==8 and det==32,'exact determinant')
need(all(row[1]==0 for row in rows),'zero start column')
for c in 'OZJ':
    need(step(unit[1],c)==unit[1],'invariant start difference')
    for i in range(9):
        need(step(unit[i],c)[1]==int(i==1),'start preserved')
all_rows=[];word_checks=0
for length in range(7):
    for w in product('OZJ',repeat=length):
        row=observation(w);all_rows.append(row)
        need(row[1]==0,'all short words annihilate start-only direction')
        need(continuation([2*x for x in unit[1]],w)==[2*x for x in unit[1]],'paired start difference continuation')
        word_checks+=1
need(echelon(all_rows)[0]==8,'full short-word observation rank')
# These synthetic deltas are algebraic hostile controls, not physical sources.
separators=[]
for i in range(9):
    if i==1:continue
    delta=unit[i][:];delta[1]=3*(delta[0]+delta[3])+2
    need(3*(delta[0]+delta[3])-delta[1]+2==0,'initial adjusted equality')
    sep=next((w for w,row in zip(words,rows) if sum(x*y for x,y in zip(row,delta))),None)
    need(sep is not None,'every non-start coordinate detected')
    separators.append([i,sep])
js=[(a,b) for a in range(9) for b in range(9) if (1+3*a)-(2+3*b)==2]
need(js==[(b+1,b) for b in range(8)],'exact admissible shifts')
print(json.dumps(dict(status='PASS',words=words,rows=rows,determinant=32,rank=8,
    all_short_word_checks=word_checks,maximum_word_length=6,algebraic_separators=separators,
    width8_shifts=js,scope='Conditional synchronization equivalence; no physical seed existence or occupancy.'),sort_keys=True,indent=2))
