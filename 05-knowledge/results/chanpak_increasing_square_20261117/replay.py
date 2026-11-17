from base_replay import validate,hams,prime
from itertools import combinations
def transitive(n):
    return tuple((1<<n)-(2<<i) for i in range(n))
def substitute(a,B):
    n=len(a)-1
    return tuple((v&((1<<n)-1))|((((1<<len(B))-1)<<n) if v>>n&1 else 0)
                 for v in a[:-1]) + tuple(a[n]|(v<<n) for v in B)
def read(a):
    n=len(a)-1
    if sum(v>>n&1 for v in a)!=1:raise ValueError("not one entry")
    validate(a)
    x,y=(hams(substitute(a,transitive(k)))[-1] for k in (1,2))
    if (y-x)%2:raise ValueError(a)
    alpha=(y-x)//2
    return alpha,x-2*alpha
R=(318,8,2,20,38,78,159,63,254)
flips=((1,6),(2,6),(2,7),(3,6),(5,7),(3,7))
states=[R]
for i,j in flips:
    a=list(states[-1]);a[i]^=1<<j;a[j]^=1<<i;states.append(tuple(a))
coeff=[read(a) for a in states]
if coeff!=[(68,5),(143,19),(216,25),(366,61),(417,57),(437,47),(476,35)]:raise ValueError(coeff)
def checks(rows):
    if any(b<=0 for a,b in rows):raise ValueError(rows)
    for (a,b),(c,d) in zip(rows,rows[1:]):
        if c<=a or 2*(c-a)+d-b<=0:raise ValueError(rows)
checks(coeff)
B=(30,204,200,160,142,87,153,33)
C=(30,204,72,160,142,87,153,37)
if tuple(x^y for x,y in zip(B,C))!=(0,0,128,0,0,0,0,4):raise ValueError("base flip")
for k,(a,h) in enumerate(zip((B,C),(485,491))):
    validate(a)
    z=hams(a);total=z[-1];P=sum(z[s]*z[len(z)-1-s] for s in range(len(z)))
    if (total,P)!=(h,6968):raise ValueError((total,P))
    if not prime(68*P+5*h):raise ValueError("prime")
    side=[x*P+y*h for x,y in coeff]
    for i in (0,3,6):
        if hams(substitute(states[i],a))[-1]!=side[i]:raise ValueError("raw replay")
    print("six-arc track",*side)
# Independent induced construction via two vertex insertions.
r=(38,20,32,7,77,26,47)
def add(a,tag):
    return tuple(row|((1-(tag>>i&1))<<len(a)) for i,row in enumerate(a))+(tag,)
def move_last(a,j):
    k=len(a)-1;order=tuple(i for i in range(len(a)) if i!=j)+(j,)
    return tuple(sum((a[i]>>x&1)<<z for z,x in enumerate(order)) for i in order)
q=add(r,3); d=add(q,1<<7)
tr=[read(move_last(a,6)) for a in (r,q,d)]
if tr!=[(34,5),(117,13),(238,35)]:raise ValueError(tr)
checks(tr)
for h in (485,491):
    print("two-insertion track",*(a*2*6968+b*h for a,b in tr))
