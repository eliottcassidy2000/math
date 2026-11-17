from itertools import combinations
from slot_profiles import end_profiles

def comp(row):
    n=len(row); N=(1<<n)-1
    reach=[(1<<i)|row[i] for i in range(n)]
    for j in range(n):
        for i in range(n):
            if reach[i]>>j&1:reach[i]|=reach[j]
    comp=sorted({sum(1<<k for k in range(n) if reach[i]>>k&1 and reach[k]>>i&1)
                 for i in range(n)})
    order=sorted(comp,key=lambda x: -reach[(x&-x).bit_length()-1].bit_count())
    return tuple(order)
for n in range(1,6):
    arcs=list(combinations(range(n),2))
    states=[]
    for mask in range(1<<len(arcs)):
        a=[0]*n
        for k,(i,j) in enumerate(arcs):
            if mask>>k&1:i,j=j,i
            a[i]|=1<<j
        H,_=end_profiles(a)
        c=comp(a)
        states.append((H[-1],tuple((x,H[x]) for x in c if x.bit_count()>1),c))
    k=0
    for m,(h,c,_) in enumerate(states):
        for j in range(len(arcs)):
            h1,c1,_=states[m^(1<<j)]
            if h1==h:
                if c1!=c:raise ValueError((n,m,j,c,c1))
                k+=1
    print(n,k)
