"""Independent SCC suffix-merge quotient of the end-state transfer."""
from rooted_growth import *
from math import prod

def components(T):
    n=len(T);reach=[]
    for i in range(n):
        mask=(1<<i)|T[i]
        while True:
            new=mask
            for j in range(n):
                if mask>>j&1:new|=T[j]
            if new==mask:break
            mask=new
        reach.append(mask)
    remaining=set(range(n));blocks=[]
    while remaining:
        i=min(remaining);block=tuple(j for j in remaining if reach[i]>>j&1 and reach[j]>>i&1)
        blocks.append(block);remaining.difference_update(block)
    snapshot=tuple(blocks)
    blocks.sort(key=lambda block:-sum(bool(T[block[0]]>>other[0]&1) for other in snapshot if other!=block))
    return blocks

def profile(sizes):
    n=sum(sizes);out=Counter({tuple(sizes)+(1,):1})
    for i,a in enumerate(sizes):
        later=sum(sizes[i+1:])
        count=((1<<a)-1)*(1<<(later-1)) if later else (1<<(a-1))-1
        if count:out[tuple(sizes[:i])+(a+later+1,)]+=count
    need(sum(out.values())==1<<(n-1),'suffix row sum')
    return out

def run():
    data=json.loads(gzip.decompress((HERE/'rooted_growth.json.gz').read_bytes()))
    weights=defaultdict(Counter);transfers=0
    for rec in data['records']:
        T=tuple(rec['rows']);n=len(T);blocks=components(T);sizes=tuple(map(len,blocks))
        need(n-1 in blocks[-1],'end root lies in terminal SCC')
        weights[n][sizes]+=rec['mass']
        if n<=7:
            actual=Counter()
            for mask in range(1<<(n-1)):
                U=cc.child(T,mask);actual[tuple(map(len,components(U)))]+=1
            need(actual==profile(sizes),('complete per-rooted-state SCC transition profile',T,sizes,actual,profile(sizes)))
            transfers+=1
    X=[1]+[1<<((n-1)*(n-2)//2) for n in range(1,31)];C=[0]*31
    for n in range(1,31):C[n]=X[n]-sum(C[i]*X[n-i] for i in range(1,n))
    for n,distribution in weights.items():
        for sizes,mass in distribution.items():need(mass==prod(C[s] for s in sizes),'complete component-word orbit weight')
        need(distribution[(n,)]==C[n],'independent strong presentation mass')
    save('component_growth',dict(scope='SCC suffix-merge profiles for every end state through7; complete size-word masses through8; formal strong-presentation recurrence through30',profiles_checked=transfers,strong_presentation_counts=C[1:],distributions=[dict(n=n,words=[dict(sizes=s,mass=m) for s,m in sorted(v.items())]) for n,v in sorted(weights.items())]))
    print('PASS',transfers,'complete profiles; strong presentation counts',C[1:11])
if __name__=='__main__':run()
