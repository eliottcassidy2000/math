from itertools import product

def outgoing(h,d):
    return d >> (1-h) & 1

def transfer(s,d,order,run):
    z=s[run]
    if z==0 and run and not outgoing(order[run-1],d):
        return None
    value=1+int(run==2 or not outgoing(order[run],d))
    return s[:run]+(value,)+s[run+1:]

def accept(s,order):
    return s[0]!=1 and (s[1]==2 or s[1]==0 and order==(0,1))

def automaton():
    ans=[]
    for order in [(0,1),(1,0)]:
        states=[(0,0,0)]
        for s in states:
            for d,r in product(range(4),range(3)):
                t=transfer(s,d,order,r)
                if t is not None and t not in states: states.append(t)
        tab=[[tuple(states.index(t) for r in range(3) if (t:=transfer(s,d,order,r)) is not None) for s in states] for d in range(4)]
        ans.append((tab, [j for j,s in enumerate(states) if accept(s,order)]))
    return ans

AUTO=automaton()
def H(w):
    ans=0
    for tab,accepts in AUTO:
        u=[1]+[0]*(len(tab[0])-1)
        for d in w:
            v=[0]*len(u)
            for paths,indices in zip(u,tab[int(d)]):
                if paths:
                    for j in indices: v[j]+=paths
            u=v
        ans+=sum(u[j] for j in accepts)
    return ans

def bar(w): return ''.join(str(3-int(d)) for d in w)
def orbit(w):
    double=w+bar(w)
    return [(double+double)[j:j+len(w)] for j in range(2*len(w))] if w else ['']

if __name__=='__main__':
    print('live and original reachable states',tuple(len(tab[0]) for tab,_ in AUTO))
    print(H('1102311'))
    for j,v in enumerate(orbit('110022311')): print(j,v,H(v))
