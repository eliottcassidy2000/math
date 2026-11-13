

def state(arcs, root):
    n=len(arcs); lim=1<<n
    d=[[0]*n for _ in range(lim)]
    for i in range(n): d[1<<i][i]=1
    for mask in range(1,lim):
        for i,v in enumerate(d[mask]):
            if v:
                for j in arcs[i]:
                    if not (mask>>j&1): d[mask|1<<j][j]+=v
    h=list(map(sum,d)); h[0]=1
    return (h[-1],state_start(arcs,root),h[lim-1-(1<<root)],sum(d[mask][root]*h[lim-1-mask] for mask in range(lim)))

def state_start(a, root):
    d={(1<<root,root):1}
    for m in range(1<<len(a)):
        for i in range(len(a)):
            v=d.get((m,i),0)
            if v:
                for j in a[i]:
                    if not m>>j&1: d[m|1<<j,j]=d.get((m|1<<j,j),0)+v
    return sum(d.get(((1<<len(a))-1,i),0) for i in range(len(a)))

def add(a, root):
    n=len(a); a=[set(s) for s in a]; a[root].add(n)
    return a+[set(range(n))-{root}]

def check():
    c=0
    for m in range(2,9):
        a=[set(range(i+1,m+1)) for i in range(m)]+[{m+1},set(range(m))]
        for l in range(m):
            r=m-1-l
            if state(a,l)!=(1+2**m,2**r,1+2**(m-1),2**(l+1)+2**(r+1)*3**l-int(r==0)):
                raise ValueError((m,l))
            c+=1
        if m>=3:
            for k in range(1,4):
                z=a
                for _ in range(k+1): z=add(z,1)
                rt=len(z)-1; z=add(z,rt)
                r=m-2; t=2**k*(3+4*2**r)+1+2**(r+1)-(2+2**r)
                b=3*t+2*(2+2**r)
                if state(z,rt)!=(b,t,t,b): raise ValueError((m,k))
    print(c, "literal rooted states; 18 explicit birth words passed")
if __name__=="__main__": check()
