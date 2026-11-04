"""Independent constructive checks; standard library only."""

def require(condition):
    if not condition: raise RuntimeError("Honeycomb verification failed")

def T(n): return n//2 if n%2==0 else 3*n+1

def v2(n):
    k=0
    while n%2==0: n//=2; k+=1
    return k

def chart(n):
    a=v2(n)
    return a,n//2**a-1

def unchart(p):
    x,y=p
    require(x>=0 and y>=0 and y%2==0)
    return 2**x*(y+1)

def arrow(p):
    x,y=p
    if x: return x-1,y
    a=v2(3*y+4)
    return a,(3*y+4)//2**a-1

def adjacent(p,q):
    x,y=p; X,Y=q
    return (y==Y and abs(x-X)==1) or (x==X and ((x+y)%2==0 and Y==y+1 or (x+y)%2==1 and Y==y-1))

def embedding(d):
    centers={'':(1,0)}; paths={}
    for t in range(d):
        for word in [w for w in centers if len(w)==t]:
            x,y=centers[word]; D=2**(d-t)
            for bit,sgn in [('0',-1),('1',1)]:
                X=x+sgn*D; bend=X-1; step=1 if bend>x else -1
                p=[(u,y) for u in range(x,bend+step,step)]
                p += [(bend,y+1),(X,y+1),(X,y+2)]
                centers[word+bit]=(X,y+2); paths[word+bit]=p
    return centers,paths

def check_embedding(d):
    centers,paths=embedding(d)
    require(len(set(centers.values()))==len(centers))
    interiors=set()
    for child,p in paths.items():
        require(p[0]==centers[child[:-1]] and p[-1]==centers[child])
        require(all(adjacent(a,b) for a,b in zip(p,p[1:])))
        require(len(p)==len(set(p)))
        internal=set(p[1:-1])
        require(not internal.intersection(centers.values()))
        require(not internal.intersection(interiors))
        interiors.update(internal)
    return len(centers),len(interiors)

def inverse_ball(root,d):
    seen={root}; front={root}
    for _ in range(d):
        nxt=set()
        for n in front:
            nxt.add(2*n)
            if n%6==4: nxt.add((n-1)//3)
        front=nxt-seen; seen|=front
    return seen

if __name__=='__main__':
    for n in range(1,100001):
        require(unchart(chart(n))==n)
        require(arrow(chart(n))==chart(T(n)))
        if n%2==0: require(adjacent(chart(n),chart(T(n))))
    for d in range(1,11): print('depth',d,'branch/interior counts',check_embedding(d))
    b=inverse_ball(5,21)
    require(len(b)==704)
    # BFS in the brick-wall honeycomb.
    seen={(0,0)}; front=seen.copy()
    for _ in range(21):
        nxt=set()
        for x,y in front: nxt.update([(x-1,y),(x+1,y),(x,y+(1 if (x+y)%2==0 else -1))])
        front=nxt-seen; seen|=front
    require(len(seen)==694)
    print('100000 chart identities passed; inverse ball',len(b),'honeycomb ball',len(seen))
