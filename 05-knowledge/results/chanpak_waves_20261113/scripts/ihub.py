# Interval-hub family: spine TT_m, hub u1 with out-set {i: x_i=1}, hub u2 with out-set I=[p,p+q), hub arc e (1: u1->u2).
NS=36
def idx(o,a,b,d): return ((o*3+a)*3+b)*2+d
T=[[None]*NS for _ in range(4)]
for let in range(4):
    xi,yi=let>>1,let&1
    for o in range(2):
        for a in range(3):
            for b in range(3):
                for d in range(2):
                    s=idx(o,a,b,d); h1=yi if o else xi; h2=xi if o else yi
                    succ=[idx(o,2 if h1 else 1,b,d)]
                    if b or h1: succ.append(idx(o,a,2 if h2 else 1,d))
                    if d or h2: succ.append(idx(o,a,b,1))
                    T[let][s]=succ
fin=[[0]*NS for _ in range(2)]
for e in range(2):
    for o in range(2):
        for a in range(3):
            for b in range(3):
                for d in range(2):
                    h12=(1-e) if o else e
                    fin[e][idx(o,a,b,d)]=0 if (a==2 or b==2 or (b==0 and not h12)) else 1
def H(x,y,e):
    v=[0]*NS; v[idx(0,0,0,0)]=1; v[idx(1,0,0,0)]=1
    for xi,yi in zip(x,y):
        let=2*xi+yi; w=[0]*NS
        for s in range(NS):
            if v[s]:
                for t in T[let][s]: w[t]+=v[s]
        v=w
    return sum(v[s] for s in range(NS) if fin[e][s])
def Hint(x,p,q,e):
    m=len(x); y=[1 if p<=i<p+q else 0 for i in range(m)]
    return H(x,y,e)
if __name__=="__main__":
    import random
    random.seed(2)
    for trial in range(2):
        m=9; x=[random.getrandbits(1) for _ in range(m)]
        X=int(''.join(map(str,x)),2)
        print("x=",''.join(map(str,x)),"X=",X,"one-hub H(x)=",H(x,[0]*m,0),H(x,[0]*m,1))
        for e in (0,1):
            print(" e=",e)
            for p in range(m+1):
                print("  p=%d:"%p, [Hint(x,p,q,e) for q in range(0,m-p+1)])
