import sympy as sp
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
M=[]
for let in range(4):
    A=sp.zeros(NS,NS)
    for s in range(NS):
        for t in T[let][s]: A[s,t]+=1     # row s -> successors; prefix row vector n' = n*A
    M.append(A)
for let in range(4):
    ev=M[let].eigenvals()
    print("letter",let,"(x,y)=",(let>>1,let&1),"eigenvalues:",dict(sorted(ev.items())))
# common left eigenvectors: v*A = lam v for all letters
for lam in (1,2,3):
    stack=sp.Matrix.vstack(*[(M[l]-lam*sp.eye(NS)).T for l in range(4)])
    ns=stack.nullspace()
    print("common left eigenvectors for lambda =",lam,": dimension",len(ns))
    for v in ns[:3]:
        print("   support:",[ (s, v[s]) for s in range(NS) if v[s]!=0 ])
