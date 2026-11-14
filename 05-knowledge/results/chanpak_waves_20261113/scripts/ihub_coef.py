from ihub import Hint
def coefs(m,p,q,e):
    base=Hint([0]*m,p,q,e)
    lin=[Hint([1 if k==i else 0 for k in range(m)],p,q,e)-base for i in range(m)]
    quad={}
    for i in range(m):
        for j in range(i+1,m):
            x=[1 if k in (i,j) else 0 for k in range(m)]
            quad[(i,j)]=Hint(x,p,q,e)-base-lin[i]-lin[j]
    return base,lin,quad
for (m,p,q,e) in [(8,0,0,0),(8,3,2,0),(8,3,2,1),(8,2,4,0)]:
    b,l,qd=coefs(m,p,q,e)
    print(f"m={m} I=[{p},{p+q}) e={e}: c0={b}, linear={l}")
    for i in range(m):
        print("   ",[qd[(i,j)] if j>i else '.' for j in range(m)])
