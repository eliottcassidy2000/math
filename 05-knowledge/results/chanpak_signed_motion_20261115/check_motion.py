from oracle import tournament, lift, path_dp
from math import prod

def TT(n): return [[int(i<j) for j in range(n)] for i in range(n)]
def dp(a):
    n=len(a); W=1<<n; d=[[0]*n for _ in range(W)]
    for u in range(n): d[1<<u][u]=1
    for S in range(1,W):
        for u,x in enumerate(d[S]):
            if x:
                for v in range(n):
                    if not S>>v&1 and a[u][v]: d[S|1<<v][v]+=x
    return d
def grades(a, mask):
    n=len(a); F=(1<<n)-1
    e=dp(a); l=dp(list(zip(*a)))
    v=[i for i in range(n) if mask>>i&1]
    h=[1]+[sum(row) for row in e[1:]]
    E=[sum(row[i] for i in v) for row in e]
    L=[sum(row[i] for i in v) for row in l]
    return h[-1],sum(h[s]*h[F^s] for s in range(F+1)),sum(E[s]*h[F^s] for s in range(F+1)),sum(h[s]*L[F^s] for s in range(F+1)),sum(E[s]*L[F^s] for s in range(F+1))
def Q(A,j):
    h=[1]+[sum(row) for row in dp(A)[1:]]
    c=[1]+[0]*(len(h)-1)
    for _ in range(j):
        k=[]
        for s in range(len(h)):
            u=s; w=0
            while True:
                w+=c[u]*h[s^u]
                if not u:break
                u=(u-1)&s
            k.append(w)
        c=k
    return c[-1]
def substitute(a, w, P):
    own=[i for i in range(len(a)) if i!=w]+[w]*len(P)
    z=len(a)-1
    return [[a[u][v] if u!=v else P[i-z][j-z] if u==w else 0
             for j,v in enumerate(own)] for i,u in enumerate(own)]
def isstrong(a):
    for t in [a,list(zip(*a))]:
        got={0}
        for _ in t: got|={v for u in got for v,w in enumerate(t[u]) if w}
        if len(got)!=len(t):return False
    return True
def real(D,mask,a,b):
    return lift([TT(a),D,TT(b)],[set(),{v for v in range(len(D)) if mask>>v&1},set()])
def run():
    three=[(4,0,9,0,3,1,16,9,9,4),(5,4,24,0,1,5,64,15,39,6),(4,63,1,0,3,1,16,1,8,0)]
    five=[(6,10,45,1,2,9,168,99,129,68),(6,10,43,1,2,9,168,95,127,66),(7,2097141,82,0,2,9,336,198,129,64),(6,10,51,1,2,9,168,91,125,62),(6,8,53,1,2,9,168,89,124,60)]
    for rows,c,d,start in [(three,128,-63,-3),(five,1344,-318,-59)]:
        for i,(n,code,v,u,w,h,q,E,L,J) in enumerate(rows):
            D=tournament(n,code)
            assert grades(D,v)==(h,q,E,L,J)
            assert (q*2**(u+w),E*2**u-L*2**w,h-J)==(c,d,start+2*i)
            W=real(D,v,1+u,1+w)
            assert isstrong(W) and path_dp(W)==4*c+2*d+start+2*i
    print("Three all-coordinate-moving and five quadratic-moving physical rows passed.")
    coefficients=[[(3,6,0),(-6,48,42),(12,14,0),(-1,6,6),(0,0,0)],
                  [(3,6,0),(-6,48,42),(6,14,0),(-4,24,12),(2,4,0)]]
    physical=[TT(m) for m in range(1,5)]+[tournament(3,2),tournament(4,4),tournament(4,19),tournament(5,76)]
    R=tournament(6,22);assert isstrong(R)
    for P in physical:
        q=[Q(P,j) for j in (1,2,3)];D=substitute(R,2,P)
        assert isstrong(D)
        for v,a,b,C,c in [(4,1,2,coefficients[0],-17),(17,2,1,coefficients[1],-15)]:
            assert grades(D,v)==tuple(sum(x*y for x,y in zip(row,q)) for row in C)
            W=real(D,v,a,b)
            assert isstrong(W) and path_dp(W)==312*q[2]+394*q[1]+c*q[0]
    print("Four Vandermonde rows and four nontransitive substitutions passed in both physical outputs.")
if __name__=="__main__":run()
