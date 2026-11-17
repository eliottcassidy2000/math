from collections import deque
from itertools import permutations
def require(ok,*x):
    if not ok:raise ArithmeticError(x)
def H(A):
    n=len(A);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for mask,r in enumerate(D):
        for i,v in enumerate(r):
            if not v:continue
            w=A[i]&~mask
            while w:
                b=w&-w;w-=b;D[mask|b][b.bit_length()-1]+=v
    return sum(D[-1])
def flip(a,i,j):
    a=list(a);a[i]^=1<<j;a[j]^=1<<i;return a
def join(a,b):
    n=len(a);m=len(b)
    return tuple(r|(((1<<m)-1)<<n) for r in a)+tuple(r<<n for r in b)
def mates(a,b):
    n=len(a)
    for i in range(n):
        for j in range(i+1,n):
            x=a[i]>>j&1;y=b[i]>>j&1
            if j==i+1:
                if (x,y)!=(1,1):return False
            elif x==y:return False
    return True
def main():
    for a,b,h,p in [
       ((14,28,24,16,1),(18,4,9,19,6),9,11),
       ((30,124,120,112,96,65,1),(98,4,9,19,39,78,30),47,91)]:
        require(mates(a,b) and (H(a),H(b))==(h,p),a,b)
    C=(2,4,1);D=(6,12,8,1)
    for A,B,w in [(C,C,17),(C,D,29),(D,C,29),(D,D,49)]:
        J=join(A,B)
        require(min(H(flip(J,i,j)) for i in range(len(A))
                    for j in range(len(A),len(J)))==w)
    cross=[(i,j) for i in range(3) for j in range(3,7)]
    J=join(C,D);W=[]
    for q in range(1<<12):
        A=list(J)
        for k,(i,j) in enumerate(cross):
            if q>>k&1:A=flip(A,i,j)
        W.append(H(A))
    def reaches(threshold):
        visited={0};Q=deque([0])
        while Q:
            q=Q.popleft()
            if q==4095:return True
            for i in range(12):
                v=q^(1<<i)
                if W[v]<=threshold and v not in visited:Q.append(v);visited.add(v)
        return False
    require(not reaches(50) and reaches(51))
    route=[(0,3),(0,4),(1,4),(2,4),(1,3),(2,3),(0,6),(1,6),(0,5),(1,5),(2,5),(2,6)]
    path=[15];s=0
    for edge in route:
        s^=1<<cross.index(edge);path.append(W[s])
    require(path==[15,29,41,51,49,51,51,51,49,51,41,29,15],path)
    print("PASS literal framed mates, four adjacent exit minima, and all 4096 frozen-internal exchange profiles.")
if __name__=='__main__':main()
