# Strong tournament classes of order <= 8 whose two-path cover count Q (= Gamma) is a power of two.
import sys
from tl import run2
S=sys.argv[1]
def decode(code,n):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1: out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
def strong(out):
    n=len(out)
    def reach(adj):
        seen={0}; st=[0]
        while st:
            x=st.pop()
            for y in range(n):
                if adj(x,y) and y not in seen: seen.add(y); st.append(y)
        return len(seen)==n
    return reach(lambda x,y: out[x]>>y&1) and reach(lambda x,y: out[y]>>x&1)
cnt=0
for n in range(3,9):
    codes=[int(l.split()[0]) for l in open(S+'/cls_%d.txt'%n)]
    tours=[decode(c,n) for c in codes]
    for c,t,(H,G,R,E) in zip(codes,tours,run2(tours)):
        if strong(t) and G&(G-1)==0:
            cnt+=1; print('n=%d code=%d H=%d Q=2^%d'%(n,c,H,G.bit_length()-1))
print('total strong classes with Q a power of two (order <= 8):',cnt)
