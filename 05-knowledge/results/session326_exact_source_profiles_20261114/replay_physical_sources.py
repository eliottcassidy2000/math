"""Independent profile/entry replay for the growing physical alphabet."""
from inverse_carries import logical,pair
from profiles_near_30030 import profile,extend,strong
from math import gcd,isqrt
from pathlib import Path
from json import loads

WORDS=((1,0,1,0,1,0,1,0,1,0,0),(1,0,1,1,1,0,1,1,1,0,0))
D=loads(Path(__file__).with_name('phase_table.json').read_text())

def word(k,s): return (1,)*k+WORDS[s]+(0,)*k
def graph(k,s):
    w=word(k,s); n=len(w)
    chain=tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
    return extend(chain,sum(b<<i for i,b in enumerate(w)))

def values(k,s):
    h=2048*4**k-1026*2**k+119+12*s
    q=354294*9**k-242340*6**k+50336*4**k-2052*2**k+2
    return h,q

def inverse(N):
    if N<1 or N%2==0: return
    V=[]
    k=0; max_q=N//(99*32)
    while values(k,0)[1]<=max_q:
        for s in (0,1):
            h,q=values(k,s)
            if gcd(h,99)==1: V.append((h,q,k,s))
        k+=1
    # q is monotone in k; both choices at k have the same q.
    def rec(i,h,q,data):
        c=N*pow(h,-1,198)%198
        L,R=divmod(N-c*h,99*q)
        if not R and L>=32 and L&(L-1)==0:
            yield (c,L.bit_length()-1,data)
        for j in range(i,len(V)):
            a,b,k,s=V[j]
            if q*b>max_q: break
            yield from rec(j,a*h,b*q,data+((k,s),))
    yield from rec(0,1,1,())

def main():
    for k in range(14):
        for s in (0,1):
            w=word(k,s); h,q=values(k,s)
            logical(pair(w)==(h,q))
            logical((gcd(h,99)==1)==(s==1 or k%10 not in (1,4)))
            if k<=2:
                T=graph(k,s); logical(strong(T) and profile(T)==(h,q))
    for k,c,N,f in ((0,17,507659921,1),(0,41,507687593,7),
                    (1,19,6120213613,1),(1,127,6120890881,7)):
        h,q=values(k,1)
        logical(D[str(c)]<=4 and N==99*32*q+c*h)
        logical((c,5,((k,1),)) in tuple(inverse(N)))
        p=N//f
        logical(N%f==0 and p>1 and all(p%d for d in range(2,isqrt(p)+1)))
    print('literal strong pairs through wrapping depth two; four inverse and prime controls')
if __name__=='__main__': main()
