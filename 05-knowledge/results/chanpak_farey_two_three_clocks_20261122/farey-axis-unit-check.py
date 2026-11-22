"""Run alongside the inherited fermat_receiver.py, with its imports available."""
from itertools import permutations
from math import gcd
from fermat_receiver import A,D,EDGE,flipped,join,receiver,tagged_data,run_machine

def advance(vec, tr, word):
    for bit in word:
        out=[0]*len(vec)
        for i,v in enumerate(vec):
            if v:
                for j in tr[i][bit]:out[j]+=v
        vec=out
    return vec

def verify_closed_form():
    # 6 direct-sum 81-state blocks. After m-1 copies of 1001, evaluate
    # the same 1010 terminal functional. Difference with an exponential
    # sum of dimension 3 obeys a monic recurrence of order at most 489.
    machines=[]
    for order in permutations("ab|"):
        tr,good=run_machine(order,None,False)
        machines.append((tr,good,[1]+[0]*(len(good)-1)))
    dimension=sum(len(vec) for _,_,vec in machines)+3
    for m in range(1,dimension+1):
        actual=0
        updated=[]
        for tr,good,vec in machines:
            actual+=sum(v for v,g in zip(advance(vec,tr,(1,0,1,0)),good) if g)
            updated.append((tr,good,advance(vec,tr,(1,0,0,1))))
        machines=updated
        assert 1170*actual==1000*81**m-441*16**m-1404
    print("receiver slope identity verified by",dimension,"exact initial values and Cayley–Hamilton")

if __name__=="__main__":
    a=join(tagged_data(A,12),tagged_data(flipped(D,*EDGE),35))
    for n in range(4,113,4):
        h,q,c,t,e,r,k=join(a,receiver(n))
        assert h==145*((1<<n)+1) and k>=h and k%2==1 and gcd(h,t)==1
    print("all 28 potentially relevant canonical axis slopes are units modulo h")
    verify_closed_form()

    for j in range(1,8):
        S=2**(4*j)+1
        t=10240*((1000*3**(4*j)-441*(S-1)-1404)//1170)+4592*S
        assert t%29==(8,18,15,17,10,16,13)[j-1]
    for n in range(4,8201,4):
        g=gcd(2**n+1,1000*3**n-963)
        assert g==(19681 if n==8200 else 1)
    p=19681
    assert all(p%d for d in range(2,141))
    for a,o,primes in [(2,3280,(2,5,41)),(3,9840,(2,3,5,41))]:
        assert pow(a,o,p)==1 and all(pow(a,o//q,p)!=1 for q in primes)
    assert pow(2,8200,p)==p-1 and 1000*pow(3,8200,p)%p==963
    for a,d,c,A0,B0,root in [(2,1,1,1745,2094,16400),(3,1000,-963,6270,7167,1828)]:
        assert (d*pow(a,8200,p*p)+c)%p==0
        assert (d*pow(a,8200,p*p)+c)//p%p==A0
        assert d*pow(a,8200,p)*((pow(a,9840,p*p)-1)//p)%p==B0
        assert (-A0*pow(B0,-1,p))%p==root
    print("first exceptional axis at n=8200, and incompatible p-square lift verified")
    chart=join(
        join(tagged_data(A,12),tagged_data(flipped(D,*EDGE),35)), receiver(8200))
    h,q,c,t,e,r,k=chart
    assert (q%p,c%p,r%p)==(16538,18124,0)
    assert gcd(gcd(h,t),gcd(c,2*q))==1
    print("the full two-port chart at n=8200 restores every odd modular phase")
