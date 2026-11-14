from itertools import product
from replay_bilateral import count, k41, k43, k65, k69

def graph(M):
    # Two singleton wings at each end; identical apex contact in each wing.
    shapes=[((0,),1)]*2+[(C.masks,C.out) for C in M]+[((0,),0)]*2
    out=[]
    N=sum(len(R) for R,U in shapes)
    apex_out=0
    for R,U in shapes:
        m=len(R); j=len(out)
        apex_out|=U<<j
        for i,a in enumerate(R):
            out.append((a<<j)|(((1<<N)-(1<<(j+m))))|(((1-(U>>i&1)))<<N))
    out.append(apex_out)
    return out

def recurrence_odd_cycles(M):
    # Test the projected formula; signed coordinates used here are the negatives
    # of positive full Hamiltonian end/start counts.
    L=[(1,2,1,0,-1,0,1)]*2
    R=[(1,2,0,1,0,-1,1)]*2
    A=L+[(C.h,C.f,C.P,C.Q,-C.S,-C.T,C.K) for C in M]+R
    from math import prod
    h,f,p,q,s,t,k=zip(*A)
    P=prod(h)
    total=sum(P//h[i]*k[i] for i in range(len(A)))-(len(A)-1)*P
    for i in range(len(A)):
        for j in range(i+1,len(A)):
            flank=prod(h[:i]+h[j+1:])
            total+=p[i]*q[j]*prod(f[i+1:j])*flank
            if j==i+1: total+=s[i]*t[j]*flank
    return total

for A,B in [((k41,k43),(k41,k43)),((k41,k43),(k65,k69)),((k65,k69),(k41,k43))]:
    P=A[0].h*B[0].h
    C0=recurrence_odd_cycles([A[0],B[0]])
    for j,l in product(range(2),repeat=2):
        a,b=A[j],B[l]
        N=count(graph([a,b]))[0][-1]
        assert N==recurrence_odd_cycles([a,b])
        assert N==C0+b.h*(a.K-A[0].K)+a.h*(b.K-B[0].K)
        if A!=B:
            z=(N-C0)//2
            r=z*pow(b.h,-1,a.h)%a.h
            s=z*pow(a.h,-1,b.h)%b.h
            assert r==(a.K-A[0].K)//2 and s==(b.K-B[0].K)//2
    print('exact baseline, switch-core counts:',C0,A[0].h,B[0].h)
for m in (2,3,5,20):
    T=[k43]+[k41]*(m-1)
    U=[k41]*(m-1)+[k43]
    assert recurrence_odd_cycles(T)==recurrence_odd_cycles(U)
    assert recurrence_odd_cycles(T[1:])+2*15**(m-2)==recurrence_odd_cycles(U[1:])
print('PASS: literal common-apex switches, CRT and neutral-word deletion.')
