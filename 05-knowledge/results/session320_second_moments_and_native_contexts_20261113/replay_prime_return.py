from math import isqrt
from replay_signed_tournaments import (strong, coeffs, inflate, hq2, info,
                                      C as C0, X as X0, Y as Y0)
R=(44,33,2,6,15,28)
X=(766,248,1018,832,264,280,944,312,3,434)
Y=(32,909,1,821,295,454,415,541,133,373)
def tour(A):
    return all((A[i]>>j&1)!=(A[j]>>i&1)
               for i in range(len(A)) for j in range(i+1,len(A)))
if not all(tour(A) and strong(A) for A in (R,X,Y)):
    raise ValueError('invalid witness')
if coeffs(R)!=[25,24,6]:raise ValueError('layout coefficients')
if info(X)!=(1915,56008,710883) or info(Y)!=(1915,55288,711603):
    raise ValueError('source profiles')
T,U=inflate(R,X),inflate(R,Y)
N=4591771
if any(N%d==0 for d in range(2,isqrt(N)+1)):raise ValueError('composite')
if not strong(T) or not strong(U) or hq2(T)!=(N,222226400) or hq2(U)!=(N,224279840):
    raise ValueError('outputs')
A,B=info(inflate(C0,X0)),info(inflate(C0,Y0))
if (A[2],B[2])!=(216103869,219387789):raise ValueError('return readiness')
print('prime exact return and post-return three-visit obstruction verified')
