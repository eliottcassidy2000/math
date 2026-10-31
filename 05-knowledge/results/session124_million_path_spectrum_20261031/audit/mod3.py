"""Matrix-product audit over F_3, without producer step-function imports."""
import itertools
# Rows act on (F,q,g,h).
Z=[[1,1,1,0],[0,2,0,0],[0,0,0,0],[0,0,0,1]]
O=[[1,0,0,0],[0,2,0,1],[0,0,0,1],[0,0,0,1]]
C=[[[3,9-3*j,3-j,0 if j in (0,3) else 2], [0,12,0,3*j], [0,0,0,j],[0,0,0,3]] for j in range(4)]
def mul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(4))%3 for j in range(4)] for i in range(4)]
def zero(A):return not any(x%3 for row in A for x in row)
def require(b,m):
 if not b:raise RuntimeError(m)
for j in (0,3):require(zero(C[j]),('constant',j))
for a,b,c in itertools.product((1,2),repeat=3):require(zero(mul(C[c],mul(C[b],C[a]))),('triple',a,b,c))
for a,b in itertools.product((1,2),repeat=2):
 for S in (Z,O):require(zero(mul(C[b],mul(S,C[a]))),('separated',a,b,S))
# After a singleton following a nonconstant triangle, rows q,g,h vanish.
for j in (1,2):
 for S in (Z,O):require(zero(mul(S,C[j])[1:]),('tail',j))
print('mod3 constant, triple, separated matrix identities PASS')
