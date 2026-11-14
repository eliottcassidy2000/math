from profiles_near_30030 import strong, profile, hp_all
from inverse_carries import logical
for A,p in [((128,129,99,71,15,27,51,124),44505),((36,97,50,7,11,216,29,95),45225)]:
    logical(strong(A) and profile(A)==(345,5640))
    h=hp_all(A); M=len(h)-1; t=0
    for S in range(M+1):
        K=M^S; R=K
        while True:
            t+=h[S]*h[R]*h[K^R]
            if not R: break
            R=(R-1)&K
    logical(t==p)
print('equal strong source H and Q2; the next slot differs by exactly 720')
