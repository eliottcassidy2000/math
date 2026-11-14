from fermat_profiles import strong, rooted, glue, submasks, hp_all, profile

def need(ok):
    if not ok: raise RuntimeError('full-slot collision replay failed')

A=(6,12,80,149,163,15,59,103)
B=(2,52,73,131,45,77,155,55)
def all_values(A):
    n=len(A); M=(1<<n)-1
    H=hp_all(A); q=[0]*(M+1);q[0]=1; values=[0]
    subs=[tuple(submasks(s)) for s in range(M+1)]
    for z in range(1,n+1):
        q=[sum(q[t]*H[s^t] for t in subs[s]) for s in range(M+1)]
        values.append(q[M])
    return values

for T,p in ((A,(349,860)),(B,(355,857))):
    need(strong(T) and list(rooted(T))[6]==p+(5536,6163,433,6))
for j in (0,1,2):
    a=all_values(glue(A,6,j) if j else A);b=all_values(glue(B,6,j) if j else B)
    need([x-y for x,y in zip(a,b)]==[f*2*z*z*(z*z-4) for z in range(len(a))
                                          for f in ((1,) if j==0 else (0,) if j==1 else (-z,))])
U,V=glue(A,6,1),glue(B,6,1)
need(strong(U) and strong(V) and profile(U)==profile(V)==(2069,58256))
need(profile(glue(U,8,1))==(6853,403552) and profile(glue(V,8,1))==(6859,403552))
p,q=349+860*2048,355+857*2048
need((p-q)%396==198 and p-q==6138)
print('full slot coalescence under the Fermat triangle and literal rerooted two-step switch verified')
