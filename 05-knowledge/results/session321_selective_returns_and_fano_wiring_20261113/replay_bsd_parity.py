import replay_fano as F
from replay_clone_profiles import C,P,X,Y,U,V,inflate,add
def phi(A):
    A=tuple(a^((a.bit_count()%2)<<i) for i,a in enumerate(A))
    out=[]
    for i,t in enumerate(A):
        s=t^(1<<i)
        for j,a in enumerate(A):
            if t>>j&1:s^=a
        out.append(s)
    return tuple(out)
def nullity(A):
    piv={}
    for s in phi(A):
        while s:
            p=s.bit_length()-1
            if p in piv:s^=piv[p]
            else:piv[p]=s;break
    N=(1<<len(A))-1
    if phi(A)!=phi(tuple(a^N^(1<<i)for i,a in enumerate(A))):
        raise ValueError('converse')
    return len(A)-len(piv)
for m,d in zip((0,1,3,7,11,13),(2,0,0,0,4,2)):
    A=tuple(sum(1<<j for j,v in enumerate(F.V)if u!=v and F.twisted(u,v,m))for u in F.V)
    if nullity(A)!=d:raise ValueError((m,d))
for r,x,y,k,d in ((C,X,Y,7167,((4,2),(2,0))),(P,U,V,6144,((2,6),(2,6)))):
    A,B=inflate(r,x),inflate(r,y)
    for m in range(8):
        if (nullity(A),nullity(B))!=d[m%2]:raise ValueError((k,m))
        A=add(A,k);B=add(B,k)
print('six Fano dimensions and both selective parity laws checked')
