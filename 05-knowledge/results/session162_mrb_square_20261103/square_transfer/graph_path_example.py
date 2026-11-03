from math import isqrt
p=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
if set(p)!=set(range(1,16)): raise AssertionError("not a permutation")
if not all(isqrt(a+b)**2==a+b for a,b in zip(p,p[1:])): raise AssertionError("non-square edge")
g=0; prev=5; P=R=1; pairs=[]
for a,b in zip(p,p[1:]):
 d=abs(a-b); s=1
 while s*s%2!=d%2 or s*s<2*(prev+1)+d: s+=1
 x=(s*s-d)//2; y=x+d; g+=x-prev-1;prev=y
 q=1 if d==1 else 1+2**(d-1)
 t=2**(d+1)+2*(3**(d-1)-1)
 P*=q;R*=t;pairs.append((d,x,y,s*s))
R*=2**g
print('gap_count',g,'largest',prev,'pairs',pairs)
print('P',P,'R',R)
print('chain',*[P+32*R,3*P+48*R,5*P+64*R])
print('jump',2*P+16*R)
