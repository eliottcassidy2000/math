import json,math

def least_prime(n):
    if n<2:raise ValueError(n)
    if n%2==0:return 2
    d=3
    while d*d<=n:
        if n%d==0:return d
        d+=2
    return n
rows=[]
for c in range(-20,21):
  if c in (-1,0,1):continue
  a=abs(c)
  for ell in (3,5,7):
    den=a+1 if c>0 else a-1
    if den%ell==0:continue
    N=(a**ell+1)//den if c>0 else (a**ell-1)//den
    p=least_prime(N)
    if N<=1 or N%2==0 or math.gcd(N,den)!=1 or pow(c,ell,p)!=p-1 or pow(c,2,p)==1 or (p-1)%(2*ell):
        raise RuntimeError((c,ell,N,p))
    rows.append({'c':c,'ell':ell,'quotient':N,'least_prime_factor':p,'order':2*ell})
print(json.dumps({'controls':len(rows),'method':'Independent trial division; exact order forced by prime ell and checked proper divisors','rows':rows},indent=2))
