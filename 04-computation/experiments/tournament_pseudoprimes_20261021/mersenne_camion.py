# Which 2^k - 1 (k <= 64) and 4^q - 1 are Camion pseudoprimes?  Decide via the criterion when every residue
# r_p = n mod (p-1) is small enough to compute D(r_p) mod p exactly (r <= 3000); else report undecided.
import subprocess, sys
def factor(n):
    out=subprocess.run(['gp','-q'],input=f'f=factor({n}); for(i=1,#f~, print(f[i,1], " ", f[i,2]))',capture_output=True,text=True).stdout
    return [tuple(map(int,l.split())) for l in out.strip().split('\n') if l.strip()]
def D_mod(r,p):
    T=[pow(2,j*(j-1)//2,p) for j in range(r+1)]; s=[0]*(r+1); row=[1]
    for k in range(1,r+1):
        row=[1]+[(row[i-1]+row[i])%p for i in range(1,k)]+[1]
        acc=0
        for i in range(1,k): acc=(acc+row[i]*s[i]*T[k-i])%p
        s[k]=(T[k]-acc)%p
    c=pow(2,(r-1)//2,p) if r%2 else 0
    return (s[r]-c)%p
res=[]
for k in range(2,65):
    n=2**k-1; fs=factor(n)
    if len(fs)==1 and fs[0][1]==1: continue  # prime
    if any(e>1 for p,e in fs if p>2): res.append((k,'squareful',[p for p,e in fs])); continue
    verdict='PASS'; rs=[]
    for p,e in fs:
        r=(n-1)%(p-1)+1; rs.append((p,r))
        if r in (1,2,3): continue
        if r<=3000:
            if D_mod(r,p)!=0: verdict='fail'; break
        else: verdict='undecided'
    res.append((k,verdict,rs))
for k,v,rs in res: print(k, v, rs if v!='squareful' else rs)
