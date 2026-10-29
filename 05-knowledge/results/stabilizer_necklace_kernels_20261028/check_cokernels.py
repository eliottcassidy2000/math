"""Integer Smith factors of the rigid-palette circulant, by Euclidean operations."""
from pathlib import Path
from math import prod
import json

CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

def smith(a):
    a=[r[:] for r in a]; n=len(a)
    for k in range(n):
        nonzero=[(abs(a[i][j]),i,j) for i in range(k,n) for j in range(k,n) if a[i][j]]
        if not nonzero: break
        _,i,j=min(nonzero); a[k],a[i]=a[i],a[k]
        for row in a: row[k],row[j]=row[j],row[k]
        while True:
            changed=False
            for i in range(k+1,n):
                if a[i][k]:
                    q=a[i][k]//a[k][k]
                    for j in range(k,n): a[i][j]-=q*a[k][j]
                    if a[i][k]: a[k],a[i]=a[i],a[k]
                    changed=True; break
            if changed: continue
            for j in range(k+1,n):
                if a[k][j]:
                    q=a[k][j]//a[k][k]
                    for i in range(k,n): a[i][j]-=q*a[i][k]
                    if a[k][j]:
                        for row in a: row[k],row[j]=row[j],row[k]
                    changed=True; break
            if changed: continue
            bad=next(((i,j) for i in range(k+1,n) for j in range(k+1,n) if a[i][j]%a[k][k]),None)
            if bad is None: break
            i,j=bad
            for j in range(k,n): a[k][j]+=a[i][j]
        if a[k][k]<0: a[k]=[-x for x in a[k]]
    need(all(a[i][j]==0 for i in range(n) for j in range(n) if i!=j),'Smith diagonal')
    diag=[a[i][i] for i in range(n)]
    need(all(diag[i+1]%diag[i]==0 for i in range(n-1)),'Smith divisibility')
    return diag

def main():
    rows=[]; d=[1]; q=[1]
    for k in range(1,7): q.append(16*q[-1]**3+3*q[-1])
    # q[0] corresponds to L_3/4, not L_1/4.
    for k in range(6):
        n=3**k; b=0 if k<2 else (3**(k-1)-1)//2
        expected=[1]+[2]*(n-1-b)
        for r in range(2,k+1): expected.extend([2*q[r-1]]*3**(k-r))
        need(len(expected)==n,'Smith dimensions')
        if k<=3:
            actual=smith([[d[(j-i)%n] for j in range(n)] for i in range(n)])
            need(actual==expected,'independent Euclidean Smith computation')
        determinant=prod(expected)
        rows.append({'depth':k,'N':n,'odd_Smith_factors':b,'Smith_factor_multiplicities':
                     {str(x):expected.count(x) for x in sorted(set(expected))},'determinant':str(determinant),
                     'integer_Smith_computed':k<=3})
        d=d+d+[-x for x in d]
    result={'status':'FINITE-EXACT independent integer checks for proved cokernel decomposition',
            'rows':rows,'checks':CHECKS}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
