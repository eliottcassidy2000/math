#!/usr/bin/env python3
"""Independent coefficient-level exhaustive audit, standard library only.
Enumerates every matrix and translation through m=4, not just cyclic maps.
"""
import hashlib,json
from pathlib import Path

def require(ok,*why):
    if not ok:raise RuntimeError(why)

def rotate(x,m):return ((x<<1)&((1<<m)-1))|(x>>(m-1))

rows=[]
for m in range(1,5):
    n=1<<m;equivariant=set();matrix_count=1<<(m*m)
    # With this bit convention T(f)=R(f)+e_0.  Hence F=Af+b is
    # equivariant iff AR=RA and A(e_0)=b+R(b).
    for matrix in range(matrix_count):
        columns=tuple((matrix>>(m*j))&(n-1) for j in range(m))
        if any(columns[(j+1)%m]!=rotate(columns[j],m) for j in range(m)):continue
        for b in range(n):
            if columns[0]==b^rotate(b,m):equivariant.add((b,columns))
    expected=set()
    for h in range(n):
        c=h^rotate(h,m);columns=[]
        for j in range(m):columns.append(c);c=rotate(c,m)
        expected.add((h,tuple(columns)))
    require(equivariant==expected,'normal form failure',m)
    rows.append({'m':m,'matrices_considered':matrix_count,'affine_parameter_pairs_covered':matrix_count*n,
                 'equivariant_affine_maps':len(equivariant)})
result={'status':'PASS','universe':'Every matrix A and every translation b for dimensions 1 through 4, using independently derived coefficient identities',
        'rows':rows,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('independent_audit.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,sort_keys=True))
