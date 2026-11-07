"""Exact bounded certificate plus all-length lifting tests. Standard library only."""
import json
from pathlib import Path
P,Q=641,6700417

def require(x, msg):
    if not x: raise RuntimeError(msg)
def prime(n):
    return n>=2 and all(n%d for d in range(2,__import__('math').isqrt(n)+1))
def move(s,c):
    D,E,U,d=s
    if c=='O': return D,0,2*U,d
    if c=='Z': return D-E+U,d,2*U-d,d
    if c=='J': return D-E+U,D-E,D+U,D
    raise ValueError(c)
def state(word):
    s=1,0,2,1
    for c in word:s=move(s,c)
    return s
def F(s):return s[0]-s[1]+s[2]
def zero_table(m):
    return [(r,s) for r in range(64) for s in range(64)
            if (2+(pow(2,r,m)-1)*(pow(2,s,m)-1))%m==0]
def zeros_independent(m):
    # Use exact powers and integer divisibility, no modular exponentiation.
    powers=[1]
    for _ in range(63): powers.append(2*powers[-1])
    return [(r,s) for r,x in enumerate(powers) for s,y in enumerate(powers)
            if divmod(2+(x-1)*(y-1),m)[1]==0]
require(P*Q==2**32+1,'F5 factorization')
for p in (P,Q):
    require(prime(p), 'trial primality')
    require(pow(2,32,p)==p-1 and pow(2,64,p)==1,'exact order 64')
pairs={str(m):zero_table(m) for m in (P,Q,P*Q)}
for m in (P,Q,P*Q):require(pairs[str(m)]==zeros_independent(m),'independent phase census')
require(pairs[str(Q)]==[(1,32),(32,1)],'Q phase certificate')
require(set(pairs[str(Q)])<=set(pairs[str(P)]),'phase containment')
require(pairs[str(P)]==[(1,32),(5,11),(11,5),(25,63),(32,1),(44,49),(49,44),(63,25)],'P certificate')
for a in range(1,33):
    for b in range(1,33):
        require(F(state('O'*a+'Z'*b))==2+(2**(a+1)-1)*(2**(b+1)-1),'two-run physical formula')
prefix='OOZJOZJJOJJOZ'
s=state(prefix)
require(s==(4909,859,6013,859),'escape profile')
rows=[]
for t in range(6):
    c=20+64*t
    h=4909+6013*2**c
    require(F(state(prefix+'O'*c))==h,'physical escape family')
    require(h%Q==0 and h%P==380,'escape residues')
    rows.append({'t':t,'source_letters':len(prefix)+c,'strong_order':len(prefix)+c+3,'h':str(h),'cofactor_Q':str(h//Q),'mod_P':h%P})
require(int(rows[0]['h'])==6305092397==941*Q,'initial escape')
require(2+(2**5-1)*(2**11-1)==99*P,'P-only seed')
# A universal Fermat factor family, including both F5 factors.
for L in (1,2,4,8,16,32):
    x=2**L
    require(2+(x-1)*(2*x*x-1)==(x+1)*(2*x*x-4*x+3),'universal shared family')
receipt={'status':'PASS','phase_universe':'all 4096 ordered pairs modulo 64, two independent arithmetic methods','zero_pairs':pairs,'prefix':prefix,'profile':s,'escape_rows':rows,'two_run_formula_checks':32**2,'bounds':'No full grammar or all tournament coverage claim.'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
