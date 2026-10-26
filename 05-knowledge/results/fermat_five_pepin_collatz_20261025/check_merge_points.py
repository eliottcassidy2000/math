# Coalescence of the orbits of 2^m+1: landing points after the valuation-2 run are 3^L c + 1 (c in {2,4}); THM-4581 pair chains.
import json
def v2(x): return (x & -x).bit_length()-1
def odd_orbit(x):
    seq=[x]
    while x!=1:
        y=3*x+1; x=y>>v2(y); seq.append(x)
    return seq
orbits={m: odd_orbit(2**m+1) for m in range(3,41)}
def merge(m1,m2):
    s1=orbits[m1]; pos={x:i for i,x in enumerate(s1)}
    for j,x in enumerate(orbits[m2]):
        if x in pos: return {'value':x,'odd_index_m1':pos[x],'odd_index_m2':j,'bits':x.bit_length()}
res={}
for m in range(3,40):
    r=merge(m,m+1); res[f'{m},{m+1}']=r
# which blocks share total stopping time and merge before 1
blocks={}
for m in range(3,41): blocks.setdefault(len(orbits[m])-1,[]).append(m)
print(json.dumps({'merge_consecutive':res,'equal_stopping_time_blocks':{str(k):v for k,v in sorted(blocks.items())}},indent=0))

# Theorem M check: for m = 2L+1, orbits of 2^m+1 and 2^(m+1)+1 first merge after exactly L+2 odd steps if L is even,
# exactly L+3 if L = 3 mod 4; for L = 1 mod 4 no merge within L+3 steps.
checks={}
for L in range(2,20):
    r=merge(2*L+1,2*L+2); i,j=r['odd_index_m1'],r['odd_index_m2']
    if L%2==0: ok=(i==j==L+2)
    elif L%4==3: ok=(i==j==L+3)
    else: ok=not (i==j and i<=L+3)
    checks[L]={'merge_index':(i,j),'ok':ok}
assert all(c['ok'] for c in checks.values()), checks
print(json.dumps({'theorem_M_checks':checks,'all_ok':True}))
