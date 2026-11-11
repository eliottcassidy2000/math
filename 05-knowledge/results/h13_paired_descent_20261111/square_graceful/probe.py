#!/usr/bin/env python3
"""Small exact controls for the symbolic source rewrite and sidecars."""
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
COMPILER = ROOT / '05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
spec = importlib.util.spec_from_file_location('compiler', COMPILER)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
g = c.g

Q = [18, 8, 3, 5, 14]
CHARTS = [('Q', Q, 3), ('S5a', [2, 4, 9, 19, 7], 20),
          ('S5b', [2, 4, 25, 3, 11], 18), ('S4', [2, 4, 9, 3], 8)]
def L(b): return 13*g(b)+32*b+6*(b&1)+18
def rev(n,k): return int(f'{n:0{k}b}'[::-1],2) if k else 0
def state(rows):
    """Independent subset paths with endpoint/root marking."""
    n=len(rows); root=n-1; full=(1<<n)-1
    h=c.hp_subsets(rows)
    dp=[[0]*n for _ in range(full+1)]
    dp[1<<root][root]=1
    for mask in range(1,full+1):
        for v,x in enumerate(dp[mask]):
            if not x: continue
            avail=rows[v]&~mask&full
            while avail:
                bit=avail&-avail;avail-=bit
                dp[mask|bit][bit.bit_length()-1]+=x
    u=sum(sum(dp[mask])*h[full^mask] for mask in range(1,full+1))
    # Every full path ending at the root is a path on its complement ending
    # at an in-neighbor. Reverse arcs and count paths starting at the root.
    rr=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
    dr=[[0]*n for _ in range(full+1)];dr[1<<root][root]=1
    for mask in range(1,full+1):
        for v,x in enumerate(dr[mask]):
            if not x: continue
            avail=rr[v]&~mask&full
            while avail:
                bit=avail&-avail;avail-=bit
                dr[mask|bit][bit.bit_length()-1]+=x
    return [h[-1],sum(dr[-1]),u,h[full^(1<<root)]]

def main():
    signatures={name:c.signature(rows,o) for name,rows,o in CHARTS}
    literal=[]
    for k in range(2,6):
        b=(1<<k)-2; h=77*(1<<(k-1))-31
        for name,rows,o in CHARTS:
            addr=4*b+3 if name=='S4' else b
            tournament=c.construct(rows,o,0,addr)
            actual=state(tournament);d=signatures[name][0]
            expected=[h,d,77*(1<<(k-1))+d,d]
            assert actual==expected,(k,name,actual,expected)
            literal.append({'k':k,'chart':name,'B':addr,'order':len(tournament),'state':actual})
    exact_zeros=[]
    for b in range(1<<16):
        expected=b>=2 and (b+2)&(b+1)==0
        for name,_,_ in CHARTS[1:]:
            addr=4*b+3 if name=='S4' else b
            same=c.count(signatures[name],0,addr)==L(b)+2
            assert same==expected,(b,name,same,expected)
        if expected: exact_zeros.append(b)
    shell_pairs=0
    for k in range(2,12):
        for ell in range(k-1):
            t=k-ell-1
            for s in range(1<<ell):
                n0=(1<<ell)+s;n1=(1<<k)-(1<<(ell+1))+s
                d0=g((1<<k)+n0)-g(n0)
                d1=g((1<<k)+n1)-g(n1)
                assert d0==d1
                bridge=102 if not ell else 96*(1<<ell)+39*((1<<ell)-1-rev(s,ell))+13*(1-(s>>(ell-1)))
                assert L((1<<k)+n0)-L(n1)==bridge
                shell_pairs+=1
    fold_cases=0
    for p in range(64):
        for n in range(16):
            u,q,r,v=[(n>>i)&1 for i in range(4)];a=u+v
            f=n+1-u-g(n)
            assert f==2*a*(a-1)+2*a*(q+r)+2*q*r
            assert L(n)==38*a+322*v+90*q+180*r+31-13*f
            assert L(16*p+n)==L(n)+13*g(p)-13+p*(707-13*rev(n,4))+13*(p&1)*(1-v)
            if u==1 and v==0:
                assert L(16*p+n+7)-L(16*p+n)==322+91*p-13*(p&1)
            if q==1 and r==0:
                assert L(16*p+n+2)-L(16*p+n)==90+26*p
            fold_cases+=1
    receipt={'status':'FINITE-EXACT controls of proved formulas','compiler_sha256':hashlib.sha256(COMPILER.read_bytes()).hexdigest(),
             'signatures':signatures,'literal_tournaments':literal,
             'zero_set_universe':'0 <= B < 65536, all three replacement charts', 'zero_set':exact_zeros,
             'shell_pairs':shell_pairs,'shell_universe':'2 <= k <= 11, every interior fibre',
             'fold_cases':fold_cases,'fold_universe':'0 <= prefix P < 64, all sixteen suffixes',
             'hostile_controls':'Every B outside the classified ray in the stated universe fails all three aligned replacements; both suffix-lift orientations are retained.',
             'coverage_claim':False}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('literal_tournaments','signatures')},indent=2))

if __name__=='__main__':main()
