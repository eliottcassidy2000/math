"""Independent exact y-frequency counts and forced-prime cover audit."""
import json,hashlib
from pathlib import Path
from itertools import combinations

def need(b,m):
    if not b:raise ArithmeticError(m)

def main():
    root=Path(__file__).resolve().parents[1];raw=root/'ranks/phase2/highrank_inputs.json'
    cert=json.loads((root/'isogeny/certificate.json').read_text());minimum=json.loads((root/'isogeny/minimum_prime_bank.json').read_text())
    need(hashlib.sha256(raw.read_bytes()).hexdigest()==cert['input_sha256'],'source pin')
    data=json.loads(raw.read_text());models={c['id']:list(map(int,c['ainvs'])) for c in data['curves']}
    ids=cert['ids'];need(len(ids)==len(models)==66,'source count');rows={}
    for ps,stored in cert['traces'].items():
        p=int(ps);need(p>=3 and all(p%d for d in range(2,int(p**0.5)+1)),'prime')
        frequency=[]
        for B in range(p):
            f=[0]*p
            for y in range(p):f[(y*y+B*y)%p]+=1
            frequency.append(f)
        row=[]
        for i in ids:
            a1,a2,a3,a4,a6=[v%p for v in models[i]]
            b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6;b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
            disc=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
            if disc%p==0:row.append(None);continue
            count=1+sum(frequency[(a1*x+a3)%p][(x**3+a2*x*x+a4*x+a6)%p] for x in range(p))
            row.append(p+1-count)
        need(row==stored,'trace bank');rows[p]=row
    bank=set(minimum['one_minimum_bank']);pairsets={}
    for a,b in combinations(range(66),2):
        witnesses={p for p,r in rows.items() if r[a] is not None and r[b] is not None and r[a]!=r[b]}
        need(bool(witnesses) and bool(bank&witnesses),'all pairs separated by named bank')
        pairsets[ids[a],ids[b]]=witnesses
    for ps,pair in minimum['forced_prime_witness_pairs'].items():
        need(pairsets[tuple(pair)]=={int(ps)},'forced unique prime')
    need(len(bank)==len(minimum['forced_prime_witness_pairs'])==11,'minimum and uniqueness')
    print(json.dumps({'curves':66,'pairwise_non_Q_isogenous_pairs':len(pairsets),
        'prime_cap':max(rows),'unique_minimum_bank_in_window':sorted(bank),
        'all_trace_rows_independently_recomputed':True},sort_keys=True))
if __name__=='__main__':main()
