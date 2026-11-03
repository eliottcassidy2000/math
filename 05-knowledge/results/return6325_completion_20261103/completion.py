"""Exact minimum one-birth observer banks on small structural-halving fibres."""
from pathlib import Path
import sys,json,subprocess,itertools
HERE=Path(__file__).resolve().parent
def need(b,s):
    if not b:raise ValueError(s)
def base(q):
    n=3*q
    return tuple(sum(1<<j for j in range(n) if (i//q==j//q and i<j) or (j//q-i//q)%3==1) for i in range(n))
def main():
    records=[]
    for q in (2,4):
        T=base(q);n=len(T);allmask=(1<<n)-1
        p=subprocess.run([sys.argv[1]],input=' '.join(map(str,[n,*T])),capture_output=True,text=True)
        need(p.returncode==0,p.stderr);values=json.loads(p.stdout)
        maps=[]
        for bits in range(8):
            perm=[]
            for b in range(3):
                local=[q-1,*range(q-1)] if bits>>b&1 else list(range(q))
                perm.extend(b*q+i for i in local)
            maps.append(perm)
        pairs=list(itertools.combinations(range(8),2));full=(1<<len(pairs))-1;bank={}
        for mask in range(1,allmask):
            response=tuple(values[sum(1<<perm[i] for i in range(n) if mask>>i&1)] for perm in maps)
            cover=sum(1<<j for j,(a,b) in enumerate(pairs) if response[a]!=response[b])
            bank.setdefault(cover,(mask,response))
        chosen=None
        for cover,row in bank.items():
            if cover==full:chosen=[row];break
        if chosen is None:
            for (c,a),(d,b) in itertools.combinations(bank.items(),2):
                if c|d==full:chosen=[a,b];break
        if chosen is None:
            for triple in itertools.combinations(bank.items(),3):
                if triple[0][0]|triple[1][0]|triple[2][0]==full:chosen=[a[1] for a in triple];break
        need(chosen is not None,'three single-block reads should recover the orientation bits')
        out=dict(block_size=q,source_rows=T,halved_rows=base(q//2),states=maps,
            H=values[0],legal_contexts=allmask-1,distinct_pair_covers=len(bank),completion_number=len(chosen),
            chosen=[dict(mask=a,responses=b) for a,b in chosen],
            maximum_single_context_outputs=max(len(set(row[1])) for row in bank.values()))
        partitions=[]
        for perm in maps:
            fibres=[]
            for b in range(3):
                ordered=sorted(range(b*q,(b+1)*q),key=lambda i:perm[i])
                fibres.extend(tuple(sorted(ordered[j:j+2])) for j in range(0,q,2))
            partitions.append(sorted(fibres))
        out['first_halving_fibres']=partitions
        out['baseline_scope']='final quotient and unordered original block fibres; intermediate pair fibres omitted for q=4'
        out['distinct_intermediate_pair_partitions']=len({tuple(map(tuple,z)) for z in partitions})
        if q==4:need(out['distinct_intermediate_pair_partitions']==8,'intermediate fibres already recover the pole phase');out['completion_with_intermediate_pairs']=0
        else:need(out['distinct_intermediate_pair_partitions']==1,'unordered pair fibres retain no pair orientation')
        records.append(out);print(q,'completion',len(chosen),'chosen',chosen,'partitions',len(bank),flush=True)
    (HERE/'completion.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
