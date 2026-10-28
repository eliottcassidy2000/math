"""Combine exact refinement with one stored aligned root row per child class.

Reuses Session88's proven transporter replay; does not materialize every
pointed child's full root records. All production canonization is on parents.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import argparse, gzip, importlib.util, json, sys

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent
sys.path.insert(0,str(RESULTS/'session88_symmetry_certificates_20261028/blockfill'))
import check_blockfill as bf
spec=importlib.util.spec_from_file_location('refine',RESULTS/'root_ports_closed_carriers_20261028/check_carriers.py')
refine=importlib.util.module_from_spec(spec);spec.loader.exec_module(refine)
old=refine.old;need=refine.need


class Carrier:
    def __init__(self,parents,n):
        start=perf_counter();self.n=n;self.engine=refine.Canon(max_order=n)
        normalized=[]
        for t in parents:
            code,_=self.engine.maps(t);rows=old.matrix(n,code)
            repeated,group=self.engine.maps(rows)
            need(repeated==code,'canonical parent idempotence')
            normalized.append((code,rows,group))
        normalized.sort();self.parents={c:t for c,t,g in normalized}
        self.states=[(c,m) for c,t,g in normalized for m in old.orbit_masks(n,g)]
        self.index={k:i for i,k in enumerate(self.states)}
        self.owner=[-1]*len(self.states);self.blocks=[];self.bank=[]
        self.port=[None]*len(self.states)
        for i,key in enumerate(self.states):
            if self.owner[i]>=0:continue
            child=self.child(i)
            row,records=bf.complete_row(self.engine,child)
            targets=tuple((self.index[k],w) for k,w in row)
            need(sum(w for j,w in targets)==n+1 and i in dict(targets),'complete root block')
            need(all(self.owner[j]<0 for j,w in targets),'block overlap')
            block=len(self.blocks)
            self.blocks.append((i,targets,records));self.bank.append(child)
            for j,w in targets:
                self.owner[j]=block
                # Retain a local anchor port, not a duplicated permutation.
                self.port[j]=next(v for v,(k,maps) in enumerate(records) if k==self.states[j])
        need(all(b>=0 for b in self.owner),'pointed coverage')
        self.stats=dict(order=n+1,classes=len(self.bank),pointed_states=len(self.states),
                        stored_anchor_root_records=(n+1)*len(self.blocks),
                        materialized_root_records=(n+1)*len(self.states),
                        classification_requests=(n+1)*len(self.blocks),
                        parent_canon_leaves=self.engine.orders_tested,
                        parent_canon_searches=self.engine.searches,
                        maximum_production_canon_order=n,seconds=perf_counter()-start)

    def child(self,index):
        c,m=self.states[index];return old.extend(self.parents[c],m)

    def row(self,index):return self.blocks[self.owner[index]][1]

    def records(self,index):
        anchor,row,records=self.blocks[self.owner[index]]
        if anchor==index:return records
        # All canonical pointed -> anchor maps at this one root work.
        v=self.port[index];to_anchor=records[v][1][0];p=bf.inverse(to_anchor)
        result=[None]*(self.n+1)
        for root,(key,maps) in enumerate(records):
            result[p[root]]=(key,tuple(sorted(bf.compose(p,g) for g in maps)))
        return tuple(result)

    def apply(self,vector):
        out=[None]*len(vector)
        for anchor,row,records in self.blocks:
            value=sum(w*vector[i] for i,w in row)
            for i,w in row:out[i]=value
        return out

    def audit(self,all_rows=False):
        N=self.n+1
        engine=refine.Canon(max_order=self.n)
        ids=list(range(len(self.states))) if all_rows else sorted({0,len(self.states)-1,*range(0,len(self.states),max(1,len(self.states)//1024))})
        for i in ids:
            row,records=bf.complete_row(engine,self.child(i))
            need(tuple((self.index[k],w) for k,w in row)==self.row(i),'lazy row exact')
            need(records==self.records(i),'ALL lazy aligned root maps exact')
        # Every retained map is checked against its actual relation, including
        # every root at order nine, not merely the sampled lazy expansions.
        maps_checked=0
        for anchor,row,records in self.blocks:
            target=self.child(anchor)
            for root,(key,maps) in enumerate(records):
                source=self.child(self.index[key])
                for p in maps:
                    need(p[-1]==root and all(((source[u]>>v)&1)==((target[p[u]]>>p[v])&1)
                         for u in range(N) for v in range(N)),'retained map is actual pointed isomorphism')
                    maps_checked+=1
        v=[i%17-8 for i in range(len(self.states))];qv=self.apply(v)
        need(self.apply(qv)==[N*x for x in qv],'factored Q squared')
        need(sum(dict(self.row(i))[i] for i in range(len(v)))==N*len(self.blocks),'exact trace')
        return dict(lazy_rows_independently_recomputed=len(ids),all_retained_maps_checked=maps_checked,
                    full_lazy_expansion=all_rows,trace=N*len(self.blocks))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-order',type=int,default=9)
    ap.add_argument('--audit-all-through',type=int,default=8);args=ap.parse_args()
    need(1<=args.max_order<=9,'declared complete universe')
    parents=[()];result=[]
    for N in range(1,args.max_order+1):
        carrier=Carrier(parents,N-1)
        audit=carrier.audit(N<=args.audit_all_through)
        parents=carrier.bank;result.append(dict(**carrier.stats,audit=audit))
        print(json.dumps(result[-1]),flush=True)
    reference=None
    if args.max_order==9:
        path=RESULTS/'relation_tower_frobenius_ports_20261028/tournament_bank.json.gz'
        raw=json.loads(gzip.decompress(path.read_bytes()))
        canonicalizer=refine.Canon()
        got={canonicalizer.maps(t)[0] for t in parents}
        # Clear the literal cache: the second bank must be computed independently.
        canonicalizer.cache.clear()
        expected={canonicalizer.maps(tuple(t))[0] for t in raw['rows']}
        need(len(got)==len(parents)==191536 and got==expected,'entire independently generated nine-bank')
        reference=dict(classes=len(got),canonical_code_sha256=sha256(json.dumps(sorted(got),separators=(',',':')).encode()).hexdigest(),
                       reference_gzip_sha256=sha256(path.read_bytes()).hexdigest())
    receipt=dict(status='FINITE-EXACT',results=result,independent_bank=reference,
                 scope='One anchor row and exact lazy root transport, no black/blue graph at order nine.')
    (HERE/'lazy-root-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
