"""Independent reconstruction and backward-recursion audit of every saved root family."""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def modules(T):
    return [M for M in range(1,1<<len(T)) if M.bit_count()>1 and all(len({bool(T[v]>>u&1) for u in range(len(T)) if M>>u&1})==1 for v in range(len(T)) if not M>>v&1)]
def xor(T,F):
    R=list(T)
    for u,v in combinations(range(len(T)),2):
        if sum(bool(M>>u&1 and M>>v&1) for M in F)%2:R[u]^=1<<v;R[v]^=1<<u
    return tuple(R)
@lru_cache(None)
def H(T):
    n=len(T);allnames=(1<<n)-1
    @lru_cache(None)
    def end_at(S,v):
        rest=S^(1<<v)
        if not rest:return 1
        return sum(end_at(rest,u) for u in range(n) if rest>>u&1 and T[u]>>v&1)
    return sum(end_at(allnames,v) for v in range(n))
def query(T,a,c):
    n=len(T);R=list(T)+[1<<(n+1),0]
    for u in range(n):
        for v,m in ((n,a),(n+1,c)):
            if m>>u&1:R[u]|=1<<v
            else:R[v]|=1<<u
    return H(tuple(R))
def main():
    p=ROOT/'covered_root_family';raw=(p/'receipt.json').read_bytes();data=json.loads(raw)
    need(raw==(p/'receipt.optimized.json').read_bytes(),'paired source receipts')
    for path,want in data['source_sha256'].items():need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==want,'source pin '+path)
    totals={'families':0,'scalar_strict':0,'scalar_equal':0,'literal_query_values':0,'full_endpoint_values':0,'deep_controls':0}
    for T,records,wantstrict,wantequal in [(tuple(data['nine_name_source']),data['family_records'],296,216),(tuple(data['unequal_transitive_source']),data['unequal_records'],56,72)]:
        full=(1<<len(T))-1;proper=[M for M in modules(T) if M!=full]
        expected={tuple([full]+[M for i,M in enumerate(proper) if bits>>i&1]) for bits in range(1<<len(proper))}
        actual={tuple(r['family']) for r in records};need(len(actual)==len(records) and actual==expected,'complete original family universe')
        h=H(T);strict=equal=0
        for r in records:
            U=xor(T,r['family']);need(list(U)==r['target_rows'],'physical XOR reconstruction')
            uh=H(U);need(uh>=h,'scalar monotonicity control');strict+=uh>h;equal+=uh==h
            a,c=r['query'];values=[query(T,a,c),query(U,a,c)];need(values==r['values'] and values[0]!=values[1],'every saved literal separator');totals['literal_query_values']+=2
            if 'endpoint_battery' in r:
                for a,c,f,g in r['endpoint_battery']:
                    need([query(T,a,c),query(U,a,c)]==[f,g],'complete endpoint battery');totals['full_endpoint_values']+=2
                first=r['endpoint_battery'][:3];last=r['endpoint_battery'][3:]
                need(sum(f-g for a,c,f,g in first)-2*sum(f-g for a,c,f,g in last)==r['battery_contrast']==192,'signed battery contrast')
        need((strict,equal)==(wantstrict,wantequal),'exact strict/equal census')
        totals['families']+=len(records);totals['scalar_strict']+=strict;totals['scalar_equal']+=equal
    old=None
    for r in data['deep_twelve_name_controls']:
        U=tuple(r['target_rows']);T=xor(U,r['family'])
        if old is None:old=T
        need(T==old and all(M in modules(T) for M in r['family']),'same deep source and original modules')
        need(all(not(A&B) or A&B in (A,B) for A,B in combinations(r['family'],2)),'deep laminar family')
        a,c=r['query'];need([query(T,a,c),query(U,a,c)]==r['values'] and r['values'][0]!=r['values'][1] and H(T)==H(U),'deep physical separator')
        totals['literal_query_values']+=2;totals['deep_controls']+=1
    return dict(status='PASS',method='No producer imports. Exhaustive original module subsets, direct physical XOR, independent backward-ending recursion.',counts=totals,source_receipt_sha256=hashlib.sha256(raw).hexdigest(),source_sha256=data['source_sha256'],checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
