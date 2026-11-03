"""Exact target matching with the old root and its unique in-port fixed."""
from itertools import permutations
import json
from support import HERE,need,save,source

def code(T):return sum(((T[i]>>j)&1)<<k for k,(i,j) in enumerate((i,j) for i in range(len(T)) for j in range(i+1,len(T))))
def main():
    U,r,b,names=source();n=len(U);incoming=[i for i in range(n) if U[i]>>r&1]
    need(len(incoming)==1,'actual unique root in-port');a=incoming[0]
    d=json.loads((HERE/'ready2107.json').read_text());base=code(U);pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    best=100;hits=[];tested=0
    for target_id,rec in enumerate(d['hits']):
        T=tuple(rec['rows']);rr=rec['root'];aa=next(i for i in range(n) if T[i]>>rr&1)
        slots=[i for i in range(n) if i not in (r,a)];other=[i for i in range(n) if i not in (rr,aa)]
        for perm in permutations(other):
            inv=[None]*n;inv[r]=rr;inv[a]=aa
            for i,j in zip(slots,perm):inv[i]=j
            z=sum(((T[inv[i]]>>inv[j])&1)<<k for k,(i,j) in enumerate(pairs));distance=(z^base).bit_count();tested+=1
            if distance<best:best=distance;hits=[];print('new minimum',best,'target',target_id,flush=True)
            if distance==best and len(hits)<80:
                Z=tuple(sum(((T[inv[i]]>>inv[j])&1)<<j for j in range(n)) for i in range(n))
                flips=[(i,j) for i,j in pairs if ((Z[i]>>j)^(U[i]>>j))&1]
                hits.append({'target_id':target_id,'rows':Z,'actual_to_target':inv,'flips':flips})
    need(best<100,'nonempty target bank')
    save('target2107.json',{'status':'FINITE_EXACT_FIXED_ROOT_AND_INPORT_TARGET_BANK','tested':tested,
      'root':r,'spare':b,'inport':a,'names':names,'minimum':best,'witnesses':hits})
    print('PASS target matching',tested,'minimum',best,flush=True)
if __name__=='__main__':main()
