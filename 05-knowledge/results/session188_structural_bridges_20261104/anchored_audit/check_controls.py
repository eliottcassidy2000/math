"""Independent complete truth-table/Mobius checks of retained larger controls."""
from check import hp,need,ROOT
import json,hashlib
from collections import Counter

def extend(a,x,c):
    n=len(a);b=list(a)+[1<<(n+1),0]
    for u in range(n):
        for mask,p in [(x,n),(c,n+1)]:
            if mask>>u&1:b[u]|=1<<p
            else:b[p]|=1<<u
    return tuple(b)

def complete(a):
    n=len(a);N=1<<(2*n);mask=(1<<n)-1
    values=[hp(extend(a,z&mask,z>>n)) for z in range(N)];coeff=list(values)
    for j in range(2*n):
        bit=1<<j
        for z in range(N):
            if z&bit:coeff[z]-=coeff[z^bit]
    return values,coeff

def main():
    receipt=json.loads((ROOT/'simultaneous_response/receipt.json').read_text());out=[];total=0
    for ex in receipt['examples']:
        T,U=tuple(ex['T']),tuple(ex['U']);n=len(T);basis=[a|(c<<n) for a,c in ex['complete_named_basis']];vs=[];cs=[]
        for a in [T,U]:
            v,c=complete(a);vs.append(v);cs.append(c);total+=len(v)
            need(all(not z or i in basis for i,z in enumerate(c)),'no dropped off-basis coefficient')
        need([cs[0][s] for s in basis]==ex['F_T_vector'],'full source polynomial')
        need([cs[1][s] for s in basis]==ex['F_U_vector'],'full target polynomial')
        diff=[a-b for a,b in zip(*cs)]
        need([diff[s] for s in basis]==ex['difference_vector'],'full difference')
        degrees=dict(sorted(Counter(i.bit_count() for i,z in enumerate(diff) if z).items()))
        need({str(k):v for k,v in degrees.items()}==ex['nonzero_difference_by_degree'],'all degree counts')
        for a,c,x,y in ex['literal_queries']:need((vs[0][a|(c<<n)],vs[1][a|(c<<n)])==(x,y),'retained literal query')
        for a,c,x,y in ex['pair_uniform_queries']:need((vs[0][a|(c<<n)],vs[1][a|(c<<n)])==(x,y),'retained uniform query')
        out.append({'n':n,'modules':ex['original_module_masks'],'full_queries_per_side':len(vs[0]),'nonzero_difference_by_degree':degrees})
    return {'status':'PASS','method':'Independent literal Hamilton counts for every full Boolean attachment, followed by full integer Mobius transform; no producer polynomial routines','literal_counts':total,'controls':out,'receipt_sha256':hashlib.sha256((ROOT/'simultaneous_response/receipt.json').read_bytes()).hexdigest()}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
