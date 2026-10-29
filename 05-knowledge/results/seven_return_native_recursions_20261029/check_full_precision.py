"""Check partial residues against independent full-prefix counts."""
import hashlib,json
import parked_counts as p
import pruned_precision as q
def verify():
    rows=[]
    for n,b in ((101,16),(201,8)):
        value=p.prefix(n)[0][-1];base=(n-1)//2
        if value%(1<<base):raise ValueError('normalization')
        residue=(value>>base)%(1<<b)
        for mode in ('completion','sum'):
            if residue!=q.count(n,b,mode=mode)[0]:raise ValueError('full count vs precision')
        rows.append(dict(n=n,bits=b,residue=residue,full_decimal_sha256=hashlib.sha256(str(value).encode()).hexdigest()))
    return dict(status='PASS',rows=rows)
if __name__=='__main__':print(json.dumps(verify(),indent=2))
