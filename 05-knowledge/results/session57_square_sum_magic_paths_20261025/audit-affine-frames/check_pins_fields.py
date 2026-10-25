import hashlib,itertools,json,pathlib
P=pathlib.Path(__file__).parent

def require(c,m):
    if not c: raise RuntimeError(m)
for pin in json.loads((P/'SOURCE_PINS.json').read_text()):
    require(hashlib.sha256((P/pin['file']).read_bytes()).hexdigest()==pin['sha256'],'source pin '+pin['file'])
    original=pathlib.Path(pin['origin']); origin=P.parent / original.relative_to('/tmp').parts[0].removeprefix('session57-') / pathlib.Path(*original.relative_to('/tmp').parts[1:])
    require(origin.exists() and hashlib.sha256(origin.read_bytes()).hexdigest()==pin['sha256'],'origin pin '+pin['file'])
r81=next(r for r in json.loads((P/'sources/results.json').read_text())['fields'] if r['q']==81)
records=[r81]+[json.loads((P/f'sources/witness-{q}.json').read_text()) for q in (729,2187)]
for r in records:
    n=r['n'];f=r['modulus']
    require(len(f)==n+1 and f[-1]==1,'monic degree')
    for d in range(1,n//2+1):
        for c in itertools.product(range(3),repeat=d):
            g=list(c)+[1]; rem=f.copy()
            for i in range(n,d-1,-1):
                t=rem[i]%3
                for j in range(d+1): rem[i-d+j]=(rem[i-d+j]-t*g[j])%3
            require(any(rem[:d]),'reducible polynomial')
    print('PASS irreducible',3**n,f)
print('PASS source hashes and original source matches')
