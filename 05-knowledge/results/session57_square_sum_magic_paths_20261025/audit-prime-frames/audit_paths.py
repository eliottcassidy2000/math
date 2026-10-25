import itertools as it, math, json

def need(x,msg):
    if not x:raise RuntimeError(msg)
out={}
for p in (3,5):
    paths=set(); B=tuple(range(p));D=tuple(range(p,2*p))
    for b in it.permutations(B):
        for d in it.permutations(D):
            path=tuple(t for pair in zip(b,d) for t in pair)
            canonical=min(path,path[::-1]);paths.add(canonical)
            partition=frozenset((frozenset(path[::2]),frozenset(path[1::2])))
            need(partition==frozenset((frozenset(B),frozenset(D))),'recover frame')
    need(len(paths)==math.factorial(p)**2,'path count')
    cp=(p+1)*p**(p-1)-p*(p+1)//2; both=p*(p+1)//2
    need(cp==both*(2*p**(p-2)-1),'line fraction')
    if p>=5:need(4*p-2<p*p,'two path obstruction')
    out[p]={'unoriented_paths':len(paths),'line_fraction_denominator':cp//both,'two_paths_maximum_labels':4*p-2,'plane_labels':p*p}
print(json.dumps(out,indent=2,sort_keys=True))
