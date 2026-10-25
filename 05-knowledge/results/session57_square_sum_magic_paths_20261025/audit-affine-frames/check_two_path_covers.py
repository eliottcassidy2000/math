import itertools as it
paths={min(p,p[::-1]) for p in it.permutations(range(6)) if all((a<3)!=(b<3) for a,b in zip(p,p[1:]))}
edges=lambda p:frozenset(frozenset((a,b)) for a,b in zip(p,p[1:]))
complete=frozenset(frozenset((a,b)) for a in range(3) for b in range(3,6))
covers=[(p,q) for p,q in it.combinations(paths,2) if edges(p)|edges(q)==complete]
if len(paths)!=36 or len(covers)!=36:raise RuntimeError('incorrect local count')
if any(len(edges(p)&edges(q))!=1 for p,q in covers):raise RuntimeError('cover overlap')
print('PASS 36 reversal classes, 36 unordered covering pairs, one overlap edge each')
