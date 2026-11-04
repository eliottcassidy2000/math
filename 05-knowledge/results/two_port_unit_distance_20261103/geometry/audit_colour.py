"""Parent's independent literal seven-vertex colour enumeration."""
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def require(ok, message):
    if not ok: raise RuntimeError(message)


def pattern(values):
    seen={}
    return tuple(seen.setdefault(v,len(seen)) for v in values)


def main():
    table=json.loads((HERE/'boundary_tables.json').read_text(encoding='utf8'))
    coords=json.loads((HERE/'coordinates.json').read_text(encoding='utf8'))
    edges=[(i,(i+1)%6) for i in range(6)]+[(i,6) for i in range(6)]
    masses=Counter()
    total=0
    for colors in itertools.product(range(6),repeat=7):
        if all(colors[u]!=colors[v] for u,v in edges):
            masses[pattern(colors[:6])]+=1;total+=1
    # Count internal colours directly, separately from the producer's formula.
    kernel={(a,b):sum(x!=a and x!=y and y!=b
                      for x in range(6) for y in range(6))
            for a in range(6) for b in range(6)}
    for record in table:
        pi=tuple(record['partition']);forest=[]
        for label in set(pi):
            block=[i for i,x in enumerate(pi) if x==label]
            forest.extend((block[0],j) for j in block[1:])
        query=0
        for beta,mass in masses.items():
            factor=1
            for u,v in forest:factor*=kernel[beta[u],beta[v]]
            query+=mass*factor
        require(record['mass']==record['recovered_mass']==masses[pi],('mass',pi))
        require(query==record['query'],('query',pi))
    require(total==24600 and len(masses)==40,'base counts')
    color=coords['colouring3']
    require(len(color)==37 and max(color)<3 and min(color)>=0,'palette')
    require(all(color[u]!=color[v] for u,v in coords['unit_edges']),'literal3color')
    # Sharp twin-neighbor control: opposite square corners have exactly two
    # common unit neighbors; the general <=2 bound has a circle proof.
    square=((0,0),(1,1),(1,0),(0,1))
    distance2=lambda a,b:sum((x-y)**2 for x,y in zip(a,b))
    require(all(distance2(square[i],square[j])==1 for i in (0,1) for j in (2,3)),'K22 witness')
    require(distance2(square[0],square[1])==distance2(square[2],square[3])==2,'square nonedges')
    require(hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==SOURCE,'source changed')
    print(json.dumps({'status':'PASS','source_sha256':SOURCE,'literal_assignments':6**7,
                      'base_colorings':total,'positive_patterns':len(masses),
                      'full_context_checks':len(table),'literal_37_vertex_coloring':True,
                      'sharp_common_neighbor_control':'K2,2 square'},indent=2,sort_keys=True))


if __name__=='__main__':main()
