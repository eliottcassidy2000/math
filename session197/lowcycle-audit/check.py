"""Verify independent census, labelled orbit mass, and all six core claims."""
import collections, itertools, json, math, pathlib
ROOT=pathlib.Path(__file__).resolve().parent

def cells(rows):
    groups=collections.defaultdict(list)
    for i,r in enumerate(rows):groups[r.bit_count()].append(i)
    return [groups[d] for d in sorted(groups)]

def canonical(rows):
    n=len(rows);best=None
    for choices in itertools.product(*(itertools.permutations(c) for c in cells(rows))):
        p=sum(choices,())
        code=sum(1<<(8*i+j) for i in range(n) for j in range(n) if rows[p[i]]&(1<<p[j]))
        if best is None or code<best:best=code
    return best

def autos(rows):
    n=len(rows);cc=cells(rows);total=0
    for pp in itertools.product(*(itertools.permutations(c) for c in cc)):
        p=list(range(n))
        for c,q in zip(cc,pp):
            for u,v in zip(c,q):p[u]=v
        total+=all(bool(rows[u]&(1<<v))==bool(rows[p[u]]&(1<<p[v])) for u in range(n) for v in range(n))
    return total

def literal(rows):
    n=len(rows);h=c=0
    for p in itertools.permutations(range(n)):
        if all(rows[u]&(1<<v) for u,v in zip(p,p[1:])):
            h+=1
            c+=p[0]==0 and bool(rows[p[-1]]&1)
    return h,c

def transitive_order(rows,vertices):
    # Topological source stripping rather than sorting induced outdegrees.
    remaining=set(vertices);order=[]
    while remaining:
        sources=[v for v in remaining if all(u==v or rows[v]&(1<<u) for u in remaining)]
        if not sources:return None
        assert len(sources)==1
        v=sources[0];order.append(v);remaining.remove(v)
    return order

def main():
    independent=json.loads((ROOT/'independent_census.json').read_text())
    source=json.loads((ROOT/'source'/'certificate.json').read_text())
    mass=[]
    for n,(actual,expected) in enumerate(zip(independent['inventory'],source['inventory']),1):
        assert actual['classes']==expected['classes']
        assert actual['C1_classes']==expected['C1_classes'] and actual['C2_classes']==expected['C2_classes']
        got={(canonical(t['rows']),t['C'],t['strong']) for t in actual['H67']}
        want={(canonical(t['rows']),t['C'],t['strong']) for t in expected['H67']}
        assert got==want
        orbit_mass=0
        for line in (ROOT/f'classes{n}.txt').read_text().splitlines():
            key,h,c,st=map(int,line.split());rows=tuple((key>>(8*i))&255 for i in range(n))
            aut=autos(rows);assert aut>0
            orbit_mass+=math.factorial(n)//aut
        assert orbit_mass==1<<(n*(n-1)//2)
        mass.append({'order':n,'unlabelled':actual['classes'],'labelled_orbit_mass':orbit_mass})
        print('Orbit mass',mass[-1],flush=True)
    keys={canonical(c['rows']):c['id'] for c in source['strong_cores']}
    details=[]
    for core in source['strong_cores']:
        rows=tuple(core['rows']);n=len(rows)
        H,C=literal(rows);assert H==67 and C==core['C']
        assert autos(rows)==1
        assert all(transitive_order(rows,[v for v in range(n) if v!=u]) is None for u in range(n))
        presentations=[]
        for a,b in itertools.permutations(range(n),2):
            if not rows[a]&(1<<b):continue
            base=transitive_order(rows,[v for v in range(n) if v not in (a,b)])
            if base is not None:
                word=[bool(rows[a]&(1<<v))+2*bool(rows[b]&(1<<v)) for v in base]
                presentations.append((a,b,tuple(base),tuple(word)))
        saved={(x['roots'][0],x['roots'][1],tuple(x['base_order']),tuple(x['word'])) for x in core['two_apex_presentations']}
        assert saved==set(presentations) and saved
        rev=tuple(sum(1<<v for v in range(n) if rows[v]&(1<<u)) for u in range(n))
        assert keys[canonical(rev)]==core['converse_id']
        details.append({'id':core['id'],'literal_H':H,'literal_C':C,'automorphisms':1,
                        'feedback_vertex_number':2,'apex_presentations':len(saved),'converse_id':core['converse_id']})
    actual8={canonical(x['rows']) for x in independent['inventory'][7]['H67']}
    padded=set()
    for core in source['strong_cores']:
        rows=tuple(core['rows'])
        # A new last-labelled source and a new last-labelled sink.
        padded.add(canonical(rows+((1<<7)-1,)))
        padded.add(canonical(tuple(r|(1<<7) for r in rows)+(0,)))
    assert actual8==padded and len(padded)==12
    result={'status':'PASS','orbit_mass_checks':mass,'strong_core_checks':details,
            'order8_targets_are_exactly_source_sink_paddings':True,
            'official_bound':'Busch 2006 Theorem 1, official PDF page 3; H^3 >= 5^(n-1)'}
    (ROOT/'AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
