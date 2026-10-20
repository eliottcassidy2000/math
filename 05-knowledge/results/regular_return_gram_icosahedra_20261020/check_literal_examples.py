"""Source-aware audit of every literal path/code example in the proof packet."""
import ast
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
from independent_audit import tournament,code_of,words,directed,complement_path,degree,directed_triangles,require

ROOT=Path(__file__).resolve().parent

def tuple_word(text):
    return ast.literal_eval(text)


def main():
    proof=(ROOT/'regular_return_theorems.md').read_text()
    gram=(ROOT/'gram_return_addendum.md').read_text()
    audit=(ROOT/'AUDIT.md').read_text()
    R=frozenset((u,(u+d)%5) for u in range(5) for d in (1,2))
    require(code_of(R,5)==947,'carousel code')
    table=[]
    for ptext,qtext in re.findall(r'^\| (\([0-9, ]+\)) \| (.*?) \|$',proof,re.M):
        p=tuple_word(ptext)
        require(len(p)==5 and set(p)==set(range(5)) and directed(R,p),'displayed P is not a carousel Hamiltonian path: '+str(p))
        S=complement_path(R,p)
        actual=sorted(q for q in words(S,5) if(q[0],q[-1])==(p[0],p[-1]))
        declared=[p]+[tuple_word(q) for q in re.findall(r'\([0-9, ]+\)',qtext)]
        require(sorted(declared)==actual,'displayed continuation list mismatch')
        table.append(dict(P=p,Q=actual,intermediate_code=code_of(S,5)))
    require(len(table)==3,'three table rows')
    require({tuple(row['P']) for row in table}=={p for p in words(R,5) if p[0]==0},'table completeness')
    bad=(0,1,3,4,2)
    missing=[uv for uv in zip(bad,bad[1:]) if uv not in R]
    require(missing==[(4,2)],'original typo hostile')
    m=re.search(r'P=(\([0-9, ]+\)),\s+Q=(\([0-9, ]+\))\s+\(11\)',proof)
    require(m is not None,'literal pair (11) not found')
    p,q=map(tuple_word,m.groups());S=complement_path(R,p)
    require(directed(S,q),'literal pair Q invalid')
    U=complement_path(S,q)
    require(degree(U,5)==(2,)*5 and len(R-U)==6,'literal pair return')
    pair11=dict(P=p,Q=q,target=code_of(U,5))
    m=re.search(r'five directed triangles are\s+([0-9, .]+)\n',proof)
    declared_tri={tuple(map(int,t)) for t in re.findall(r'\b\d{3}\b',m.group(1))}
    actual_tri={tuple(sorted(t)) for t in directed_triangles(R,5)}
    require(declared_tri==actual_tri,'literal triangle list')
    m=re.search(r'Exactly five pairs meet in one vertex:\s+([^\n]+)',proof)
    pairs={tuple(sorted((tuple(map(int,a)),tuple(map(int,b))))) for a,b in re.findall(r'\((\d{3}),(\d{3})\)',m.group(1))}
    expected={tuple(sorted((a,b))) for a,b in combinations(actual_tri,2) if len(set(a)&set(b))==1}
    require(pairs==expected and len(pairs)==5,'literal triangle pair list')
    m=re.search(r'N=(\d+), S=(\d+),\s+U=(\([0-9, ]+\)),\s+L=(\([0-9, ]+\))',proof)
    north,south=int(m[1]),int(m[2]);upper,lower=tuple_word(m[3]),tuple_word(m[4])
    cert={north:set(upper),south:set(lower)}
    for i in range(5):
        cert[upper[i]]={north,upper[(i-1)%5],upper[(i+1)%5],lower[i],lower[(i+1)%5]}
        cert[lower[i]]={south,lower[(i-1)%5],lower[(i+1)%5],upper[i],upper[(i-1)%5]}
    require(len(cert)==12,'certificate distinct vertices')
    for t in cert:
        require(degree(tournament(5,t),5)==(2,)*5,'certificate regularity')
        require(cert[t]=={u for u in cert if(t^u).bit_count()==6},'certificate row')
        require(t^1023 in cert,'converse antipode in sheet')
    require(north==947 and south==(north^1023),'literal poles')
    rotated=set()
    for k in range(5):
        rotated.add(code_of(frozenset(((u+k)%5,(v+k)%5) for u,v in U),5))
    require(rotated==set(upper) and {u^1023 for u in upper}==set(lower),'rotation/converse code provenance')
    m=re.search(r'outgoing masks (\([0-9, ]+\)), code(\d+)',gram)
    masks=tuple_word(m[1]);sid=int(m[2]);a=tournament(5,sid)
    require(masks==tuple(sum(1<<v for u,v in a if u==i) for i in range(5)),'literal intermediate masks/code')
    m=re.search(r'endpoints are s=(\d+),t=(\d+), its interior cycle is ([0-9>-]+)',gram)
    start,end=int(m[1]),int(m[2]);cycle=tuple(map(int,m[3].split('->')))
    require(degree(a,5)[start]==3 and degree(a,5)[end]==1 and directed(a,cycle),'literal endpoints/cycle')
    m=re.search(r'neighbours have codes([0-9,]+), the face',gram)
    neighbors={int(x) for x in m[1].split(',')}
    actual={code_of(complement_path(a,q),5) for q in words(a,5) if(q[0],q[-1])==(start,end)}
    require(neighbors==actual=={north,upper[0],upper[1]},'literal face neighbors')
    m=re.search(r'P=(\([0-9, ]+\)), Q=(\([0-9, ]+\))\.',audit)
    p7,q7=map(tuple_word,m.groups());R7=frozenset((u,(u+d)%7) for u in range(7) for d in (1,2,3))
    U7=complement_path(complement_path(R7,p7),q7)
    require(degree(U7,7)==(3,)*7 and len(R7-U7)==6,'literal audit seven-vertex return')
    cycle7=(1,3,6,2,4,5,1)
    require(directed(R7-U7,cycle7),'literal audit six-cycle')

    svg=ET.fromstring((ROOT/'regular_return_icosahedra.svg').read_bytes())
    ns={'s':'http://www.w3.org/2000/svg'}
    labels=[t for t in svg.findall('s:text',ns) if (t.text or '').isdigit()]
    ellipses=svg.findall('s:ellipse',ns)
    require(len(labels)==len(ellipses)==24,'SVG 24 labels/nodes')
    position={}
    for e in ellipses:
        x,y=float(e.attrib['cx']),float(e.attrib['cy'])
        hits=[t for t in labels if abs(float(t.attrib['x'])-x)<1e-6 and abs(float(t.attrib['y'])-y-4)<1e-6]
        require(len(hits)==1,'SVG label alignment')
        position[(x,y)]=int(hits[0].text)
    regular={c for c in range(1024) if degree(tournament(5,c),5)==(2,)*5}
    require(set(position.values())==regular,'SVG all 24 regular codes')
    edges=[]
    for line in svg.findall('s:line',ns):
        a=position[(float(line.attrib['x1']),float(line.attrib['y1']))]
        b=position[(float(line.attrib['x2']),float(line.attrib['y2']))]
        require((a^b).bit_count()==6,'SVG edge six-arc condition')
        edges.append(tuple(sorted((a,b))))
    require(len(edges)==len(set(edges))==60,'SVG exactly 60 distinct edges')
    require(set(edges)=={(a,b) for a,b in combinations(sorted(regular),2) if(a^b).bit_count()==6},'SVG complete actual graph')
    require(set(Counter(v for e in edges for v in e).values())=={5},'SVG degrees')
    # Supplementary finite classification requested after the literal repair.
    def intermediate_metrics(arcs):
        hp=words(arcs,5)
        aut=sum(frozenset((pi[u],pi[v]) for u,v in arcs)==arcs for pi in permutations(range(5)))
        blue=0
        for p in hp:
            rev={p[i]:p[4-i] for i in range(5)}
            blue+=all((rev[v],rev[u]) in arcs for u,v in arcs)
        require(blue%aut==0,'blue path orbit integrality')
        return len(hp),aut,blue//aut
    supports={}
    for t in regular:
        source=tournament(5,t)
        for p in words(source,5):
            sid2=code_of(complement_path(source,p),5)
            supports.setdefault(sid2,set()).add(t)
    kinds=Counter((len(support),)+intermediate_metrics(tournament(5,c)) for c,support in supports.items())
    require(kinds=={(1,11,1,3):120,(1,13,1,3):120,(3,15,3,1):40},'intermediate type refinement')
    for row in table:
        h,aut,blue=intermediate_metrics(tournament(5,row['intermediate_code']))
        row.update(intermediate_H=h,intermediate_aut_order=aut,intermediate_blue_fibre=blue)
    result=dict(status='PASS',proof_table=table,
                intermediate_type_control=[dict(column_degree=c,H=h,aut_order=a,blue_fibre=b,labelled_columns=count) for (c,h,a,b),count in sorted(kinds.items())],original_typo_hostile=dict(P=bad,missing_arcs=missing),
                literal_pair_11=pair11,
                directed_triangles=sorted(actual_tri),triangle_pairs=sorted(pairs),
                certificate_vertices=sorted(cert),literal_intermediate=dict(code=sid,rows=masks,neighbors=sorted(neighbors)),
                audit_seven_vertex_example=dict(P=p7,Q=q7,target=code_of(U7,7)),
                figure=dict(labels=24,distinct_edges=60,all_regular_codes=True,all_edges_exact=True,
                            visual_review='Final-reviewer visual inspection; this checker validates the literal labels and edges only.'),
                checked_file_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in ['regular_return_theorems.md','gram_return_addendum.md','AUDIT.md','regular_return_icosahedra.svg']},
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
