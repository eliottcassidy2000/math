#!/usr/bin/env python3
"""Exact Parts seam classification, SAT witness verification, optional native SAT/proof runs.
Reuses hn_udg.py and type_m_probe.py. Requires their numpy/scipy/python-sat dependencies.
Default --verify checks SAT witnesses and exact geometry, NOT the UNSAT proof.
"""
import argparse,gzip,hashlib,json,math,os
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
os.environ['HN_PRIMES']='3,5,11'
import hn_udg as h
from type_m_probe import build
from pysat.solvers import Solver

ROOT=Path(__file__).resolve().parents[3]
RES=ROOT/'05-knowledge/results'
CERT=RES/'hn_two_tag_seam_certificates_20261014.json'
CNF=RES/'hn_two_tag_seam_minus_plus_20261014.cnf'

def require(x,msg):
    if not x:raise RuntimeError(msg)

def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def encode(p):return [[str(c) for c in a] for a in p]
def encode_tag(t):return tuple(tuple(str(c) for c in a) for a in t)

def source_check(cert):
    for name,want in cert['source_blobs'].items():
        p=Path(h.__file__).parent/name;b=p.read_bytes()
        got=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
        require(got==want,'toolkit source blob mismatch: '+name)

def geometry():
    L,S,G=build('rho');RS=h.rotate(S,h.RHO)
    require((len(L),len(S),len(G))==(727,361,1087),'piece cardinality mismatch')
    require(h.is_unit(h.RHO),'rotation is not unit')
    pts,xy,E=h.unit_graph(G);idx={p:i for i,p in enumerate(pts)}
    li={idx[p] for p in L};si={idx[p] for p in RS};internal=[];tagged={}
    norms=[h.cnorm(p) for p in pts]
    for u,v in E:
        if (u in li and v in li) or (u in si and v in si):internal.append((u,v))
        else:
            if u not in li:u,v=v,u
            require(u in li-si and v in si-li,'cross-edge membership error')
            tagged.setdefault(encode_tag((norms[u],norms[v])),[]).append((u,v))
    tags=[{'tag':list(t),'edges':tagged[t]} for t in sorted(tagged)]
    # Completeness check extends the toolkit candidate search with an exact
    # full-pair rational-component rejection, using the same field table.
    den=math.lcm(*(q.denominator for p in pts for a in p for q in a))
    ip=[[[int(q*den) for q in a] for a in p] for p in pts]
    weights=[h.MULC[i][i] for i in range(h.DIM)];exact=[];cand=0
    for u,v in combinations(range(len(pts)),2):
        tr=sum(weights[i]*((ip[u][0][i]-ip[v][0][i])**2+(ip[u][1][i]-ip[v][1][i])**2) for i in range(h.DIM))
        if tr!=den*den:continue
        cand+=1
        if h.is_unit(h.csub(pts[u],pts[v])):exact.append((u,v))
    require(exact==E,'exact all-pair check disagrees with toolkit')
    require((len(E),len(internal),[len(t['edges']) for t in tags])==(5892,5862,[6,6,18]),'edge/tag census mismatch')
    require(li&si=={idx[h.C0]},'overlap is not exactly the origin')
    return pts,idx,E,internal,tags,cand

def geometric_symmetries(pts,idx,E,tags):
    edge_set={tuple(sorted(e)) for e in E}
    # Real Galois conjugation, not a claimed Euclidean isometry.
    sigma=lambda p:tuple(tuple(-c if k&4 else c for k,c in enumerate(a)) for a in p)
    gp=[idx[sigma(p)] for p in pts]
    rp=[idx[h.cmul(p,h.HEX[1])] for p in pts]
    for p in [gp,rp]:
        require(len(set(p))==len(pts),'not a permutation')
        require({tuple(sorted([p[u],p[v]])) for u,v in E}==edge_set,'not a graph automorphism')
    for i,j in [(0,1),(1,0),(2,2)]:
        require({(gp[u],gp[v]) for u,v in tags[i]['edges']}=={tuple(e) for e in tags[j]['edges']},'wrong Galois tag swap')
    zeta=(h.el(r3=F(1,2)),h.el(r1=F(1,2)))
    rm=(h.el(r11=F(1,2),r3=-F(1,6)),h.ZERO)
    rr=(h.el(r11=F(1,2),r3=F(1,6)),h.ZERO)
    rays=[h.cmul(zeta,z) for z in h.HEX]
    for t,(a,b) in [(0,(rm,rr)),(1,(rr,rm))]:
        predicted={(idx[h.cmul(a,z)],idx[h.cmul(h.RHO,h.cmul(b,z))]) for z in rays}
        require(predicted=={tuple(e) for e in tags[t]['edges']},'six-ray formula mismatch')
    circle=h.rotate(h.H_rot([-1,0,1]),(h.el(r1=2),h.ZERO))-{h.C0}
    require({(idx[p],idx[h.cmul(h.RHO,p)]) for p in circle}=={tuple(e) for e in tags[2]['edges']},'radius-two formula mismatch')
    return gp,rp

def selected_edges(internal,tags,active):return internal+[e for i in active for e in tags[i]['edges']]
def colouring_check(cols,n,edges):
    require(len(cols)==n and all(c in range(4) for c in cols),'invalid colour vector')
    require(all(cols[u]!=cols[v] for u,v in edges),'bad SAT colour certificate')

def clauses_for(n,edges,tri):
    x=lambda v,c:4*v+c+1
    return [[x(v,c) for c in range(4)] for v in range(n)]+[[-x(u,c),-x(v,c)] for u,v in edges for c in range(4)]+[[x(v,c)] for c,v in enumerate(tri)]

def cnf_bytes(n,clauses):
    return (f'p cnf {4*n} {len(clauses)}\n'+''.join(' '.join(map(str,c))+' 0\n' for c in clauses)).encode()

def verify(cert,pts,idx,E,internal,tags,cand):
    require(digest([encode(p) for p in pts])==cert['points_sha256'],'point ordering/digest mismatch')
    require(digest(E)==cert['edges_sha256'],'edge digest mismatch')
    require(digest(internal)==cert['internal_edges_sha256'],'internal edge digest mismatch')
    require(json.loads(json.dumps(tags))==cert['tags'],'tag data mismatch')
    gp,rp=geometric_symmetries(pts,idx,E,tags);n=len(pts)
    tri=cert['fixed_triangle'];iset={tuple(sorted(e)) for e in internal}
    require(all(tuple(sorted(e)) in iset for e in combinations(tri,2)),'fixed triangle missing')
    sat=[]
    for q in cert['queries']:
        if q['status']=='SAT':
            colouring_check(q['colouring'],n,selected_edges(internal,tags,q['active_tags']));sat.append(q['name'])
        else:require(q['status']=='UNSAT_SOLVER_VERDICT' and not q['proof_checked'],'incorrect UNSAT status')
    anchor=selected_edges(internal,tags,[2]);six=tags[0]['edges']
    seen=set()
    for r in cert['edge_deletion_certificates']:
        e=tuple(r['deleted_edge']);seen.add(e)
        colouring_check(r['colouring'],n,anchor+[s for s in six if tuple(s)!=e])
        require(r['colouring'][e[0]]==r['colouring'][e[1]],'deleted edge unexpectedly satisfied')
        # Conjugate witness gives the corresponding opposite-tag deletion.
        transformed=[None]*n
        for v,w in enumerate(gp):transformed[w]=r['colouring'][v]
        ge=(gp[e[0]],gp[e[1]])
        colouring_check(transformed,n,anchor+[s for s in tags[1]['edges'] if tuple(s)!=ge])
    require(seen=={tuple(e) for e in six},'incomplete edge-deletion coverage')
    cls=clauses_for(n,selected_edges(internal,tags,[2,0]),tri);b=cnf_bytes(n,cls)
    require(b==CNF.read_bytes(),'committed standalone CNF mismatch')
    require(hashlib.sha256(b).hexdigest()==cert['native_proof']['files']['two_tag_minus_plus.cnf']['sha256'],'CNF digest mismatch')
    out={'scope':'Exact geometry and SAT witnesses verified; UNSAT remains solver-reported, native trace unchecked.',
      'vertices':n,'edges':len(E),'internal_edges':len(internal),'cross_edges':30,'tag_counts':[6,6,18],
      'full_pair_count':n*(n-1)//2,'rational_component_candidates':cand,'exact_all_pair_agreement':True,
      'verified_sat_cases':sat,'verified_single_edge_deletions_per_additional_tag':6,
      'galois_tag_swap_checked':True,'six_ray_formulas_checked':True,
      'two_tag_verdicts':{'minus_plus':'UNSAT_SOLVER_VERDICT','plus_minus':'UNSAT_SOLVER_VERDICT'},
      'native_proof_checked':False,'compressed_trace_in_git':False,
      'standalone_cnf_variables':4*n,'standalone_cnf_clauses':len(cls)}
    print(json.dumps(out,indent=2,sort_keys=True))

def solve(cert,pts,idx,E,internal,tags,cand):
    n=len(pts);col=h.Colourer(n,internal,k=4,solver='cadical195',fixed_triangle=cert['fixed_triangle'])
    sels=[5*n+i+1 for i in range(3)]
    for i,t in enumerate(tags):
        for u,v in t['edges']:
            for c in range(4):col.solver.add_clause([-sels[i],-col.x(u,c),-col.x(v,c)])
    for name,active in [('internal_only',[]),('radius_four_control',[2]),('two_tag_minus_plus',[2,0]),('two_tag_plus_minus',[2,1]),('all_tags',[0,1,2])]:
        ok=col.solver.solve(assumptions=col.sel+[s if i in active else -s for i,s in enumerate(sels)])
        out={'name':name,'active_tags':active,'status':'SAT' if ok else 'UNSAT_SOLVER_VERDICT'}
        if ok:
            model=set(col.solver.get_model());cs=[next(c for c in range(4) if col.x(v,c) in model) for v in range(n)]
            colouring_check(cs,n,selected_edges(internal,tags,active));out['colouring']=cs
        print(json.dumps(out,sort_keys=True),flush=True)
    col.close()

def emit_proof(cert,pts,idx,E,internal,tags,cand,outdir):
    outdir.mkdir(parents=True,exist_ok=True);n=len(pts)
    cls=clauses_for(n,selected_edges(internal,tags,[2,0]),cert['fixed_triangle'])
    (outdir/'two_tag_minus_plus.cnf').write_bytes(cnf_bytes(n,cls))
    with Solver(name='cadical195',bootstrap_with=cls,with_proof=True) as solver:
        require(not solver.solve(),'standalone unexpectedly SAT')
        proof=solver.get_proof();require(proof is not None,'solver returned no native trace')
        raw=('\n'.join(proof)+'\n').encode();(outdir/'two_tag_minus_plus.drat').write_bytes(raw)
        compressed=gzip.compress(raw,compresslevel=9,mtime=0);(outdir/'two_tag_minus_plus.drat.gz').write_bytes(compressed)
        print(json.dumps({'status':'UNSAT_SOLVER_VERDICT','proof_checked':False,'raw_bytes':len(raw),'compressed_bytes':len(compressed),'raw_sha256':hashlib.sha256(raw).hexdigest(),'compressed_sha256':hashlib.sha256(compressed).hexdigest()},indent=2,sort_keys=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group()
    group.add_argument('--verify',action='store_true');group.add_argument('--solve',action='store_true');group.add_argument('--emit-proof',action='store_true')
    parser.add_argument('--output-dir',type=Path,default=Path('hn-two-tag-local-proof'))
    args=parser.parse_args();cert=json.loads(CERT.read_text());source_check(cert);data=geometry()
    if args.solve:solve(cert,*data)
    elif args.emit_proof:emit_proof(cert,*data,args.output_dir)
    else:verify(cert,*data)
