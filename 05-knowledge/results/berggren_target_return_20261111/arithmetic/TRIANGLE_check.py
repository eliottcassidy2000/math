"""Exact C2 coefficient gcd plus seeded hostile k3 tests of triangle chains."""
from pathlib import Path
from math import gcd
from hashlib import sha256
from itertools import product
import argparse,json,random,subprocess,tempfile

HERE=Path(__file__).resolve().parent

def need(ok,why):
    if not ok:raise RuntimeError(why)

def join(blocks):
    n=sum(map(len,blocks));out=[];shift=0
    for block in blocks:
        later=((1<<n)-1)^((1<<(shift+len(block)))-1)
        out.extend((row<<shift)|later for row in block);shift+=len(block)
    return out

def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]

def chain(r,gaps):
    blocks=[tt(gaps[0])]
    for j in range(r):blocks.extend(([2,4,1],tt(gaps[j+1])))
    return join(blocks)

def extend(rows,mask):
    n=len(rows)
    return [row|((1<<n) if not(mask>>i&1) else 0) for i,row in enumerate(rows)]+[mask]

def matching(gaps):
    total=0;last=False
    for gap in gaps[1:-1]:
        take=gap==0 and not last
        total+=take;last=take
    return total

def graph_input(rows):return f'{len(rows)} 0\n'+' '.join(map(str,rows))+'\n'

def run():
    rng=random.Random(2026111103)
    exact=[];evaluations=[];inputs=[];probes=[]
    with tempfile.TemporaryDirectory(prefix='triangle_exact_') as td:
        exact_exe=Path(td)/'exact.exe';audit_exe=Path(td)/'audit.exe'
        for source,out in [('TRIANGLE_exact.cpp',exact_exe),('audit_state.cpp',audit_exe)]:
            subprocess.run(['g++','-O3','-std=c++17',str(HERE/source),'-o',str(out)],check=True,capture_output=True,text=True)
        cases=[(2,(0,0,0)),(2,(0,1,0)),(3,(0,0,0,0)),(3,(0,0,1,0)),
               (3,(0,1,0,0)),(3,(0,1,1,0)),(3,(0,2,1,0)),(3,(1,1,1,1))]
        cases.extend((4,(0,*middle,0)) for middle in product((0,1),repeat=3))
        for r,gaps in cases:
            data=f'{r} '+' '.join(map(str,gaps))+'\n'
            result=subprocess.run([str(exact_exe)],input=data,check=True,capture_output=True,text=True)
            lines=result.stdout.splitlines();n,states,covers,terms,content=map(int,lines[0].split())
            polynomial={int(a):int(b) for a,b in (line.split() for line in lines[1:])}
            need(covers==33**r*3**sum(gaps),'complete ordered-cover universe')
            need(polynomial.get(0)==3**r,'all-sink polynomial constant')
            expected=3**max(r-2-min(2,matching(gaps)),0)
            owner={};at=gaps[0]
            for block in range(r):
                for vertex in range(at,at+3):owner[vertex]=block
                at+=3+gaps[block+1]
            group_gcd={}
            for monomial,coefficient in polynomial.items():
                touched={owner[j % n] for j in range(2*n) if monomial>>j&1 and j % n in owner}
                label=str(len(touched))
                group_gcd[label]=gcd(group_gcd.get(label,0),coefficient)
            # Do not abort on a conjecture mismatch: preserve the exact refutation.
            witness=next((dict(monomial=m,coefficient=c) for m,c in polynomial.items() if c%(3*content)),None)
            exact.append(dict(r=r,gaps=gaps,order=n,cover_endpoint_states=states,total_ordered_three_path_covers=covers,
                              nonzero_coefficients=terms,C2=content,matching_prediction=expected,
                              agrees=content==expected,minimal_3adic_coefficient=witness,
                              coefficient_gcd_by_marked_triangle_count=group_gcd,
                              coefficient_sha256=sha256(result.stdout.encode()).hexdigest()))
            old=chain(r,gaps)
            for trial in range(6):
                assignment=rng.randrange(1<<(2*n+1))
                m0=assignment&((1<<n)-1);m1=(assignment>>n)&((1<<n)-1)
                # y=1 means new0 -> new1, so new1 does not beat new0.
                m1|=(1-((assignment>>(2*n))&1))<<n
                graph=extend(extend(old,m0),m1)
                value=sum(c for m,c in polynomial.items() if m&assignment==m)
                evaluations.append(dict(case=len(exact)-1,assignment=assignment,polynomial_value=value));inputs.append(graph_input(graph))
        for r,gaps,trials in [(4,(0,1,1,1,0),64),(5,(0,1,1,1,1,0),12)]:
            old=chain(r,gaps)
            for trial in range(trials):
                graph=old;masks=[]
                for _ in range(3):
                    mask=rng.randrange(1<<len(graph));masks.append(mask);graph=extend(graph,mask)
                probes.append(dict(r=r,gaps=gaps,trial=trial,masks=masks,rows=graph));inputs.append(graph_input(graph))
        counted=subprocess.run([str(audit_exe)],input=''.join(inputs),check=True,capture_output=True,text=True)
        values=[int(line.split()[0]) for line in counted.stdout.splitlines()]
        need(len(values)==len(inputs),'graph output length')
        for test,value in zip(evaluations,values):need(test['polynomial_value']==value,('independent literal polynomial control',test,value))
        for probe,value in zip(probes,values[len(evaluations):]):probe['H']=value
    summary=[]
    for r in (4,5):
        group=[p for p in probes if p['r']==r];content=3**r
        for p in group:content=gcd(content,p['H'])
        prediction=3**max(r-3,0)
        summary.append(dict(r=r,k=3,all_internal_gaps=1,trials=len(group),sample_gcd=content,
                            matching_prediction=prediction,refuted=content%prediction!=0,
                            first_refutation=next((p for p in group if p['H']%prediction),None)))
    return dict(status='FINITE-EXACT complete C2 polynomials for named sources; seeded k3 hostile probes only; all-k law proved separately in TRIANGLE_ALL_K_PROOF.md',
                exact_C2=exact,literal_polynomial_controls=evaluations,k3_summary=summary,k3_probes=probes,seed=2026111103)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    report=run();data=json.dumps(report,indent=2,sort_keys=True)+'\n';path=HERE/'TRIANGLE_receipt.json'
    if args.write:path.write_text(data,encoding='utf-8')
    else:need(path.read_text(encoding='utf-8')==data,'receipt mismatch')
    print(json.dumps(dict(exact_C2=[{k:v for k,v in x.items() if k in ('r','gaps','C2','matching_prediction','agrees','nonzero_coefficients')} for x in report['exact_C2']],k3=report['k3_summary']),indent=2))

if __name__=='__main__':main()
