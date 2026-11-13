"""Reproduce the independent complete ready-target audit without leaving a binary."""
from pathlib import Path
from tempfile import TemporaryDirectory
from subprocess import run
import json
import intertwiner_probe as i

HERE=Path(__file__).resolve().parent


def advance(z,w):
    return tuple(sum(a*b for a,b in zip(row,z)) for row in i.matrix(w))


def unpruned_small(N):
    cap=(2*N).bit_length()-1
    bank={1:[(1,1,1,2,1,1)]}
    for n in range(3,cap+1):
        bank[n]=[advance(z,w) for m,rows in list(bank.items())
                 for w in i.words(n-m-1) if len(w)==n-m-1
                 for z in rows]
    counts={n:sum(z[3]<=2*N for z in rows) for n,rows in bank.items()}
    hits=0;tests=0
    for n,rows in bank.items():
        for h,e,s,p,f,_ in rows:
            if p>2*N:continue
            G=s+p;d=G-e-f-1
            for x in range(int(n==1),max(-1,(N-h-d)//G)+1):
                word=tuple(map(int,bin(x)[2:][::-1])) if x else ()
                g=i.coeff(word,1)[2]  # Direct word automaton, not binary-shell inverse.
                tests+=1
                hits+=h*g+G*x+e*(x&1)+d==N
    actual=json.loads((HERE/f'complete_audit_{N}.json').read_text())
    assert all(counts.get(r['order'],0)==r['retained'] for r in actual['layers'])
    assert tests==actual['direct_outer_tests'] and hits==len(actual['hits'])
    return dict(target=N,unfiltered_cores=sum(map(len,bank.values())),filtered_counts=counts,
                direct_outer_tests=tests,hits=hits)


def main():
    with TemporaryDirectory(prefix='.audit_build_',dir=HERE) as temp:
        exe=Path(temp)/'audit'
        run(['c++','-O3','-std=c++17',str(HERE/'complete_audit.cpp'),'-o',str(exe)],check=True)
        for N in (881,1997,53084231,159252697):
            run([str(exe),str(N),str(HERE/f'complete_audit_{N}.json')],check=True)
    comparisons=[]
    for N,path in [(53084231,HERE.parent/'complete_53084231.json'),
                   (159252697,HERE.parent/'constructors/ready_after_defect_two.json')]:
        other=json.loads(path.read_text())
        own=json.loads((HERE/f'complete_audit_{N}.json').read_text())
        by_order={r['order']:r['retained'] for r in own['layers']}
        assert all(by_order[r['order']]==r['retained_cores'] for r in other['layers'])
        assert not own['found']
        comparisons.append(dict(target=N,full_core_layers_agree=True,total_cores=own['total_cores'],
                                independent_direct_outer_tests=own['direct_outer_tests']))
    b=json.loads((HERE/'complete_audit_53084231.json').read_text())
    small=[unpruned_small(N) for N in (881,1997)]
    receipt=dict(status='PASS',full_core_layers_agree=True,total_cores=b['total_cores'],
                 independent_direct_outer_tests=b['direct_outer_tests'],small_unpruned_controls=small,
                 complete_target_comparisons=comparisons)
    (HERE/'complete_audit_controls.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
