#!/usr/bin/env python3
"""Independent proof-family and saved-witness audit, no researcher imports."""
import csv, hashlib, itertools, json
from functools import lru_cache
from pathlib import Path
from audit_binary_cycles import binary_matrix,cycles,brute_cycles,start_paths,strong,masks
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
LEVELS=next(p for p in (BASE/"cycle-levels",BASE/"session212-cycle-levels") if p.exists())
SPECTRUM=next(p for p in (BASE/"cycle-spectrum",BASE/"session212-cycle-spectrum") if p.exists())

def ck(ok,msg='check failed'):
    if not ok:raise RuntimeError(msg)

def matrix(rows):
    n=len(rows)
    a=[[bool(rows[i]>>j&1) for j in range(n)] for i in range(n)]
    ck(all(not a[i][i] for i in range(n)),'loop')
    ck(all(a[i][j]!=a[j][i] for i in range(n) for j in range(i)),'incomplete orientation')
    return a

def gatefree(a):
    n=len(a);outs=list(map(sum,a));ins=[sum(a[i][j] for i in range(n)) for j in range(n)]
    return not any(a[i][j] and outs[i]==ins[j]==1 for i in range(n) for j in range(n))

def comb(m,q,tau,direct):
    a=binary_matrix(q+1,m)
    z=m+2;x=m+3
    a[1][z]=bool(tau);a[z][1]=not tau
    a[x][z]=bool(direct);a[z][x]=not direct
    return a

def near_transitive(m):
    a=[[i<j for j in range(m)] for i in range(m)]
    a[0][-1]=False;a[-1][0]=True
    return a

def endpoint_gate(m,b):
    a=[r+[False, bool(b>>i&1)] for i,r in enumerate(near_transitive(m))]
    a+=[[True]*m+[False,False], [not bool(b>>j&1) for j in range(m)]+[True,False]]
    return a

def close_end(a,root):
    n=len(a)
    return [r+[i==root] for i,r in enumerate(a)]+[[i!=root for i in range(n)]+[False]]

def pad_last(a,extra):
    n=len(a)-1;k=extra+1
    return [r[:n]+[r[-1]]*k for r in a[:n]] + [a[-1][:n]+[i<j for j in range(k)] for i in range(k)]


def main():
    out={'status':'PASS','checks':{}}
    cnt=0
    for m in range(1,9):
        for q in range(2**m):
            for tau in (0,1):
                for direct in (0,1):
                    a=comb(m,q,tau,direct); expected=q+tau+direct
                    ck(cycles(a)==expected,('comb',m,q,tau,direct))
                    ck(strong(a)==bool(expected),('strong comb',m,q,tau,direct))
                    if not direct:ck(gatefree(a),'reversed comb gate')
                    else:ck(gatefree(a)!=(tau==1 and q==2**m-1),'forward gate boundary')
                    if m<=4:ck(brute_cycles(a)==expected,'brute comb')
                    cnt+=1
    out['checks']['comb_both_directions_m_1_to_8']=cnt
    out['checks']['comb_permutation_oracle_m_1_to_4']=120
    # Record the honest m=0 direct-orientation degeneracy.
    for direct,tau in itertools.product((0,1),repeat=2):
        ck(cycles(comb(0,0,tau,direct))==tau,'m0 degeneration')
    out['checks']['comb_m0_boundary']=4
    cnt=0
    for m in range(3,9):
        t=near_transitive(m);rev=list(map(list,zip(*t)))
        expected=[1]+[2**(j-1) for j in range(1,m-1)]+[1]
        ck([start_paths(rev,j) for j in range(m)]==expected,'endpoint profile')
        observed=set()
        for b in range(2**m):
            a=endpoint_gate(m,b); val=sum(e for j,e in enumerate(expected) if b>>j&1)
            ck(cycles(a)==val,'endpoint gate count');observed.add(val)
            ck(strong(a)==bool(val),'endpoint gate strong')
            # Delete v and old0: remainder is the source z followed by a chain.
            idx=[i for i in range(m+1) if i!=0]
            ck(all(not (a[i][j] and a[j][k] and a[k][i]) for i,j,k in itertools.permutations(idx,3)),'FVS upper bound')
            cnt+=1
        ck(observed==set(range(2**(m-2)+2)),'endpoint subset interval')
    out['checks']['endpoint_gate_all_subsets_m_3_to_8']=cnt
    # Existing research output is data only; verify each physical witness afresh.
    C={};E={};saved=0;closures=0
    for n in range(1,10):
        C[n]=set();E[n]=set()
        for kind in ('cycles','endpoints'):
            for row in csv.DictReader((LEVELS/f'n{n}_{kind}.tsv').open(),delimiter='\t'):
                val=int(row['value']);a=matrix(list(map(int,row['rows'].split(','))))
                if kind=='cycles':
                    ck(cycles(a)==val,(kind,n,val));C[n].add(val)
                else:
                    root=int(row['root']);rev=list(map(list,zip(*a)))
                    ck(start_paths(rev,root)==val,(kind,n,val));E[n].add(val)
                    if n>=2:
                        b=close_end(a,root)
                        ck(cycles(b)==val,'end closure')
                        ck(cycles(pad_last(b,2))==val,'degree-one padding')
                        closures+=1
                saved+=1
    ck(all(C[n]<=E[n] for n in range(3,10)),'cycle support lacks next-order tail')
    ck(set(range(594))<=E[9] and 594 not in E[9],'endpoint prefix boundary')
    out['checks']['saved_cycle_and_endpoint_witnesses']=saved
    out['checks']['indegree_one_closure_and_length3_padding']=closures
    out['checks']['cycle_support_in_same_order_endpoint_support_n_3_to_9']=True
    out['checks']['endpoint_n9_prefix_inclusive']=593
    inverse={str(k):next((n for n in range(3,10) if k in C[n]),10) for k in range(1,594)}
    out['inverse_counts']={str(n):sum(int(v)==n for v in inverse.values()) for n in range(3,11)}
    # Compare two independent author enumerators, rather than conflate cycles and presentations.
    for n in range(3,9):
        data=json.loads((SPECTRUM/f'fixedcycle{n}.json').read_text())
        ck(set(data['spectrum'])|{0}==C[n],'fixed-cycle mismatch')
    data=json.loads((SPECTRUM/'spectrum9.json').read_text())
    ck(set(data['spectrum'])==C[9],'n9 second counter mismatch')
    out['checks']['two_independent_author_enumerators_agree_n_3_to_9']=True
    out['checks']['n9_gaps']=sorted(set(range(max(C[9])+1))-C[9])
    out['input_sha256']={label:hashlib.sha256(p.read_bytes()).hexdigest() for label,p in [('cycle-levels/build_inputs.py',LEVELS/'build_inputs.py'),('cycle-levels/census.cpp',LEVELS/'census.cpp'),('cycle-spectrum/census.cpp',SPECTRUM/'census.cpp')]}
    out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'packet_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
