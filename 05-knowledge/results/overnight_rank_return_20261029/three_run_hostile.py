"""Full twelve-vertex three-run controls, with independent exact verification."""
import hashlib
from itertools import permutations
import json
from pathlib import Path
from random import Random
import shutil
import subprocess
import tempfile

import numpy as np
from atom_precision import hostile, automorphisms
from complement_kernel import HERE, need, rank_mod
from defect_transfer import defects, query, subset_query
from two_run import subset_paths


def packed_rank3(matrix):
    """Row-leading elimination using two Python-integer bitplanes over F3."""
    basis={};mask=(1<<len(matrix[0]))-1
    for row in matrix:
        one=two=0
        for j,x in enumerate(row):
            if int(x)%3==1:one|=1<<j
            elif int(x)%3==2:two|=1<<j
        while one or two:
            bit=(one|two)&-(one|two);col=bit.bit_length()-1
            if col not in basis:
                basis[col]=(one,two) if one&bit else (two,one)
                break
            a,b=basis[col]
            # Add minus coefficient times the normalized pivot.
            if one&bit:a,b=b,a
            z=mask^(one|two);w=mask^(a|b)
            one,two=(z&a)|(one&w)|(two&b),(z&b)|(two&w)|(one&a)
    return len(basis),sorted(basis)


def mass_and_leaves(rows,table):
    n=len(rows);full=(1<<n)-1
    h=[sum(map(sum,value)) for value in table]
    pairs=[0]*(full+1)
    for mask in range(1,full+1):
        subset=(mask-1)&mask
        while subset:
            pairs[mask]+=h[subset]*h[mask^subset]
            subset=(subset-1)&mask
    mass=sum(h[mask]*pairs[full^mask] for mask in range(1,full))
    dp=[[[0]*3 for _ in range(n)] for _ in range(full+1)]
    for v in range(n):dp[1<<v][v][0]=1
    for mask in range(1,full+1):
        for u,counts in enumerate(dp[mask]):
            if not any(counts):continue
            for v in range(n):
                if mask>>v&1:continue
                bad=int(not(rows[u]>>v&1))
                for d in range(3-bad):dp[mask|1<<v][v][d+bad]+=counts[d]
    return mass,sum(sum(row) for row in dp[full])


def run():
    compiler=shutil.which('g++')
    need(compiler is not None,'g++ C++17 compiler required for exact full census')
    records=[];rng=Random(202610291230)
    with tempfile.TemporaryDirectory(prefix='math-three-run-') as tmp:
        folder=Path(tmp);exe=folder/'census.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'three_run_hostile.cpp'),'-o',str(exe)],check=True)
        ports=list(permutations(range(12),3))
        for flip in (False,True):
            output=folder/f'{int(flip)}.bin'
            result=subprocess.run([str(exe),str(int(flip)),str(output)],capture_output=True,text=True,check=True)
            lines=result.stdout.splitlines();fields=lines[0].split()
            need(list(map(int,lines[1].split()[1:]))==hostile(flip),'independent C++ literal source')
            raw=output.read_bytes()
            need(len(raw)==1320*1320*8,'full response dimensions')
            M=np.frombuffer(raw,dtype='<u8').reshape(1320,1320)
            rows=hostile(flip);table=subset_paths(rows)
            mass,leaves=mass_and_leaves(rows,table)
            need(mass==int(M.sum())==int(fields[1])==80614764,'independent full cover mass')
            need(leaves==int(fields[3])==20338857,'independent budgeted permutation DP')
            selected=[];D=defects(rows)
            nonzero=list(zip(*np.nonzero(M)))
            rng.shuffle(nonzero)
            largest=tuple(map(int,np.unravel_index(M.argmax(),M.shape)))
            pairs=[largest]+nonzero[:48]+[(rng.randrange(1320),rng.randrange(1320)) for _ in range(48)]
            for i,j in pairs:
                s,t=ports[i],ports[j]
                value,stats=query(12,D,s,t)
                other=subset_query(table,s,t)
                need(value==other==int(M[i,j]),'frontier, subset and permutation/cut entry agreement')
                selected.append(dict(starts=s,ends=t,value=value,stats=stats))
            rank,pivots=packed_rank3(M)
            need(rank==rank_mod(M,3)==(435 if flip else 430),'independent exact F3 eliminations')
            # Simultaneously permuting run slots must preserve every entry.
            ix={p:i for i,p in enumerate(ports)}
            group=automorphisms(rows)
            need(len(group)==(27 if flip else 81),'full hostile automorphism orders')
            for g in group:
                action=[ix[tuple(g[v] for v in p)] for p in ports]
                need(np.array_equal(M,M[np.ix_(action,action)]),'full vertex transport in a common frame')
            for perm in ((1,2,0),(1,0,2)):
                action=[ix[tuple(p[q] for q in perm)] for p in ports]
                need(np.array_equal(M,M[np.ix_(action,action)]),'full run-slot transport')
            records.append(dict(label='27-control' if flip else '81-control',source=rows,
                                dimension=1320,mass=mass,permutations=leaves,max_entry=int(M.max()),
                                max_ports=[ports[largest[0]],ports[largest[1]]],automorphism_order=len(group),
                                rank_mod3=rank,sha256=hashlib.sha256(raw).hexdigest(),
                                pivot_columns=pivots,queries=selected))
            print('full three-run',records[-1]['label'],'mass',mass,'rank3',rank,flush=True)
    return dict(objects=records)


if __name__=='__main__':
    import sys
    receipt=json.loads(json.dumps(run()))
    path=HERE/'three_run_hostile.json'
    if '--write' in sys.argv:path.write_text(json.dumps(receipt,separators=(',',':'))+'\n')
    else:need(receipt==json.loads(path.read_text()),'unchanged exact full-matrix receipt')
    print('PASS full 81/27 three-run matrices and ternary ranks',flush=True)
