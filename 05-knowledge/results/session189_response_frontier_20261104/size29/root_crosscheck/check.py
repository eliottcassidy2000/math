from pathlib import Path
import csv,json,subprocess,math
from itertools import permutations
P=Path(__file__).resolve().parent
B=P.parent
# Re-enumerate every small physical child profile used in this check.
physical=json.loads((P/'physical_profiles.json').read_text())
for label,item in physical.items():
    graph=item['rows']; n=len(graph); q=[0]*n
    for perm in permutations(range(n)):
        for cut in range(1<<(n-1)):
            if all((cut>>i&1) or (graph[perm[i]]>>perm[i+1]&1) for i in range(n-1)):
                q[cut.bit_count()]+=1
    if q!=item['q']:raise ValueError(('Physical profile mismatch',label,q))
lines=(B/'root_bank_input.txt').read_text().splitlines()
Q={int(v[0]):list(map(int,v[3:])) for line in lines[1:] if (v:=line.split())}
rows=list(csv.DictReader((B/'root_residues.tsv').open(),delimiter='\t'))
extra=list(csv.DictReader((B/'root_selected_residues.tsv').open(),delimiter='\t'))
pairs=[(i,j) for i in range(9) for j in range(i+1,9)]
samples=[]
for n,m in enumerate(json.loads((B/'root_initial_misses.json').read_text())):
    conf=m['configuration']
    row=next(r for r in extra if int(r['node'])==m['node'] and (int(r['a']),int(r['c'])) in [(0,2),(1,17)] and r[f'S{conf}']!=r[f'T{conf}'])
    samples.append((row,conf,'extra',n%3))
for k,conf in enumerate(list(range(1,10))+[10,13,17,20,24,28,32,36,40,43,45]):
    offset=(k*1600)%len(rows)
    row=next(r for r in rows[offset:]+rows[:offset] if r[f'S{conf}']!=r[f'T{conf}'])
    samples.append((row,conf,'initial',k%3))
# A five-child with normalized C3 profile must use baseline, not a five column.
for k in [0,17000]:
    row=next(r for r in rows[k:]+rows[:k] if r['S0']!=r['T0'])
    samples.append((row,0,'inherited_C3_five',k%9))
inputs=[]; manifest=[]
for idx,(row,conf,kind,variant) in enumerate(samples):
    profiles=[[3,6,6] for _ in range(9)]+[[1],[1]]
    if conf==0:profiles[variant]=[9,66,186,240,120]
    elif conf<10:profiles[conf-1]=[3,42,162,240,120]
    else:
        i,j=pairs[conf-10]
        profiles[i]=[3,18,36,24] if variant!=2 else [5,22,36,24]
        profiles[j]=[3,18,36,24] if variant==0 else [5,22,36,24]
    if sum(map(len,profiles))!=31:raise ValueError('Not 29 old names plus ports')
    factor=math.prod(p[0] for p in profiles)%8
    expected=[factor*int(row[f'{s}{conf}'])%8 for s in ['S','T']]
    inp=[idx,int(row['a']),int(row['c'])]+Q[int(row['node'])]
    for p in profiles:inp += [len(p)]+p
    inputs.append(' '.join(map(str,inp)))
    manifest.append(dict(sample=idx,node=int(row['node']),a=int(row['a']),c=int(row['c']),configuration=conf,kind=kind,profiles=profiles,physical_factor=factor,expected=expected))
(P/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(P/'input.txt').write_text(str(len(inputs))+'\n'+'\n'.join(inputs)+'\n')
subprocess.run(['g++','-O3','-std=c++17',str(P/'full_occupancy.cpp'),'-o',str(P/'full_occupancy')],check=True)
with (P/'input.txt').open() as i,(P/'output.txt').open('w') as o:
    subprocess.run([str(P/'full_occupancy')],stdin=i,stdout=o,check=True)
output=(P/'output.txt').read_text().splitlines()
if len(output)!=len(manifest):raise ValueError('Incomplete result')
for row,m in zip(output,manifest):
    idx,a,c,s,t=map(int,row.split())
    if (idx,a,c)!=(m['sample'],m['a'],m['c']) or [s,t]!=m['expected']:raise ValueError((m,row))
    if s==t:raise ValueError('Not a separating sample')
print(f'PASS: {len(manifest)} complete 31-name source/target physical residue pairs, including all 11 initial-bank missed configurations at selected extra queries, 20 initial witnesses, and 2 inherited-C3 five-child cases. No visit truncation or normalized transition factors.')
(P/'RESULT.txt').write_text(f'PASS: {len(manifest)} pairs ({2*len(manifest)} full-occupancy counts). All 11 initial misses, 20 initial witnesses, 2 C3-profile five cases. Bounded cross-check only.\n')
