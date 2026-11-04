#!/usr/bin/env python3
"""Deterministically materialize fixed physical masks and selected witnesses."""
import csv,json,pathlib
P=pathlib.Path(__file__).resolve().parent

def expand(m):
    return sum(7<<(3*i) for i in range(9) if m>>i&1)
def inflate(q):
    r=[0]*27
    for i in range(9):
        for v in range(3):
            r[3*i+v]=(1<<(3*i+(v+1)%3))|expand(q[i])
    return r
def augment(r,a,c):
    n=len(r);f=(1<<n)-1
    return [x|((a>>i&1)<<n)|((c>>i&1)<<(n+1)) for i,x in enumerate(r)]+[(f^a)|(1<<(n+1)),f^c]
def write(name,value):
    (P/name).write_text(json.dumps(value,separators=(',',':'))+'\n')
def main():
    battery=json.loads((P/'battery.json').read_text())
    candidates=json.loads((P/'candidates.json').read_text())['candidates']
    bank=list(csv.DictReader((P/'residues.tsv').open(),delimiter='\t'))
    if len(candidates)!=len(bank): raise ValueError('bank size')
    witnesses=[]
    for q,values in zip(candidates,bank):
        if q['node']!=int(values['node']): raise ValueError('bank ordering')
        k=next(i for i in range(12) if values[f'S{i}']!=values[f'T{i}'])
        b=battery[k];a,c=b['physical_a'],b['physical_c']
        rows=q['rows'];op=[511^x^(1<<i) for i,x in enumerate(rows)]
        witnesses.append(dict(node=q['node'],code=q['code'],query_index=k,
            coarse_a=b['a'],coarse_c=b['c'],physical_a=a,physical_c=c,
            source_residue=int(values[f'S{k}']),target_residue=int(values[f'T{k}']),
            augmented_source_rows=augment(inflate(rows),a,c),
            augmented_target_rows=augment(inflate(op),a,c)))
    write('physical_witnesses.json',witnesses)
    anchors=[]
    with (P/'filter.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['p_mask']=='-' or row['p_mask']=='0':continue
            p=int(row['p_mask']);i=(p&-p).bit_length()-1
            anchors.append(dict(node=int(row['node']),p_mask=p,anchor=i,coarse_a=1<<i,
                coarse_c=0,physical_a=expand(1<<i),physical_c=0,contrast_mod4=2))
    with (P/'endpoint_witnesses.tsv').open('w') as f:
        cols=['node','p_mask','anchor','coarse_a','coarse_c','physical_a','physical_c','contrast_mod4']
        f.write('\t'.join(cols)+'\n')
        for x in anchors:f.write('\t'.join(str(x[c]) for c in cols)+'\n')
    print(json.dumps({'physical_residual_witnesses':len(witnesses),'endpoint_witnesses':len(anchors)}))
if __name__=='__main__':main()
