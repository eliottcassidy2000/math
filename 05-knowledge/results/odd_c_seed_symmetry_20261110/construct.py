"""Emit a literal tournament for any certified ray, without searching.

python3 -B construct.py --c 245 --exponent 8
Exponents below the stored ray threshold are not handled by this constructor.
"""
import argparse
import json
from pathlib import Path
from verify import closure, decode, induced

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def construct(c, exponent=None):
    if c<1 or c>4097 or not c%2:raise ValueError('c must be odd and between 1 and 4097')
    new={x['c']:x for x in json.loads((HERE/'verification.json').read_text())['new_seeds']}
    if c in new:
        x=new[c];B=x['rows'];O=x['O'];M=x['M'];s0=x['s0'];kind='tail'
    else:
        bank=ROOT/'05-knowledge/results/kps_s21_two_adic_lines_20261110/line_seeds.txt'
        fields=next(l.split() for l in bank.read_text().splitlines() if not l.startswith('#') and int(l.split()[0])==c)
        _,s0,kind,n,code,v,*_=fields;s0=int(s0);n=int(n);v=int(v)
        K=decode(n,int(code))
        if kind=='tail':
            keep=[i for i in range(n) if i!=v];B=induced(K,keep)
            O=sum(((K[v]>>i)&1)<<j for j,i in enumerate(keep));M=c*(1<<(s0-1))
    s=s0 if exponent is None else exponent
    if s<s0:raise ValueError(f'this certified ray starts at exponent {s0}; a below-ray witness requires the separate finite construction')
    if kind=='tail':
        t=s-s0+1;rows=closure(B,O,t)
    else:
        copies=s-s0+1
        owners=[(i,j) for i in range(n) for j in range(copies if i==v else 1)]
        rows=[]
        for i,a in owners:
            rows.append(sum(1<<z for z,(j,b) in enumerate(owners) if (i==j and a<b) or (i!=j and K[i]>>j&1)))
    result=dict(c=c,exponent=s,start_exponent=s0,mechanism=kind,expected_H=1+c*(1<<s),order=len(rows),rows=rows)
    if kind=='tail':result['seed']=dict(base_rows=B,selected_mask=O,offset=1,slope=M)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--c',type=int,required=True);p.add_argument('--exponent',type=int)
    a=p.parse_args();print(json.dumps(construct(a.c,a.exponent),indent=2))
