"""Construct the d=35 grid at exponent one or any exponent at least six."""
from pathlib import Path
import argparse,json,sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_coverage_compiler_20261111'))
import construct as old


def construct(parameter,exponent,data=None):
    if parameter<1 or exponent<1 or 2<=exponent<=5:
        raise ValueError('Require parameter >=1 and exponent 1 or >=6; exponents 2..5 remain open here.')
    if data is None:data=json.loads((HERE/'residual_grids.json').read_text())
    if exponent==1:
        q=next(q for q in data['head_queries'] if q['exponent']==1 and q['shift']==0)
        w=q['witness']
        if w['exterior_H']!=1:raise RuntimeError('Unexpected exterior factor')
        _,rows=old.ib.inflate(w['core_rows'],w['core_port'],parameter)
    else:
        item=next(x for x in data['binary_records'] if x['d']==35)
        labels,base=old.ib.inflate(item['rows'],item['port'],parameter+5)
        selected=[i for i,v in enumerate(labels) if item['selected']>>v&1]
        rows=old.tail(base,selected,exponent-5)
    return dict(d=35,parameter=parameter,exponent=exponent,c=35*2**parameter+1,
                expected_H=1+(35*2**parameter+1)*2**exponent,order=len(rows),rows=rows)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--parameter',type=int,required=True)
    parser.add_argument('--exponent',type=int,required=True);args=parser.parse_args()
    print(json.dumps(construct(args.parameter,args.exponent),indent=2))


if __name__=='__main__':main()
