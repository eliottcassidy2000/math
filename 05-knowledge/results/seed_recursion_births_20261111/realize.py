"""Literal tournaments for the completed three/twentyone coefficient lines."""
from pathlib import Path
import argparse
import json
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_coverage_compiler_20261111'))
import construct as inherited


def head(name,parameter):
    data=json.loads((HERE/'head_completion.json').read_text())['completed_heads'][name]
    if parameter<data['start']:
        return next(x['rows'] for x in data['prefix'] if x['parameter']==parameter),'finite prefix'
    w=data['witness']
    _,rows=inherited.ib.inflate(w['rows'],w['port'],parameter-data['start']+1)
    return rows,'certified Hamilton kernel'


def realize(family,parameter,exponent):
    if parameter<1 or exponent<1:
        raise ValueError('Both parameters must be positive integers.')
    if family in ('three','twentyone'):
        multiplier,valuation=(3,4) if family=='three' else (21,2)
        c=multiplier*2**parameter+1
        size=parameter+1
    elif family=='mixed11':
        valuation=2
        c=2*3**parameter+3*2**parameter-1
        size=parameter
        if exponent==2:
            raise ValueError('The mixed11 exponent-two head is not uniformly closed by this packet.')
    else:
        raise ValueError('Unknown family.')
    if exponent>valuation:
        rows=inherited.construct(family,size,exponent-valuation)['rows']
        mechanism='inherited unit seed and tail'
    else:
        rows,mechanism=head(f'{family}_s{exponent}',parameter)
    return dict(family=family,parameter=parameter,c=c,exponent=exponent,
                expected_H=1+c*2**exponent,order=len(rows),rows=rows,mechanism=mechanism)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--family',required=True,choices=['three','twentyone','mixed11'])
    parser.add_argument('--parameter',required=True,type=int)
    parser.add_argument('--exponent',required=True,type=int)
    args=parser.parse_args()
    print(json.dumps(realize(args.family,args.parameter,args.exponent),indent=2))


if __name__=='__main__':main()
