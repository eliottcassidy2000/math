"""Literal adjacency constructors for certified unbounded seed families."""
from pathlib import Path
import argparse
import json
import inflation_bank as ib
from ray_probe import oddpart

HERE = Path(__file__).resolve().parent
# (ordinary order, canonical code, selected mask, inflated vertex)
FAMILIES = {
    'seven': (5,258,3,3), 'three': (6,321,12,4),
    'twentyone': (6,12866,17,5), 'binary29': (7,524673,5,5),
    'ternary30': (6,385,12,0), 'ternary50': (7,1281,12,0),
    'ternary78': (7,1537,20,0), 'ternary90': (7,1537,24,0),
    'mixed': (5,258,9,1), 'mixed11': (5,258,24,1),
}


def tail(rows, selected, length):
    if length < 1:
        raise ValueError('The proved tail law requires length >= 1.')
    n = len(rows)
    out = list(rows) + [0]*(length+1)
    root = n+length
    for old in range(n):
        for new in range(n,n+length):
            out[old] |= 1 << new
        out[root if old in selected else old] |= 1 << (old if old in selected else root)
    for new in range(n,n+length):
        for later in range(new+1,root+1):
            out[new] |= 1 << later
    return out


def seed(family, size):
    if size < 1:
        raise ValueError('The inflated block must be nonempty.')
    key = FAMILIES[family]
    bank = json.loads((HERE / 'inflation_bank_7.json').read_text())['records']
    item = next(r for r in bank if (r['n'],r['code'],r['selected'],r['inflated_vertex']) == key)
    labels, rows = ib.inflate(item['rows'],item['inflated_vertex'],size)
    selected = [i for i,v in enumerate(labels) if item['selected'] >> v & 1]
    slope = sum(u*ib.onto(size,j) for j,u in enumerate(item['slope_cover_coefficients'],1))
    c, j = oddpart(slope)
    return rows, selected, dict(family=family,block_size=size,offset=1,U=slope,c=c,first_line_exponent=j+1)


def construct(family, size=1, length=1):
    if family == 'ray99':
        rows, selected = [40,17,3,6,13,30], [0,3]
        metadata = dict(family=family,offset=23,U=99)
    else:
        rows, selected, metadata = seed(family,size)
    out = tail(rows,selected,length)
    return dict(**metadata,tail_length=length,base_rows=rows,selected_vertices=selected,
                rows=out,order=len(out),expected_H=metadata['offset']+(1 << length)*metadata['U'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--family', required=True, choices=['ray99']+list(FAMILIES))
    parser.add_argument('--size',type=int,default=1)
    parser.add_argument('--tail',type=int,default=1)
    args = parser.parse_args()
    print(json.dumps(construct(args.family,args.size,args.tail),indent=2))


if __name__ == '__main__':
    main()
