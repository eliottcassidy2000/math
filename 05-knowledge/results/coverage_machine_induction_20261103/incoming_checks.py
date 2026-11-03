#!/usr/bin/env python3
"""The single declared colour hostile in incoming priority tilings."""
import json
from itertools import combinations


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def main():
    denominators = (143,85,19,17,69,29,23,62,37,31,287,43,41)
    status = tuple(int(41 % d == 0) for d in denominators)
    need(status == (0,)*12+(1,), 'declared arithmetic status')
    chords = tuple((i,j) for i,j in combinations(range(15),2) if j>i+1)
    def encode(e):
        special = {(i,i+2):e[i] for i in range(13)}
        return {pair:special.get(pair,e[0]) for pair in chords}
    t = encode(status)
    partner = encode(tuple(1-b for b in status))
    need(len(t)==91 and all(partner[e]==1-t[e] for e in chords),'literal native complement')
    defect = [pair for pair in chords if t[pair]!=t[(14-pair[1],14-pair[0])]]
    need(defect==[(0,2),(12,14)],'reflection defect and black line')
    print(json.dumps(dict(status='FINITE-EXACT local correction control',input=41,
                          guard_bits=status,free_chords=91,reflection_defect=defect,
                          line_colour='black',native_subcube_identity=True),indent=2,sort_keys=True))


if __name__=='__main__':
    main()
