#!/usr/bin/env python3
"""Reproduce pairing-obstruction.svg from the independent n=5 census.

Run from any directory. Standard library only; no repository producer imports.
"""
from hashlib import sha256
from html import escape
from itertools import combinations, permutations
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'audit/small.normal.json'
NOTE = ROOT / 'PAIRING_WIRING.md'
OUT = ROOT / 'pairing-obstruction.svg'


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def ordinary_converse(code, n=5):
    pairs = list(combinations(range(n), 2))
    index = {p: k for k, p in enumerate(pairs)}
    converse = code ^ ((1 << len(pairs))-1)
    def arc(a, b):
        return bool((converse >> index[tuple(sorted((a, b)))]) & 1) == (a < b)
    return min(sum(int(arc(p[a], p[b])) << k for k, (a, b) in enumerate(pairs))
               for p in permutations(range(n)))


def main():
    row = next(x for x in json.loads(DATA.read_text())['rows'] if x['n'] == 5)
    nodes = {x['code']: x for x in row['node_metadata']}
    need(row['ordinary'] == 12 and row['self_converse'] == 8, 'wrong census universe')
    need(min(map(int, row['pairing_guard']['abstract_fixed_histogram'])) == 4,
         'abstract minimum changed')
    need(min(map(int, row['pairing_guard']['wired_fixed_histogram'])) == 8,
         'wiring minimum changed')
    need(len(row['weighted_vertex_automorphisms']) == 1, 'merged vertices not rigid')
    chiral = {x['H']: (x['code'], ordinary_converse(x['code']))
              for x in nodes.values() if not x['self_converse']}
    need(chiral == {3: (2, 11), 5: (4, 9)}, 'canonical converse pairs changed')
    self_nodes = {h: sorted(x['code'] for x in nodes.values()
                           if x['H'] == h and x['self_converse'])
                  for h in (1, 3, 5, 9, 11, 13, 15)}
    h9 = self_nodes[9]
    h15 = self_nodes[15]
    tt = self_nodes[1][0]
    need(h9 == [8, 10] and h15 == [40, 76] and tt == 0, 'canonical labels changed')
    degree9 = [nodes[c]['degree'] for c in h9]
    incidence9 = [nodes[c]['neighbors'].get(str(tt), 0) for c in h9]
    degree15 = [nodes[c]['degree'] for c in h15]
    need(degree9 == [9, 9] and incidence9 == [0, 1] and degree15 == [5, 3],
         'wiring obstruction changed')

    metadata = {'source_sha256': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                                 for p in (DATA, NOTE)},
                'ordinary_chiral_pairs': chiral,
                'class_codes': sorted([x['code'] for x in nodes.values()] +
                                      [pair[1] for pair in chiral.values()]),
                'scope': 'All 12 ordinary tournament classes at n=5; full merged dart wiring.'}
    need(len(set(metadata['class_codes'])) == 12, 'classes not distinct')
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1440" height="1040" viewBox="0 0 1440 1040" role="img" aria-labelledby="title description">',
           '<title id="title">Order-five pairing obstruction: four abstract fixed classes, eight with wiring</title>',
           '<desc id="description">Two panels show all twelve ordinary tournament classes grouped by Hamilton-path count. The abstract optimal involution swaps the two self-converse classes of count nine and the two of count fifteen, leaving four fixed. Preserving the merged dart diagram blocks both swaps and forces eight fixed classes. The two chiral converse pairs remain exchangeable.</desc>',
           '<metadata>' + escape(json.dumps(metadata, sort_keys=True)) + '</metadata>',
           '<defs><style>text{font-family:Arial,Helvetica,sans-serif;fill:#172b3a} .small{font-size:17px} .muted{fill:#526573} .rowlabel{font-size:19px;font-weight:700} .code{font-size:19px;font-weight:700;text-anchor:middle;dominant-baseline:central} .paneltitle{font-size:23px;font-weight:700}</style>',
           '<marker id="arrow" markerWidth="7" markerHeight="7" refX="5" refY="3.5" orient="auto-start-reverse"><path d="M0,0 L6,3.5 L0,7" fill="none" stroke="#167d82" stroke-width="1.5"/></marker>',
           '</defs><rect width="1440" height="1040" fill="#ffffff"/>']

    def text(x, y, value, size=20, weight=None, fill=None, anchor=None):
        attrs = f'x="{x}" y="{y}" font-size="{size}"'
        if weight: attrs += f' font-weight="{weight}"'
        if fill: attrs += f' fill="{fill}"'
        if anchor: attrs += f' text-anchor="{anchor}"'
        svg.append(f'<text {attrs}>{escape(str(value))}</text>')

    def line(x1, y1, x2, y2, color='#dbe2e7', width=1, extra=''):
        svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}" {extra}/>')

    def node(x, y, code, exchanged=False, is_chiral=False):
        fill, border = ('#e3f3f1', '#167d82') if exchanged else ('#f0f3f6', '#677b8b')
        if is_chiral:
            svg.append(f'<rect x="{x-24}" y="{y-22}" width="48" height="44" rx="7" fill="{fill}" stroke="{border}" stroke-width="2"/>')
        else:
            svg.append(f'<circle cx="{x}" cy="{y}" r="23" fill="{fill}" stroke="{border}" stroke-width="2"/>')
        svg.append(f'<text x="{x}" y="{y}" class="code">{code}</text>')

    def exchange(x1, x2, y):
        line(x1+31, y, x2-31, y, '#167d82', 2.4,
             'marker-start="url(#arrow)" marker-end="url(#arrow)"')

    def blocked(x1, x2, y):
        line(x1+28, y, x2-28, y, '#b55f30', 2, 'stroke-dasharray="4 4"')
        middle = (x1+x2)/2
        line(middle-6, y-7, middle+6, y+7, '#b55f30', 2.5)
        line(middle-6, y+7, middle+6, y-7, '#b55f30', 2.5)

    text(48, 55, 'Pairings at order five', 34, '700')
    text(48, 89, '12 ordinary classes • preserve H and commute with converse R', 21, fill='#526573')
    for x, wired in ((48, False), (736, True)):
        svg.append(f'<rect x="{x}" y="124" width="656" height="580" rx="13" fill="#ffffff" stroke="#cdd8df" stroke-width="1.5"/>')
        svg.append(f'<path d="M{x+13},124 H{x+643} Q{x+656},124 {x+656},137 V213 H{x} V137 Q{x},124 {x+13},124" fill="{"#f7f0e8" if wired else "#eef5f6"}"/>')
        text(x+26, 160, 'B  Preserve the dart wiring' if wired else 'A  Abstract pairing', 23, '700')
        text(x+26, 194, '8 fixed classes · 2 exchanged pairs' if wired else '4 fixed classes · 4 exchanged pairs', 21, '700', '#9b4d24' if wired else '#167d82')
        text(x+26, 237, 'H', 18, '700', '#526573')
        text(x+127, 237, 'ordinary class codes', 16, fill='#526573')
        if wired: text(x+337, 237, 'why this exchange fails', 16, fill='#526573')
        for index, h in enumerate((1, 3, 5, 9, 11, 13, 15)):
            y = 275 + 60*index
            if h in (9, 15) and wired:
                svg.append(f'<rect x="{x+13}" y="{y-28}" width="630" height="56" rx="6" fill="#fff5eb"/>')
            text(x+26, y+7, h, 20, '700')
            single = self_nodes[h]
            if len(single) == 1:
                node(x+145, y, single[0])
            elif len(single) == 2:
                node(x+145, y, single[0], exchanged=not wired)
                node(x+255, y, single[1], exchanged=not wired)
                (blocked if wired else exchange)(x+145, x+255, y)
            if h in chiral:
                a,b = chiral[h]
                node(x+255, y, a, exchanged=True, is_chiral=True)
                node(x+355, y, b, exchanged=True, is_chiral=True)
                exchange(x+255, x+355, y)
                if wired: text(x+425, y+6, 'still allowed', 17, fill='#167d82')
            if wired and h == 9:
                text(x+337, y-4, 'dart degree: 9 = 9', 18, '700', '#9b4d24')
                text(x+337, y+19, 'TT incidence: 0 ≠ 1', 18, fill='#9b4d24')
            if wired and h == 15:
                text(x+337, y+6, 'dart degree: 5 ≠ 3', 18, '700', '#9b4d24')
        text(x+26, 683, 'Four self-converse classes can be paired.' if not wired
             else 'All eight self-converse classes must stay fixed.', 17, fill='#526573')

    svg.append('<rect x="48" y="734" width="656" height="144" rx="10" fill="#f5f8fa"/>')
    text(72, 768, 'H = 9: degree alone misses the obstruction', 21, '700')
    text(72, 802, 'Codes 8 and 10 each have dart degree 9.', 19)
    text(72, 833, 'Only code 10 meets the transitive class TT₅ (code 0).', 19)
    text(72, 858, 'Exchanging them changes that incidence.', 18, fill='#526573')
    svg.append('<rect x="736" y="734" width="656" height="144" rx="10" fill="#f5f8fa"/>')
    text(760, 768, 'H = 15: the degrees already disagree', 21, '700')
    text(760, 802, 'Code 40 has dart degree 5; code 76 has degree 3.', 19)
    text(760, 833, 'Their exchange cannot preserve the merged diagram.', 19)
    text(760, 858, 'The two chiral pairs remain independently flippable.', 18, fill='#526573')

    node(72, 925, '·')
    text(108, 932, 'self-converse class', 17)
    node(343, 925, '·', is_chiral=True)
    text(379, 932, 'chiral class', 17)
    exchange(532, 652, 925)
    text(671, 932, 'exchanged by J', 17)
    blocked(913, 1023, 925)
    text(1044, 932, 'blocked by wiring', 17)
    text(48, 982, 'Fixed classes are gray; exchanged classes are teal. Wiring means the full merged dart multigraph, including multiplicity.', 17, fill='#526573')
    text(48, 1011, 'FINITE-EXACT  •  Codes are minimum orientation bitmasks over vertex relabelings; the diagram displays ordinary classes.', 16, fill='#526573')
    svg.append('</svg>')
    OUT.write_text('\n'.join(svg)+'\n')
    print(json.dumps({'svg':str(OUT), **metadata}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
