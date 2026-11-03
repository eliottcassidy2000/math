#!/usr/bin/env python3
"""Exact reflection/physical-port toothpick model; no external dependencies.

Run from any directory. Writes toothpick.out and toothpick.svg alongside itself.
Optional Pillow, when installed, also writes toothpick.png from the same drawing
primitives. Universe: the unique classical origin-seeded process through 128
integer stages, not an enumeration of arbitrary seed configurations.
"""
from collections import Counter
from functools import lru_cache
from pathlib import Path
import json
import math

HERE = Path(__file__).resolve().parent

def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)

def ends(stick):
    x, y, o = stick
    return ((x, y-1), (x, y+1)) if o == 0 else ((x-1, y), (x+1, y))

def reflect(stick, sx, sy):
    x, y, o = stick
    return (sx*x, sy*y, o)

def exposed(sticks):
    endpoint_inc = Counter(p for s in sticks for p in ends(s))
    centers = {s[:2] for s in sticks}
    return {p for p, weight in endpoint_inc.items()
            if weight == 1 and p not in centers}

def advance(sticks):
    ports = exposed(sticks)
    newborn = {}
    for s in sticks:
        for p in ends(s):
            if p in ports:
                child = (*p, 1-s[2])
                require(child not in newborn, 'duplicate exposed source port')
                require(child not in sticks, 'exposed port already occupied')
                newborn[child] = s
    return set(newborn), newborn

@lru_cache(None)
def births_rec(n):
    if n <= 1:
        return n
    power = 1 << (n.bit_length()-1)
    if n == power:
        return power
    i = n-power
    return 2*births_rec(i)+births_rec(i+1)

def model(limit):
    states = {0:set()}
    sticks = {(0,0,0)}
    born = {(0,0,0):1}
    parents = {}
    rows = [{'n':0, 'total':0, 'births':0, 'exposed':0}]
    for n in range(1,limit+1):
        if n > 1:
            newborn, p = advance(sticks)
            sticks |= newborn
            born.update({s:n for s in newborn})
            parents.update(p)
        else:
            newborn = set(sticks)
        for sx,sy in ((1,-1),(-1,1),(-1,-1)):
            require({reflect(s,sx,sy) for s in sticks} == sticks,
                    f'configuration reflection failed at {n}')
            for s in newborn:
                require(born[reflect(s,sx,sy)] == n,'birth clock lost')
                if n > 1:
                    require(parents[reflect(s,sx,sy)] == reflect(parents[s],sx,sy),
                            'physical source transport failed')
        axes = {s for s in sticks if s[0] == 0 or s[1] == 0}
        expected_axes = {(0,0,0)} if n == 1 else {(0,0,0),(0,1,1),(0,-1,1)}
        require(axes == expected_axes,'axis fixed locus changed')
        require(len(newborn) == births_rec(n),f'independent recurrence failed {n}')
        require(len(sticks) == sum(births_rec(j) for j in range(n+1)),
                'independent cumulative recurrence failed')
        if n&(n-1) == 0:
            require(len(sticks) == (2*n*n+1)//3,'dyadic total failed')
        if n >= 3:
            require(len(newborn)%4 == 0,'newborn four-orbits failed')
        rows.append({'n':n, 'total':len(sticks),'births':len(newborn),
                     'exposed':len(exposed(sticks))})
        if n in (2,3,4,8,16,32,64,128):
            states[n] = set(sticks)
    return rows, states

def drawing(states):
    # Drawing primitives are shared by SVG and optional PNG renderers.
    from html import escape
    width,height=1300,760
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="#f8fafc"/>']
    try:
        from PIL import Image, ImageDraw, ImageFont
        pic=Image.new('RGB',(width,height),'#f8fafc'); d=ImageDraw.Draw(pic)
        font_path='/System/Library/Fonts/Supplemental/Arial.ttf'
        def font(size):
            try:return ImageFont.truetype(font_path,size)
            except OSError:return ImageFont.load_default()
    except ImportError:
        pic=d=None
    def line(x1,y1,x2,y2,color='#334155',w=2):
        parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{w}" stroke-linecap="round"/>')
        if d: d.line((x1,y1,x2,y2),fill=color,width=w)
    def txt(x,y,s,size=20,color='#0f172a'):
        parts.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{size}" fill="{color}">{escape(s)}</text>')
        if d: d.text((x,y-size),s,fill=color,font=font(size))
    def circle(x,y,r,color):
        parts.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color}"/>')
        if d: d.ellipse((x-r,y-r,x+r,y+r),fill=color)
    txt(35,42,'Reflection gives a pairing only after its fixed points and physical ports are retained',25)
    txt(35,79,'Classical toothpicks: exact integer coordinates, synchronous birth clocks, unique exposed parent',19,'#475569')
    line(30,105,1270,105,'#cbd5e1')
    txt(35,145,'Stage 8: 43 sticks = 1 + 2 + 10 x 4',22)
    scale=32; cx,cy=260,330
    line(cx-160,cy,cx+160,cy,'#cbd5e1'); line(cx,cy-165,cx,cy+165,'#cbd5e1')
    for x,y,o in sorted(states[8]):
        color=('#dc2626' if (x,y)==(0,0) else '#b7791f' if x==0 else '#0f766e')
        p,q=ends((x,y,o));line(cx+p[0]*scale,cy-p[1]*scale,cx+q[0]*scale,cy-q[1]*scale,color,4)
    txt(35,530,'Red: the fixed seed. Gold: one orbit of size 2.',18)
    txt(35,558,'All later sticks have free four-reflection orbits.',18)
    txt(35,586,'Hence T(n) = 3 mod 4 for every n >= 2.',18)
    txt(35,614,'Half-turn alone fixes only the seed: T(n) is odd.',18)
    line(530,125,530,680,'#cbd5e1')
    txt(565,145,'Stage 3: why orbit counts lose the growth rule',22)
    cx,cy,scale=730,300,45
    for x,y,o in sorted(states[3]):
        p,q=ends((x,y,o));line(cx+p[0]*scale,cy-p[1]*scale,cx+q[0]*scale,cy-q[1]*scale,'#0f766e' if x else '#64748b',5)
    for x in (-1,1):circle(cx+x*scale,cy,7,'#dc2626')
    txt(875,245,'4 reflected sticks',18)
    txt(875,275,'2 physical collision ports',18)
    txt(875,305,'4 endpoint incidences / 2 ports',18)
    txt(875,335,'= 2 incidences at each port',18,'#dc2626')
    txt(565,435,'The ports are blocked. One adjacent orbit does not mean one source.',18)
    txt(565,485,'Tournament insertion has the analogous root-fiber loss:',20)
    txt(565,525,'one order-2 class  ->  four rooted order-3 classes',19)
    txt(565,560,'three roots of the transitive class  +  one root orbit of the cycle',18)
    txt(565,595,'Forgetting the root gives only two ordinary classes.',18)
    txt(565,630,'Both are self-converse; H = 1 and H = 3.',18,'#dc2626')
    txt(35,718,'Exact bridge: weighted orbit incidence. No toothpick-to-tournament bijection or shared dyadic recursion is claimed.',18,'#475569')
    parts.append('</svg>')
    (HERE/'toothpick.svg').write_text('\n'.join(parts)+'\n')
    if pic:pic.save(HERE/'toothpick.png')

def main():
    rows,states=model(128)
    # Hostile quotient: one stick orbit meets two ports, twice per port.
    orbit={(-1,-1,0),(-1,1,0),(1,-1,0),(1,1,0)}
    ports={(-1,0),(1,0)}
    incidence=sum(p in ports for s in orbit for p in ends(s))
    require(incidence==4 and len(ports)==2,'hostile incidence changed')
    require(ports.isdisjoint(exposed(states[3])),'collision ports incorrectly exposed')
    true_new,_=advance(states[3]); require(len(true_new)==4,'stage4 hostile failed')
    # Rooted insertion of a named third vertex into 0 -> 1.
    rooted=[]
    from itertools import permutations
    for mask in range(4):
        arcs={(0,1)}|{(2,v) if mask>>v&1 else (v,2) for v in range(2)}
        h=sum(all((p[i],p[i+1]) in arcs for i in range(2)) for p in permutations(range(3)))
        rooted.append(h)
    require(Counter(rooted)==Counter({1:3,3:1}),'rooted insertion fiber failed')
    # Explicitly demonstrate the old missing -1 recurrence.
    require(rows[3]['total']==7 and rows[2]['total']+2*rows[1]['total']+rows[2]['total']==8,
            'missing minus one witness changed')
    result={'scope':'FINITE-EXACT: classical origin seed, stages 0..128; all reflections/parents/clocks checked',
            'selected_stages':[r for r in rows if r['n']<=16 or r['n'] in (32,64,128)],
            'independent_control':'All newborn and total counts equal cited dyadic recurrence through 128.',
            'hostile_ports':{'stage':3,'port_orbit':sorted(ports),'stick_orbit_size':len(orbit),
                             'incidences':incidence,'per_port':incidence//len(ports),
                             'naive_neighbor_orbit_count':1,'true_newborn_at_stage4':len(true_new),
                             'false_added_by_support_only_port_rule':2},
            'rooted_order3_insertion_H_fibers':dict(sorted(Counter(rooted).items())),
            'normal_or_optimized':'Checks use explicit exceptions and remain active under python -O.'}
    out=json.dumps(result,indent=2,sort_keys=True)+'\n'
    (HERE/'toothpick.out').write_text(out);print(out,end='')
    drawing(states)

if __name__=='__main__':main()
