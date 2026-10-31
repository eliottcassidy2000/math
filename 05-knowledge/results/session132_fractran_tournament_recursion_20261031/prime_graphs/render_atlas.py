#!/usr/bin/env python3
"""Render the two exact 25-graph contact sheets as code-native SVG."""
from pathlib import Path
from math import sin,cos,pi,isqrt
from html import escape
import json,sys
ROOT=Path(__file__).resolve().parent

def main():
 output=Path(sys.argv[sys.argv.index('--output-dir')+1]) if '--output-dir' in sys.argv else ROOT
 output.mkdir(parents=True,exist_ok=True)
 data=json.loads((ROOT/'graphs.json').read_text());rows={r['n']:r for r in data['all_graphs_1_through_97']};ps=data['first_25_primes'];groups=data['groups_as_supplied']
 for name,orders in [('prime_order',ps),('rank_order',list(range(1,26)))]:
  title='Prime-order reading: p_i labels G_(p_i)' if name=='prime_order' else 'Rank reading: p_i labels G_i'
  pieces=['<svg xmlns="http://www.w3.org/2000/svg" width="1280" height="1410" viewBox="0 0 1280 1410">','<rect width="1280" height="1410" fill="#f7f8fb"/>',f'<text x="28" y="36" font-family="sans-serif" font-size="25" font-weight="bold" fill="#17243a">{title}</text>','<text x="28" y="63" font-family="sans-serif" font-size="14" fill="#39465d">Exact square-sum graphs • 71 is listed in C and D; 67 is missing • No silent correction</text>','<line x1="29" y1="84" x2="54" y2="84" stroke="#267b70" stroke-width="2"/><text x="62" y="89" font-family="sans-serif" font-size="13">Odd square-sum root</text>','<line x1="240" y1="84" x2="265" y2="84" stroke="#5266b3" stroke-width="2"/><text x="273" y="89" font-family="sans-serif" font-size="13">Even square-sum root</text>']
  for i,(p,n) in enumerate(zip(ps,orders)):
   x=20+(i%5)*252;y=110+(i//5)*252;cx=x+120;cy=y+139;r=rows[n];gs=','.join(g for g,vals in groups.items() if p in vals) or 'MISSING'
   pieces.append(f'<g data-graph-n="{n}" data-prime-label="{p}" data-rank="{i+1}">')
   pieces.append(f'<rect x="{x}" y="{y}" width="240" height="240" rx="10" fill="white" stroke="#d8deea"/>')
   pieces.append(f'<text x="{x+12}" y="{y+23}" font-family="sans-serif" font-size="15" font-weight="bold" fill="#17243a">i={i+1} · p={p} · G_{n}</text>')
   pieces.append(f'<text x="{x+12}" y="{y+42}" font-family="sans-serif" font-size="11" fill="#536178">Group {gs} · {r["edges"]} edges · {r["components"]} components</text>')
   pos={v:(cx+79*cos(-pi/2+2*pi*(v-1)/n),cy+79*sin(-pi/2+2*pi*(v-1)/n)) for v in range(1,n+1)}
   if n==1:pos[1]=(cx,cy)
   for u,v in r['edge_list']:
    a,b=pos[u];c,d=pos[v];k=isqrt(u+v);color='#267b70' if k%2 else '#5266b3'
    pieces.append(f'<line data-u="{u}" data-v="{v}" x1="{a:.4f}" y1="{b:.4f}" x2="{c:.4f}" y2="{d:.4f}" stroke="{color}" stroke-width="0.7" opacity="0.72"><title>{u}+{v}={k}²</title></line>')
   fs=8 if n<=32 else 6 if n<=61 else 4.8
   for v,(a,b) in pos.items():
    pieces.append(f'<circle data-vertex="{v}" cx="{a:.4f}" cy="{b:.4f}" r="{2.6 if n<=32 else 1.7}" fill="#19273c"/>')
    angle=-pi/2+2*pi*(v-1)/n;tx=cx+90*cos(angle);ty=cy+90*sin(angle)+fs/3
    if n==1:tx=cx;ty=cy-9
    pieces.append(f'<text x="{tx:.4f}" y="{ty:.4f}" text-anchor="middle" font-family="sans-serif" font-size="{fs}" fill="#19273c">{v}</text>')
   pieces.append('</g>')
  pieces.append('<text x="28" y="1390" font-family="sans-serif" font-size="12" fill="#536178">Vertices are labelled 1..n; every displayed edge is checked against the exact edge list. Vector labels can be zoomed.</text></svg>')
  (output/(name+'_first25.svg')).write_text('\n'.join(pieces)+'\n')
 print('Rendered two 25-graph SVG contact sheets.')
if __name__=='__main__':main()
