"""Dependency-free scientific SVG; exact geometry is checked by verify.py."""
from pathlib import Path
import json,xml.etree.ElementTree as ET
HERE=Path(__file__).resolve().parent
pts=json.loads((HERE/'verification.json').read_text())['geometry']['upper_vertices_numeric']
parts=['<svg xmlns="http://www.w3.org/2000/svg" width="960" height="520" viewBox="0 0 960 520">',
       '<rect width="960" height="520" fill="white"/>',
       '<g font-family="Arial,sans-serif" fill="#263238" text-anchor="middle">']
def text(x,y,s,size=16):parts.append(f'<text x="{x}" y="{y}" font-size="{size}">{s}</text>')
def edge(a,b,shared=False):
    colour='#bb4430' if shared else '#245b85'
    parts.append(f'<line x1="{a[0]:.6f}" y1="{a[1]:.6f}" x2="{b[0]:.6f}" y2="{b[1]:.6f}" stroke="{colour}" stroke-width="3"/>')
def dot(p):parts.append(f'<circle cx="{p[0]:.6f}" cy="{p[1]:.6f}" r="3"/>')
text(480,31,'Shared boundaries leave the pseudoforest class',23)
text(250,78,'Capped total area',20)
text(700,78,'Product of areas — no unit-area cap',18)
grid=lambda x,y:(155+95*x,140+95*y)
for x in range(3):
    for y in range(2):edge(grid(x,y),grid(x,y+1),x==1)
for y in range(3):
    for x in range(2):edge(grid(x,y),grid(x+1,y),y==1)
for x in range(3):
    for y in range(3):dot(grid(x,y))
for x in [.5,1.5]:
    for y in [.5,1.5]:a,b=grid(x,y);text(a,b+6,'1',21)
text(250,395,'12 unit fences · total area 4',17)
text(250,422,'Cycle rank 4; pseudoforest maximum 3',15)
for sign in [-1,1]:
    points=[(700+90*x,255-90*sign*y) for x,y in pts]
    poly=' '.join(f'{x:.6f},{y:.6f}' for x,y in points)
    parts.append(f'<polygon points="{poly}" fill="#245b85" fill-opacity="0.08"/>')
    for a,b in zip(points,points[1:]):edge(a,b)
    for p in points:dot(p)
edge((610,255),(700,255),True);edge((700,255),(790,255),True);dot((700,255))
text(700,188,'A',22);text(700,332,'A',22)
text(700,455,'12 unit edges · A² = (89 + 26√13)/16',16)
text(700,482,'Cycle rank 2',16)
parts.append('</g></svg>')
svg='\n'.join(parts)+'\n';ET.fromstring(svg)
(HERE/'shared_boundaries.svg').write_text(svg,encoding='utf-8')
print('Wrote checked SVG with exact-source coordinates.')
