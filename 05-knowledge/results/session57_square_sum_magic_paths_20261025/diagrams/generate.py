#!/usr/bin/env python3
"""Deterministic text-only SVG figures. No third-party Python dependencies."""
from pathlib import Path
from html import escape
import json, hashlib
D=Path(__file__).resolve().parent
INK='#14263d'; MUT='#53657a'; BLUE='#1469ad'; ORANGE='#b85d13'; PURPLE='#7c3aad'; BG='#f5f7fb'
class SVG:
 def __init__(self,w,h): self.s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',f'<rect width="{w}" height="{h}" fill="{BG}"/>','<g font-family="DejaVu Sans, sans-serif">']
 def rect(self,x,y,w,h,fill='white',stroke='none',r=16): self.s.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}"/>')
 def text(self,x,y,t,size=24,fill=INK,weight='normal',anchor='start'): self.s.append(f'<text x="{x}" y="{y}" font-size="{size}" fill="{fill}" font-weight="{weight}" text-anchor="{anchor}">{escape(str(t))}</text>')
 def line(self,x,y,X,Y,color=INK,width=3): self.s.append(f'<path d="M{x},{y} L{X},{Y}" stroke="{color}" stroke-width="{width}" fill="none"/>')
 def save(self,name): (D/name).write_text('\n'.join(self.s+['</g></svg>'])+'\n')
P=['A0','B0','A1','B1','A2','B2']; Q=['B1','A0','B2','A1','B0','A2']
def edges(path): return {tuple(sorted(e)) for e in zip(path,path[1:])}
p,q=edges(P),edges(Q)
def color(a,b): e=tuple(sorted((a,b))); return PURPLE if e in p&q else BLUE if e in p else ORANGE
vals=dict(zip(['A0','A1','A2','B0','B1','B2'],[48,0,11,46,1,15]))
s=SVG(1500,1160)
s.text(55,62,'Repair the representation, not the magic square',37,weight='bold')
s.text(55,104,'The same nine square values in F₅₉: vertex adjacency loses information; AP edge labels retain it.',23,fill=MUT)
s.rect(45,140,465,560); s.rect(555,140,900,560)
s.text(70,185,'1  ENTRY VALUES AS VERTICES',23,weight='bold')
s.text(70,223,'Original magic grid • every line sums to 3',19,fill=MUT)
grid=[[15,35,12],[57,1,4],[49,26,46]]
for i,row in enumerate(grid):
 for j,v in enumerate(row):
  s.rect(95+j*120,250+i*87,108,76,'#fff1ea' if v==57 else '#edf2f8',r=8);s.text(149+j*120,300+i*87,v,34,ORANGE if v==57 else INK,anchor='middle')
s.text(70,556,'57 is isolated in the entry-value graph.',19,weight='bold')
s.text(70,592,'For every other entry v, 57 + v is a',20,fill=MUT)
s.text(70,622,'nonsquare modulo 59. No Hamilton path.',20,fill=MUT)
s.text(70,665,'The grid is magic. Its adjacency graph fails.',18,fill=MUT)
s.text(580,185,'2  THE SAME VALUES AS EDGE LABELS',23,weight='bold')
s.text(580,223,'Retain AP coordinates; move values to edge labels →',23,fill=BLUE,weight='bold')
s.text(580,265,'A = (48, 0, 11), step 11     B = (46, 1, 15), step 14',22)
s.text(580,307,'Cross-edge label Lᵢⱼ = Aᵢ + Bⱼ  (mod 59)',22,fill=MUT)
for j,b in enumerate([46,1,15]): s.text(935+j*180,357,f'B{j} = {b}',22,anchor='middle',weight='bold')
for i,a in enumerate([48,0,11]):
 s.text(610,415+i*85,f'A{i} = {a}',24,weight='bold')
 for j,b in enumerate([46,1,15]):
  c=color(f'A{i}',f'B{j}');s.rect(860+j*180,378+i*85,150,68,'#f6f3fa',r=7);s.text(935+j*180,423+i*85,(a+b)%59,34,c,weight='bold',anchor='middle')
s.text(580,665,'Cross-skeleton K₃,₃ only: ambient within-part edges may exist.',21,fill=MUT)
s.text(55,745,'Two Hamilton paths cover all nine labeled edges',28,weight='bold')
s.text(845,745,'P only',22,BLUE); s.text(985,745,'Q only',22,ORANGE); s.text(1125,745,'shared: label 46',22,PURPLE)
for path,y,title in [(P,832,'P'),(Q,960,'Q')]:
 s.text(55,y+8,title,29,weight='bold')
 for k,(a,b) in enumerate(zip(path,path[1:])):
  x=155+k*240; X=x+240; c=color(a,b);s.line(x,y,X,y,c,6);s.rect((x+X)/2-31,y-24,62,46,'white',r=7);s.text((x+X)/2,y+9,(vals[a]+vals[b])%59,25,c,weight='bold',anchor='middle')
 for k,v in enumerate(path):
  x=155+k*240;s.rect(x-51,y-25,102,50,INK,r=25);s.text(x,y+8,v,23,'white',weight='bold',anchor='middle');s.text(x,y+55,f'value {vals[v]}',17,MUT,anchor='middle')
s.text(55,1060,'Certificate structure: 4 additive cycle constraints + 2 AP guards = 6 independent linear constraints.',23,weight='bold')
s.text(55,1103,'All arithmetic above is modulo 59. This finite-field witness is not an integer solution or a prize claim.',22,fill=MUT)
s.save('edge-label-repair.svg')
s=SVG(1400,940)
s.text(55,65,'Square line totals are not the same as equal line totals',34,weight='bold')
s.text(55,110,'A literal integer example: nine distinct positive square entries, eight square line totals.',23,fill=MUT)
s.rect(45,150,760,590);s.rect(840,150,515,590)
s.text(75,195,'Each displayed cell is a square',25,weight='bold')
roots=[[2655,2976,320],[2624,2520,1665],[1440,895,3624]]
for i,row in enumerate(roots):
 for j,v in enumerate(row):
  s.rect(75+j*190,230+i*120,175,100,'#edf2f8',r=9);s.text(162+j*190,291+i*120,str(v)+'²',34,anchor='middle')
 s.text(675,290+i*120,'4001²',26,BLUE,weight='bold')
for j in range(3): s.text(162+j*190,640,'4001²',26,BLUE,weight='bold',anchor='middle')
s.text(75,700,'All 3 rows and all 3 columns: 16,008,001 = 4001²',24,BLUE)
s.text(870,202,'The diagonals disagree',27,weight='bold')
s.text(870,267,'Main diagonal  ↘',24,PURPLE,weight='bold')
s.text(870,307,'2655² + 2520² + 3624²',23)
s.text(870,350,'= 26,532,801 = 5151²',25,PURPLE)
s.text(870,423,'Other diagonal  ↙',24,ORANGE,weight='bold')
s.text(870,463,'320² + 2520² + 1440²',23)
s.text(870,506,'= 8,526,400 = 2920²',25,ORANGE)
s.rect(870,556,455,132,'#fff1ea',r=12)
s.text(895,605,'NOT A MAGIC SQUARE',28,ORANGE,weight='bold')
s.text(895,651,'4001² ≠ 5151² ≠ 2920²',26,ORANGE)
s.text(55,793,'A separate obstruction when equality AND square totals are required',27,weight='bold')
s.text(55,836,'Every 3 × 3 magic array has total T = 3 × center. If its nonzero center is e², then T = 3e²',23)
s.text(55,875,'cannot be an integer square. This does not settle the usual problem with unrestricted T.',23)
s.text(55,917,'Exact integer arithmetic • no distinct-square magic-square solution or prize claim',19,fill=MUT)
s.save('square-totals-not-magic.svg')
print('Generated two deterministic SVG figures.')
