#!/usr/bin/env python3
"""Deterministic, dependency-free, exact SVG scientific figures."""
from pathlib import Path
from html import escape
import itertools,json
P=Path(__file__).resolve().parent
W=[0,9,1,16,15,14,29,21,30]
PAIRS=list(itertools.combinations(range(9),2))
def mask(a):return sum(1<<k for k,(i,j) in enumerate(PAIRS) if a[i]==a[j])
D={'schema':'session55-exact-diagrams-v1','crt':{'modulus':35,'entries':W,'mod5':[x%5 for x in W],'mod7':[x%7 for x in W],'masks':[hex(mask([x%n for x in W])) for n in [5,7]],'intersection':0,'ordinary_row_sums':[sum(W[i:i+3]) for i in [0,3,6]],'modular_line_sum':10,'square_roots':[next(r for r in range(35) if r*r%35==x) for x in W]},'field_ring_comparison':[{'q':27,'field_parker':True,'ring_parker':False},{'q':81,'field_parker':False,'ring_parker':False},{'q':243,'field_parker':True,'ring_parker':False}],'doubling':[{'n':50,'parker':False},{'n':100,'parker':True},{'n':200,'parker':False}],'thresholds':{'field_characteristic_2':'all positive exponents Parker','ring_2_power_nonparker_min_exponent':7,'ring_3_power_nonparker_min_exponent':3}}
BG='#f8fafc'; INK='#102a43'; MUT='#486581'; BLUE='#d9eaff'; GOLD='#ffedba'; MINT='#d5f1e5'; COLORS=[BLUE,GOLD,MINT]
class SVG:
 def __init__(self,w,h,title,desc):
  self.o=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-labelledby="title desc">',f'<title id="title">{escape(title)}</title><desc id="desc">{escape(desc)}</desc>',f'<rect width="{w}" height="{h}" fill="{BG}"/>']
 def text(self,x,y,s,size=20,color=INK,weight=400,anchor='start'):
  self.o.append(f'<text x="{x}" y="{y}" font-family="DejaVu Sans, sans-serif" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{color}">{escape(s)}</text>')
 def rect(self,x,y,w,h,color='#ffffff',stroke='#bcccdc',rx=8):self.o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{color}" stroke="{stroke}"/>')
 def line(self,x,y,x2,y2): self.o.append(f'<path d="M{x},{y} L{x2},{y2}" stroke="{MUT}" stroke-width="2" fill="none"/>')
 def save(self,name): (P/name).write_text('\n'.join(self.o+['</svg>'])+'\n')
def grid(s,x,y,vals,n=None):
 palette={v:COLORS[k] for k,v in enumerate(sorted(set(vals)))} if n else {}
 for i,v in enumerate(vals):
  xx=x+90*(i%3);yy=y+82*(i//3)
  s.rect(xx,yy,84,76,palette.get(v,'#ffffff'))
  s.text(xx+9,yy+18,f'cell {i}',11,MUT)
  s.text(xx+42,yy+52,str(v),29,weight=600,anchor='middle')
s=SVG(1180,790,'CRT resolves local collisions, not integer squarehood','The same nine cell labels in two local square-residue grids modulo 5 and 7 combine to nine distinct square residues modulo 35. Their ordinary row sums differ.')
s.text(36,49,'Shared cell labels are the missing CRT information',29,weight=700)
s.text(36,82,'Parker = no magic grid with nine distinct square entries. Zero is allowed.',17,MUT)
s.text(36,128,'LOCAL: modulo 5',21,weight=700);s.text(414,128,'LOCAL: modulo 7',21,weight=700);s.text(792,128,'PRODUCT: modulo 35',21,weight=700)
s.text(36,157,'Parker; 3 square values here',16,MUT);s.text(414,157,'Parker; 3 square values here',16,MUT);s.text(792,157,'Non-Parker; 9 distinct values',16,MUT)
grid(s,36,181,D['crt']['mod5'],5);grid(s,414,181,D['crt']['mod7'],7);grid(s,792,181,W)
s.text(350,307,'×',33,anchor='middle');s.text(728,307,'→',33,anchor='middle')
s.text(36,459,'A: {0,4,8}  {2,3,7}  {1,5,6}',17,weight=600)
s.text(414,459,'B: {2,4,6}  {0,5,7}  {1,3,8}',17,weight=600)
s.text(792,459,'Equal modulo 35 ⇔ equal locally',16,weight=600)
s.text(36,490,'Colour = equal residues within each local grid; cell numbers stay fixed.',17,MUT)
s.rect(36,515,1108,79,'#eaf2f8')
s.text(58,547,'C₃₅ = C₅ ∩ C₇ = ∅',26,weight=700)
s.text(465,545,'Each A-block meets each B-block in one cell.',19)
s.text(465,574,'No pair of different cells remains equal in both factors.',17,MUT)
s.rect(36,616,1108,128,'#fff4e8',stroke='#dca45d')
s.text(58,650,'Modular magic is not integer magic.',23,weight=700)
s.text(58,684,'Ordinary row sums: 10, 45, 80. They differ, but all are 10 modulo 35.',20)
s.text(58,716,'Every entry is a square residue modulo 35; that does not make it an integer square.',18,MUT)
s.text(36,773,'Source: CRT collision-wiring theorem, §§2, 4–5. Exact witness and collision masks are replay-checked.',13,MUT)
s.save('crt-shared-cells.svg')
s=SVG(1180,805,'Fields and quotient rings have different Parker behavior','Characteristic-two fields are always Parker, whereas the rings modulo 2 to the k are non-Parker exactly at k at least seven. Orders 27, 81, 243 distinguish field and ring behavior. Ring orders 50, 100, 200 alternate non-Parker, Parker, non-Parker.')
s.text(36,49,'Same cardinality does not mean the same arithmetic',28,weight=700)
s.text(36,82,'P = Parker (no distinct square-entry magic grid); non-P = such a grid exists.',17,MUT)
s.rect(36,110,1108,110,'#eaf2f8')
s.text(58,148,'Fields of order 2ʳ: P for every r ≥ 1',21,weight=700)
s.text(58,185,'Characteristic 2 forces opposite cells to coincide.',18,MUT)
s.text(655,148,'Rings ℤ/2ᵏℤ: non-P iff k ≥ 7',24,weight=700)
s.text(655,185,'First non-P ring in this family: 128.',18,MUT)
s.text(58,266,'Order q',20,weight=700);s.text(314,266,'Field of order q',20,weight=700);s.text(650,266,'Ring ℤ/qℤ',20,weight=700)
for i,row in enumerate(D['field_ring_comparison']):
 y=286+64*i;s.rect(36,y,1108,56,'#ffffff');s.text(58,y+36,str(row['q']),25,weight=600)
 for xx,key in [(314,'field_parker'),(650,'ring_parker')]:
  p=row[key];s.text(xx,y+36,'P' if p else 'non-P',24,'#854d0e' if p else '#16624b',600)
s.text(36,512,'F₈₁ is non-P, but F₂₄₃ is P: field embeddings follow degree divisibility, not size.',18,MUT)
s.line(36,537,1144,537)
s.text(36,575,'Smallest doubling re-entry into the Parker rings',23,weight=700)
for x,n,p in [(36,50,False),(444,100,True),(852,200,False)]:
 s.rect(x,601,292,91,'#fff0cd' if p else '#e3f3eb')
 s.text(x+146,637,f'ℤ/{n}ℤ',25,weight=700,anchor='middle')
 s.text(x+146,672,'P' if p else 'non-P',23,anchor='middle')
for x in [386,794]:s.text(x,650,'×2 →',22,anchor='middle')
s.text(36,730,'Distinctness can collapse under reduction; a witness need not lift through doubling.',19)
s.text(36,761,'These are exact field/ring statements, not integer-square existence claims.',18,MUT)
s.text(36,790,'Sources: finite-field classification; 2-adic/3-adic thresholds; dyadic re-entry certificate.',13,MUT)
s.save('fields-rings-reentry.svg')
(P/'countdata.json').write_text(json.dumps(D,indent=2,sort_keys=True)+'\n')
