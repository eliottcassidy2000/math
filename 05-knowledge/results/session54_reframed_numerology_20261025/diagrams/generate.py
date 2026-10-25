#!/usr/bin/env python3
"""Deterministic, standard-library SVG figures from proved Session54 packets."""
from pathlib import Path
from html import escape
import json
P=Path(__file__).resolve().parent
class Figure:
 def __init__(self,title,desc,h):
  self.h=h
  self.s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="{h}" viewBox="0 0 1200 {h}" role="img" aria-labelledby="title desc">',f'<title id="title">{escape(title)}</title><desc id="desc">{escape(desc)}</desc>', '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8" fill="#30465c"/></marker></defs>',f'<rect width="1200" height="{h}" fill="#ffffff"/>']
 def text(self,x,y,s,size=21,weight='normal',color='#152b40',anchor='start'):
  self.s.append(f'<text x="{x}" y="{y}" font-family="DejaVu Sans, sans-serif" font-size="{size}" font-weight="{weight}" fill="{color}" text-anchor="{anchor}">{escape(str(s))}</text>')
 def box(self,x,y,w,h,title,lines,fill='#eef4f8'):
  self.s.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{fill}" stroke="#8b9ba8"/>')
  self.text(x+20,y+34,title,23,'bold')
  for i,line in enumerate(lines): self.text(x+20,y+67+i*29,line,20)
 def arrow(self,x,y,x2,y2,label=None):
  self.s.append(f'<path d="M{x},{y} L{x2},{y2}" fill="none" stroke="#30465c" stroke-width="2" marker-end="url(#arrow)"/>')
  if label: self.text((x+x2)/2,(y+y2)/2-12,label,17,anchor='middle')
 def save(self,name): (P/name).write_text('\n'.join(self.s+['</svg>'])+'\n')

def generate():
 f=Figure('Finite digit budgets: exact counts, different limits','Base-nine digit-budget state transfer and exact residue counts. Fixed finite budget gives completed 3-adic dimension one half, not finiteness of ordinary exponents.',1000)
 f.text(40,47,'Finite digit budgets: exact counts, different limits',30,'bold')
 f.text(40,82,'Lowest r base-9 digits of 2^k • r ≥ 1 • fixed integer budget b ≥ 0',21)
 f.box(40,110,1120,113,'Digit rule and the correct initial level',['Normal digits {0,1,2}; debt digits {3,4,5}; digits {6,7,8} forbidden.', 'At r = 1: unit digits {1,2} give debt 0; {4,5} give debt 1.'])
 f.text(40,262,'State j = number of debt digits already used',23,'bold')
 f.box(40,291,300,111,'State j',['0 ≤ j < b', 'Count v_r(j)'])
 f.box(650,270,510,85,'State j',['3 lifts: append 0, 1, or 2'])
 f.box(650,381,510,85,'State j + 1',['3 lifts: append 3, 4, or 5'])
 f.arrow(340,323,650,313,'no new debt')
 f.arrow(340,365,650,424,'spend one debt')
 f.text(40,445,'At j = b: only the 3 no-debt lifts survive.',20)
 f.text(40,477,'Debt overflow is discarded; it never wraps around.',20)
 f.box(40,505,1120,158,'Exact transfer and residue-class counts',[
  'v_(r+1) = 3(I + N)v_r;  N e_j = e_(j+1),  N e_b = 0.',
  'v_1 = (2) for b = 0;  v_1 = (2,2,0,…,0) for b ≥ 1.',
  'A_(r,b) = 2 · 3^(r−1) Σ_{j=0}^b binom(r,j);  M_r = 6 · 9^(r−1).'])
 f.box(40,687,540,121,'Finite-level check: r = 10',['b = 0: 39,366 classes', 'b = 1: 433,026 classes'])
 f.box(610,687,550,121,'Fixed b: polynomial × exponential',['One Jordan block: size b + 1; eigenvalue 3.', 'Extra debt changes r^b, not the base 3.'])
 f.box(40,831,1120,115,'Keep the two questions separate',[
  'Admissible set in (Z/2Z) × Z₃: dimension 1/2 and Haar measure zero.',
  'Ordinary k ≥ 0, full expansion of 2^k: finiteness or infinitude is NOT settled.'], '#fff4df')
 f.text(40,978,'Sources: sources/cross-base.md; sources/budget-transfer.md. b ≥ 2 is combinatorial, not a coefficient congruence.',16)
 f.save('digit-budget.svg')
 f=Figure('Native boundary quotient and a separately defined incidence lift','A commuting square identifies 63 native blockers modulo their boundary fibers with 21 projective points. A new typed trace incidence has 42 objects and 105 flags; its threefold vertex blow-up has 126 vertices, 945 edges, and degree 15.',1140)
 f.text(40,47,'Native boundary quotient → a separately defined incidence',29,'bold')
 f.text(40,82,'Solid arrows below are the proved quotient maps and bijections, not graph-cover claims.',20)
 f.box(40,115,420,100,'63 native Paley-7 blockers',['B: exact deleted-arc subsets'])
 f.box(740,115,420,100,'63 nonzero field vectors',['Ψ(B) = g^e ∈ F₆₄× (adjusted coding)'])
 f.arrow(460,161,740,161,'bijection Ψ')
 f.box(40,330,420,100,'21 native boundary images',['∂B: endpoint XOR of deleted arcs'])
 f.box(740,330,420,100,'21 points of PG(2,4)',['[g^e] = g^e F₄×'])
 f.arrow(250,215,250,330)
 f.text(270,271,'∂',22)
 f.arrow(950,215,950,330)
 f.text(970,271,'quotient',18)
 f.arrow(460,378,740,378,'ρ⁻¹ identification')
 f.text(40,476,'Commuting square: ∂B = ρ(e mod 21). Each fiber has exactly 3 blockers.',22,'bold')
 f.text(40,511,'Scalar fiber: {g^e, g^(e+21), g^(e+42)} represents ONE point, not a line.',21)
 f.text(40,546,'A projective line has 5 points. Native arc intersection is not the relation below.',21)
 f.s.append('<path d="M40,579 L1160,579" stroke="#8b9ba8" stroke-width="2"/>')
 f.text(40,618,'NEW relation: add a point type and a line type',25,'bold')
 f.box(40,645,540,158,'Typed quotient incidence',[
  '21 points e + 21 lines ℓ = 42 objects.',
  'Incident exactly when Tr(g^(e−ℓ)) = 0.',
  '105 flags; every object has degree 5.'])
 f.box(610,645,550,158,'Lift through actual boundary fibers',[
  '63 point-type + 63 line-type blockers = 126.',
  'Each flag lifts to a complete bipartite K₃,₃.',
  '105 × 9 = 945 edges; degree 5 × 3 = 15.'])
 f.box(40,828,1120,117,'The precise lift: a threefold vertex blow-up',[
  'Each quotient object has 3 blockers above it; every base flag becomes all 9 pairs.',
  'Not a locally bijective graph covering: one base neighbor lifts to 3 neighbors.'], '#fff4df')
 f.text(40,987,'Converse with a type switch: (e, point) ↔ (−e, line) is trace polarity.',22,'bold')
 f.text(40,1022,'Point-only inversion [x] → [x⁻¹] is NOT a collineation of this PG(2,4).',21)
 f.text(40,1062,'Tr(z) = z + z⁴ + z¹⁶. New trace incidence is not native arc-intersection incidence.',20)
 f.text(40,1113,'Sources: sources/boundary.md; sources/two-sheet.md. The figures introduce no additional relation or theorem.',16)
 f.save('boundary-incidence.svg')
 data={'r10':{'b0':39366,'b1':433026},'native':{'blockers':63,'boundary_images':21,'fiber_size':3},'typed':{'objects':42,'flags':105,'degree':5,'lift_vertices':126,'lift_edges':945,'lift_degree':15}}
 (P/'counts.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
if __name__=='__main__': generate()
