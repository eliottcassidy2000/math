"""Deterministic, standard-library SVG source. No raster assets."""
from pathlib import Path
from html import escape
OUT=Path(__file__).resolve().parent
INK='#172b40'; MUTED='#4e6175'; GRAY='#e2e8ee'; BLUE='#d7eafe'; ORANGE='#ffe5c8'
class SVG:
 def __init__(self,w,h,title,desc):
  self.parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-labelledby="title desc">',f'<title id="title">{escape(title)}</title><desc id="desc">{escape(desc)}</desc>', '<rect width="100%" height="100%" fill="#ffffff"/>']; self.w=w; self.h=h
 def rect(self,x,y,w,h,fill='#f7f9fb',stroke='none',r=12):self.parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}"/>')
 def text(self,x,y,s,size=20,fill=INK,weight=400,anchor='start'):
  self.parts.append(f'<text x="{x}" y="{y}" font-family="DejaVu Sans, sans-serif" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{escape(s)}</text>')
 def line(self,x,y,xx,yy,color='#a5b2bf',width=2):self.parts.append(f'<path d="M{x} {y} L{xx} {yy}" fill="none" stroke="{color}" stroke-width="{width}"/>')
 def save(self,name): (OUT/name).write_text('\n'.join(self.parts+['</svg>'])+'\n')
def swatch(s,x,y,color,label):
 s.rect(x,y-16,20,20,color,r=4);s.text(x+30,y,label,17,MUTED)
def cell(s,x,y,name,value,color,w=150,h=90):
 s.rect(x,y,w,h,color,'#ffffff',7);s.text(x+12,y+23,name,15,MUTED,600);s.text(x+w/2,y+61,value,23,INK,500,'middle')
def main():
 s=SVG(1440,900,'Two hourglasses, one labelled carrier','The horizontal and vertical seven-cell hourglasses share four corners and center. In the additive K3,3 edge-value grid these five shared cells form a spanning tree; each hourglass contributes two chords.')
 s.text(50,58,'Two hourglasses, one labelled carrier',34,weight=700)
 s.text(50,94,'Same nine values, different placements. The arithmetic travels with the labels.',21,MUTED)
 swatch(s,50,133,GRAY,'Shared: 5 cells');swatch(s,330,133,BLUE,'Horizontal only: 2');swatch(s,660,133,ORANGE,'Vertical only: 2')
 s.rect(36,163,650,482);s.rect(710,163,694,482)
 s.text(60,202,'MAGIC-SQUARE PLACEMENT',19,weight=700)
 s.text(60,235,'Horizontal = top + center + bottom; vertical = rotation.',17,MUTED)
 vals=[('A','E + x',GRAY),('B','E − x − z',BLUE),('C','E + z',GRAY),('D','E − x + z',ORANGE),('E','E',GRAY),('F','E + x − z',ORANGE),('G','E − z',GRAY),('H','E + x + z',BLUE),('I','E − x',GRAY)]
 for k,(name,v,col) in enumerate(vals): cell(s,100+(k%3)*174,260+(k//3)*106,name,v,col,164,96)
 s.text(60,607,'7 + 7 − 5 = 9 cells; 5 + 5 − 2 = 8 magic lines.',19,weight=500)
 s.text(734,202,'AP–K₃,₃ EDGE-VALUE GRID',19,weight=700)
 s.text(734,235,'Each table entry labels an edge, not a graph vertex.',17,MUTED)
 for j,b in enumerate(['B₋₁ = E − z','B₀ = E','B₁ = E + z']):s.text(910+j*158,284,b,17,weight=600,anchor='middle')
 grid=[('B','E − x − z',BLUE),('I','E − x',GRAY),('D','E − x + z',ORANGE),('G','E − z',GRAY),('E','E',GRAY),('C','E + z',GRAY),('F','E + x − z',ORANGE),('A','E + x',GRAY),('H','E + x + z',BLUE)]
 for i,a in enumerate(['A₋₁ = −x','A₀ = 0','A₁ = x']):s.text(823,353+i*96,a,17,weight=600,anchor='end')
 for k,(name,v,col) in enumerate(grid):cell(s,835+(k%3)*158,300+(k//3)*96,name,v,col,150,88)
 s.text(734,607,'Vertex potentials: Aᵢ = ix, Bⱼ = E + jz. Edge: Lᵢⱼ = Aᵢ + Bⱼ.',17)
 s.text(50,687,'Five shared edges form a spanning tree',25,weight=700)
 s.text(50,725,'Gray: i = 0 or j = 0. Blue: (−1,−1), (1,1). Orange: (−1,1), (1,−1).',20)
 s.text(50,760,'The two central AP equations reduce the five tree labels to E, x, z.',20)
 s.text(50,795,'Each colored pair adds two chords; together they close all four fundamental cycles.',20)
 s.line(50,822,1390,822)
 s.text(50,854,'Square gate: E and E ± x, E ± z, E ± (x + z), E ± (x − z) must be squares.',19,weight=600)
 s.text(50,881,'Structural equivalence only: no integer-square existence claim. Indexed potential vertices stay distinct even if numbers coincide.',16,MUTED)
 s.save('hourglass-carrier.svg')
 s=SVG(1440,780,'Distinct profiles need not be linearly complete','At n=6, 34 distinct tournament converse-class profiles have linear rank 32. Adding two squared coordinates from observations conditional on the same latent class gives rank 34. This is not an individual-profile collision or an all-n quadratic claim.')
 s.text(50,58,'Distinct profiles ≠ linear completeness',34,weight=700)
 s.text(50,95,'Order-averaged graph-pair observer • exact finite result at n = 6',22,MUTED)
 cards=[(50,'34','converse classes','34 distinct individual profiles',GRAY),(510,'32','first-moment rank','Two mixture directions are lost',BLUE),(970,'34','augmented rank','Two squared coordinates repair it',ORANGE)]
 for x,num,head,sub,col in cards:
  s.rect(x,135,420,211,col);s.text(x+25,214,num,64,weight=700);s.text(x+25,259,head,24,weight=600);s.text(x+25,310,sub,18)
 s.text(480,247,'→',26,anchor='middle');s.text(940,247,'→',26,anchor='middle')
 s.text(50,390,'Keep 32 independent linear coordinates; append Aₜ(3)² and Aₜ(15)².',24,weight=600)
 s.text(50,427,'q = 3: K₁,₂ + 3 isolated vertices.   q = 15: K₁,₄ + 1 isolated vertex.   Each paired with its complement.',18,MUTED)
 s.rect(50,458,1340,166,'#f7f9fb')
 s.text(75,496,'SAME LATENT CLASS, TWO CONDITIONALLY INDEPENDENT ORDERS',18,weight=700)
 s.text(75,538,'Draw T once from μ. Then observe twice:  P(q, q) = Σₜ μ(T) [Aₜ(q) / 6!]².',23)
 s.text(75,581,'Fresh independent draws of T give the square of a marginal, not this second moment.',21,MUTED)
 s.text(50,665,'Exact certificate: a 34 × 34 integer minor has determinant',20)
 s.text(50,701,'211991592542375584355844096000000 ≠ 0.',25,weight=600)
 s.text(50,750,'No individual-profile collision at n = 6. No claim that two squared coordinates suffice for all n.',20,weight=600)
 s.save('observer-profile.svg')
if __name__=='__main__':main()
