#!/usr/bin/env python3
"""Deterministic SVGs; independent finite recomputation; Python standard library only."""
from pathlib import Path
from itertools import product,permutations,combinations
import sys,json,hashlib,math,html
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P/'sources'))
import session55_verify as old

def need(v,m):
 if not v:raise RuntimeError(m)
def facts():
 states=[]
 for q in (2,5,3):
  Q=old.roots(q)
  states.append(sorted({old.magic(q,c,a,b) for c,a,b in product(Q,repeat=3) if all(v in Q for v in old.magic(q,c,a,b))}))
 masks=[[old.mask(a) for a in S] for S in states]
 left=list(product(range(len(states[0])),range(len(states[1]))));right=list(product(range(len(states[1])),range(len(states[2]))))
 bad=sum((masks[0][a]&masks[1][b]&masks[1][c]&masks[2][d])==0 for (a,b),(c,d) in product(left,right))
 good=[(a,b,d) for (a,b),(c,d) in product(left,right) if b==c]
 distinct=sum((masks[0][a]&masks[1][b]&masks[2][d])==0 for a,b,d in good)
 overlap=dict(local_state_counts=list(map(len,states)),left_rows=len(left),right_rows=len(right),independent_rows=len(left)*len(right),fabricated_distinct_count=bad,repaired_rows=len(good),repaired_distinct_count=distinct)
 source=json.loads((P/'sources/gluing-certificate.json').read_text())['hostile']
 need(all(source[k]==v for k,v in overlap.items()),'overlap source mismatch')
 edges=list(combinations(range(5),2));locked={(i,i+1) for i in range(4)};free=[e for e in edges if e not in locked]
 fibre=[]
 for bits in product((0,1),repeat=len(free)):
  arcs=locked|{e if b else e[::-1] for e,b in zip(free,bits)}
  score=[sum(a==i for a,b in arcs) for i in range(5)]
  if score!=[3,2,1,2,2]:continue
  H=sum(all((p[i],p[i+1]) in arcs for i in range(4)) for p in permutations(range(5)))
  tri=[p for p in permutations(range(5),3) if p[0]==min(p) and all((p[i],p[(i+1)%3]) in arcs and tuple(sorted((p[i],p[(i+1)%3]))) in free for i in range(3))]
  fibre.append(dict(arcs=sorted(map(list,arcs)),Hamilton_paths=H,directed_triangle_moves=len(tri)))
 fibre.sort(key=lambda a:a['Hamilton_paths'])
 source=json.loads((P/'sources/cycle-certificate.json').read_text())['sharpness']
 need(len(fibre)==2 and [x['Hamilton_paths'] for x in fibre]==[x['Hamilton_paths'] for x in source['entire_fibre']],'cycle source mismatch')
 cycle=[0,4,1,3];A=set(map(tuple,fibre[0]['arcs']));B=set(map(tuple,fibre[1]['arcs']));C={(cycle[i],cycle[(i+1)%4]) for i in range(4)}
 need(C<=A and (A-C)|{(b,a) for a,b in C}==B,'wrong cycle reversal')
 need(all(x['directed_triangle_moves']==0 for x in fibre),'triangle present')
 for x,y in zip(fibre,source['entire_fibre']):need(sorted([list(e) for e in edges if e not in map(tuple,x['arcs'])])==y['reversed_forward_edges'],'wrong arcs')
 return dict(overlap=overlap,locked_path=[0,1,2,3,4],ordered_score=[3,2,1,2,2],fibre=fibre,cycle_reversal=cycle+[0])

BG='#f6f8fc';INK='#172b43';MUT='#50637b';TEAL='#007f79';RED='#b34836';BLUE='#326ba8'
def start(title,desc,h=800):
 return [f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="{h}" viewBox="0 0 1200 {h}" role="img" aria-labelledby="title desc"><title id="title">{html.escape(title)}</title><desc id="desc">{html.escape(desc)}</desc><defs>'+''.join(f'<marker id="{name}" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0 0L8 4L0 8Z" fill="{col}"/></marker>' for name,col in [('teal',TEAL),('mut',MUT),('blue',BLUE),('red',RED)])+f'</defs><rect width="1200" height="{h}" fill="{BG}"/>']
def text(s,x,y,t,size=22,col=INK,weight=400,anchor='start'):
 s.append(f'<text x="{x}" y="{y}" fill="{col}" font-family="DejaVu Sans, sans-serif" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}">{html.escape(t)}</text>')
def box(s,x,y,w,h,stroke='#d4dfe9',fill='white'):
 s.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="16" fill="{fill}" stroke="{stroke}" stroke-width="2"/>')
def arrow(s,x1,y1,x2,y2,c='teal',dash=False,width=3):
 col=dict(teal=TEAL,mut=MUT,blue=BLUE,red=RED)[c]
 s.append(f'<path d="M{x1:.3f},{y1:.3f} L{x2:.3f},{y2:.3f}" stroke="{col}" stroke-width="{width}" fill="none" marker-end="url(#{c})"'+(' stroke-dasharray="7 6"' if dash else '')+'/>')
def overlap_svg(F):
 s=start('Shared choice, not independent resampling','Exact overlap counts for CRT blocks [2,5] and [5,3]. Independent pairing fabricates 128 distinct outcomes; the same F5 state join recovers 112 true arrays and zero distinct arrays.')
 text(s,45,55,'01  /  SHARED CHOICE ≠ INDEPENDENT RESAMPLING',29,weight=700)
 text(s,45,92,'Finite-exact CRT witness: local state counts 8 × 7 × 2',22,MUT)
 box(s,45,125,1110,150)
 text(s,85,169,'Block [2,5]',24,weight=700);text(s,85,215,'8 × 7 = 56 rows',27)
 text(s,600,166,'ONE physical F₅ choice',22,TEAL,700,'middle')
 text(s,600,206,'7 actual states',22,TEAL,400,'middle')
 text(s,880,169,'Block [5,3]',24,weight=700);text(s,880,215,'7 × 2 = 14 rows',27)
 text(s,600,250,'The shared state ID must agree in both occurrences.',19,MUT,anchor='middle')
 box(s,45,315,535,240,RED);box(s,620,315,535,240,TEAL)
 arrow(s,300,278,300,310,'red');arrow(s,890,278,890,310)
 text(s,75,355,'WRONG COUPLING',18,RED,700);text(s,650,355,'EXACT COUPLING',18,TEAL,700)
 text(s,75,396,'Resample F₅ independently',24,weight=700);text(s,650,396,'Join on the same F₅ state',24,weight=700)
 text(s,75,440,'56 × 14 = 784 pairs',28,RED,700);text(s,650,440,'8 × 7 × 2 = 112 true arrays',27,TEAL,700)
 text(s,75,482,'128 fabricated all-distinct outcomes',21,RED);text(s,650,482,'0 all-distinct arrays',23,TEAL)
 text(s,75,522,'Here, independently paired zero masks are fake.',18,MUT);text(s,650,522,'Compatible union first; count each choice once.',18,MUT)
 box(s,45,590,1110,154)
 text(s,70,627,'FIXED MASK',18,MUT,700);text(s,300,627,'m AND m = m',23,weight=700)
 text(s,70,665,'CORRELATED COPY',18,TEAL,700);text(s,300,665,'Choose s once, then use (s,s): meet returns the same law.',21)
 text(s,70,703,'RESAMPLING',18,RED,700);text(s,300,703,'Choose s and t separately: a different coupling, not idempotence.',21)
 text(s,45,779,'Scope: this exact product-family repair does not reconstruct arbitrary global laws from marginals.',17,MUT)
 return '\n'.join(s+['</svg>'])+'\n'
def cycle_svg(F):
 s=start('A locked path and equal scores do not fix Hamilton count','Two full five-vertex tournaments with fixed path 0 to 1 to 2 to 3 to 4 and score vector (3,2,1,2,2). Their only fibre transition is a directed quadrilateral reversal; Hamilton path counts are 11 and 13.',900)
 text(s,45,55,'02  /  SAME LOCKED PATH + SCORES, DIFFERENT H',29,weight=700)
 text(s,45,94,'H counts Hamilton paths. This fibre has two states and no legal residual triangle move.',22,MUT)
 box(s,45,123,500,535);box(s,655,123,500,535)
 for off,state,label in [(0,F['fibre'][0],'A'),(610,F['fibre'][1],'B')]:
  text(s,80+off,163,f'STATE {label}     H = {state["Hamilton_paths"]}',28,weight=700)
  pos={0:(125+off,235),4:(465+off,235),1:(465+off,450),3:(125+off,450),2:(295+off,565)}
  cyc={frozenset(e) for e in [(0,4),(4,1),(1,3),(3,0)]}
  for a,b in state['arcs']:
   lock=b==a+1;typ='mut' if lock else ('teal' if frozenset((a,b)) in cyc else 'blue')
   x,y=pos[a];xx,yy=pos[b];d=math.hypot(xx-x,yy-y);u,v=(xx-x)/d,(yy-y)/d
   arrow(s,x+24*u,y+24*v,xx-27*u,yy-27*v,typ,lock,2.5 if lock else 3.5)
  for i,(x,y) in pos.items():
   s.append(f'<circle cx="{x}" cy="{y}" r="23" fill="white" stroke="{INK}" stroke-width="2"/>');text(s,x,y+8,str(i),24,weight=700,anchor='middle')
  text(s,295+off,622,'Ordered score: (3, 2, 1, 2, 2)',22,anchor='middle')
 arrow(s,560,368,640,368);arrow(s,640,389,560,389)
 text(s,600,340,'flip',18,TEAL,700,'middle')
 text(s,45,703,'Locked path (dashed): 0 → 1 → 2 → 3 → 4',22,MUT)
 text(s,45,741,'Teal cycle in A: 0 → 4 → 1 → 3 → 0. Reverse all four arcs to get B.',22,TEAL)
 text(s,45,779,'Blue arcs stay fixed: 0 → 2 and 4 → 2. All ten tournament arcs are shown.',21,BLUE)
 text(s,45,828,'All n: directed 3-/4-cycle moves on free edges connect every fixed-score, forest-locked fibre.',20,weight=700)
 text(s,45,863,'They preserve the locked forest and scores, but need not preserve H mod 8 or full native wiring.',20,RED)
 return '\n'.join(s+['</svg>'])+'\n'
def products():
 F=facts();return {'counts.json':json.dumps(F,indent=2,sort_keys=True)+'\n','shared-choice-overlap.svg':overlap_svg(F),'path-locked-score-fibre.svg':cycle_svg(F)}
def main():
 out=products()
 for name,data in out.items():
  if '--verify' in sys.argv:need((P/name).read_text()==data,'generated mismatch: '+name)
  else:(P/name).write_text(data)
 print('PASS: recomputed exact counts, all tournament arcs, entire fibre, H counts, no triangle move, and deterministic SVGs.')
if __name__=='__main__':main()
