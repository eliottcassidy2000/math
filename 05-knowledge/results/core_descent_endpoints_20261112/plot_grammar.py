"""Small exportable diagram of the certified grammar, not the full metagraph."""
from pathlib import Path
import gzip,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
HERE=Path(__file__).resolve().parent
d=json.loads(gzip.decompress((HERE/'descent_grammar.json.gz').read_bytes()))
fig,ax=plt.subplots(figsize=(12,6.2));fig.patch.set_facecolor('#fafaf7');ax.set_facecolor('#fafaf7')
positions={'even:0':(0,2.6),'odd:0':(-1,1),'odd:1':(1,1),'even:5':(4,2.6),'odd:15':(3,1),'odd:16':(5,1),'odd:2':(8,2.8),'odd:3':(8,1.6),'odd:17':(8,.4)}
charts=json.loads(gzip.decompress((HERE/'synchronized_descent.json.gz').read_bytes()))['records']
bases={r['family']:r['anchor']+2*sum(b*c for b,c in r['B_odd']) for r in charts}
for child,parent,text in [('odd:0','even:0',r'$(5c+24)/13$'),('odd:1','even:0',r'$(c+6)/3$'),('odd:15','even:5',r'$(5c+24)/13$'),('odd:16','even:5',r'$(c+6)/3$')]:
    x,y=positions[child];u,v=positions[parent]
    ax.add_patch(FancyArrowPatch((x,y+.34),(u,v-.34),arrowstyle='-|>',mutation_scale=16,color='#37728a',linewidth=1.8))
    ax.text((x+u)/2+(-.36 if x<u else .36),(y+v)/2,text,ha='center',va='center',fontsize=11,color='#153d51',bbox={'facecolor':'#fafaf7','edgecolor':'none','pad':1.5})
for name,(x,y) in positions.items():
    terminal=name in d['terminal_families'];color='#dcecf2' if terminal else '#edf0e8'
    ax.add_patch(FancyBboxPatch((x-.85,y-.32),1.7,.64,boxstyle='round,pad=.06,rounding_size=.08',facecolor=color,edgecolor='#37728a' if terminal else '#879582',linewidth=1.3))
    ax.text(x,y+.1,name,ha='center',va='center',fontsize=12,weight='bold',color='#183342')
    ax.text(x,y-.16,f'base coefficient {bases[name]:,}',ha='center',va='center',fontsize=9.5,color='#384e56')
ax.text(-1.9,4,'Nine seed charts, five terminal charts',fontsize=21,weight='bold',color='#183342')
ax.text(-1.9,3.55,'Four certified inverse changes of core; blue boxes are terminal chart types.',fontsize=12,color='#4e6064')
ax.text(-1.9,-.15,'Every chart also has decreasing coordinate moves:',fontsize=12,weight='bold',color='#183342')
ax.text(-1.9,-.57,r'$k>1:\ (m,k)\mapsto(m,k-1)$     |     $m>1:\ (m,k)\mapsto(m-1,k)$',fontsize=13,color='#183342')
ax.text(-1.9,-1.03,'Retain the selected gate and ray threshold. This is a sparse subgrammar; Chan–Pak remains open.',fontsize=10.5,color='#5d686b')
ax.set_xlim(-2.1,9.2);ax.set_ylim(-1.3,4.35);ax.axis('off');fig.tight_layout(pad=.5)
fig.savefig(HERE/'grammar.png',dpi=180,facecolor=fig.get_facecolor())
fig.savefig(HERE/'grammar.svg',facecolor=fig.get_facecolor())
svg=HERE/'grammar.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf8').splitlines())+'\n',encoding='utf8',newline='\n')
