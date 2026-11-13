"""Static figures drawn only from the checked mathematical receipts."""
from pathlib import Path
import json,gzip,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
def read(name):return json.loads(gzip.decompress((HERE/(name+'.json.gz')).read_bytes()))

def run():
    a=read('square_covers')['records'];b=read('square_cut_audit')['records']
    x=[r['n'] for r in a];width=[r['minimum_path_cover'] for r in a];cut=[r['maximum_components_minus_deleted'] for r in b]+[1]
    fig,ax=plt.subplots(figsize=(11,4.1));fig.set_facecolor('#fafaf8');ax.set_facecolor('#fafaf8')
    ax.plot(x,width,'o-',color='#1d527a',lw=2.1,ms=5,label='Minimum spanning path-cover width')
    ax.plot(x,cut,'s--',color='#d49b28',lw=1.6,ms=4,label='Best bound from every vertex cut')
    ax.plot([24,24],[1,2],color='#b94143',lw=5,alpha=.5)
    ax.annotate('24: cuts miss the obstruction',xy=(24,1.5),xytext=(13.5,2.6),arrowprops={'arrowstyle':'->','color':'#b94143'},color='#9e3032')
    ax.annotate('25: repair through ports 11, 24',xy=(25,1),xytext=(13.5,.65),arrowprops={'arrowstyle':'->','color':'#1d527a'},color='#1d527a')
    ax.set(xlim=(.5,25.5),ylim=(.45,3.4),xticks=range(1,26),yticks=[1,2,3],xlabel='Vertices 1 through n, joined when their sum is a square',ylabel='Number of paths')
    ax.set_title('The first 25 square-sum graphs',loc='left',weight='bold',pad=15)
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.legend(loc='upper right',frameon=False)
    fig.tight_layout();fig.savefig(HERE/'square_widths.png',dpi=170);plt.close(fig)
    data=read('graceful_padding')['one_leaf_bridge'];old=json.loads((HERE.parent/'graceful_fixed_color_counterexample_20261019/witness.json').read_text(encoding='utf8'))
    F=[tuple(e) for e in old['fixed_forest'].values()];colors=data['free_colors'];states=data['edge_start_tuples'];palette=['#20857f','#d38d21','#9356a0']
    # Keep vertex0 away from the diameter edges of difference6.
    pos={0:(-.12,.28),13:(1.5,1.25)}
    pos.update({v:(math.cos(2*math.pi*(v-1)/12),math.sin(2*math.pi*(v-1)/12)) for v in range(1,13)})
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.7));fig.set_facecolor('#fafaf8')
    for ax,state,title in zip(axes,states,['Old component A','New bridge state','Old component B']):
        ax.set_facecolor('#fafaf8')
        for u,v in F:ax.plot([pos[u][0],pos[v][0]],[pos[u][1],pos[v][1]],color='#aab0b4',lw=1.2,zorder=1)
        ax.plot([pos[0][0],pos[13][0]],[pos[0][1],pos[13][1]],color='#444',lw=1.3,ls='--',zorder=1)
        for c,start,color in zip(colors,state,palette):
            u,v=start,start+c;ax.plot([pos[u][0],pos[v][0]],[pos[u][1],pos[v][1]],color=color,lw=2.7,zorder=2)
        for v,(px,py) in pos.items():
            ax.scatter(px,py,s=260,c='#fff' if v!=13 else '#ffe4a0',edgecolors='#404953',lw=1,zorder=3)
            ax.text(px,py,str(v),ha='center',va='center',fontsize=8.5,zorder=4)
        ax.set_title(title+'\nstarts '+str(tuple(state)),fontsize=11,pad=12)
        ax.set(xlim=(-1.18,1.7),ylim=(-1.2,1.45),aspect='equal');ax.axis('off')
    fig.suptitle('One new frozen star leaf repairs a disconnected graceful fiber',x=.04,ha='left',weight='bold',fontsize=13)
    fig.text(.04,.025,'Gray: original frozen edges    Dashed: new frozen edge (0,13)    Teal / gold / purple: free differences 6 / 9 / 11',fontsize=10)
    fig.tight_layout(rect=[0,.07,1,.92]);fig.savefig(HERE/'graceful_bridge.png',dpi=170);plt.close(fig)
    print('Rendered square_widths.png and graceful_bridge.png')
if __name__=='__main__':run()
