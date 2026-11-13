"""Publication-style view of the smallest informative metagraph and growth states."""
from pathlib import Path
import gzip,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch,Arc
HERE=Path(__file__).resolve().parent

def main():
    d=json.loads(gzip.decompress((HERE/'rooted_growth.json.gz').read_bytes()))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(13,5.6),gridspec_kw={'width_ratios':[1.05,1]})
    ax=axes[0];ax.set_xlim(-1.8,1.8);ax.set_ylim(-1.65,2.1);ax.axis('off')
    pos={0:(-1.05,-.65),2:(1.05,-.65),4:(0,.65)};blue='#2875c7';black='#27313e'
    ax.plot([-1.05,0],[-.65,.65],color=blue,lw=3,zorder=1)
    for rad in (-.15,.15):
        ax.add_patch(FancyArrowPatch(pos[2],pos[4],arrowstyle='-',connectionstyle=f'arc3,rad={rad}',lw=2.4,color=black,shrinkA=28,shrinkB=30))
    ax.add_patch(Arc((0,.94),.6,.65,theta1=-30,theta2=210,color=blue,lw=3))
    for code,(x,y) in pos.items():
        ax.add_patch(Circle((x,y),.29,facecolor='white',edgecolor=black,lw=1.5,zorder=3))
        ax.text(x,y,str({0:1,2:2,4:5}[code]),ha='center',va='center',fontsize=23,zorder=4)
    ax.text(-1.05,-1.14,'Transitive\nH = 1',ha='center',va='top')
    ax.text(1.05,-1.14,'Converse pair\nH = 3, |Aut| = 3',ha='center',va='top')
    ax.text(-1.13,.86,'Strong\nself-converse\nH = 5',ha='center',va='top',fontsize=10.5)
    ax.text(-1.65,2.02,'Order 4: three merged nodes',fontsize=15,fontweight='bold',va='top')
    ax.text(-1.65,1.75,'Circle numbers are tiling weights.',color='#4c5968',va='top')
    ax.text(-.88,-.03,'blue',color=blue,rotation=52,ha='center')
    ax.text(.95,.18,'two black\nmatches',ha='center',fontsize=10)
    ax=axes[1];n=[r['n'] for r in d['levels']]
    ax.semilogy(n,[r['tilings'] for r in d['levels']],'-o',color='#b77b28',label='Full path presentations')
    ax.semilogy(n,[r['end_rooted_types'] for r in d['levels']],'-o',color=blue,label='End states: exact forward growth')
    ax.semilogy(n,[r['ordinary_types'] for r in d['levels']],'-o',color=black,label='Ordinary nodes: A000568')
    ax.set_xticks(n);ax.set_xlabel('Tournament order');ax.set_ylabel('Number of states (log scale)');ax.grid(alpha=.15)
    ax.set_title('The information needed depends on the task',loc='left',fontsize=13,fontweight='bold',pad=16)
    ax.legend(loc='upper left',frameon=False,fontsize=9)
    fig.text(.52,.02,'End states retain exact forward counts; complement matches still require internal path order.',ha='center',fontsize=11,color='#364454')
    fig.subplots_adjust(left=.03,right=.98,bottom=.15,top=.88,wspace=.19)
    fig.savefig(HERE/'metagraph_levels.png',dpi=180,facecolor='white')
    fig.savefig(HERE/'metagraph_levels.svg',facecolor='white')
if __name__=='__main__':main()
