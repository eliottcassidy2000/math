"""A static scientific map: three actual reductions, no numerical identity."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
HERE=Path(__file__).resolve().parent
def run():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12})
    fig,ax=plt.subplots(figsize=(16,9));fig.patch.set_facecolor('#f8fafc');ax.set_facecolor('#f8fafc')
    ax.set_xlim(0,16);ax.set_ylim(0,9);ax.axis('off')
    ax.text(.4,8.65,'The recurring five: keep the operator and its observer',fontsize=22,fontweight='bold',color='#132238')
    ax.text(.4,8.13,'Three concrete constructions. Equal numbers do not identify their state spaces.',fontsize=13,color='#44546a')
    blue='#e2eefb';orange='#fff0da';gray='#e9edf1'
    rows=[
      (6.8,'ELASTICITY',[
       ('Fixed isothermal net','regular holomorphic Z',blue),
       ('5 real parameters','4 amplitudes + pressure',blue),
       ('Commutator signs','sample sites; exclude ties',blue),
       ('Locally transitive T','exact curve-family counts',blue),
       ('Direct image stops','H = 13 is absent',gray)]),
      (4.4,'KING-STRIP PAPER',[
       ('19 physical states','local candidate automaton',orange),
       ('5 persistent modes','ranks: 9, 7, 5, 5',blue),
       ('Reflection split','3 even + 2 odd modes',blue),
       ('3 observed modes','exact rational generating law',blue),
       ('Global metric test','local legality is insufficient',gray)]),
      (2.0,'TOURNAMENT SEEDS',[
       ('5-state marking DFA','zero or one internal gap',blue),
       ('Weighted endpoints','retain defect and slope U',blue),
       ('Stable inflation','full mixed coefficients',blue),
       ('Contract clone twins','161 states through order 8',orange),
       ('Every odd c','construction still OPEN',gray)])]
    xs=[.4,3.55,6.7,9.85,13.0];width=2.6;height=1.12
    for y,label,boxes in rows:
        ax.text(.4,y+.85,label,fontsize=10,fontweight='bold',color='#5a6577')
        for i,(title,sub,col) in enumerate(boxes):
            x=xs[i];ax.add_patch(FancyBboxPatch((x,y-.45),width,height,boxstyle='round,pad=0.12',facecolor=col,edgecolor='#8797ab',linewidth=.9))
            ax.text(x+width/2,y+.26,title,ha='center',va='center',fontsize=12.3,fontweight='bold',color='#132238')
            ax.text(x+width/2,y-.12,sub,ha='center',va='center',fontsize=10.1,color='#36465c')
            if i<4:ax.add_patch(FancyArrowPatch((x+width+.14,y+.1),(xs[i+1]-.15,y+.1),arrowstyle='-|>',mutation_scale=14,color='#56697e'))
    ax.text(.4,5.75,'Required sidecars: global periods, displacement compatibility, magnitudes and sensor placement.',fontsize=11,color='#44546a')
    ax.text(.4,3.35,'The symmetric observer hides two modes; it does not remove the global resolving-set condition.',fontsize=11,color='#44546a')
    ax.text(.4,.95,'A unit defect gives a ray only after its finite head. Joint stability needs more than separate axis tests.',fontsize=11,color='#44546a')
    ax.text(.4,.34,'Blue: proved reduction / exact law     Amber: supplied or finite census carrier     Grey: boundary or open target',fontsize=11,color='#44546a')
    fig.savefig(HERE/'connection_map.png',dpi=160,bbox_inches='tight',facecolor=fig.get_facecolor())
    fig.savefig(HERE/'connection_map.svg',bbox_inches='tight',facecolor=fig.get_facecolor())
    svg=HERE/'connection_map.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf8').splitlines())+'\n',encoding='utf8')
    plt.close(fig)
if __name__=='__main__':run()
