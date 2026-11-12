"""Scientific illustration of the two proved unit loci."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent

def run():
    fig,axes=plt.subplots(1,2,figsize=(11,5.3),layout='constrained')
    for ax,center,title,formula in zip(axes,[1,2],['Stable ports: joint growth fails','Contractible blocks: ancestry branches'],[r'$\delta=1-(2^a-2)(2^b-2)$',r'$\delta=1+(2^a-4)(2^b-4)$']):
        for a in range(1,6):
            for b in range(1,6):
                unit=a==center or b==center
                ax.scatter(a,b,s=90 if unit else 28,color='#176ca4' if unit else '#cfd6dc',zorder=3)
        ax.axvline(center,color='#176ca4',alpha=.45,lw=2);ax.axhline(center,color='#176ca4',alpha=.45,lw=2)
        ax.set(xlim=(.55,5.45),ylim=(.4,5.5),xticks=range(1,6),yticks=range(1,6),xlabel='First block length a',ylabel='Second block length b',aspect='equal')
        ax.set_title(title+'\n'+formula,fontsize=12,pad=12);ax.grid(alpha=.15)
        for spine in ax.spines.values():spine.set_color('#aab5be')
    axes[0].scatter(1,1,s=240,marker='*',color='#bb5127',zorder=5)
    axes[0].annotate('Birth: order 9',xy=(1,1),xytext=(1.7,.65),fontsize=10,color='#90401f')
    axes[0].scatter(2,2,marker='x',s=100,lw=2,color='#bb5127',zorder=5)
    axes[0].annotate(r'$\delta(2,2)=-3$',xy=(2,2),xytext=(2.4,2.25),fontsize=11,color='#90401f')
    axes[1].scatter(2,2,s=180,facecolors='white',edgecolors='#bb5127',lw=2,zorder=5)
    for a,b in [(1,2),(2,1)]:
        axes[1].scatter(a,b,s=240,marker='*',color='#bb5127',zorder=5)
        axes[1].annotate('',xy=(a,b),xytext=(2,2),arrowprops={'arrowstyle':'->','color':'#bb5127','lw':2,'shrinkA':10,'shrinkB':10})
    axes[1].annotate('Order 11\nH = 2729',xy=(2,2),xytext=(2.35,2.4),fontsize=10,color='#90401f')
    axes[1].text(.65,1.45,'Birth A\nU = 8158',fontsize=9,color='#90401f')
    axes[1].text(2.3,.65,'Birth B\nU = 10718',fontsize=9,color='#90401f')
    fig.suptitle('The blue lattice points are exactly the unit seeds',fontsize=15)
    fig.savefig(HERE/'crossings.png',dpi=180,facecolor='white')
    plt.close(fig)
if __name__=='__main__':run()
