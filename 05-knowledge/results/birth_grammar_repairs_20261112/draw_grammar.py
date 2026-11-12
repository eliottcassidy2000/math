"""Export the declared construction graph, with its scope printed on the figure."""
from pathlib import Path
import gzip,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
HERE=Path(__file__).resolve().parent

def run():
    d=json.loads(gzip.decompress((HERE/'birth_inventory.json.gz').read_bytes()))
    if d['retained_births_through14']!=4042:raise RuntimeError('figure census changed')
    fig,ax=plt.subplots(figsize=(15,8),dpi=160)
    fig.patch.set_facecolor('#f8fafc');ax.set_facecolor('#f8fafc')
    ax.set_xlim(0,15);ax.set_ylim(0,8);ax.axis('off')
    ax.text(.3,7.65,'A finite grammar with infinitely many stable-clone births',fontsize=19,weight='bold',color='#10233f')
    ax.text(.3,7.18,'Complete old bank: 1,150 through order 10. Retained union: 4,042 through order 14.',fontsize=12,color='#34465b')
    ax.text(.3,6.83,'After order 10, counts are complete only inside the declared construction routes below.',fontsize=11,color='#34465b')
    def box(x,y,w,h,s,color='#dfeafa'):
        ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.08,rounding_size=.12',facecolor=color,edgecolor='#8da5c2',linewidth=1.1))
        ax.text(x,y,s,ha='center',va='center',fontsize=11,color='#10233f',linespacing=1.45)
    def arrow(x1,y1,x2,y2,label=None,ly=.16):
        ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops={'arrowstyle':'->','color':'#476481','lw':1.6})
        if label:ax.text((x1+x2)/2,(y1+y2)/2+ly,label,ha='center',va='bottom',fontsize=9,color='#34465b')
    for x,s in [(1.35,'Order 10\n850 births'),(4.15,'Order 11\n464 births'),(6.95,'Order 12\n834 births'),(9.75,'Order 13\n1,414 births')]:box(x,5.8,2.05,.9,s)
    for x in [1.35,4.15,6.95]:arrow(x+1.1,5.8,x+1.7,5.8)
    ax.text(5.6,6.43,'All one-vertex orientations; both markings; strong births only',ha='center',fontsize=10,color='#34465b')
    box(13,5.8,2.3,.9,'Order 14\n44 births')
    arrow(10.85,5.8,11.77,5.8)
    ax.text(11.38,6.37,'Only the 8 parents\nwith five stable ports',ha='center',fontsize=9,color='#34465b')
    box(1.35,3.65,2.05,.95,'12 old births\n24 coupled charts','#e0f0eb')
    arrow(1.35,5.28,1.35,4.2,'subset',.02)
    box(4.55,3.65,2.2,.95,'Order 12\n24 new births','#e0f0eb')
    arrow(2.47,3.65,3.35,3.65,'synchronize')
    box(7.85,3.65,2.2,.95,'Order 13\n88 new births','#e0f0eb')
    arrow(5.75,3.65,6.65,3.65,'extend')
    box(12,3.65,3.6,.95,'46 distinct coupled charts\nwith a free stable port','#e0f0eb')
    arrow(9.05,3.65,10.1,3.65,'certify')
    box(4.55,1.8,3.4,.95,'24 births at every\neven order 12, 14, 16, …','#cdebe1')
    arrow(4.55,3.08,4.55,2.37)
    box(12,1.8,3.6,.95,'46 births at every\nodd order 13, 15, 17, …','#cdebe1')
    arrow(12,3.08,12,2.37)
    ax.text(.3,.62,'Infinite subgrammar:  (24 z¹² + 46 z¹³) / (1 − z²)',fontsize=15,weight='bold',color='#155b48')
    ax.text(.3,.15,'These are marked contraction births, not all tournaments or all newly attained path counts. Chan–Pak remains open.',fontsize=10,color='#34465b')
    fig.subplots_adjust(left=.025,right=.99,top=.98,bottom=.03)
    fig.savefig(HERE/'birth-grammar.png',facecolor=fig.get_facecolor())
    fig.savefig(HERE/'birth-grammar.svg',facecolor=fig.get_facecolor(),metadata={'Date':None})
    plt.close(fig)

if __name__=='__main__':run()
