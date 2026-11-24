import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.4,'figure.facecolor':'white'})
fig,ax=plt.subplots(figsize=(7.2,4.4)); ax.set_xlim(0,100); ax.set_ylim(0,64); ax.axis('off')
BL='#2a78d6'; OR='#eb6834'; GY='#48545F'; LG='#F2F5F8'
def box(x,y,w,h,txt,fc=LG,ec='#C9D2DA',fs=7.4,bold=False,tc='#151A20'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0,rounding_size=1.2",
        fc=fc,ec=ec,lw=.9,zorder=2))
    ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=fs,zorder=3,
            color=tc,fontweight='bold' if bold else 'normal',linespacing=1.5)
def arrow(x1,y1,x2,y2,c=GY):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=8,
        lw=.9,color=c,zorder=1,shrinkA=0,shrinkB=0))
def lab(x,y,t,c=GY,fs=6.6,ha='center',w='bold'):
    ax.text(x,y,t,ha=ha,va='center',fontsize=fs,color=c,fontweight=w,zorder=4)

lab(3,61.5,'A',fs=10,ha='left',c='#151A20')
box(9,56,38,7,'GSE28497  ·  CuMiDa\n281 samples · 7 B-ALL subtypes',fc='#E8F1FB',ec='#B9D5F0')
box(53,56,38,7,'GSE71935  ·  CuMiDa\n46 samples · JMML vs normal',fc='#FCEDE6',ec='#F3CDB9')
box(9,47,82,6.5,'Restrict to the 22,215 probe sets shared by GPL96 and GPL570  ·  remove AFFX controls  ·  verify log$_2$ scale',fs=6.9)
arrow(28,56,28,53.5); arrow(72,56,72,53.5)

lab(3,43.5,'B',fs=10,ha='left',c='#151A20')
box(9,36,82,6.5,'30 stratified subsamples at 80%  ·  unsupervised variance prefilter to 2,000 probes inside every training split',fs=6.9)
arrow(50,47,50,42.5)

lab(3,32.5,'C',fs=10,ha='left',c='#151A20')
xs=[9,25.6,42.2,58.8,75.4]
names=['ANOVA $F$','mRMR','ReliefF','LASSO','Random\nforest']
kinds=['filter','filter','filter','embedded','embedded']
cols=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4']
for x,n,k,c in zip(xs,names,kinds,cols):
    box(x,23,15.6,7.5,n,fc='white',ec=c,fs=7.4)
    ax.add_patch(plt.Rectangle((x,23),15.6,.9,fc=c,ec='none',zorder=3))
    lab(x+7.8,21.4,k.upper(),c=c,fs=6.0)
arrow(50,36,50,30.5)
for x in xs: arrow(50,32.2,x+7.8,30.6)

lab(3,18,'D',fs=10,ha='left',c='#151A20')
box(9,10.5,25,7,'Within-method\nstability\nKuncheva · Jaccard',fc='#F2F5F8')
box(37.5,10.5,25,7,'Between-method\nagreement\npairwise Jaccard · multiplicity',fc='#F2F5F8',fs=7.0)
box(66,10.5,25,7,'Consequences\nper-subtype stability ·\nnested-CV accuracy',fc='#F2F5F8',fs=7.0)
ax.plot([16.8,83.2],[18.9,18.9],color=GY,lw=.9,zorder=1)
for x in xs: ax.plot([x+7.8,x+7.8],[23,18.9],color=GY,lw=.9,zorder=1)
for x in [21.5,50,78.5]: arrow(x,18.9,x,17.7)

box(9,1.5,82,6,'Primary question  —  do the five algorithms select the same probe sets?',fc='#EDF1F5',ec='#9FB0BE',bold=True,fs=8)
arrow(50,10.5,50,7.6)
fig.savefig('figs/Figure1.png',dpi=400,bbox_inches='tight')
fig.savefig('figs/Figure1.tif',dpi=400,bbox_inches='tight')
print('fig1 done')
