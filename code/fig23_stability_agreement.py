import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.8,
 'xtick.major.width':.8,'ytick.major.width':.8,'axes.edgecolor':'#333333',
 'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
A=json.load(open('results_A.json'))['A']; B=json.load(open('results_B.json'))['B']
METH=['ANOVA','mRMR','ReliefF','LASSO','RF']
LBL={'ANOVA':'ANOVA F','mRMR':'mRMR','ReliefF':'ReliefF','LASSO':'LASSO','RF':'RF'}
COL={'ANOVA':'#2a78d6','mRMR':'#eb6834','ReliefF':'#1baf7a','LASSO':'#eda100','RF':'#e87ba4'}
MK={'ANOVA':'o','mRMR':'s','ReliefF':'^','LASSO':'D','RF':'v'}
KS=[10,25,50,100]

# ---------- Figure 2 ----------
fig,axes=plt.subplots(1,2,figsize=(7.2,3.0),sharey=True)
for ax,(R,nm,sub) in zip(axes,[(A,'GSE28497','7-class B-ALL, n = 281'),(B,'GSE71935','JMML vs normal, n = 46')]):
    ends=[]
    for m in METH:
        v=[R['within'][m][str(k)]['kuncheva'] for k in KS]
        ax.plot(range(len(KS)),v,color=COL[m],marker=MK[m],ms=4.2,lw=1.8,label=LBL[m],
                markeredgecolor='white',markeredgewidth=.7,clip_on=False,zorder=3)
        ends.append([v[-1],m])
    ends.sort(key=lambda z:z[0]); MIN=.052
    for i in range(1,len(ends)):
        if ends[i][0]-ends[i-1][0]<MIN: ends[i][0]=ends[i-1][0]+MIN
    off=(ends[-1][0]-1.0)
    if off>0:
        for e in ends: e[0]-=off
    for yy,m in ends:
        ax.annotate(LBL[m],(len(KS)-1,yy),xytext=(6,0),textcoords='offset points',
                    fontsize=7,color=COL[m],va='center',fontweight='bold',annotation_clip=False)
    ax.set_xticks(range(len(KS))); ax.set_xticklabels(KS)
    ax.set_xlabel('Signature size $k$ (probes)')
    ax.set_ylim(0,1.0); ax.set_xlim(-.15,len(KS)-1+.02)
    ax.grid(axis='y',color='#E6E6E6',lw=.7,zorder=0); ax.set_axisbelow(True)
    for s in ['top','right']: ax.spines[s].set_visible(False)
    ax.set_title(f'{nm}\n{sub}',fontsize=8.5,pad=6)
axes[0].set_ylabel('Kuncheva stability index')
for ax in axes: ax.margins(x=.30)
fig.tight_layout(); fig.subplots_adjust(right=.88)
fig.savefig('figs/Figure2.png',dpi=400,bbox_inches='tight')
fig.savefig('figs/Figure2.tif',dpi=400,bbox_inches='tight'); plt.close(fig)

# ---------- Figure 3 ----------
seq=LinearSegmentedColormap.from_list('bl',['#F4F8FC','#c3dcf2','#7fb3e5','#3d86d3','#14487f'])
fig=plt.figure(figsize=(7.4,3.15))
gs=fig.add_gridspec(2,3,width_ratios=[1,1,1.12],height_ratios=[1,.055],wspace=.62,hspace=.06)
for j,(R,nm) in enumerate([(A,'GSE28497'),(B,'GSE71935')]):
    ax=fig.add_subplot(gs[0,j])
    M=np.full((5,5),np.nan)
    for i1,m1 in enumerate(METH):
        for i2,m2 in enumerate(METH):
            if i1==i2: M[i1,i2]=np.nan; continue
            key=f"{m1}|{m2}" if f"{m1}|{m2}" in R['cross']['50'] else f"{m2}|{m1}"
            M[i1,i2]=R['cross']['50'][key]
    im=ax.imshow(M,cmap=seq,vmin=0,vmax=.8)
    ax.set_xticks(range(5)); ax.set_yticks(range(5))
    ax.set_xticklabels([LBL[m] for m in METH],rotation=45,ha='right',fontsize=7)
    ax.set_yticklabels([LBL[m] for m in METH],fontsize=7)
    for i1 in range(5):
        for i2 in range(5):
            if np.isnan(M[i1,i2]):
                ax.add_patch(plt.Rectangle((i2-.5,i1-.5),1,1,facecolor='#EFEFEF',edgecolor='white',lw=1.4))
                continue
            ax.text(i2,i1,f"{M[i1,i2]:.2f}",ha='center',va='center',fontsize=6.6,
                    color='white' if M[i1,i2]>.45 else '#1a1a1a')
    ax.set_xticks(np.arange(-.5,5,1),minor=True); ax.set_yticks(np.arange(-.5,5,1),minor=True)
    ax.grid(which='minor',color='white',lw=1.4); ax.tick_params(which='minor',length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(f'({chr(97+j)}) {nm}',fontsize=8.5,pad=6)
cax=fig.add_subplot(gs[1,0:2]); cb=fig.colorbar(im,cax=cax,orientation='horizontal')
cb.set_label('Jaccard index, $k$ = 50',fontsize=7,labelpad=2); cb.ax.tick_params(labelsize=6.5,length=2.5)
ax=fig.add_subplot(gs[0:2,2])
w=.38; x=np.arange(1,6)
for j,(R,nm,c) in enumerate([(A,'GSE28497','#2a78d6'),(B,'GSE71935','#eb6834')]):
    tot=sum(R['multiplicity']['50'].values())
    v=[100*R['multiplicity']['50'][str(n)]/tot for n in range(1,6)]
    bars=ax.bar(x+(j-.5)*w,v,w*.92,color=c,label=nm,zorder=3)
    for b,val in zip(bars,v):
        ax.text(b.get_x()+b.get_width()/2,val+1.2,f"{val:.0f}",ha='center',fontsize=6.2,color='#333')
ax.set_xticks(x); ax.set_xlabel('Selected by $n$ of 5 methods'); ax.set_ylabel('% of union set')
ax.grid(axis='y',color='#E6E6E6',lw=.7,zorder=0); ax.set_axisbelow(True)
for s in ['top','right']: ax.spines[s].set_visible(False)
ax.legend(frameon=False,fontsize=7,loc='upper right')
ax.set_title('(c) Selection multiplicity',fontsize=8.5,pad=6)
fig.savefig('figs/Figure3.png',dpi=400,bbox_inches='tight')
fig.savefig('figs/Figure3.tif',dpi=400,bbox_inches='tight'); plt.close(fig)
print("figs 2,3 done")
