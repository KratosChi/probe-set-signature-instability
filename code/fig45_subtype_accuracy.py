import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr, wilcoxon
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.8,
 'axes.edgecolor':'#333333','figure.facecolor':'white','savefig.facecolor':'white'})
A=json.load(open('results_A.json')); probes=np.array(A['probes']); A=A['A']
PC=json.load(open('results_perclass.json')); perclass=PC['perclass']; corr=PC['corr']
PF=json.load(open('results_perf.json'))
METH=['ANOVA','mRMR','ReliefF','LASSO','RF']
LBL={'ANOVA':'ANOVA F','mRMR':'mRMR','ReliefF':'ReliefF','LASSO':'LASSO','RF':'RF'}
COL={'ANOVA':'#2a78d6','mRMR':'#eb6834','ReliefF':'#1baf7a','LASSO':'#eda100','RF':'#e87ba4'}
MK={'ANOVA':'o','mRMR':'s','ReliefF':'^','LASSO':'D','RF':'v'}

stab=[A['within'][m]['50']['kuncheva'] for m in METH]
bacc=[PF[m]['bacc_mean'] for m in METH]
f1=[PF[m]['f1_mean'] for m in METH]
r,p=spearmanr(stab,bacc); r2,p2=spearmanr(stab,f1)
print(f"Spearman stability vs bAcc: rho={r:.3f} p={p:.3f}")
print(f"Spearman stability vs macroF1: rho={r2:.3f} p={p2:.3f}")
best=max(METH,key=lambda m:PF[m]['bacc_mean']); worst=min(METH,key=lambda m:PF[m]['bacc_mean'])
w=wilcoxon(PF[best]['bacc_all'],PF[worst]['bacc_all'])
print(f"Wilcoxon {best} vs {worst} bAcc across 15 folds: W={w.statistic:.1f} p={w.pvalue:.4f}")
mostS='ReliefF'; leastS='RF'
w2=wilcoxon(PF[mostS]['bacc_all'],PF[leastS]['bacc_all'])
print(f"Wilcoxon most-stable({mostS}) vs least-stable({leastS}): W={w2.statistic:.1f} p={w2.pvalue:.4f}")
rng=max(bacc)-min(bacc); rngs=max(stab)-min(stab)
print(f"bAcc range {min(bacc):.3f}-{max(bacc):.3f} (spread {rng:.3f}); stability spread {rngs:.3f}")
json.dump({'spearman_stab_bacc':[float(r),float(p)],'spearman_stab_f1':[float(r2),float(p2)],
 'wilcoxon_best_worst':[float(w.statistic),float(w.pvalue),best,worst],
 'wilcoxon_moststable_leaststable':[float(w2.statistic),float(w2.pvalue)],
 'bacc_spread':float(rng),'stab_spread':float(rngs)},open('results_stats.json','w'),indent=1)

# ---------- Figure 4 ----------
CL=[c for c in perclass]
SHORT={'B-CELL_ALL':'Unclassified','B-CELL_ALL_ETV6-RUNX1':'ETV6-RUNX1','B-CELL_ALL_HYPERDIP':'Hyperdiploid',
 'B-CELL_ALL_T-ALL':'T-ALL','B-CELL_ALL_TCF3-PBX1':'TCF3-PBX1','B-CELL_ALL_HYPO':'Hypodiploid','B-CELL_ALL_MLL':'MLL'}
fig,axes=plt.subplots(1,2,figsize=(7.2,3.1),sharey=True)
for ax,(xkey,xlab) in zip(axes,[('n','Samples in subtype'),('silhouette','Mean silhouette width (PC1–50)')]):
    for m in ['ANOVA','RF','LASSO']:
        xs=[perclass[c][xkey] for c in CL]; ys=[perclass[c][m] for c in CL]
        ax.scatter(xs,ys,s=30,facecolor=COL[m],edgecolor='white',lw=.7,marker=MK[m],label=LBL[m],zorder=3)
    ax.set_xlabel(xlab); ax.grid(color='#EAEAEA',lw=.7,zorder=0); ax.set_axisbelow(True)
    for s in ['top','right']: ax.spines[s].set_visible(False)
    ax.set_ylim(0,1.16); ax.set_yticks([0,.2,.4,.6,.8,1.0])
def place(ax,xkey,thr):
    pts=sorted(CL,key=lambda c: perclass[c][xkey])
    placed=[]   # (x, y, dy)
    for c in pts:
        x=perclass[c][xkey]; y=perclass[c]['ANOVA']
        for dy in (9,-14,20,-25,31):
            ok=all(not(abs(x-px)<thr and abs(dy-pdy)<9 and abs(y-py)<0.13)
                   for px,py,pdy in placed)
            if ok: break
        ax.annotate(SHORT[c],(x,y),xytext=(0,dy),textcoords='offset points',
                    ha='center',fontsize=6.1,color='#48545F')
        placed.append((x,y,dy))
place(axes[0],'n',15); place(axes[1],'silhouette',0.022)
axes[0].set_ylabel('Kuncheva stability index ($k$ = 50)')
axes[0].legend(frameon=False,fontsize=7,loc='lower right')
axes[0].set_title(f"(a) Subtype size   $\\rho$ = {corr['ANOVA']['size_rho']:.2f}, $P$ = {corr['ANOVA']['size_p']:.2f}",fontsize=8,pad=6)
axes[1].set_title(f"(b) Subtype separability   $\\rho$ = {corr['ANOVA']['sil_rho']:.2f}, $P$ = {corr['ANOVA']['sil_p']:.2f}",fontsize=8,pad=6)
fig.tight_layout(); fig.savefig('figs/Figure4.png',dpi=400,bbox_inches='tight')
fig.savefig('figs/Figure4.tif',dpi=400,bbox_inches='tight'); plt.close(fig)

# ---------- Figure 5 ----------
Xa=np.load('Xa_shared.npy'); ya=np.load('ya.npy',allow_pickle=True)
cons=A['consensus']['50']; cnt={}
for m in METH:
    for g in cons[m]: cnt[g]=cnt.get(g,0)+1
top=[g for g,v in sorted(cnt.items(),key=lambda x:(-x[1],x[0])) if v>=3]
mult={g:cnt[g] for g in top}
order=sorted(top,key=lambda g:-cnt[g])
classes=['B-CELL_ALL_ETV6-RUNX1','B-CELL_ALL_TCF3-PBX1','B-CELL_ALL_MLL','B-CELL_ALL_T-ALL',
         'B-CELL_ALL_HYPERDIP','B-CELL_ALL_HYPO','B-CELL_ALL']
Z=(Xa-Xa.mean(0))/ (Xa.std(0)+1e-9)
M=np.array([[Z[ya==c][:,int(g)].mean() for c in classes] for g in order])
fig=plt.figure(figsize=(7.2,4.6))
gs=fig.add_gridspec(1,2,width_ratios=[1,1.18],wspace=.42)
ax=fig.add_subplot(gs[0,0])
for m in METH:
    ax.scatter(A['within'][m]['50']['kuncheva'],PF[m]['bacc_mean'],s=55,facecolor=COL[m],
               edgecolor='white',lw=.9,marker=MK[m],zorder=3)
    ax.errorbar(A['within'][m]['50']['kuncheva'],PF[m]['bacc_mean'],yerr=PF[m]['bacc_sd'],
                color=COL[m],lw=1,capsize=2.2,zorder=2)
    OFF={'ANOVA':(-4,12,'right'),'ReliefF':(6,-4,'left'),'mRMR':(-6,-2,'right'),
         'LASSO':(7,4,'left'),'RF':(-7,4,'right')}
    dx,dy,ha_=OFF[m]
    ax.annotate(LBL[m],(A['within'][m]['50']['kuncheva'],PF[m]['bacc_mean']),xytext=(dx,dy),
                textcoords='offset points',ha=ha_,fontsize=6.8,color=COL[m],fontweight='bold')
ax.set_xlabel('Kuncheva stability index ($k$ = 50)'); ax.set_ylabel('Balanced accuracy, nested CV')
ax.set_xlim(.45,1.0); ax.set_ylim(.68,.90)
ax.grid(color='#EAEAEA',lw=.7,zorder=0); ax.set_axisbelow(True)
for s in ['top','right']: ax.spines[s].set_visible(False)
ax.set_title(f"(a) Stability does not predict accuracy\n$\\rho$ = {r:.2f}, $P$ = {p:.2f}",fontsize=8,pad=6)
ax2=fig.add_subplot(gs[0,1])
div=LinearSegmentedColormap.from_list('d',['#14487f','#7fb3e5','#F5F5F3','#f0a07f','#a8341a'])
lim=np.abs(M).max()
im=ax2.imshow(M,cmap=div,vmin=-lim,vmax=lim,aspect='auto')
ax2.set_xticks(range(len(classes)))
ax2.set_xticklabels([SHORT[c] for c in classes],rotation=45,ha='right',fontsize=6.8)
ax2.set_yticks(range(len(order)))
ax2.set_yticklabels([f"{probes[int(g)]}" for g in order],fontsize=5.6,fontfamily='DejaVu Sans Mono')
for i,g in enumerate(order):
    ax2.text(len(classes)-.34,i,f"{cnt[g]}",fontsize=5.6,va='center',color='#48545F',ha='left')
ax2.text(len(classes)-.34,-0.75,'methods',fontsize=5.8,va='center',ha='left',color='#48545F')
for s in ax2.spines.values(): s.set_visible(False)
ax2.tick_params(length=0)
cb=fig.colorbar(im,ax=ax2,fraction=.035,pad=.17); cb.set_label('Mean expression ($z$)',fontsize=7)
cb.ax.tick_params(labelsize=6.5)
ax2.set_title('(b) Probe sets selected by $\\geq$3 of 5 methods',fontsize=8,pad=6)
fig.savefig('figs/Figure5.png',dpi=400,bbox_inches='tight')
fig.savefig('figs/Figure5.tif',dpi=400,bbox_inches='tight')
print("n consensus probes >=3 methods:",len(order))
print("figs 4,5 done")
