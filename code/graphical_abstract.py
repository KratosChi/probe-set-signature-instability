import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch
plt.rcParams.update({'font.family':'DejaVu Sans','figure.facecolor':'white','savefig.facecolor':'white'})

A=json.load(open('results_A.json'))['A']
METH=['ANOVA','mRMR','ReliefF','LASSO','RF']
NAME={'ANOVA':'ANOVA $F$-test','mRMR':'mRMR','ReliefF':'ReliefF','LASSO':'LASSO','RF':'Random forest'}
SHORT={'ANOVA':'ANOVA','mRMR':'mRMR','ReliefF':'ReliefF','LASSO':'LASSO','RF':'RF'}
FAM={'ANOVA':'univariate filter','mRMR':'multivariate filter','ReliefF':'instance-based filter',
     'LASSO':'embedded, $L$1 penalty','RF':'embedded, ensemble'}
COL={'ANOVA':'#2a78d6','mRMR':'#eb6834','ReliefF':'#1baf7a','LASSO':'#eda100','RF':'#e87ba4'}
cons={m:set(A['consensus']['50'][m]) for m in METH}
union=sorted(set().union(*cons.values()))
mult={g:sum(g in cons[m] for m in METH) for g in union}
order=sorted(union,key=lambda g:(-mult[g],g)); NU=len(union)
n1=sum(1 for g in union if mult[g]==1); n5=sum(1 for g in union if mult[g]==5)
nsh=NU-n1

PXW,PXH=2656,1062; DPI=200
fig=plt.figure(figsize=(PXW/DPI,PXH/DPI),dpi=DPI)
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,PXW); ax.set_ylim(0,PXH); ax.axis('off')
INK='#151C24'; SOFT='#4A5560'; FAINT='#7A8792'; RED='#A94C36'; NAVY='#1B4C7E'
HDR='#26384A'; LINE='#C3CDD6'
boxes=[]
def T(x,y,s,size=11.5,color=SOFT,bold=False,ha='left',va='center',zone=None,ls=1.4,italic=False):
    t=ax.text(x,y,s,fontsize=size,color=color,ha=ha,va=va,linespacing=ls,
              fontweight='bold' if bold else 'normal',style='italic' if italic else 'normal')
    boxes.append((t,zone,s)); return t
def panel(x0,x1,y0,y1,num,title,zone):
    ax.add_patch(FancyBboxPatch((x0,y0),x1-x0,y1-y0,boxstyle='round,pad=0,rounding_size=10',
        fc='white',ec=LINE,lw=2.0,zorder=1))
    ax.add_patch(FancyBboxPatch((x0,y1-54),x1-x0,54,boxstyle='round,pad=0,rounding_size=10',
        fc=HDR,ec=HDR,lw=0,zorder=2))
    ax.add_patch(Rectangle((x0,y1-54),x1-x0,18,fc=HDR,ec='none',zorder=2))
    T(x0+22,y1-27,f'{num}',13.5,'#8FB6D9',bold=True,zone=zone).set_zorder(3)
    T(x0+60,y1-27,title,13,'white',bold=True,zone=zone).set_zorder(3)
def sub(x0,x1,y0,y1,fc='#F4F7FA',ec='#D6DEE6',lw=1.4,r=7):
    ax.add_patch(FancyBboxPatch((x0,y0),x1-x0,y1-y0,boxstyle=f'round,pad=0,rounding_size={r}',
        fc=fc,ec=ec,lw=lw,zorder=2))
def down(x,y0,y1):
    ax.add_patch(FancyArrowPatch((x,y0),(x,y1),arrowstyle='-|>',mutation_scale=13,lw=1.7,
        color='#8896A3',zorder=3))
def right(x0,x1,y):
    ax.add_patch(FancyArrowPatch((x0,y),(x1,y),arrowstyle='-|>',mutation_scale=26,lw=4.2,
        color='#5E7183',zorder=4))

ax.add_patch(FancyBboxPatch((14,14),PXW-28,PXH-28,boxstyle='round,pad=0,rounding_size=14',
    fc='none',ec='#AEBAC5',lw=2.4,zorder=0))
PY0,PY1=64,996

# ============ PANEL 1 ============
panel(52,592,PY0,PY1,'1','INPUT DATA','P1')
sub(74,570,822,930)
T(94,900,'GSE28497',13,INK,bold=True,zone='P1')
T(94,868,'Pediatric B-lineage ALL',11,SOFT,zone='P1')
T(94,842,'281 samples  ·  7 subtypes',11,SOFT,zone='P1')
sub(74,570,690,798)
T(94,768,'GSE71935',13,INK,bold=True,zone='P1')
T(94,736,'Juvenile myelomonocytic leukemia',9.5,SOFT,zone='P1')
T(94,710,'46 samples  ·  JMML vs normal',11,SOFT,zone='P1')
down(322,682,646)
sub(74,570,532,636,fc='#EAF1F8',ec='#BFD4E8')
T(94,606,'22,215 shared probe sets',11.5,INK,bold=True,zone='P1')
T(94,576,'GPL96 ∩ GPL570  ·  AFFX removed',10,SOFT,zone='P1')
T(94,552,'log$_2$ verified  ·  no missing values',10,SOFT,zone='P1')
down(322,524,488)
sub(74,570,374,478)
T(94,448,'30 stratified subsamples',12.5,INK,bold=True,zone='P1')
T(94,418,'80% of each cohort, balance kept',10,SOFT,zone='P1')
T(94,394,'variance prefilter → 2,000 probes',10,SOFT,zone='P1')
T(76,296,'Curated by CuMiDa from NCBI GEO,',10,FAINT,zone='P1')
T(76,268,'both series in the public domain.',10,FAINT,zone='P1')
right(604,644,530)

# ============ PANEL 2 ============
panel(652,1122,PY0,PY1,'2','FIVE ALGORITHMS','P2')
T(674,900,'Identical splits and prefilter;',9.5,FAINT,zone='P2')
T(674,874,'top 50 probes from each fit.',9.5,FAINT,zone='P2')
for i,m in enumerate(METH):
    yb=754-i*112
    sub(674,1100,yb,yb+92,fc='white',ec=COL[m],lw=1.8)
    ax.add_patch(Rectangle((674,yb),13,92,fc=COL[m],ec='none',zorder=3))
    T(706,yb+58,NAME[m],13,INK,bold=True,zone='P2')
    T(706,yb+30,FAM[m],10,SOFT,zone='P2')
T(674,150,'Signature sizes $k$ = 10–100',10,FAINT,zone='P2')
right(1134,1174,530)

# ============ PANEL 3 ============
panel(1182,1922,PY0,PY1,'3','DO THEY AGREE?','P3')
T(1204,900,'Union of the five signatures  ·  142 probe sets',10.5,SOFT,zone='P3')
MX0,MX1,MYT,MYB=1324,1900,860,644
cw=(MX1-MX0)/NU; rh=(MYT-MYB)/5
seq=['#DCE8F5','#A5C7E9','#6BA3DB','#3A7CC8','#14487F']
for r,m in enumerate(METH):
    y=MYT-(r+1)*rh
    T(MX0-14,y+rh/2,SHORT[m],10.5,SOFT,ha='right',zone='P3')
    for c,g in enumerate(order):
        fc=seq[mult[g]-1] if g in cons[m] else '#EFF2F4'
        ax.add_patch(Rectangle((MX0+c*cw,y+2),max(cw*0.88,1.5),rh-4,fc=fc,ec='none',zorder=3))
ax.add_patch(Rectangle((MX0,MYB),MX1-MX0,MYT-MYB,fc='none',ec='#CBD4DC',lw=1.2,zorder=4))
xb=MX0+nsh*cw
ax.plot([MX0,xb-3],[MYB-16,MYB-16],color='#98A5B1',lw=1.6,zorder=3)
ax.plot([xb+3,MX1],[MYB-16,MYB-16],color='#98A5B1',lw=1.6,zorder=3)
for xx in (MX0,xb-3,xb+3,MX1): ax.plot([xx,xx],[MYB-16,MYB-8],color='#98A5B1',lw=1.6,zorder=3)
T((MX0+xb)/2,MYB-40,'shared by 2–5',10.5,SOFT,ha='center',zone='P3')
T((xb+MX1)/2,MYB-40,'ONE method only',10.5,RED,ha='center',bold=True,zone='P3')
sub(1204,1900,398,556,fc='#FBF2EF',ec='#E3C3B8')
T(1226,477,f'{100*n1/NU:.0f}%',30,RED,bold=True,zone='P3')
T(1452,494,'of the pooled signature',11,INK,zone='P3')
T(1452,462,'came from ONE algorithm',11,INK,zone='P3')
sub(1204,1900,228,376,fc='#EEF3F9',ec='#BFD4E8')
T(1226,302,f'{100*n5/NU:.1f}%',30,NAVY,bold=True,zone='P3')
T(1478,318,'was chosen by all five',11,INK,zone='P3')
T(1478,286,'algorithms',11,INK,zone='P3')
for i in range(5):
    ax.add_patch(Rectangle((1204+i*38,140),30,22,fc=seq[i],ec='none',zorder=3))
T(1408,151,'selected by 1 … 5 of 5 methods',10,FAINT,zone='P3')
right(1934,1974,530)

# ============ PANEL 4 ============
panel(1982,2604,PY0,PY1,'4','WHAT IT MEANS','P4')
def spread(y,lo,hi,lab,val,color):
    x0,x1=2004,2360; rng=(0.5,1.0)
    ax.plot([x0,x1],[y,y],color='#E4E9ED',lw=13,solid_capstyle='round',zorder=2)
    a=x0+(lo-rng[0])/(rng[1]-rng[0])*(x1-x0); b=x0+(hi-rng[0])/(rng[1]-rng[0])*(x1-x0)
    ax.plot([a,b],[y,y],color=color,lw=13,solid_capstyle='round',zorder=3)
    T(2004,y+46,lab,11,INK,bold=True,zone='P4')
    T(2004,y-38,f'{lo:.3f} – {hi:.3f}',10.5,FAINT,zone='P4')
    T(2382,y,val,17,color,bold=True,zone='P4')
T(2004,900,'Range across the five algorithms,',10,FAINT,zone='P4')
T(2004,874,'plotted on one common scale',10,FAINT,zone='P4')
spread(778,0.558,0.898,'Signature stability (Kuncheva)','0.341',NAVY)
spread(618,0.775,0.798,'Balanced accuracy, nested CV','0.024',RED)
T(2004,508,'No pairwise accuracy difference',10,FAINT,zone='P4')
T(2004,482,'was significant (all $P$ > 0.2).',10,FAINT,zone='P4')
sub(2004,2582,142,414,fc='#F2F5F8',ec='#C6D0D9',lw=1.6,r=10)
T(2026,364,'Stability varies 14× more',13,INK,bold=True,zone='P4')
T(2026,330,'than accuracy',13,INK,bold=True,zone='P4')
T(2026,232,'Signature composition reflects\nthe choice of algorithm as much\nas the biology — and accuracy\ncannot warn you.',11,SOFT,zone='P4',ls=1.5)

Z={'P1':(52,592),'P2':(652,1122),'P3':(1182,1922),'P4':(1982,2604)}
fig.canvas.draw(); r=fig.canvas.get_renderer(); inv=ax.transData.inverted()
bad=0
for t,z,s in boxes:
    if not z: continue
    bb=t.get_window_extent(renderer=r)
    (x0,y0)=inv.transform((bb.x0,bb.y0)); (x1,y1)=inv.transform((bb.x1,bb.y1))
    zx0,zx1=Z[z]
    if x0<zx0+8 or x1>zx1-8 or y0<PY0+4 or y1>PY1-4:
        bad+=1; print(f"OVERFLOW[{z}] x{x0:.0f}-{x1:.0f} zone{zx0}-{zx1} y{y0:.0f}-{y1:.0f} :: {s[:44]!r}")
print("overflows:",bad)
fig.savefig('figs/Graphical_Abstract.png',dpi=DPI,pad_inches=0)
fig.savefig('figs/Graphical_Abstract.tif',dpi=DPI,pad_inches=0)
