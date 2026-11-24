"""Full stability analysis: 5 feature selection algorithms x 2 CuMiDa leukemia cohorts."""
import numpy as np, pandas as pd, json, time, sys, warnings
warnings.filterwarnings("ignore")
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import StratifiedShuffleSplit, RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.metrics import balanced_accuracy_score, f1_score
from skrebate import ReliefF
from mrmr import mrmr_classif

D="/mnt/user-data/uploads/Project Gene (Shilpy mam)/"
PREFILTER=2000; N_RESAMP=30; KS=[10,25,50,100]; SEED=42
log=lambda *a: (print(*a, flush=True))

def load(f):
    d=pd.read_csv(D+f)
    y=d['type'].values
    X=d.drop(columns=['samples','type'])
    X=X.drop(columns=[c for c in X.columns if c.startswith('AFFX')])
    return X, y

Xa,ya=load("Leukemia_GSE28497.csv")
Xb,yb=load("Leukemia_GSE71935.csv")
_sb=set(Xb.columns)
shared=[c for c in Xa.columns if c in _sb]
log(f"shared probes: {len(shared)}")
Xa=Xa[shared].values.astype(np.float32); Xb=Xb[shared].values.astype(np.float32)
probes=np.array(shared)
log("A",Xa.shape,"B",Xb.shape)

def rank_all(Xt,yt):
    """Return dict method -> ranked feature indices (into the prefiltered space) + the prefilter index."""
    v=Xt.var(axis=0); pre=np.argsort(-v)[:PREFILTER]
    Xp=Xt[:,pre]
    out={}
    F,_=f_classif(Xp,yt); F=np.nan_to_num(F)
    out['ANOVA']=np.argsort(-F)
    Z=StandardScaler().fit_transform(Xp)
    base=LogisticRegression(penalty='l1',solver='liblinear',C=0.1,max_iter=400)
    if len(np.unique(yt))>2:
        m=OneVsRestClassifier(base).fit(Z,yt)
        imp=np.abs(np.vstack([e.coef_[0] for e in m.estimators_])).max(axis=0)
    else:
        imp=np.abs(base.fit(Z,yt).coef_[0])
    out['LASSO']=np.argsort(-imp)
    rf=RandomForestClassifier(n_estimators=500,n_jobs=-1,random_state=1).fit(Xp,yt)
    out['RF']=np.argsort(-rf.feature_importances_)
    rl=ReliefF(n_neighbors=10,n_jobs=-1).fit(Xp,yt)
    out['ReliefF']=np.argsort(-rl.feature_importances_)
    sel=mrmr_classif(X=pd.DataFrame(Xp),y=pd.Series(yt),K=max(KS),show_progress=False)
    rest=[i for i in range(Xp.shape[1]) if i not in set(sel)]
    out['mRMR']=np.array(list(sel)+rest)
    return {m: pre[r] for m,r in out.items()}   # map back to global probe index

def kuncheva(s1,s2,k,P):
    r=len(s1&s2); e=k*k/P
    return (r-e)/(k-e) if k!=e else 0.0

def run_cohort(X,y,name,n_resamp=N_RESAMP,methods=None,seed=SEED):
    P=X.shape[1]
    sss=StratifiedShuffleSplit(n_splits=n_resamp,train_size=0.8,random_state=seed)
    store={m:[] for m in ['ANOVA','LASSO','RF','ReliefF','mRMR']}
    t0=time.time()
    for i,(tr,_) in enumerate(sss.split(X,y)):
        rk=rank_all(X[tr],y[tr])
        for m,r in rk.items(): store[m].append(r)
        log(f"  [{name}] resample {i+1}/{n_resamp}  {time.time()-t0:.0f}s")
    res={'within':{}, 'cross':{}, 'consensus':{}}
    for m,R in store.items():
        res['within'][m]={}
        for k in KS:
            S=[set(r[:k]) for r in R]
            ku=[];ja=[]
            for i in range(len(S)):
                for j in range(i+1,len(S)):
                    inter=len(S[i]&S[j]); uni=len(S[i]|S[j])
                    ja.append(inter/uni); ku.append(kuncheva(S[i],S[j],k,P))
            res['within'][m][k]={'kuncheva':float(np.mean(ku)),'kuncheva_sd':float(np.std(ku)),
                                 'jaccard':float(np.mean(ja))}
    for k in KS:
        cons={}
        for m,R in store.items():
            c={}
            for r in R:
                for g in r[:k]: c[g]=c.get(g,0)+1
            cons[m]=set(sorted(c,key=lambda x:-c[x])[:k])
            res['consensus'].setdefault(str(k),{})[m]=[int(x) for x in sorted(cons[m])]
        ms=list(cons); cr={}
        for i in range(len(ms)):
            for j in range(i+1,len(ms)):
                A,B=cons[ms[i]],cons[ms[j]]
                cr[f"{ms[i]}|{ms[j]}"]=len(A&B)/len(A|B)
        res['cross'][str(k)]=cr
        allsel={}
        for m in ms:
            for g in cons[m]: allsel[g]=allsel.get(g,0)+1
        res.setdefault('multiplicity',{})[str(k)]={str(n):int(sum(1 for v in allsel.values() if v==n)) for n in range(1,6)}
    res['probe_freq_k50']={}
    for m,R in store.items():
        c={}
        for r in R:
            for g in r[:50]: c[g]=c.get(g,0)+1
        res['probe_freq_k50'][m]={str(int(g)):int(v) for g,v in sorted(c.items(),key=lambda x:-x[1])[:80]}
    return res

log("=== COHORT A: GSE28497 ===")
RA=run_cohort(Xa,ya,"GSE28497")
json.dump({'probes':list(map(str,probes)),'A':RA}, open('results_A.json','w'))
log("=== COHORT B: GSE71935 ===")
RB=run_cohort(Xb,yb,"GSE71935")
json.dump({'B':RB}, open('results_B.json','w'))
log("DONE stability")
