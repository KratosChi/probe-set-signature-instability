"""Per-class stability, classification performance, and expression-level outputs."""
import numpy as np, pandas as pd, json, time, warnings
warnings.filterwarnings("ignore")
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import StratifiedShuffleSplit, RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_samples
from scipy.stats import spearmanr
from skrebate import ReliefF
from mrmr import mrmr_classif
log=lambda *a: print(*a, flush=True)
D="/mnt/user-data/uploads/Project Gene (Shilpy mam)/"
PREFILTER=2000; SEED=42

def load(f):
    d=pd.read_csv(D+f); y=d['type'].values
    X=d.drop(columns=['samples','type'])
    X=X.drop(columns=[c for c in X.columns if c.startswith('AFFX')])
    return X,y
Xa_df,ya=load("Leukemia_GSE28497.csv"); Xb_df,yb=load("Leukemia_GSE71935.csv")
_sb=set(Xb_df.columns); shared=[c for c in Xa_df.columns if c in _sb]
Xa=Xa_df[shared].values.astype(np.float32); probes=np.array(shared)
del Xa_df, Xb_df
log("loaded", Xa.shape)

def rank_fast(Xt,yt):
    v=Xt.var(axis=0); pre=np.argsort(-v)[:PREFILTER]; Xp=Xt[:,pre]; out={}
    F,_=f_classif(Xp,yt); out['ANOVA']=pre[np.argsort(-np.nan_to_num(F))]
    Z=StandardScaler().fit_transform(Xp)
    base=LogisticRegression(penalty='l1',solver='liblinear',C=0.1,max_iter=400)
    if len(np.unique(yt))>2:
        m=OneVsRestClassifier(base).fit(Z,yt)
        imp=np.abs(np.vstack([e.coef_[0] for e in m.estimators_])).max(axis=0)
    else: imp=np.abs(base.fit(Z,yt).coef_[0])
    out['LASSO']=pre[np.argsort(-imp)]
    rf=RandomForestClassifier(n_estimators=500,n_jobs=-1,random_state=1).fit(Xp,yt)
    out['RF']=pre[np.argsort(-rf.feature_importances_)]
    return out

def rank_all(Xt,yt):
    out=rank_fast(Xt,yt)
    v=Xt.var(axis=0); pre=np.argsort(-v)[:PREFILTER]; Xp=Xt[:,pre]
    rl=ReliefF(n_neighbors=10,n_jobs=-1).fit(Xp,yt); out['ReliefF']=pre[np.argsort(-rl.feature_importances_)]
    sel=mrmr_classif(X=pd.DataFrame(Xp),y=pd.Series(yt),K=50,show_progress=False)
    rest=[i for i in range(Xp.shape[1]) if i not in set(sel)]
    out['mRMR']=pre[np.array(list(sel)+rest)]
    return out

def kuncheva(s1,s2,k,P):
    r=len(s1&s2); e=k*k/P
    return (r-e)/(k-e) if k!=e else 0.0

# ---- 1. per-class one-vs-rest stability
log("=== per-class stability ===")
P=Xa.shape[1]; K=50; NR=25
perclass={}
for cl in pd.Series(ya).value_counts().index:
    yy=(ya==cl).astype(int); t=time.time()
    sss=StratifiedShuffleSplit(n_splits=NR,train_size=0.8,random_state=SEED)
    store={'ANOVA':[],'LASSO':[],'RF':[]}
    for tr,_ in sss.split(Xa,yy):
        r=rank_fast(Xa[tr],yy[tr])
        for m in store: store[m].append(set(r[m][:K]))
    d={}
    for m,S in store.items():
        ku=[kuncheva(S[i],S[j],K,P) for i in range(len(S)) for j in range(i+1,len(S))]
        d[m]=float(np.mean(ku))
    perclass[cl]={'n':int(yy.sum()),**d}
    log(f"  {cl} n={int(yy.sum())} {d} [{time.time()-t:.0f}s]")

# class separability (silhouette on PCA-50)
Z=StandardScaler().fit_transform(Xa); pcs=PCA(n_components=50,random_state=0).fit_transform(Z)
sil=silhouette_samples(pcs,ya)
for cl in perclass: perclass[cl]['silhouette']=float(np.mean(sil[ya==cl]))
CLS=list(perclass.keys())
ns=[perclass[c]['n'] for c in CLS]; ss=[perclass[c]['silhouette'] for c in CLS]
corr={}
for m in ['ANOVA','LASSO','RF']:
    vs=[perclass[c][m] for c in CLS]
    r1,p1=spearmanr(ns,vs); r2,p2=spearmanr(ss,vs)
    corr[m]={'size_rho':float(r1),'size_p':float(p1),'sil_rho':float(r2),'sil_p':float(p2)}
    log(f"  corr {m}: size rho={r1:.3f} p={p1:.3f} | silhouette rho={r2:.3f} p={p2:.3f}")
json.dump({'perclass':perclass,'corr':corr},open('results_perclass.json','w'),indent=1)

# ---- 2. classification performance with FS inside the fold
log("=== nested classification ===")
rskf=RepeatedStratifiedKFold(n_splits=5,n_repeats=3,random_state=SEED)
perf={m:{'bacc':[],'f1':[]} for m in ['ANOVA','LASSO','RF','ReliefF','mRMR']}
t0=time.time()
for fi,(tr,te) in enumerate(rskf.split(Xa,ya)):
    rk=rank_all(Xa[tr],ya[tr])
    for m,r in rk.items():
        idx=r[:50]
        sc=StandardScaler().fit(Xa[tr][:,idx])
        clf=LinearSVC(C=1.0,max_iter=5000).fit(sc.transform(Xa[tr][:,idx]),ya[tr])
        pred=clf.predict(sc.transform(Xa[te][:,idx]))
        perf[m]['bacc'].append(float(balanced_accuracy_score(ya[te],pred)))
        perf[m]['f1'].append(float(f1_score(ya[te],pred,average='macro')))
    log(f"  fold {fi+1}/15 [{time.time()-t0:.0f}s]")
summ={m:{'bacc_mean':float(np.mean(v['bacc'])),'bacc_sd':float(np.std(v['bacc'])),
         'f1_mean':float(np.mean(v['f1'])),'f1_sd':float(np.std(v['f1'])),
         'bacc_all':v['bacc']} for m,v in perf.items()}
for m,v in summ.items(): log(f"  {m}: bAcc={v['bacc_mean']:.3f}+-{v['bacc_sd']:.3f} macroF1={v['f1_mean']:.3f}")
json.dump(summ,open('results_perf.json','w'),indent=1)

# ---- 3. expression matrix for consensus probes (filled by figure script)
np.save('Xa_shared.npy',Xa); np.save('ya.npy',ya); np.save('probes.npy',probes)
log("DONE analysis2")
