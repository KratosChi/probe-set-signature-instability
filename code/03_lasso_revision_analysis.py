"""Reviewer response: LASSO sparsity + cross-validated regularisation tuning."""
import numpy as np, pandas as pd, json, time, warnings, itertools
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold
from sklearn.preprocessing import StandardScaler
log=lambda *a: print(*a, flush=True)
D="/mnt/user-data/uploads/Project Gene (Shilpy mam)/"
PRE=2000; NR=30; SEED=42; KS=[10,25,50,100]
CGRID=[0.003,0.01,0.03,0.1,0.3,1.0,3.0]

def load(f):
    d=pd.read_csv(D+f); y=d['type'].values
    X=d.drop(columns=['samples','type'])
    X=X.drop(columns=[c for c in X.columns if c.startswith('AFFX')])
    return X,y
Xa_df,ya=load("Leukemia_GSE28497.csv"); Xb_df,yb=load("Leukemia_GSE71935.csv")
_sb=set(Xb_df.columns); shared=[c for c in Xa_df.columns if c in _sb]
Xa=Xa_df[shared].values.astype(np.float32); Xb=Xb_df[shared].values.astype(np.float32)
del Xa_df,Xb_df
log("loaded",Xa.shape,Xb.shape)

def coef_rank(Z,yt,C):
    base=LogisticRegression(penalty='l1',solver='liblinear',C=C,max_iter=400)
    if len(np.unique(yt))>2:
        m=OneVsRestClassifier(base).fit(Z,yt)
        W=np.vstack([e.coef_[0] for e in m.estimators_])
        imp=np.abs(W).max(axis=0)
        nz_per=[int((np.abs(e.coef_[0])>0).sum()) for e in m.estimators_]
        nz_union=int((imp>0).sum())
    else:
        m=base.fit(Z,yt); imp=np.abs(m.coef_[0])
        nz_per=[int((imp>0).sum())]; nz_union=int((imp>0).sum())
    return imp,nz_union,nz_per

def tune_C(Z,yt,seed):
    """inner 5-fold CV over the C grid, balanced accuracy."""
    from sklearn.metrics import balanced_accuracy_score
    skf=StratifiedKFold(n_splits=5,shuffle=True,random_state=seed)
    best=(None,-1)
    for C in CGRID:
        sc=[]
        for tr,te in skf.split(Z,yt):
            base=LogisticRegression(penalty='l1',solver='liblinear',C=C,max_iter=400)
            mdl=OneVsRestClassifier(base) if len(np.unique(yt))>2 else base
            mdl.fit(Z[tr],yt[tr])
            sc.append(balanced_accuracy_score(yt[te],mdl.predict(Z[te])))
        s=float(np.mean(sc))
        if s>best[1]: best=(C,s)
    return best[0]

def kuncheva(s1,s2,k,P):
    r=len(s1&s2); e=k*k/P
    return (r-e)/(k-e)

def run(X,y,name):
    P=X.shape[1]
    sss=StratifiedShuffleSplit(n_splits=NR,train_size=0.8,random_state=SEED)
    fixed=[]; tuned=[]; nz_fix=[]; nz_tun=[]; Cs=[]; ties=[]
    t0=time.time()
    for i,(tr,_) in enumerate(sss.split(X,y)):
        Xt,yt=X[tr],y[tr]
        pre=np.argsort(-Xt.var(axis=0))[:PRE]
        Z=StandardScaler().fit_transform(Xt[:,pre])
        imp,nzu,nzp=coef_rank(Z,yt,0.1)
        fixed.append(pre[np.argsort(-imp)]); nz_fix.append(nzu)
        ties.append({str(k):int(max(0,k-nzu)) for k in KS})
        C=tune_C(Z,yt,SEED+i); Cs.append(C)
        imp2,nzu2,_=coef_rank(Z,yt,C)
        tuned.append(pre[np.argsort(-imp2)]); nz_tun.append(nzu2)
        log(f"  [{name}] {i+1}/{NR} nz(C=0.1)={nzu} C*={C} nz(C*)={nzu2}  {time.time()-t0:.0f}s")
    out={'name':name,'nz_fixed':nz_fix,'nz_tuned':nz_tun,'C_selected':Cs,'ties':ties}
    for lab,store in [('fixed',fixed),('tuned',tuned)]:
        out[lab]={}
        for k in KS:
            S=[set(r[:k]) for r in store]
            ku=[kuncheva(S[a],S[b],k,P) for a,b in itertools.combinations(range(len(S)),2)]
            ja=[len(S[a]&S[b])/len(S[a]|S[b]) for a,b in itertools.combinations(range(len(S)),2)]
            out[lab][str(k)]={'kuncheva':float(np.mean(ku)),'kuncheva_sd':float(np.std(ku)),
                              'jaccard':float(np.mean(ja))}
    return out

RA=run(Xa,ya,'GSE28497'); json.dump(RA,open('lasso_rev_A.json','w'),indent=1)
RB=run(Xb,yb,'GSE71935'); json.dump(RB,open('lasso_rev_B.json','w'),indent=1)
log("DONE")
