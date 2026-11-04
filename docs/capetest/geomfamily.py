import numpy as np
LAM=0.1; d0=LAM/2; N=24
def centre(p): p=np.asarray(p,float); return p-(p.max()+p.min())/2
def rescale(p,D): p=centre(p); return p*(D/(p.max()-p.min()))
def nested(N):
    N1=N//2;N2=N-N1
    return np.array(sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1) for i in range(N2)])))[:N]
def endload(N,fm,D):
    nm=int(round(N*fm)); ne=(N-nm)//2; nm=N-2*ne
    return np.sort(np.concatenate([-D/2+d0*np.arange(ne), d0*(np.arange(nm)-(nm-1)/2), D/2-d0*np.arange(ne)[::-1]]))
def hier(N,ncore,D,dcore=d0):
    """R1: dense unambiguous core + sparse outriggers."""
    core=dcore*(np.arange(ncore)-(ncore-1)/2)
    no=N-ncore; half=no//2
    out=np.concatenate([np.linspace(-D/2,-D/6,half), np.linspace(D/6,D/2,no-half)])
    return np.sort(np.concatenate([core,out]))
def powerlaw(N,D,p):
    """R2a: symmetric power-law (aperiodic) spacing."""
    M=N//2; t=(np.arange(1,M+1)/M)**p
    h=t*D/2
    return np.sort(np.concatenate([-h,h]))[:N]
def logarr(N,D):
    M=N//2; t=np.logspace(0,1,M); t=(t-t.min())/(t.max()-t.min())
    h=t*D/2+d0
    return np.sort(np.concatenate([-h,h]))[:N]
def golombish(N,D,seed=0):
    """R2b: greedy Sidon-like (near-Golomb) ruler -> aperiodic, min redundancy."""
    rng=np.random.default_rng(seed); best=None;bs=-1
    for trial in range(400):
        pos=[0]; diffs=set()
        cand=list(range(1,400)); rng.shuffle(cand)
        for c in cand:
            if len(pos)>=N: break
            nd={abs(c-q) for q in pos}
            if nd & diffs or len(nd)<len(pos): continue
            pos.append(c); diffs|=nd
        if len(pos)==N:
            sc=-max(pos)
            if sc>bs: bs=sc; best=sorted(pos)
    return rescale(np.array(best,float)*d0,D) if best else None
def psl(pos):
    s=np.linspace(-1,1,40001)
    B=np.abs(np.exp(1j*2*np.pi*np.outer(s,pos)/LAM).sum(1))/len(pos)
    i0=np.argmax(B); i=i0
    while i+1<len(B) and B[i+1]<B[i]: i+=1
    j=i0
    while j-1>=0 and B[j-1]<B[j]: j-=1
    m=B.copy(); m[j:i+1]=0; return m.max()
def crb(pos,r=100.0):
    x=np.asarray(pos,float); A=np.vstack([np.ones_like(x),x,x*x]).T
    Mi=np.linalg.inv(A.T@A)
    sa=np.sqrt(Mi[1,1]/2)*LAM/(2*np.pi)
    sr=np.sqrt(Mi[2,2]/2)*LAM*r*r/np.pi
    return sa,sr
G={}
G["dense ULA D=1.15"]=rescale(np.arange(N)*d0,1.15)
G["nested D=7.75"]=rescale(nested(N)*d0,7.75)
G["end-load D=2.0"]=endload(N,0.5,2.0)
for D in (3.0,6.0):
    G[f"R1 hier core16 D={D}"]=hier(N,16,D)
    G[f"R2 power p=2.0 D={D}"]=powerlaw(N,D,2.0)
    G[f"R2 log D={D}"]=logarr(N,D)
    g=golombish(N,D)
    if g is not None: G[f"R2 golomb D={D}"]=g
print(f"N={N}  (sa=angle CRB in sin units, sr=range CRB at r=100 m)\n")
print(f"{'geometry':26s} {'aper[m]':>8s} {'PSL':>7s} {'sa':>10s} {'sr[m]':>8s}")
for k,p in G.items():
    ap=p.max()-p.min(); sa,sr=crb(p)
    print(f"{k:26s} {ap:8.2f} {psl(p):7.3f} {sa:10.2e} {sr:8.2f}")
np.save("geoms.npy",G,allow_pickle=True)
