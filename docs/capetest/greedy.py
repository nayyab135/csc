import numpy as np
def runlen(X):
    X=sorted(set(int(v) for v in X)); best=cur=1
    for i in range(1,len(X)):
        cur=cur+1 if X[i]==X[i-1]+1 else 1; best=max(best,cur)
    return best
def udof(half):
    n=sorted(set([-v for v in half]+list(half)))
    return runlen({a-b for a in n for b in n}), n
def greedy(N0,bound_mult=5):
    ncore=2*max(1,N0//3); ncore=min(ncore,max(2,N0-2)); ncore-=ncore%2
    c=ncore//2; nout=max(1,(N0-ncore)//2)
    half=list(range(1,c+1)); bound=bound_mult*N0
    for k in range(nout):
        best=None;bs=-1
        for p in range(max(half)+1,bound+1):
            r,_=udof(half+[p])
            if r>bs: bs=r; best=p
        half.append(best)
    r,n=udof(half)
    return half,r,n
for N0 in (24,8,16,32):
    half,r,n=greedy(N0)
    ap=(max(n)-min(n))*0.05
    print(f"N0={N0:3d}  core=1..{N0//3 if True else 0}  outriggers={half[len([h for h in half if h<=max(1,N0//3)]):]}")
    print(f"        half={half}  N={len(n)}  uDOFdiff={r}  aperture={ap:.2f} m")
