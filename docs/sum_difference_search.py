import random
random.seed(7)
def runlen(X):
    X=sorted(set(X)); best=cur=1
    for i in range(1,len(X)):
        cur=cur+1 if X[i]==X[i-1]+1 else 1; best=max(best,cur)
    return best
def ud(S): return runlen({a-b for a in S for b in S})
def us(S): return runlen({a+b for a in S for b in S})
def nested(N0):
    N1=N0//2; N2=N0-N1
    return sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1) for i in range(N2)]))[:N0]
def anneal(N0,span,iters=160000,wd=1.0,ws=1.0):
    cur=sorted(random.sample(range(0,span+1),N0))
    if cur[0]!=0: cur[0]=0; cur=sorted(set(cur))
    while len(cur)<N0: cur=sorted(set(cur+[random.randint(0,span)]))
    def score(S): return wd*ud(S)+ws*us(S)
    cs=score(cur); best=cur[:]; bs=cs
    T=6.0
    for it in range(iters):
        T=max(0.02,6.0*(1-it/iters))
        cand=cur[:]
        i=random.randrange(1,N0)
        cand[i]=random.randint(0,span)
        cand=sorted(set(cand))
        if len(cand)!=N0: continue
        s2=score(cand)
        if s2>=cs or random.random()<pow(2.718,(s2-cs)/T):
            cur,cs=cand,s2
            if s2>bs: best,bs=cand[:],s2
    return best
for N0,span in ((8,24),(24,170)):
    nb=nested(N0)
    print(f"=== N0={N0} ===  nested: uDOFdiff={ud(nb)} uDOFsum={us(nb)}")
    for wd,ws,tag in ((1,1,'balanced'),(1,2,'sum-weighted'),(2,1,'diff-weighted')):
        B=anneal(N0,span,wd=wd,ws=ws)
        print(f"  {tag:13s}: uDOFdiff={ud(B):3d} uDOFsum={us(B):3d}  S={B}")
    print()
