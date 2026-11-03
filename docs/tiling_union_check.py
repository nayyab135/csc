def dca(n): return {a-b for a in n for b in n}
def runlen(S):
    S=sorted(S); best=cur=1
    for i in range(1,len(S)):
        cur=cur+1 if S[i]==S[i-1]+1 else 1; best=max(best,cur)
    return best
def stretched_nested(N0,s):
    N1=N0//2; N2=N0-N1
    return sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1)*s for i in range(N2)]))[:N0]
LAM=0.1
for N0 in (24,8):
    print(f"=== N0={N0}/AP, L=8 ===")
    base=stretched_nested(N0,1); db=dca(base)
    print(f" identical nested: DOF={len(db)} contig={runlen(db)} aper={(max(base)-min(base))*0.05:.2f}m dR={2*((max(base)-min(base))*0.05)**2/LAM:.0f}m")
    U=set(); print(" heterogeneous (stretch factors 1,1,2,2,3,3,4,4):")
    for l,s in enumerate([1,1,2,2,3,3,4,4],1):
        n=stretched_nested(N0,s); d=dca(n); U|=d
        ap=(max(n)-min(n))*0.05
        print(f"   AP{l} s={s}: aper={ap:6.2f}m dR={2*ap*ap/LAM:7.0f}m  DOF={len(d):4d} contig={runlen(d):4d}")
    print(f" >> union DOF={len(U)} contig={runlen(U)} virtual-aperture=+-{max(U)} (vs +-{max(db)})")
    print(f"    total-DOF gain {len(U)/len(db):.2f}x ; contiguous gain {runlen(U)/runlen(db):.2f}x\n")
