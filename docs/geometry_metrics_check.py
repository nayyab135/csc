LAM=0.1; U=LAM/4   # positions are integers in units of lambda/4
def run_len(S):
    S=sorted(set(S))
    if not S: return 0
    best=cur=1
    for i in range(1,len(S)):
        cur = cur+1 if S[i]==S[i-1]+1 else 1
        best=max(best,cur)
    return best
def metrics(pos,name):
    pos=sorted(set(pos)); N=len(pos)
    D=[a-b for a in pos for b in pos]
    S=[a+b for a in pos for b in pos]
    # standard uDOF on the lambda/2 lattice: even lags, rescaled
    ud=run_len([x//2 for x in D if x%2==0])
    us=run_len([x//2 for x in S if x%2==0])
    clean=sorted({2*x for x in pos if -x in pos})       # clean angle lags (mirror pairs)
    uc=run_len(clean)                                    # on lambda/4 grid
    ap=(max(pos)-min(pos))*U; dR=2*ap*ap/LAM
    w={k:sum(1 for a in pos for b in pos if a-b==2*k)//1 for k in (1,2,3)}  # k*lambda/2
    print(f"{name:28s} N={N:3d} ap={ap:6.2f}m dR={dR:7.0f}m uDOFdiff={ud:4d} uDOFsum={us:4d} cleanAng={uc:3d} w1={w[1]:3d} w2={w[2]:3d} w3={w[3]:3d}")

print("positions in lambda/4 units; uDOF measured on the lambda/2 lattice; cleanAng on lambda/4\n")
h=sorted(set(list(range(1,7))+[i*7 for i in range(1,7)]))
metrics([2*x for x in ([-y for y in h]+h)],"CURRENT symnested N=24")
nst=sorted(set(list(range(1,13))+[i*13 for i in range(1,13)]))
metrics([2*x for x in nst],"nested N=24")
metrics([2*i for i in range(24)],"dense ULA N=24")
print()
def tsmisc(nc,wings,quarter=True):
    step = 1 if quarter else 2          # core spacing: lambda/4 or lambda/2
    core=[]
    for k in range(1,nc+1):
        core += [k*step,-k*step]
    w=[]
    for g in wings: w += [2*g,-2*g]     # wings on lambda/2 lattice
    return core+w
for nc,wg,q in [(4,[6,13,21,30,40,51,63,76],True),
                (4,[6,13,21,30,40,51,63,76],False),
                (6,[8,17,27,38,50,63],True),
                (4,[7,15,24,34,45,57,70,84],True),
                (3,[5,11,18,26,35,45,56,68,81],True)]:
    metrics(tsmisc(nc,wg,q),f"TS-MISC core={nc}{'(l/4)' if q else '(l/2)'} wings={len(wg)}")
