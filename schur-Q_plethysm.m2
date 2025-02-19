-- computes plethysm of Schur's Q-functions
-- converts between q_n and p_n bases
-- computes plethysm recurrence from [GJ24b]

restart
n = 30 --at most ~30
R = QQ[q_(-n)..q_n]--/(ideal(join({q_0-1},for k from -n to -1 list q_k)))
S = QQ[p_(-n)..p_n]/(ideal(join({p_0-1},for k from -n to -1 list p_k)))

--mod out by relations from 8.2' [Macdonald p.251]
R = R/ideal(join(for m from 1 to n//2 list (-q_(2*m)+(1/2)*(-1)^(m-1)*q_m^2 + sum for r from 1 to (m-1) list (-1)^(r-1)*q_r*q_(2*m-r)),{q_0-1},for i from (-n) to (-1) list q_i))

---------- conversions between p_i's and q_i's

--list of odd partitions of numbers 0..m
oddPartitionList = (
    ans := {{{}}};
    oddNums := for i from 0 to (n-1)//2 list (2*i+1);
    
    for k from 1 to n do (
        currOddNums := for i from 0 to (k-1)//2 list oddNums#i;
        currPartitions := {};
        
        for i from 0 to #currOddNums-1 do (
            --for j from 0 to k-currOddNums#i do (
                for l from 0 to #ans#(k-currOddNums#i)-1 do (
                    currPartitions = append(currPartitions,append(ans#(k-currOddNums#i)#l,currOddNums#i));
                    );
                --);
            );
        currPartitions = for theList in currPartitions list rsort theList;
        currPartitions = toList set currPartitions;
        
        ans = append(ans,currPartitions);
        );
    
    Bag ans
    )

--computes p_lam
plam = lam -> (
    product for i in lam list p_i
    )

--computes z_lam(-1) [Macdonald p.255, p.24]
zlam = lam -> (
    product for i in unique lam list (
        mi := number(lam,k -> k == i);
        i^mi*mi!
        )
    )

--list used for writing q_m in p_k basis [Macdonald p.260]
qTOpList = Bag for m from 0 to n list (sum for lam in oddPartitionList#m list (zlam lam)^(-1)*2^(#lam)*(plam lam))

--maps q_i -> p_j's
qTOpMap = map(S,R,toList(n:0)|toList(qTOpList));

--list used for writing p_m in q_k basis [Macdonald p.260]
pTOqList = (
    ans := {1};
    
    for m from 1 to n do (
        pm := 0;
        if odd m then (
            pm = (zlam({m})/2)*(q_m - sum for lam in delete({m},oddPartitionList#m) list 
                ((zlam lam)^(-1)*2^(#lam)*(product for thePart in lam list ans#thePart)))
            );
        ans = append(ans,pm);
        );
    
    Bag ans
    )

--maps p_i -> q_j's
pTOqMap = map(R,S,toList(n:0)|toList(pTOqList));

--maps any f -> q_i's
TOq = f -> (
    if (ring f) === (ring q_1) then return(f);
    if (ring f) === ZZ or (ring f) === QQ then return(sub(f,ring q_1));
    if pmEven(f) then print("TOq warning: p_m with m even");
    pTOqMap f
    )

--maps any f -> p_i's
TOp = f -> (
    if (ring f) === (ring p_1) then return(f);
    if (ring f) === ZZ or (ring f) === QQ then return(sub(f,ring p_1));
    qTOpMap f
    )

--returns true if f has a term with p_m, m even, false otherwise
pmEven = f -> (
    if (ring f) === (ring q_1) then return(false);
    
    theList := listForm f;
    
    listToLam := termList -> (
        newAns := {};
        for i from 0 to #termList-1 do (
            for j from 1 to termList#i do (
                newAns = append(newAns,i-n);
                );
            );
        newAns
        );
    
    for theTerm in theList do (
        if any(listToLam(theTerm#0),even) then return(true);
        );
    
    false
    )

--if polyn=2p_1^2p_5-p_2^3, returns {({5,1,1},2),({2,2,2},-1)}
fToLamList = polyn -> (
    if polyn == 0 then return({({0},0)});

    ans := {};
    theList := listForm polyn;

    listToLam := termList -> (
        newAns := {};
        for i from 0 to #termList-1 do (
            for j from 1 to termList#i do (
                newAns = append(newAns,i-n);
                );
            );
        newAns
        );

    for theTerm in theList do (
        ans = append(ans,(rsort listToLam(theTerm#0),theTerm#1));
        );
    
    ans
    )

--computes plethysm f(g) [Loehr-Remmel]
pleth = (f,g) -> (
    fInS := ((ring f) === (ring p_1));
    gInS := ((ring g) === (ring p_1));
    if (fInS xor gInS) and (pmEven(f) or pmEven(g)) then print("warning: p_m with m even");
    
    fList := fToLamList TOp f;
    gList := fToLamList TOp g;
    
    fMaxWeight := max for fTerm in fList list sum fTerm#0;
    gMaxWeight := max for gTerm in gList list sum gTerm#0;
    
    if gMaxWeight*fMaxWeight > n then print("plethysm warning: weight might be too high");
    
    ans := sum for fTerm in fList list ((fTerm#1)*(
            product for i from 0 to #fTerm#0-1 list (
                sum for gTerm in gList list ((gTerm#1)*(
                        product for j from 0 to #gTerm#0-1 list(
                            p_((fTerm#0#i)*(gTerm#0#j))))))));
    
    if fInS or gInS then return(ans);
    TOq ans
    )

---------- Pfaffian definition of Schur's Q-functions

--appends 0's to the end of lam
appendZeros = (lam,num) -> (
    for i from 0 to #lam+num-1 list (
        if i < #lam then (lam#i) else (0)
        )
    )

--computes Q_(r,s)
Qtwo = (r,s) -> (
    q_r*q_s + 2*(sum for i from 1 to s list (-1)^i*q_(r+i)*q_(s-i))
    )

--computes M(lam), where Q_lam=Pf(M(lam))
M = lam -> (
    theLam := lam;
    if #lam % 2 == 1 then (
        theLam = appendZeros(lam,1);
        );
    
    map(R^(# theLam),R^(# theLam),(i,j)->
        (if i > j then (
                -Qtwo(theLam_j,theLam_i)
                ) else if i < j then (
                Qtwo(theLam_i,theLam_j)
                ) else (0)))
    )

--computes the Pfaffian of the matrix
pfaff = mat -> (
    if numRows mat == 0 then (
        1
        ) else (
        sum for j from 1 to (numRows mat -1) list (-1)^(j+1)*mat_(0,j)*pfaff(submatrix'(mat,{0,j},{0,j}))
        )
    )

--computes Q_lam
Q = lam -> (pfaff M lam)

--computes M(lam,mu), where Q_{lam/mu}=Pf(M(lam,mu))
skewM = (lam,mu) -> (
    theLam := lam;
    if (#lam + #mu) % 2 == 1 then (
        theLam = appendZeros(lam,1);
        );
    
    Mlam = map(R^(# theLam),R^(# theLam),(i,j)->
        (if i > j then (
                -Qtwo(theLam_j,theLam_i)
                ) else if i < j then (
                Qtwo(theLam_i,theLam_j)
                ) else (0)));
    Nmu = map(R^(# theLam),R^(# mu),(i,j)->q_(theLam_i-mu_((# mu) - j - 1)));
    
    (Mlam|Nmu)||((-transpose(Nmu))|map(R^(# mu),R^(# mu),(i,j)->q_(-1)))
    )

--computes Q_{lam/mu}
skewQ = (lam,mu) -> (
    pfaff skewM(lam,mu)
    )

---------- plethysm recurrence formulas

--computes A.B
dotProd = (A,B) -> (
    if #A != #B then (print("dot product error")) else (
        sum for k from 0 to #A-1 list ((A#k)*(B#k))
        )
    )

--adds 1 to last entry of seq, and rolls over if result > theMax
seqAdd = (seq,theMax) -> (
    ans = new MutableList from seq;
    for antiIndx from 1 to #seq do (
        ans#(#ans-antiIndx) = ans#(#ans-antiIndx)+1;
        if ans#(#ans-antiIndx) > theMax then (
            ans#(#ans-antiIndx) = 0;
            if antiIndx == #seq then (break (toSequence for i from 1 to #seq list 0))
            ) else (
            break toSequence ans
            )
        )
    )

--lists all indices J needed for BK formula
Jlist = (k,n,m) -> (
    -- J.d-k >= 0
    -- n-|J| >= 0
    ans := {};
    del := (1..m);
    ones := toSequence for i from 0 to m-1 list 1;
    curr := toSequence for i from 0 to m-1 list 0;
    for i from 1 to (n+1)^m-1 list (
        curr = seqAdd(curr,n);
        if (dotProd(curr,del) < k) or (dotProd(curr,ones) > n) then continue;
        curr
        )
    )

--computes D_k(q_n(q_m))
BK = (k,n,m) -> (
    indexList := Jlist(k,n,m);
    del := (1..m);
    ones := toSequence for i from 0 to m-1 list 1;
    
    sum for J in indexList list (
        (-1)^(dotProd(J,del)-k)*q_(dotProd(J,del)-k)*(pleth(q_(n-dotProd(J,ones)),q_m))*
        product for s from 1 to m list pleth(q_(J#(s-1)),q_(m-s))
        )
    )
