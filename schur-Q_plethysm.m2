-- computes plethysm of Schur's Q-functions and power sums
-- converts between q_n and p_n bases
-- computes plethysm recurrence from [GJ24b]

restart
largestIndex = 30 --at most ~30
R = QQ[q_(-largestIndex)..q_largestIndex]
S = QQ[p_(-largestIndex)..p_largestIndex]/(ideal(join({p_0-1},for k from -largestIndex to -1 list p_k)))

--mod out by relations from 8.2' [Macdonald p.251]
R = R/ideal(join(for m from 1 to largestIndex//2 list (-q_(2*m)+(1/2)*(-1)^(m-1)*q_m^2 + sum for r from 1 to (m-1) list (-1)^(r-1)*q_r*q_(2*m-r)),{q_0-1},for i from (-largestIndex) to (-1) list q_i))

protect symbol q
protect symbol p
protect symbol R
protect symbol S
protect symbol largestIndex

---------- conversions between p_i's and q_i's

--list of odd partitions of numbers 0..m
oddPartitionList = (
    ans := {{{}}};
    oddNums := for i from 0 to (largestIndex-1)//2 list (2*i+1);
    
    for k from 1 to largestIndex do (
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
qTOpList = Bag for m from 0 to largestIndex list (sum for lam in oddPartitionList#m list (zlam lam)^(-1)*2^(#lam)*(plam lam))

--maps q_i -> p_j's
qTOpMap = map(S,R,toList(largestIndex:0)|toList(qTOpList));

--list used for writing p_m in q_k basis [Macdonald p.260]
pTOqList = (
    ans := {1};
    
    for m from 1 to largestIndex do (
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
pTOqMap = map(R,S,toList(largestIndex:0)|toList(pTOqList));

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
                newAns = append(newAns,i-largestIndex);
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
                newAns = append(newAns,i-largestIndex);
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
    
    if gMaxWeight*fMaxWeight > largestIndex then print("plethysm warning: weight might be too high");
    
    ans := sum for fTerm in fList list ((fTerm#1)*(
            product for i from 0 to #fTerm#0-1 list (
                sum for gTerm in gList list ((gTerm#1)*(
                        product for j from 0 to #gTerm#0-1 list(
                            p_((fTerm#0#i)*(gTerm#0#j))))))));
    
    if fInS or gInS then return(ans);
    TOq ans
    )

-- decomposes Q_{lam/mu} as a linear combination of Q_nu
decomposeSkew = {doPrint => true} >> o -> (lam,mu) -> (
    if lam == rsort lam and lam#-1 >= 0 and lam == mu then return(sub(1,S));
    
    ans := {};
    numVars := #lam;
    currPolyn := skewQ(lam,mu,numVars);
    
    while currPolyn != 0 do (
        leadPart := delete(0,(listForm leadTerm currPolyn)#0#0);
        leadCoeff := leadCoefficient currPolyn;
        
        basisFunction := Qlam(leadPart,numVars);
        basisCoeff := leadCoefficient basisFunction;
        
        theCoeff := (sub(leadCoeff,R))//(sub(basisCoeff,R));
        ans = append(ans,(leadPart,theCoeff));
        
        currPolyn = currPolyn - theCoeff * (basisFunction);
        );
    
    ans = rsort ans;
    
    if o.doPrint then (
        -- prints decomposition into latex code
        print("\n");
        print(concatenate("Q_{("|toString(lam)|")/("|toString(mu)|")}&=",decompToTex(ans)));
        
        -- prints decomposition into m2 functions
        print("\n");
        print(decompToM2(ans));
        
        if #lam <= 3 then (
            return(decompToPretty(ans));
            );
        );
    
    ans
    )

-- computes the inner product (F,G)
innerProd = (F,G) -> (
    decompF := decomposeQ(F,doPrint => false);
    decompG := decomposeQ(G,doPrint => false);
    
    ans := 0;
    
    for termF in decompF do (
        theInd := positions(decompG,i -> i#0 == termF#0);
        if #theInd > 0 then (
            ans = ans + (termF#1)*(decompG#(theInd#0)#1)*2^(#(termF#0));
            );
        );
    
    ans
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

--decompose f into linear combination of Q_nu's
decomposeQ = {doPrint => true} >> o -> f -> (
    fList := fToLamList TOq f;
    
    ans := {};
    currPolyn := f;
    
    while currPolyn != 0 do (
        leadLam := positions((listForm leadTerm currPolyn)#0#0,i -> i>0);
        leadLam = rsort(leadLam - toList((#leadLam):largestIndex));
        leadCoeff := leadCoefficient currPolyn;
        
        basisFunction := Q leadLam;
        basisCoeff := leadCoefficient sub(basisFunction,R);
        
        theCoeff := (sub(leadCoeff,R))//(sub(basisCoeff,R));
        ans = append(ans,(leadLam,theCoeff));
        
        currPolyn = currPolyn - theCoeff * (basisFunction);
        );
    
    ans = rsort ans;
    
    if o.doPrint then (
        -- prints decomposition into latex code
        print("\n");
        print(decompToTex(ans));
        
        -- prints decomposition into m2 functions
        print("\n");
        print(decompToM2(ans));
        );
    
    ans
    )

-- returns decomp as string of LaTeX
decompToTex = (theDecomp) -> (
    if theDecomp == {} then return (toString(0));
    replace("\\+-","-",replace("\\*","",replace("\\+$","",concatenate(for theTerm in theDecomp list (
                    if theTerm#0 == {} then (toString(theTerm#1)|"+")
                        else if theTerm#1 != 1 then (toString(theTerm#1)|"Q_{("|toString(theTerm#0)|")}+")
                        else ("Q_{("|toString(theTerm#0)|")}+")
                        )))))
    )

-- returns decomp as string of M2 code
decompToM2 = (theDecomp) -> (
    replace("\\+$","",concatenate(for theTerm in theDecomp list ("("|toString(theTerm#1)|")*Q("|toString(theTerm#0)|")+")))
    )

---------- Plethysm recurrence formulas

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

--computes list of I,J for BK formula
IJlist = (k,n,m) -> (
    ans := {(toList(m:0),toList(m:0))};
    
    curr := (2*m):0;
    while true do (
        curr = seqAdd(curr,n);
        if curr == (2*m):0 then break;
        I := toList take(curr,{0,m-1});
        J := toList take(curr,{m,2*m-1});
        ans = append(ans, (I,J));
        );
    
    ans
    )

--computes RHS of D_k(q_n(q_m)) formula
BKformula = (k,n,m) -> (
    indexList := IJlist(k,n,m);
    del := (1..m);
    
    sum for theInd in indexList list (
        I := toList theInd#0;
        J := toList theInd#1;
        
        (-1)^(dotProd(I+J,del)-k)*q_(dotProd(I+J,del)-k)*pleth(
            q_(n-(sum I)-(sum J)),q_m)*product for s from 0 to m-1 list (
            pleth(q_(I#s),q_(m-s-1))*pleth(q_(J#s),q_(m-s-1))
            )
        )
    )

--computes D_k(F)
Dk = (k,F) -> (
    theDecomp := decomposeQ(F,doPrint => false);
    
    ansDecomp := {};
    
    for theTerm in theDecomp do (
        if any(theTerm#0,i -> i == k) then (
            ansDecomp = append(ansDecomp,(delete(k,theTerm#0),(-1)^(position(theTerm#0,i -> i == k))*2*theTerm#1))
            );
        );
    
    sum for theTerm in ansDecomp list ((theTerm#1)*Q(theTerm#0))
    )
