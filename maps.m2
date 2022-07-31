if vSymbolic then (
    s_1 = transpose C.dd_1;
    s_2 = transpose C.dd_2;
    s_3 = transpose C.dd_3;
    (f0,f1,f2,f3) = toSequence frmt; ic=1;
    dv={b,c,d};
    --defect variables
    seqToNum = L -> (
    	sum(#L, i -> 10^i*(reverse L)_i)
    	);
    inds = L -> (
    	apply(L, l -> toSequence(seqToNum\l))
    	);
    --the subscripts on b match the basis order for \bigwedge^2 F_1^* \otimes F_3
    --c for wedge 2 f3
    varlist = (
    	(i -> (dv_0)_i)\inds(subsets(ic..(ic+f1-1),2)**subsets(ic..(ic+f3-1),1))|
    	(i -> (dv_1)_i)\inds(subsets(ic..(ic+f1-1),4)**subsets(ic..(ic+f3-1),2))
    	);
    deglist = (
    	((degrees exteriorPower(2,C_1))**(degrees C_3))|
    	((degrees exteriorPower(4,C_1))**(degrees exteriorPower(2,C_3)))
    	);
    --for the third set of defect vars, depends if format is
    --(1,5,*,*) or (1,*,*,2)
    if (f1>=5 and f3>2) then (
    	varlist = varlist | (
	    (i -> (dv_2)_i)\inds((p -> {p_0_0,p_0_1,p_1})\(subsets(ic..(ic+f1-1),5)**subsets(ic..(ic+f1-1),1)**subsets(ic..(ic+f3-1),3)))
            );
    	deglist = deglist | (
	    ((degrees (exteriorPower(5,C_1)**C_1))**(degrees exteriorPower(3,C_3)))
	    );
    	);
    -*
    if (f1>5 and f3==2) then (
    	varlist = varlist | (
	    (i -> (dv_2)_i)\inds(subsets(ic..(ic+f1-1),6)**subsets(ic..(ic+f3-1),1))
	    );
    	deglist = deglist | (
	    ((degrees exteriorPower(6,C_1))**(degrees (exteriorPower(2,C_3)**C_3)))
	    );
    	);
    *-
    Sdef = S[varlist,Degrees => deglist/(L -> L_0-L_1),Join => false];
    --generic wedge^2 F1 -> F3
    usedvars = 0; C = C**Sdef;
    C.cache.StructureMapCache = new MutableHashTable;
    --v^(3)_1
    f = (
    	(C.dd_1**id_(C_1))
    	*wedgeDiag(1,1,C_1)
    	);
    if (f1 >= 2) then (
    	defmatrix1 = genericMatrix(Sdef,Sdef_usedvars,f3,binomial(f1,2));
    	usedvars = usedvars + f3*binomial(f1,2);
    	C.cache.StructureMapCache#(3,1) = (
            s_2*f
	    +(C.dd_3)*defmatrix1
	    );
    	);
    --v^(3)_2
    f = (1/2)*(
    	symProd(1,1,C_2)
    	*(structureMap(C,3,1)**structureMap(C,3,1))
    	*wedgeDiag(2,2,C_1)
    	);
    if (f3 >= 2 and f1 >= 4) then (
    	defmatrix2 = genericMatrix(Sdef,Sdef_usedvars,binomial(f3,2),binomial(f1,4));
    	usedvars = usedvars + binomial(f3,2)*binomial(f1,4);
    	C.cache.StructureMapCache#(3,2) = (
            (s_3**(s_2*C.dd_2 + (1/2)*C.dd_3*s_3))*(dual divProd(1,1,C_2))*f
	    +(id_(C_3)**C.dd_3)*(dual wedgeProduct(1,1,C_3))*defmatrix2
	    );
    	);
    --v^(3)_3,1
    f1' = (
    	(id_(C_3) ** symProd(1,1,C_2))
    	*(structureMap(C,3,2)**structureMap(C,3,1))
    	*(id_(exteriorPower(4,C_1))**wedgeProduct(1,1,C_1))
    	*((transpose wedgeProduct(4,1,C_1))**id_(C_1))
    	);
    f2' = (
    	(id_(C_3) ** symProd(1,1,C_2))
    	*(structureMap(C,3,2)**structureMap(C,3,1))
    	*flip(exteriorPower(2,C_1),exteriorPower(4,C_1))
    	*(id_(exteriorPower(2,C_1))**wedgeProduct(3,1,C_1))
    	*((transpose wedgeProduct(2,3,C_1))**id_(C_1))
    	);
    f = f1' - (1/2)*f2';
    sectiontest = (
    	(wedgeProduct(1,1,C_3)**id_(C_2))
    	*(id_(C_3)**((s_3**((1/2)**s_2*C.dd_2+(1/3)**C.dd_3*s_3)))*(dual divProd(1,1,C_2)))
    	);
    if (f3 >= 3 and f1 >= 5) then (
    	defmatrix3 = genericMatrix(Sdef,Sdef_usedvars,binomial(f3,3),binomial(f1,5)*f1);
    	usedvars = usedvars + binomial(f3,3)*(binomial(f1,4)*f1);
    	C.cache.StructureMapCache#(3,3) = (
            sectiontest*f
	    +(id_(exteriorPower(2,C_3))**C.dd_3)*(dual wedgeProduct(2,1,C_3))*defmatrix3
	    );
    	);
    )

---------------------------------
-- v^(3)_1,2,3,4 by random lifts --
---------------------------------
if not C.cache.?StructureMapCache then C.cache.StructureMapCache = new MutableHashTable;

--v^(3)_1
if not C.cache.StructureMapCache#?(3,1) then (
    f := (
    	(C.dd_1**id_(C_1))
    	*wedgeDiag(1,1,C_1)
    	);
    g := C.dd_2;
    C.cache.StructureMapCache#(3,1) = randomLift(f,g)
    )

--v^(3)_2
if not C.cache.StructureMapCache#?(3,2) then (
    f = (1/2)*(
    	symProd(1,1,C_2)
    	*(structureMap(C,3,1)**structureMap(C,3,1))
    	*wedgeDiag(2,2,C_1)
    	);
    g = (
    	symProd(1,1,C_2)
    	*(C.dd_3**id_(C_2))
    	);
    C.cache.StructureMapCache#(3,2) = randomLift(f,g);
    )

--v^(3)_3
if not C.cache.StructureMapCache#?(3,3) then (
    f1 := (
    	(id_(C_3) ** symProd(1,1,C_2))
    	*(structureMap(C,3,2)**structureMap(C,3,1))
    	*(id_(exteriorPower(4,C_1))**wedgeProduct(1,1,C_1))
    	*((transpose wedgeProduct(4,1,C_1))**id_(C_1))
    	);
    f2 := (
    	(id_(C_3) ** symProd(1,1,C_2))
    	*(structureMap(C,3,2)**structureMap(C,3,1))
    	*flip(exteriorPower(2,C_1),exteriorPower(4,C_1))
    	*(id_(exteriorPower(2,C_1))**wedgeProduct(3,1,C_1))
    	*((transpose wedgeProduct(2,3,C_1))**id_(C_1))
    	);
    f = f1 - (1/2)*f2; --works
    g = (
    	(id_(C_3) ** symProd(1,1,C_2))
    	*(id_(C_3) ** C.dd_3 ** id_(C_2))
    	*((transpose wedgeProduct(1,1,C_3))**id_(C_2))
    	);
    C.cache.StructureMapCache#(3,3) = (2/3)*randomLift(f,g);
    )

if not C.cache.StructureMapCache#?(3,4) then (
    f1 = (
    	(wedgeProduct(1,1,C_3)**symProd(1,1,C_2))
    	*(id_(C_3)**flip(C_2,C_3)**id_(C_2))
    	*(structureMap(C,3,2)**structureMap(C,3,2))
    	*(id_(exteriorPower(4,C_1))**wedgeProduct(1,3,C_1))
    	*((dual wedgeProduct(4,1,C_1))**id_(exteriorPower(3,C_1)))
    	);
    f2 = (
    	(id_(exteriorPower(2,C_3))**symProd(1,1,C_2))
    	*(structureMap(C,3,3)**structureMap(C,3,1))
    	*(id_(exteriorPower(5,C_1))**(dual wedgeProduct(1,2,C_1)))
    	);
    f3 = (
    	(id_(exteriorPower(2,C_3))**symProd(1,1,C_2))
    	*flip(C_2,exteriorPower(2,C_3)**C_2)
    	*(structureMap(C,3,1)**structureMap(C,3,3))
    	*(id_(exteriorPower(2,C_1))**wedgeProduct(3,2,C_1)**id_(C_1))
    	*((dual wedgeProduct(2,3,C_1))**(dual wedgeProduct(2,1,C_1)))
    	);
    g = (
    	(id_(exteriorPower(2,C_3)) ** symProd(1,1,C_2))
    	*(id_(exteriorPower(2,C_3)) ** C.dd_3 ** id_(C_2))
    	*((transpose wedgeProduct(2,1,C_3))**id_(C_2))
    	);
    f = (-1)*f1 + (3/2)*f2 + (1/2)*f3;
    if vSymbolic then (
	C.cache.StructureMapCache#(3,4) = genericLift(f,g,e)
    	) else (
	C.cache.StructureMapCache#(3,4) = randomLift(f,g)
	) 
    )

---------------
-- v^(i)_1,2 --
---------------

--v^(2)_1
f = (	 
    (C.dd_1 ** id_(C_2)) - (structureMap(C,3,1))*(wedgeProduct(1,1,C_1))*(id_(C_1) ** C.dd_2)
    );
g = C.dd_3;
C.cache.StructureMapCache#(2,1) = randomLift(f,g);

--v^(1)_1
f = (C.dd_1 ** structureMap(C,3,1))*(transpose wedgeProduct(1,2,C_1));
g = C.dd_3
C.cache.StructureMapCache#(1,1) = randomLift(f,g);

--v^(2)_2
f = (
    (structureMap(C,1,1)**id_(C_2))
    -
    flip(C_2,C_3)*(structureMap(C,3,1)**structureMap(C,2,1))*(wedgeDiag(2,1,C_1)**id_(C_2))
    -
    (structureMap(C,3,2)*wedgeProduct(3,1,C_1)*(id_(exteriorPower(3,C_1))**C.dd_2))
    );
g = (id_(C_3)**C.dd_3)*(wedgeDiag(1,1,C_3));
C.cache.StructureMapCache#(2,2) = randomLift(f,g);

--v^(1)_2
g = (
    (id_(C_3)**C.dd_3)
    *wedgeDiag(1,1,C_3)
    );
f1 = (
    (id_(C_3)**structureMap(C,3,1))
    *(structureMap(C,1,1)**wedgeProduct(1,1,C_1))
    *(wedgeDiag(3,1,C_1)**id_(C_1))
    );
f2 = (
    ((structureMap(C,1,1)*wedgeProduct(1,2,C_1))**structureMap(C,3,1))
    *(id_(C_1)**wedgeDiag(2,2,C_1))
    *flip(exteriorPower(4,C_1),C_1)
    );
f3 = (
    (structureMap(C,3,2)**C.dd_1)
    );
f4 = (
    (C.dd_1**(structureMap(C,3,2)*wedgeProduct(3,1,C_1)))
    *(wedgeDiag(1,3,C_1)**id_(C_1))
    );
f' = (-1)*f1 + f4;
C.cache.StructureMapCache#(1,2) = randomLift(f',g); --A

--v^(2)_3
f1' = (
    structureMap(C,1,2)**id_(C_2)
    );
f2' = (
    flip(C_2,exteriorPower(2,C_3))
    *(structureMap(C,3,1)**id_(exteriorPower(2,C_3)))
    *(wedgeProduct(1,1,C_1)**structureMap(C,2,2))
    *(id_(C_1)**wedgeDiag(1,3,C_1)**id_(C_2))
    *(flip(exteriorPower(4,C_1),C_1)**id_(C_2))
    );
f3' = (
    (wedgeProduct(1,1,C_3)**id_(C_2))
    *(id_(C_3)**flip(C_2,C_3))
    *(structureMap(C,3,2)**id_(C_3))
    *(wedgeProduct(1,3,C_1)**structureMap(C,2,1))
    *(id_(C_1)**wedgeDiag(3,1,C_1)**id_(C_2))
    *(flip(exteriorPower(4,C_1),C_1)**id_(C_2))
    );
f4' = (
    structureMap(C,3,3)
    *(wedgeProduct(4,1,C_1)**C.dd_2)
    );
f4'' = (
    structureMap(C,3,3)
    *(wedgeProduct(4,1,C_1)**id_(C_1))
    *(id_(exteriorPower(4,C_1))**(flip(C_1,C_1)*(id_(C_1)**C.dd_2)))
    );

g = (
    (id_(exteriorPower(2,C_3))**C.dd_3)
    *wedgeDiag(2,1,C_3)
    );

C.cache.StructureMapCache#(2,3) = randomLift(f1' + f2' - f3' - f4' - f4'',g);

--v^(1)_3
--maps wedge^4 F1 ** wedge^3 F1 -> wedge^3 F3
f1' = (
    (structureMap(C,1,2)**structureMap(C,3,1))
    *(id_(exteriorPower(4,C_1))**wedgeDiag(1,2,C_1))
    );
f2' = (
    (wedgeProduct(1,1,C_3)**id_(C_2))
    *(id_(C_3)**flip(C_2,C_3))
    *(structureMap(C,3,2)**structureMap(C,1,1))
    *(wedgeProduct(3,1,C_1)**id_(exteriorPower(3,C_1)))
    *(id_(exteriorPower(3,C_1))**wedgeDiag(1,3,C_1))
    *flip(exteriorPower(4,C_1),exteriorPower(3,C_1))
    );
f3' = (
    structureMap(C,3,3)
    *(wedgeProduct(3,2,C_1)**C.dd_1**id_(C_1))
    *(id_(exteriorPower(3,C_1))**flip(C_1,exteriorPower(2,C_1))**id_(C_1))
    *(wedgeDiag(3,1,C_1)**wedgeDiag(2,1,C_1))
    );
g = (
    (id_(exteriorPower(2,C_3))**C.dd_3)
    *wedgeDiag(2,1,C_3)
    );
C.cache.StructureMapCache#(1,3) = randomLift(f1' - f2' + f3',g);

--v^(2)_4
--try to define q^(2)_4, a map wedge^4 F1 ** wedge^3 F1 ** F2 -> wedge^4 F3
--We will only define it on S_2221 F1 ** F2.

f1' = (
    structureMap(C,3,4)
    *(wedgeProduct(4,1,C_1)**id_(exteriorPower(3,C_1)))
    *(id_(exteriorPower(4,C_1))**((C.dd_2**id_(exteriorPower(3,C_1)))*flip(exteriorPower(3,C_1),C_2)))
    );
f2' = (
    (structureMap(C,1,3)**id_(C_2))
    );
g = (
    (id_(exteriorPower(3,C_3))**C.dd_3)
    *wedgeDiag(3,1,C_3)
    );
f3' = (
    (structureMap(C,2,3)**structureMap(C,3,1))
    *(id_(exteriorPower(4,C_1)**C_1)**flip(exteriorPower(2,C_1),C_2))
    *(id_(exteriorPower(4,C_1))**wedgeDiag(1,2,C_1)**id_(C_2))
    );
f4' = (
    (wedgeProduct(1,2,C_3)**id_(C_2))
    *(id_(C_3)**flip(C_2,exteriorPower(2,C_3)))
    *(structureMap(C,3,2)**structureMap(C,2,2))
    );

needsPackage "SchurFunctors"
p0 = wedgeProduct(4,3,C_1);
p1 = (
    ((schurModule({2,1,1,1,1,1},C_1)).cache#"Schur"#0)
    *(wedgeProduct(4,2,C_1)**id_(exteriorPower(1,C_1)))
    *(id_(exteriorPower(4,C_1))**wedgeDiag(2,1,C_1))
    );
p2 = (
    ((schurModule({2,2,1,1,1},C_1)).cache#"Schur"#0)
    *(wedgeProduct(4,1,C_1)**id_(exteriorPower(2,C_1)))
    *(id_(exteriorPower(4,C_1))**wedgeDiag(1,2,C_1))
    );
p3 = (schurModule({2,2,2,1},C_1)).cache#"Schur"#0;

M0 = schurModule({1,1,1,1,1,1,1},C_1)
M1 = schurModule({2,1,1,1,1,1},C_1)
M2 = schurModule({2,2,1,1,1},C_1)
M3 = schurModule({2,2,2,1},C_1)

M = directSum(M0,M1,M2,M3)

decomp = (
    M_[0]*p0
    +M_[1]*p1
    +M_[2]*p2
    +M_[3]*p3
    );

inccomp = inverse decomp;

i0=inccomp*M_[0];
i1=inccomp*M_[1];
i2=inccomp*M_[2];
i3=inccomp*M_[3];

f = ((1/2)*f1' - f2' + f3' + f4')*(i3**id_(C_2));

g = (
    (id_(exteriorPower(3,C_3))**C.dd_3)
    *wedgeDiag(3,1,C_3)
    );

C.cache.StructureMapCache#(2,4) = randomLift(f,g)*(p3**id_(C_2));
