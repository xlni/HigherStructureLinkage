-------------
-- (W32,3) --
-------------

--wedge^3 F1 ** F1 ** S2 F2 -> 
--LHS (up to sign...)

u0 = (
    structureMap(C,2,3)
    *(wedgeProduct(1,3,C_1)**id_(C_1**C_2))
    *(flip(exteriorPower(3,C_1)**C_1,C_1)**id_(C_2))
    *(id_(exteriorPower(3,C_1)**C_1)**(
	    (C.dd_2**id_(C_2))
	    *(dual divProd(1,1,C_2))
	    ))
    );

--RHS (up to sign...)
u1 = (
    wedgeProduct(2,1,C_3)
    *(structureMap(C,2,2)**structureMap(C,2,1))
    *(id_(exteriorPower(3,C_1))**flip(C_1,C_2)**id_(C_2))
    *(id_(exteriorPower(3,C_1)**C_1)**(
	    (dual divProd(1,1,C_2))
	    ))
    );

M211 = schurModule({2,1,1},C_1);
M = directSum(211=>M211,1111=>exteriorPower(4,C_1))

p211 = M211.cache#"Schur"#0;

decomp = p211||wedgeProduct(3,1,C_1);
inccomp = inverse decomp;

i211=inccomp*M_[211];
i1111=inccomp*M_[1111];

<< zero ((u0 + u1)*(i211**id_(schurModule({2},C_2)))) << endl
