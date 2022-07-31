-------------
-- (W32,4) --
-------------

--wedge^3 F1 ** wedge^3 F1 ** S2 F2 -> wedge^4 F3
--LHS (up to sign...)

u0 = (
    structureMap(C,2,4)
    *(wedgeProduct(1,3,C_1)**id_(exteriorPower(3,C_1)**C_2))
    *(flip(exteriorPower(3,C_1)**exteriorPower(3,C_1),C_1)**id_(C_2))
    *(id_(exteriorPower(3,C_1)**exteriorPower(3,C_1))**(
	    (C.dd_2**id_(C_2))
	    *(dual divProd(1,1,C_2))
	    ))
    );

--RHS (up to sign...)
u1 = (
    (1/2)*wedgeProduct(2,2,C_3)
    *(structureMap(C,2,2)**structureMap(C,2,2))
    *(id_(exteriorPower(3,C_1))**flip(exteriorPower(3,C_1),C_2)**id_(C_2))
    *(id_(exteriorPower(3,C_1)**exteriorPower(3,C_1))**(
	    (dual divProd(1,1,C_2))
	    ))
    );

--need to restrict to just S222 F1
p0 = wedgeProduct(3,3,C_1);
p1 = (
    ((schurModule({2,1,1,1,1},C_1)).cache#"Schur"#0)
    *(wedgeProduct(3,2,C_1)**id_(exteriorPower(1,C_1)))
    *(id_(exteriorPower(3,C_1))**wedgeDiag(2,1,C_1))
    );
p2 = (
    ((schurModule({2,2,1,1},C_1)).cache#"Schur"#0)
    *(wedgeProduct(3,1,C_1)**id_(exteriorPower(2,C_1)))
    *(id_(exteriorPower(3,C_1))**wedgeDiag(1,2,C_1))
    );
p3 = (schurModule({2,2,2},C_1)).cache#"Schur"#0;

M0 = schurModule({1,1,1,1,1,1},C_1)
M1 = schurModule({2,1,1,1,1},C_1)
M2 = schurModule({2,2,1,1},C_1)
M3 = schurModule({2,2,2},C_1)

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

<< zero ((u0 - u1)*(i3**id_(schurModule({2},C_2)))) << endl
