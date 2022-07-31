-------------
-- (W22,7) --
-------------

u0 = (
    structureMap(C,3,4)
    *(wedgeProduct(4,1,C_1)**id_(exteriorPower(3,C_1)))
    *(id_(exteriorPower(4,C_1))**((C.dd_2**id_(exteriorPower(3,C_1)))*flip(exteriorPower(3,C_1),C_2)))
    );
u1 = (
    (structureMap(C,1,3)**id_(C_2))
    );
u2 = (
    (id_(exteriorPower(3,C_3))**C.dd_3)
    *wedgeDiag(3,1,C_3)
    *structureMap(C,2,4)
    );
u3 = (
    (structureMap(C,2,3)**structureMap(C,3,1))
    *(id_(exteriorPower(4,C_1)**C_1)**flip(exteriorPower(2,C_1),C_2))
    *(id_(exteriorPower(4,C_1))**wedgeDiag(1,2,C_1)**id_(C_2))
    );
u4 = (
    (wedgeProduct(1,2,C_3)**id_(C_2))
    *(id_(C_3)**flip(C_2,exteriorPower(2,C_3)))
    *(structureMap(C,3,2)**structureMap(C,2,2))
    );

--need to verify on S2221 F1
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

<< zero (((1/2)*u0 - (u1 + u2 - u3 - u4))*(i3**id_(C_2))) << endl
