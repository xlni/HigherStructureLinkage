-------------
-- (W22,5) --
-------------
u0 = (
    structureMap(C,3,3)
    *(wedgeProduct(4,1,C_1)**id_(C_1))
    *(id_(exteriorPower(4,C_1))**((C.dd_2**id_(C_1))*flip(C_1,C_2)))
    );
u1 = (
    (structureMap(C,1,2)**id_(C_2))
    );
--zero when rank F3 <= 2
u2 = (
    (id_(exteriorPower(2,C_3))**C.dd_3)
    *wedgeDiag(2,1,C_3)
    *structureMap(C,2,3)
    );

u3 = (
    (structureMap(C,2,2)**structureMap(C,3,1))
    *(id_(exteriorPower(3,C_1))**(flip(exteriorPower(2,C_1),C_2)*(wedgeProduct(1,1,C_1)**id_(C_2))))
    *(wedgeDiag(3,1,C_1)**id_(C_1**C_2))
    );
u4 = (
    (wedgeProduct(1,1,C_3)**id_(C_2))
    *(id_(C_3)**flip(C_2,C_3))
    *(structureMap(C,3,2)**structureMap(C,2,1))
    );

p11111 = wedgeProduct(4,1,C_1);
M2111 = schurModule({2,1,1,1},C_1);

M = directSum(2111=>M2111,11111=>exteriorPower(5,C_1))

p2111 = M2111.cache#"Schur"#0;

decomp = p2111||p11111;
inccomp = inverse decomp;

i2111=inccomp*M_[2111];
i11111=inccomp*M_[11111];

<< zero ((u0 - (u1 - u2 + u3 + u4))*(i2111**id_(C_2))) << endl
