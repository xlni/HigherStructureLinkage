-------------
-- (W11,3) --
-------------
-- wedge^3 F1 ** wedge^3 F1 ** F2 -> wedge^3 F3

--LHS
u0 = (
    structureMap(C,1,3)
    *(wedgeProduct(1,3,C_1)**id_(exteriorPower(3,C_1)))
    *(C.dd_2**id_(exteriorPower(3,C_1)**exteriorPower(3,C_1)))
    *flip(exteriorPower(3,C_1)**exteriorPower(3,C_1),C_2)
    );
--RHS
-*
u1 = (
    wedgeProduct(2,1,C_3)
    *(structureMap(C,1,2)**id_(C_3))
    *(wedgeProduct(3,1,C_1)**id_(C_1**C_3))
    *(id_(exteriorPower(3,C_1))**wedgeDiag(1,1,C_1)**structureMap(C,2,1))
    *(id_(exteriorPower(3,C_1))**wedgeDiag(2,1,C_1)**id_(C_2))
    );
*-
u2 = (
    wedgeProduct(1,2,C_3)
    *(structureMap(C,1,1)**structureMap(C,2,2))
    );
u3 = (
    structureMap(C,2,3)
    *(wedgeProduct(3,1,C_1)**id_(C_1**C_2))
    *(id_(exteriorPower(3,C_1))**wedgeDiag(1,1,C_1)**id_(C_2))
    *(C.dd_1**flip(exteriorPower(2,C_1),exteriorPower(3,C_1))**id_(C_2))
    *(wedgeDiag(1,2,C_1)**id_(exteriorPower(3,C_1)**C_2))
    );

--need this to be true on subrep:
p111111 = wedgeProduct(3,3,C_1);
p21111 = (
    ((schurModule({2,1,1,1,1},C_1)).cache#"Schur"#0)
    *(wedgeProduct(3,2,C_1)**id_(exteriorPower(1,C_1)))
    *(id_(exteriorPower(3,C_1))**wedgeDiag(2,1,C_1))
    );
p2211 = (
    ((schurModule({2,2,1,1},C_1)).cache#"Schur"#0)
    *(wedgeProduct(3,1,C_1)**id_(exteriorPower(2,C_1)))
    *(id_(exteriorPower(3,C_1))**wedgeDiag(1,2,C_1))
    );
p222 = (schurModule({2,2,2},C_1)).cache#"Schur"#0;


M111111 = schurModule({1,1,1,1,1,1},C_1)
M21111 = schurModule({2,1,1,1,1},C_1)
M2211 = schurModule({2,2,1,1},C_1)
M222 = schurModule({2,2,2},C_1)

M = directSum(M222,M2211,M21111,M111111)

decomp = (
    M_[0]*p222
    +M_[1]*p2211
    +M_[2]*p21111
    +M_[3]*p111111
    );

--decomp = p1111111||p211111||p22111||p2221;

inccomp = inverse decomp;

i222=inccomp*M_[0];
i2211=inccomp*M_[1];
i21111=inccomp*M_[2];
i111111=inccomp*M_[3];

<< (
    (zero ((u0 - (-u2 + u3))*(i2211**id_(C_2))))
    and
    (zero ((u0 - (-u2 + u3))*(i222**id_(C_2))))
    ) << endl
