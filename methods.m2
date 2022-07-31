mapInfo = method()
mapInfo Matrix := Net => M -> (
    return (
	(net "degrees: "|netList apply((entries M),rrow->rrow/(entry->if not zero entry then return degree entry else return null)))||
	(net "# terms: "|netList apply((entries M),rrow->rrow/(entry->#terms entry)))
	)
    )

numterms = method()
numterms Matrix := List => M -> apply((entries M),rrow->rrow/(entry->(
	    if not zero (#terms entry) then return (#terms entry) else return "  "
	    )
	    ))

tr = method()
tr Module := Matrix => M -> (
    return adjoint'(id_M,dual M,(ring M)^1)
    )

wedgeDiag = method();
wedgeDiag (ZZ,ZZ,Module) := Matrix => (a,b,M) -> (
    return dual wedgeProduct(a,b,dual M)
    )

genericLift = method();
genericLift (Matrix,Matrix,Thing) := Matrix => (f,g,b) -> (
    v := f//map(target f,,g);
    if not zero (g*v - f) then error "lifting failed";
    if (class b === Nothing) then return v;
    R := ring f;
    def := gens trim ker g;
    m := rank source def; n := rank source f;
    S := R[b_1..b_(m*n)];
    return v+def*genericMatrix(S,m,n)
    )

randomLift = method();
randomLift (Matrix,Matrix) := Matrix => (f,g) -> (
    v := f//map(target f,,g);
    if not zero (g*v - f) then error "lifting failed";
    R := ring f;
    def := gens trim ker g;
    m := rank source def; n := rank source f;
    return v+def*random(R^m,R^n)
    )

testCoeffs = method();
testCoeffs (Matrix,Matrix) := Matrix => (f,g) -> (
    g' := dual gens trim ker dual g;
    return trim ideal (g'*f)
    )

symProd = method()
symProd (ZZ,ZZ,Module) := Matrix => (a,b,M) -> (
    R := ring M;
    d := rank M;
    if a==0 or b==0 then return id_((ring M)^(binomial(d-1+a+b,a+b)));
    LA := LB := L1 := apply(d,i->{i});
    for i from 2 to a do (
	LA = select((LA**L1)/flatten, j->j_(i-1)>=j_(i-2))
	);
    for i from 2 to b do (
	LB = select((LB**L1)/flatten, j->j_(i-1)>=j_(i-2))
	);
    --find the basis element of sym that corresponds to a
    --weakly increasing list of indices
    symLookup := p -> (
	binomial(d+#p-1,#p)-1-sum(#p, i -> binomial(d-(p_i)-2+#p-i,#p-i))
	);
    return map(
	R^(binomial(d-1+a+b,a+b)),R^(binomial(d-1+a,a)*binomial(d-1+b,b)),
	apply(apply(#LA,i->i)**apply(#LB,i->i),
	    i -> (symLookup sort (LA_(i_0)|LB_(i_1)), (i_0)*#LB+i_1)=>1)
	);
    )

--multiplication D^a \otimes D^b \to D^{a+b}
divProd = method()
divProd (ZZ,ZZ,Module) := Matrix => (a,b,M) -> (
    R := ring M;
    d := rank M;
    if a==0 or b==0 then return id_((ring M)^(binomial(d-1+a+b,a+b)));
    LA := LB := L1 := apply(d,i->{i});
    for i from 2 to a do (
	LA = select((LA**L1)/flatten, j->j_(i-1)>=j_(i-2))
	);
    for i from 2 to b do (
	LB = select((LB**L1)/flatten, j->j_(i-1)>=j_(i-2))
	);
    --find the basis element of sym that corresponds to a
    --weakly increasing list of indices
    symLookup := p -> (
	binomial(d+#p-1,#p)-1-sum(#p, i -> binomial(d-(p_i)-2+#p-i,#p-i))
	);
    coeff := (L1,L2) -> (
	L := L1|L2;
	return product(unique L, a -> binomial(number(L,i->i==a),number(L1,i->i==a)))
	);
    return map(
	R^(binomial(d-1+a+b,a+b)),R^(binomial(d-1+a,a)*binomial(d-1+b,b)),
	apply(apply(#LA,i->i)**apply(#LB,i->i),
	    i -> (symLookup sort (LA_(i_0)|LB_(i_1)), (i_0)*#LB+i_1)=>coeff(LA_(i_0),LB_(i_1)))
	);
    )

wedgeDiag = method();
wedgeDiag (ZZ,ZZ,Module) := Matrix => (a,b,M) -> (
    return dual wedgeProduct(a,b,dual M)
    )

structureMap = method()
structureMap (ChainComplex,ZZ,ZZ) := Matrix => (C,i,j) -> (
    return C.cache.StructureMapCache#(i,j)
    )


