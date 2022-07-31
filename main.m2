restart
needsPackage "SchurFunctors"

kk=ZZ/101; S=kk --base ring

--format for the complex
frmt={1,8,11,4}

--the complex
C = (
    chainComplex(-id_(S^(frmt_1-frmt_0)))[-1]
    ++chainComplex(id_(S^(frmt_0)))
    ++chainComplex(id_(S^(frmt_3)))[-2]
    )

vSymbolic = false --use random lifts
vSymbolic = true --use generic lifts (a lot slower)

--methods used in the subsequent files
load "methods.m2"

--compute maps v^(i)_j
load "maps.m2"

--verify relations in the paper
load "./relations/W11-3.m2"
load "./relations/W22-5.m2"
load "./relations/W22-7.m2"
load "./relations/W32-3.m2"
load "./relations/W32-4.m2"
