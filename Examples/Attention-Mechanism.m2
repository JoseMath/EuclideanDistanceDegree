needsPackage "EuclideanDistanceDegree"

ds = ((1,1), (1,2), (2,1), (2,2), (1,3), (2,3), (3,1), (3,2), (3,3));
t = 2;
for d in ds do (
    l = #d - 1;

    -- Parameter space: weights
    R = QQ[A_0..A_(d_0 * d_0)] ** QQ[V_0..V_(d_1 * d_0)];
    M = genericMatrix(R, A_0, d_0, d_0);
    N = genericMatrix(R, V_0, d_1, d_0);

    -- Function space
    T = R[x_0..x_(d_0 * t)];
    X = genericMatrix(T, x_0, d_0, t);  -- input
    Phi = N * X * transpose X * M * X; -- output

    -- Create the ambient space (space of symmetric tensors)
    mons = flatten entries basis(3^l, T);  -- homogeneous monomials of degree 3^l
    S = QQ[c_0..c_(d_1 * t * #mons - 1)];

    -- Get image of the parameterization map
    PhiFlat = flatten entries Phi;
    im = flatten apply(d_1 * t, i -> (
        f = PhiFlat_i;
        apply(mons, m -> coefficient(m, f))
    ));

    -- Get the defining ideal
    paramMap := map(R, S, im);
    I = kernel paramMap;
    F = flatten entries gens I;
    print(toString d | " " | codim I);
)

-- Prune variety down to complete variety (for left kernel)
I = getI((1,1,2), 2);
c = codim I;
G = apply(c, i -> sum apply(#F, j -> (random(QQ) * F_j)));
GED_left = leftKernelGenericEDDegree(G)
GED_num = numericGenericEDDegree(F, F)
GED_Left === GED_num