-- needsPackage "EuclideanDistanceDegree"

getAttentionIdeal = method()
getAttentionIdeal(Sequence, ZZ) := (d, t) -> (
    l := #d - 1;

    -- Parameter space: weights
    varList := flatten apply(l, i -> (apply(d_(i+1) * d_i, j -> w_(i,j))));
    varList = varList | flatten apply(l, i -> (apply(d_i * d_i, j -> v_(i,j))));
    R := QQ[varList];
    Ws := apply(l, i -> genericMatrix(R, w_(i,0), d_(i+1), d_i));
    As := apply(l, i -> genericMatrix(R, v_(i,0), d_i, d_i));

    -- Function space
    T := R[x_0..x_(d_0 * t)];
    X := genericMatrix(T, x_0, d_0, t);  -- input
    apply(l, i -> (
        X = Ws_i * X * transpose X * As_i * X;
    ));
    Phi := flatten entries X;

    -- Create the ambient space (space of symmetric tensors)
    mons := flatten entries basis(3^l, T);  -- homogeneous monomials of degree 3^l
    S := QQ[c_0..c_(d_l * t * #mons - 1)];

    -- Get image of the parameterization map
    im := flatten apply(d_l * t, i -> (
        f := Phi_i;
        apply(mons, m -> coefficient(m, f))
    ));

    -- Get the defining ideal
    paramMap := map(R, S, im);
    I := kernel paramMap;
    F := flatten entries gens I;

    -- Prune variety down to a complete intersection
    c := codim I;
    G := F;
    if #G != c then G = apply(c, i -> sum apply(#F, j -> (random(QQ) * F_j)));
    (I, F, G)
)

end

I = getAttentionIdeal((1,1,2), 2);
c = codim I;
G = apply(c, i -> sum apply(#F, j -> (random(QQ) * F_j)));
GED_left = leftKernelGenericEDDegree(G)
GED_num = numericGenericEDDegree(F, F)
GED_Left === GED_num