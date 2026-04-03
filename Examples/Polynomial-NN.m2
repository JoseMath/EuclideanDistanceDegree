needsPackage "EuclideanDistanceDegree"

-- d is a vector defining the layer widths, r is the activation degree
-- Returns ideal I, list of generators F, and pruned list (complete intersection) G
getPNNIdeal = method()
getPNNIdeal(Sequence, ZZ) := (d, r) -> (
    l := #d - 1;
    
    -- Parameter space: weights
    varList := flatten apply(l, i -> (apply(d_(i+1) * d_i, j -> w_(i,j))));
    R := QQ[varList];
    Ws := apply(l, i -> genericMatrix(R, w_(i,0), d_(i+1), d_i));

    -- Function space
    T := R[x_0..x_(d_0 - 1)];
    X := transpose matrix{apply(d_0, i -> x_i)};  -- input
    apply(l-1, i -> (
        X = Ws_i * X;
        X = matrix table(d_(i+1), 1, (j, k) -> (X_(j, k))^r);
    ));
    Phi := Ws_(l-1) * X;  -- output function

    -- Create the ambient space (space of symmetric tensors)
    mons = flatten entries basis(r^(l-1), T);  -- homogeneous monomials of degree r^{l-1}
    S := QQ[c_0..c_(d_l * #mons - 1)];

    -- Get image of the parameterization map
    im := flatten apply(d_l, i -> (
        f := Phi_(i,0);
        apply(mons, m -> coefficient(m, f))
    ));

    -- Get the defining ideal
    paramMap := map(R, S, im);
    I := kernel paramMap;
    F = flatten entries gens I;

    -- Prune variety down to a complete intersection (for left kernel)
    c := codim I;
    G := F;
    if #G != c then G := apply(c, i -> sum apply(#F, j -> (random(QQ) * F_j)));
    (I, G)
)

(I, F) = getPNNIdeal((3,1,1), 2);
GED_left = leftKernelGenericEDDegree(F)
GED_num = numericGenericEDDegree(F, F)
GED_Left === GED_num