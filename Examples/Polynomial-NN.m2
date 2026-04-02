needsPackage "EuclideanDistanceDegree"

-- Define the architecture
d = (3,1,1);  -- size of layers
r = 2;  -- degree of activation function

-- Parameter space: weights
R = QQ[W_0..W_(d_1 * d_0), V_0..V_(d_2 * d_1)];
W = genericMatrix(R, W_0, d_1, d_0);
V = genericMatrix(R, V_0, d_2, d_1);

-- Function space
T = R[x_0..x_(d_0 - 1)];
X = transpose matrix{apply(d_0, i -> x_i)};  -- input
Z = W * X;
A = matrix table(d_1, 1, (i, j) -> (Z_(i,j))^r);
Phi = V * A;  -- output function

-- Create the ambient space (space of symmetric tensors)
mons = flatten entries basis(r, T);  -- homogeneous monomials of degree r
S = QQ[c_0..c_(d_2 * #mons - 1)];

-- Get image of the parameterization map
im = flatten apply(d_2, i -> (
    f = Phi_(i,0);
    apply(mons, m -> coefficient(m, f))
));

-- Get the defining ideal
paramMap = map(R, S, im);
I = kernel paramMap;
F = flatten entries gens I;

-- Prune variety down to a complete intersection (for left kernel)
c = codim I;
G = apply(c, i -> sum apply(#F, j -> (random(QQ) * F_j)));
GED_left = leftKernelGenericEDDegree(G)
GED_num = numericGenericEDDegree(F, F)
GED_Left === GED_num