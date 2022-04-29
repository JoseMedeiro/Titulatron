
function D = Formula_AcidoFraco(A, H3OInicial, RInicial, RHInicial, k1)
    D = (A(1) + RInicial)*(A(1) + A(2) + H3OInicial)./(k1*(RHInicial - A(1))) - 1;
end