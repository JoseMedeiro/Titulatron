
function D = Formula_AcidoForte(A, H3OInicial, OHInicial, KH2O)
    D = (A + H3OInicial).*(A + OHInicial)/KH2O - 1;
end