%   Função de Ácido Forte
%%  INPUT
%
%   OHInicial   - OH inicial                        (mol/l)
%   H3OInicial  - H3O inicial                       (mol/l)
%   KH2O        - Constante de ionização da água
%
%%  OPERAÇÃO
%
%   Realiza o método da bisseção à equação de equilíbrio da água para
%   encontrar o incremento y
%
%%  DEPENDÊNCIAS
%
%   Nenhuma
%
%%  OUTPUT
%   
%   D - valor do incremento y                       (mol/l)
%
%%
function D = AcidoForteEX(OHInicial, H3OInicial, KH2O)

d = 0;
%   Fase Inicial
if OHInicial < H3OInicial
    y = [-OHInicial 0];
else
    y = [-H3OInicial 0];
end

ERR = (y + H3OInicial).*(y + OHInicial)/KH2O - 1;
d = 0;

%   Ciclo para encontrar uma solução
while abs(y(1)-y(2)) > 10^(-14)
    
    yMeio = (y(1)+y(2))/2;
    ERRMeio = (yMeio + H3OInicial)*(yMeio + OHInicial)/KH2O - 1;
    
    if ERRMeio < 0
        if ERR(1) < 0
            y(1)   = yMeio;
            ERR(1) = ERRMeio;
            FLAG   = 1;
        else
            y(2)   = yMeio;
            ERR(2) = ERRMeio;
            FLAG   = 2;
        end
    else
        if ERR(1) > 0
            y(1)   = yMeio;
            ERR(1) = ERRMeio;
            FLAG   = 1;
        else
            y(2)   = yMeio;
            ERR(2) = ERRMeio;
            FLAG   = 2;
        end
    end
    d = d+1;
    
end

D = y(FLAG);

end
