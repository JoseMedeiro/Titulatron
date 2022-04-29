%   Função de Ácido Fraco
%%  INPUT
%
%   OHInicial   - OH inicial                        (mol/l)
%   H3OInicial  - H3O inicial                       (mol/l)
%   RHInicial   - Ácido não ionizado inicial        (mol/l)
%   RInicial    - Ácido ionizado inicial            (mol/l)
%   k1          - Constante de ionização do ácido
%   KH2O        - Constante de ionização da água
%
%%  OPERAÇÃO
%
%   Realiza o método da bisseção em R2:
%   Realiza o método da bisseção à equação de equilíbrio do ácido fraco com
%   o incremento y obtido através do método da bisseção para os valores w
%   estipulados (através da função AcidoForteEX)
%
%%  DEPENDÊNCIAS
%
%   AcidoForteEX
%
%%  OUTPUT
%   
%   D - vetor linha com o valor do incremento y e do incremento w,
%   respetivamente                                  (mol/l)
%
%%
function D = AcidoFracoEX(OHInicial, H3OInicial, RHInicial, RInicial, k1, KH2O)

d = 0;
%   Fase Inicial

w = [-RInicial+2^(-26) RHInicial-2^(-26)];

y(1) = AcidoForteEX(OHInicial, H3OInicial+w(1), KH2O);
y(2) = AcidoForteEX(OHInicial, H3OInicial+w(2), KH2O);

ERR = (w+RInicial).*(w + y + H3OInicial)./(k1*(RHInicial-w)) - 1;
d = 0;
FLAG = 2;

%   Ciclo para encontrar uma solução
while ((abs(y(1)-y(2)) > 10^(-14))&&(abs((y(1)-y(2))/(y(1)-y(2))*2) > 10^(-6)))||...
      ((abs(w(1)-w(2)) > 10^( -8))&&(abs((w(1)-w(2))/(w(1)-w(2))*2) > 10^(-6)))
    
    wMeio = (w(1)+w(2))/2;
    yMeio = AcidoForteEX(OHInicial, H3OInicial+wMeio, KH2O);
    ERRMeio = (wMeio+RInicial)*(wMeio + yMeio + H3OInicial)./(k1*(RHInicial-wMeio)) - 1;
    
    if ERRMeio < 0
        if ERR(1) < 0
            w(1)   = wMeio;
            y(1)   = yMeio;
            ERR(1) = ERRMeio;
            FLAG   = 1;
        else
            w(2)   = wMeio;
            y(2)   = yMeio;
            ERR(2) = ERRMeio;
            FLAG   = 2;
        end
    else
        if ERR(1) > 0
            w(1)   = wMeio;
            y(1)   = yMeio;
            ERR(1) = ERRMeio;
            FLAG   = 1;
        else
            w(2)   = wMeio;
            y(2)   = yMeio;
            ERR(2) = ERRMeio;
            FLAG   = 2;
        end
    end
    d = d+1;
    
end

D = [y(FLAG) w(FLAG)];

end
