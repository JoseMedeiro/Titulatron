%   Função de Equações
%%  INPUT
%
%   REAGENTES   - OH inicial                        (mol/l)
%   EXPOENTES   - H3O inicial                       (mol/l)
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
function D = EQUACOES_APP(REAGENTES, EXPOENTES)
            
    D = 1;
    for c=1:size(REAGENTES,1)
        D = D * REAGENTES(c)^EXPOENTES(c);
    end
	D = D - 1;

end