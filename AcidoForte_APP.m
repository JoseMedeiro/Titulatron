%   Funcao de Acido Forte
%%  INPUT
%
%   app         - dados da aplicacao
%   OHInicial   - OH inicial                        (mol/l)
%   H3OInicial  - H3O inicial                       (mol/l)
%   KH2O        - Constante de ionização da água
%
%%  OPERACAO
%
%   Resgata de app
%       KH2O        - Constante de ionizacao da agua
%   Realiza o metodo da bissecao a equacao de equilibrio da agua para
%   encontrar o incremento y
%
%%  DEPENDENCIAS
%
%   Cancelar
%
%%  OUTPUT
%   
%   D - valor do incremento y                       (mol/l)
%
%%  NOTAS
%
%   QUIM
%       QUIM (1, :) - Incremento de OH = incremento de H3O
%       QUIM (2, :) - Erro da equacao com os valores supramencionados
%
%%
function D = AcidoForte_APP(app, REAGENTES_INICIAIS)

    % Preparacao dos dados
    KH2O = app.KH2O;
    QUIM        = zeros(2,2);
    QUIM_Meio   = zeros(2,1);
    COUNTER     = 0;
    %   Fase Inicial
    if REAGENTES_INICIAIS(1) < REAGENTES_INICIAIS(2)
        QUIM(1,:) = [-REAGENTES_INICIAIS(1), 0];
    else
        QUIM(1,:) = [-REAGENTES_INICIAIS(2), 0];
    end
    
    REAGENTES_EQUACAO   = [ QUIM(1,:) + REAGENTES_INICIAIS(2)   ;...
                            QUIM(1,:) + REAGENTES_INICIAIS(1)   ;...
                            KH2O, KH2O                          ];
    EXPOENTES           = [ 1                                   ;...
                            1                                   ;...
                            -1                                  ];
    
    QUIM(2,1) = EQUACOES_APP(REAGENTES_EQUACAO(:,1), EXPOENTES);
    QUIM(2,2) = EQUACOES_APP(REAGENTES_EQUACAO(:,2), EXPOENTES);
    %   Ciclo para encontrar uma solução
    while (DIST_ABS(QUIM(1,1),QUIM(1,2)) > 10^(-14))&&...
           COUNTER < 400

        QUIM_Meio(1) = (QUIM(1,1)+QUIM(1,2))/2;
        
        REAGENTES_EQUACAO   = [ QUIM_Meio(1) + REAGENTES_INICIAIS(2)    ;...
                                QUIM_Meio(1) + REAGENTES_INICIAIS(1)    ;...
                                KH2O                                    ];
                            
        QUIM_Meio(2) = EQUACOES_APP(REAGENTES_EQUACAO, EXPOENTES);
        
        QUIM = Bissecao_APP(QUIM, QUIM_Meio);
        %   LOOP COUNTER
        COUNTER = COUNTER+1;
    end
    D.ERROS       = DIST_ABS(QUIM(1,1),QUIM(1,2));
    D.INCREMENTOS = QUIM(1:(size(QUIM,1)-1));

end
