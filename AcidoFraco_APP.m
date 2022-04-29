%   Funcao de Acido Fraco
%%  INPUT
%
%   app         - dados da aplicacao
%   OHInicial   - OH inicial                        (mol/l)
%   H3OInicial  - H3O inicial                       (mol/l)
%   RHInicial   - Ácido não ionizado inicial        (mol/l)
%   RInicial    - Ácido ionizado inicial            (mol/l)
%
%%  OPERACAO
%
%   Resgata de app
%       k1          - Constante de ionizacao do acido
%   Realiza o metodo da bissecao em R2:
%   Realiza o metodo da bissecao a equacao de equilibrio do acido fraco com
%   o incremento y obtido atraves do metodo da bissecao para os valores w
%   estipulados (atraves da função AcidoForteEX)
%
%%  DEPENDENCIAS
%
%   AcidoForte_APP
%   Cancelar
%
%%  OUTPUT
%   
%   D - vetor linha com o valor do incremento y e do incremento w,
%   respetivamente                                  (mol/l)
%
%%  NOTAS
%
%   QUIM
%       QUIM (1, :) - Incremento em R (eq.1)
%       QUIM (2, :) - Incremento em H3O (eq.2)
%       QUIM (3, :) - Erro da equacao com os valores supramencionados
%
%%
function D = AcidoFraco_APP(app, REAGENTES_INICIAIS)
    
    % Preparacao dos dados
    k1 = app.K11.Value;
    QUIM        = zeros(3,2);
    QUIM_Meio   = zeros(3,1);
    COUNTER     = 0;
    %   Fase Inicial
    QUIM(1,:)   = [-REAGENTES_INICIAIS(2)+2^(-26), REAGENTES_INICIAIS(1)-2^(-26)];
    REAGENTES   = [REAGENTES_INICIAIS(3), REAGENTES_INICIAIS(3); ...
                   REAGENTES_INICIAIS(4) + QUIM(1,:)];
    HOLDER      = AcidoForte_APP(app, REAGENTES(:,1));
    QUIM(2,1)   = HOLDER.INCREMENTOS;
    HOLDER      = AcidoForte_APP(app, REAGENTES(:,2));
    QUIM(2,2)   = HOLDER.INCREMENTOS;
    
    REAGENTES_EQUACAO   = [ QUIM(1,:) + REAGENTES_INICIAIS(2)               ;...
                            QUIM(1,:) + QUIM(2,:) + REAGENTES_INICIAIS(4)   ;...
                            REAGENTES_INICIAIS(1) -  QUIM(1,:)            	;...
                            k1, k1                                      	];
    EXPOENTES           = [ 1                                            	;...
                            1                                               ;...
                            -1                                             	;...
                            -1                                             	];
    
    QUIM(3,1) = EQUACOES_APP(REAGENTES_EQUACAO(:,1), EXPOENTES);
    QUIM(3,2) = EQUACOES_APP(REAGENTES_EQUACAO(:,2), EXPOENTES);

    %   Ciclo para encontrar uma solução
    while (((DIST_ABS(QUIM(2,1),QUIM(2,2)) > 10^(-14))&&(DIST_REL(QUIM(2,1),QUIM(2,2)) > 10^(-6)))  ||...
           ((DIST_ABS(QUIM(1,1),QUIM(1,2)) > 10^( -8))&&(DIST_REL(QUIM(1,1),QUIM(1,2)) > 10^(-6)))) &&...
             COUNTER < 400

        QUIM_Meio(1) = (QUIM(1,1)+QUIM(1,2))/2;
        REAGENTES    = [REAGENTES_INICIAIS(3)                   ;...
                        REAGENTES_INICIAIS(4) + QUIM_Meio(1)    ];
        HOLDER       = AcidoForte_APP(app, REAGENTES);
        QUIM_Meio(2) = HOLDER.INCREMENTOS;
        
        REAGENTES_EQUACAO   = [ QUIM_Meio(1) + REAGENTES_INICIAIS(2)                ;...
                                QUIM_Meio(1) + QUIM_Meio(2) + REAGENTES_INICIAIS(4) ;...
                                REAGENTES_INICIAIS(1) -  QUIM_Meio(1)               ;...
                                k1                                                  ];
        QUIM_Meio(3)        = EQUACOES_APP(REAGENTES_EQUACAO, EXPOENTES);
        
        QUIM = Bissecao_APP(QUIM, QUIM_Meio);
        %   LOOP COUNTER
        COUNTER = COUNTER + 1;
    end
    
    D.ERROS         = [ DIST_ABS(QUIM(1,1),QUIM(1,2))   ;...
                        DIST_ABS(QUIM(2,1),QUIM(2,2))   ];
    D.INCREMENTOS   = QUIM(1:(size(QUIM,1)-1), 1);

end