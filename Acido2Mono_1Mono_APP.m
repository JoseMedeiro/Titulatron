%   Funcao de Acido Diprotico
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
%       k2          - Constante de ionizacao do acido
%   Realiza o metodo da bisseção em R2:
%   Realiza o metodo da bissecao a equacao de equilibrio do acido fraco com
%   o incremento y obtido atraves do metodo da bissecao para os valores w
%   estipulados (atraves da função AcidoFraco_APP)
%
%%  DEPENDENCIAS
%
%   AcidoFraco_APP
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
%       QUIM (1, :) - Incremento em RH                              (eq.1)
%       QUIM (2, :) - Incremento em R                               (eq.2)
%       QUIM (3, :) - Incremento em H3O                             (eq.3)
%       QUIM (4, :) - Erro da equacao com os valores supramencionados
%
%%
function D = Acido2Mono_1Mono_APP(app, REAGENTES_INICIAIS, REACAO_2)
    
    %Retirada de variaveis
    P1  = REAGENTES_INICIAIS(1);
    P0  = REAGENTES_INICIAIS(2);
    R1  = REAGENTES_INICIAIS(3);
    R0  = REAGENTES_INICIAIS(4);
    OH  = REAGENTES_INICIAIS(5);
    H3O = REAGENTES_INICIAIS(6);
    
    % Preparacao dos dados
    k1 = app.K21.Value;
    QUIM        = zeros(4,2);
    QUIM_Meio   = zeros(4,1);
    COUNTER     = 0;
    %   Fase Inicial
    QUIM(1,:)   = [-P0 + 2^(-26), P1 - 2^(-26)];
    % Para resolver as equações abaixo
    if REACAO_2 == 1
        REAGENTES   = [ R1  , R1        ;...
                        R0  , R0        ;...
                        OH  , OH        ;...
                        H3O + QUIM(1,:)	];
    else
        REAGENTES   = [ R1  , R1        ;...
                        R0  , R0        ;...
                        OH  + QUIM(1,:) ;...
                        H3O , H3O       ];     
    end
    % Resolver as equações abaixo
    HOLDER      = AcidoFraco_APP(app, REAGENTES(:,1));
    QUIM(2:3,1) = HOLDER.INCREMENTOS;
    HOLDER      = AcidoFraco_APP(app, REAGENTES(:,2));
    QUIM(2:3,2) = HOLDER.INCREMENTOS;
    
    REAGENTES_EQUACAO   = [ QUIM(1,:) + P0                          ;...
                            QUIM(1,:) + QUIM(2,:) + QUIM(3,:) + H3O ;...
                            QUIM(1,:) + QUIM(3,:) + OH              ;...
                            P1 - QUIM(1,:)                          ;...
                            k1, k1                                  ];
    if REACAO_2 == 1
        EXPOENTES           = [ 1   ;...
                                1 	;...
                                0   ;...
                                -1	;...
                                -1	];
    else
        EXPOENTES           = [ 1   ;...
                                0   ;...
                                1 	;...
                                -1	;...
                                -1	];
    end
    
    QUIM(4,1)   = EQUACOES_APP(REAGENTES_EQUACAO(:,1), EXPOENTES);
    QUIM(4,2)   = EQUACOES_APP(REAGENTES_EQUACAO(:,2), EXPOENTES);

    %   Ciclo para encontrar uma solução
    while (((DIST_ABS(QUIM(3,1),QUIM(3,2)) > 10^(-14))&&(DIST_REL(QUIM(3,1),QUIM(3,2)) > 10^(-5))) ||...
           ((DIST_ABS(QUIM(2,1),QUIM(2,2)) > 10^(-14))&&(DIST_REL(QUIM(2,1),QUIM(2,2)) > 10^(-5))) ||...
           ((DIST_ABS(QUIM(1,1),QUIM(1,2)) > 10^(-14))&&(DIST_REL(QUIM(1,1),QUIM(1,2)) > 10^(-5))))&&...
            COUNTER < 400

        QUIM_Meio(1)	= (QUIM(1,1)+QUIM(1,2))/2;
        if REACAO_2 == 1
            REAGENTES       = [ R1                  ;...
                                R0                  ;...
                                OH                  ;...
                                H3O + QUIM_Meio(1)  ];
        else
            REAGENTES       = [ R1                  ;...
                                R0                  ;...
                                OH  + QUIM_Meio(1)  ;...
                                H3O                 ];
        end
        HOLDER          = AcidoFraco_APP(app, REAGENTES);
        QUIM_Meio(2:3)  = HOLDER.INCREMENTOS;

        REAGENTES_EQUACAO 	= [ QUIM_Meio(1) + P0                               ;...
                                QUIM_Meio(1) + QUIM_Meio(2) + QUIM_Meio(3) + H3O;...
                                QUIM_Meio(1) + QUIM_Meio(3) + OH                ;...
                                P1 - QUIM_Meio(1)                               ;...
                                k1                                              ];

        QUIM_Meio(4)    = EQUACOES_APP(REAGENTES_EQUACAO, EXPOENTES);
            
        QUIM = Bissecao_APP(QUIM, QUIM_Meio);
        %   LOOP COUNTER
        COUNTER = COUNTER + 1;
    end

    D.ERROS         = [ DIST_ABS(QUIM(1,1),QUIM(1,2))   ;...
                        DIST_ABS(QUIM(2,1),QUIM(2,2))   ;...
                        DIST_ABS(QUIM(3,1),QUIM(3,2))   ];
    D.INCREMENTOS   = QUIM(1:(size(QUIM,1)-1), 1);

end