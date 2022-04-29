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
function D = AcidoDiprotico_APP(app, REAGENTES_INICIAIS)
    
    %Retirada de variaveis
    R2  = REAGENTES_INICIAIS(1);
    R1  = REAGENTES_INICIAIS(2);
    R0  = REAGENTES_INICIAIS(3);
    OH  = REAGENTES_INICIAIS(4);
    H3O = REAGENTES_INICIAIS(5);
    
    % Preparacao dos dados
    k2 = app.K12.Value;
    QUIM        = zeros(4,2);
    QUIM_Meio   = zeros(4,1);
    COUNTER     = 0;
    %   Fase Inicial
    QUIM(1,:)   = [-R0, R2 + R1 - 2^(-25)];
    % Para resolver as equações abaixo
   	REAGENTES   = [ R2  , R2        ;...
           	        R1  - QUIM(1,:)	;...
              	    OH  , OH        ;...
                   	H3O + QUIM(1,:)	];
    % Resolver as equações abaixo
    HOLDER      = AcidoFraco_APP(app, REAGENTES(:,1));
    QUIM(2:3,1) = HOLDER.INCREMENTOS;
    HOLDER      = AcidoFraco_APP(app, REAGENTES(:,2));
    QUIM(2:3,2) = HOLDER.INCREMENTOS;
    
    REAGENTES_EQUACAO   = [ QUIM(1,:) + R0                          ;...
                            QUIM(1,:) + QUIM(2,:) + QUIM(3,:) + H3O	;...
                            R1 + QUIM(2,:) - QUIM(1,:)              ;...
                            k2, k2                                  ];
    EXPOENTES           = [ 1   ;...
                            1 	;...
                            -1	;...
                            -1	];
    
    QUIM(4,1)   = EQUACOES_APP(REAGENTES_EQUACAO(:,1), EXPOENTES);
    QUIM(4,2)   = EQUACOES_APP(REAGENTES_EQUACAO(:,2), EXPOENTES);

    %   Ciclo para encontrar uma solução
    while (((DIST_ABS(QUIM(3,1),QUIM(3,2)) > 10^(-14))&&(DIST_REL(QUIM(3,1),QUIM(3,2)) > 10^(-5))) ||...
           ((DIST_ABS(QUIM(2,1),QUIM(2,2)) > 10^(-14))&&(DIST_REL(QUIM(2,1),QUIM(2,2)) > 10^(-5))) ||...
           ((DIST_ABS(QUIM(1,1),QUIM(1,2)) > 10^(-14))&&(DIST_REL(QUIM(1,1),QUIM(1,2)) > 10^(-5))))&&...
            COUNTER < 400

        QUIM_Meio(1)	= (QUIM(1,1)+QUIM(1,2))/2;
        
        REAGENTES       = [ R2                  ;...
                            R1 - QUIM_Meio(1)   ;...
                            OH                  ;...
                            H3O + QUIM_Meio(1)  ];
        HOLDER          = AcidoFraco_APP(app, REAGENTES);
        QUIM_Meio(2:3)  = HOLDER.INCREMENTOS;
        
        REAGENTES_EQUACAO   = [ QUIM_Meio(1) + R0                                   ;...
                                QUIM_Meio(1) + QUIM_Meio(2) + QUIM_Meio(3) + H3O 	;...
                                R1 + QUIM_Meio(2) - QUIM_Meio(1)                    ;...
                                k2                                              	];

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