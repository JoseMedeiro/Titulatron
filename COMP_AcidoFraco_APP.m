%   Função de Ácido Fraco - Compilação
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
function D = COMP_AcidoFraco_APP(app, QuantidadeInicial)
            
    ponto_titulacao = Ponto_Equilibrio_Simples(app);
    a = size(ponto_titulacao,2);

    z  = 0:(ponto_titulacao(a)*(1+1/a))/(10^(4)):ponto_titulacao(a)*(1+1/a);
    Volume = Volume_Inicial(app) + z;
    Titulante = z*Concentracao_Titulante(app);

    ph = zeros(size(z));
    RH = zeros(size(z));
    R  = zeros(size(z));
    
    INCREMENTOS = zeros(2,size(z,2));
    ERROS       = zeros(2,size(z,2));
    OHInicial = 10^(-7);
    H3OInicial = 10^(-7);
            
    for c = 1:size(z,2)
        if ~Cancelar(app)
            break;
        end
        OH          = OHInicial + Titulante(c)/Volume(c);
        H3O         = H3OInicial;
        RHInicial   = QuantidadeInicial/Volume(c);
        RInicial    = 0;
        
        REAGENTES   = [ RHInicial;  ...
                        RInicial;   ...
                        OH;         ...
                        H3O         ];
                
        y = AcidoFraco_APP(app, REAGENTES);
        
        ph(c) = -log(H3O+y.INCREMENTOS(1)+y.INCREMENTOS(2))/log(10);
        RH(c) = RHInicial-y.INCREMENTOS(1);
        R(c)  = y.INCREMENTOS(1);
        INCREMENTOS(:,c)  = y.INCREMENTOS;
        ERROS(:,c)        = y.ERROS;
        
        app.Percentagem.Value = c/size(z,2)*100;
        if mod(c/size(z,2)*100,5)<mod((c-1)/size(z,2)*100,5)
            Tempo_Restante(app, c/size(z,2));
        end
        drawnow;
        
    end
    if strcmp(app.TitulanteSwitch.Value,'Base')
        D.pH.DADOS              = ph;
        D.CONCENTRACAO.LEGENDAS = ["R"; "RH"];
    else
        D.pH.DADOS              = 14 - ph;
        D.CONCENTRACAO.LEGENDAS = ["R"; "R(OH)"];
    end
    D.ABSISSA               = z;
    D.pH.LEGENDAS          	= 'pH';
  	D.CONCENTRACAO.DADOS    = [R; RH];
    D.INCREMENTOS.DADOS     = INCREMENTOS;
    D.INCREMENTOS.LEGENDAS 	= ["x_1"; "x_2"];
    D.ERROS.DADOS           = ERROS;
    D.ERROS.LEGENDAS        = D.INCREMENTOS.LEGENDAS;
    
end