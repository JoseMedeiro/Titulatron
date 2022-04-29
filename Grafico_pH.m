%   Funcao do Grafico - 1 Acido Diprotico
%%  INPUT
%
%   P               - Vetor com os Pontos (P; M; Q)
%   Arco            - linha com todos os pontos equidistantes
%   RESOLUCAO_TERRA - Resolucao da Longitude = 2*Resolucao da Latitude
%   UIAxes          - Eixos para o gráfico
%
%%  OPERACAO
%
%   Cria um Globo;
%   
%   Plota:
%       O Globo;
%       Os Pontos e a Linha Equidistante, por ordem (LE1, P1, M1, Q1,...),
%       criando simultaneamente as legendas para as Trajetorias e Pontos
%       (A CONSIDERAR SE QUEREMOS FAZER MAIS DO QUE UMA LINHA)
%       Condensa todas os plots a legendar para um vetor linha e as 
%       legendas correspondentes para uma celula linha;
%       Define a vista inicial do plot;
%   Acaba o Plot.
%
%%  DEPENDENCIAS
%
%   Grafico_Globo
%   Grafico_Equador
%   Grafico_Meridiano
%
%%  OUTPUT
%   
%   Nenhum (ignorar)
%
%%
function grafequi = Grafico_pH(app, UIAxes, Values)
    
    % R Plot
    Plot_Legenda = plot(UIAxes, Values.ABSISSA(1,:)*10^(3), Values.pH.DADOS);
    hold(UIAxes,'on');
    
    %Plot legends
    title(UIAxes, app.title);
    xlabel(UIAxes, 'Titulante adicionado (mL)');
    ylabel(UIAxes, 'pH');
    % Legenda
    Celula_Texto =  cellstr(Values.pH.LEGENDAS); 
    legend(UIAxes,Plot_Legenda, Celula_Texto', 'Location', 'northeast');
        
    hold(UIAxes,'off');
    
    grafequi = 1;
end