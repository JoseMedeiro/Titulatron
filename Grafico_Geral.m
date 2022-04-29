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
function grafequi = Grafico_Geral(app, UIAxes, Values)
    
    switch app.Grafico_Group.SelectedObject
        case app.Grafico_pH
            Grafico_pH(app, UIAxes, Values);
        case app.Grafico_Concentracao
            Grafico_Concentracao(app, UIAxes, Values);
        case app.Grafico_Reacoes
            Grafico_Incrementos(app, UIAxes, Values);
        case app.Grafico_ErroAbsoluto
            Grafico_Erros(app, UIAxes, Values);
    end
    
    
    grafequi = 1;
end