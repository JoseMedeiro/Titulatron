%   K H2O
KH2O = 10^(-14);

%   Acido Forte
ACIDOINICIAL = 10^(-2);

OH  = 10^(-7);
H3O = 10^(-7)+ ACIDOINICIAL;

x = -(H3O + OH)/2 + sqrt((H3O + OH)^2-4*(H3O*OH-KH2O))/2;

H3OInicial = H3O + x;
OHInicial  = OH  + x;

%   Agora vamos introduzindo o titulante (sempre forte)

z  = 0:0.1*10^(-3):20*10^(-3);

%   Formula deduzida

x = -(H3OInicial + (OHInicial+z))/2 + sqrt((H3OInicial + (OHInicial+z)).^2-4*(H3OInicial*(OHInicial+z)-KH2O))/2;

pH = - log(H3OInicial+x)/log(10);

%Experimental

a  = size(z);
b  = a(1,2);

ph = zeros(2:a);
%   Método da bisseção
for c = 1:b
    
    y = AcidoForteEX(z(c)+OHInicial, H3OInicial, KH2O);
    
    %   Troco o OH pelo H3O para evitar ao máximo cancelamento subtrativo
    %   (mas não afetava muito)
    ph(c) = -log(H3OInicial+y)/log(10);
    TEXT = [ num2str(c) ' de 201 feitos'];
    disp(TEXT)
    
end

plot(z,pH,'c',z,ph,'r')

disp("Pressione uma tecla para continuar")
pause;
disp(" ")

%   Vamos tentar fazer um acido fraco
ACIDOINICIAL = 5*10^(-3);
k1 = 10^(-3);

OH  = 10^(-7);
H3O = 10^(-7);

ph = zeros(2:a);
CONC = zeros(2:a);

H3OInicial = H3O;
OHInicial  = OH;
yy = zeros(1:4);
y = zeros(1:2);
w = zeros(1:2);

for c = 1:b
    
    d = 0;
    %   Fase Inicial
    yy = AcidoFracoEX_Teste(OHInicial+z(c), H3OInicial, ACIDOINICIAL, 0, k1, KH2O);
    y(1) = yy(1);
    y(2) = yy(2);
    w(1) = yy(3);
    w(2) = yy(4);
    
    ph(:,c) = 14+log(OHInicial+z(c)+y)/log(10);
    CONC(:,c) = w;
    TEXT = [ num2str(c) ' de 201 feitos'];
    disp(TEXT)
    
end

plot(z,ph(1,:)-ph(2,:))
disp("Pressione uma tecla para continuar Next")
pause;
plot(z,CONC(2,:),'r')


disp("Pressione uma tecla para continuar AA")
pause;
disp(" ")

%   Vamos tentar fazer um ácido dipróticooooooo
ACIDOINICIAL = 10^(-2)/2;
k1 = 10^(-3);
k2 = 10^(-9);

OH  = 10^(-7);
H3O = 10^(-7);

H3OInicial = H3O;
OHInicial  = OH;

ph2 = zeros(2:a);
x = zeros(1:2);
y = zeros(1:2);
w = zeros(1:2);
holder = zeros(1:2);

RH2Equilibrio = ph2;
RHEquilibrio = ph2;
REquilibrio = ph2;
H3OEquilibrio = ph2;
%   y é para a reação agua
%   x é para a reação 1
%   w é para a reação 2

for c = 1:b
    
    %   Fase Inicial
    
    w = [0 ACIDOINICIAL-2^(-25)];
    
    holder = AcidoFracoEX(OHInicial+z(c), H3OInicial+w(1), ACIDOINICIAL, -w(1), k1, KH2O);
    y(1) = holder(1);
    x(1) = holder(2);
    holder = AcidoFracoEX(OHInicial+z(c), H3OInicial+w(2), ACIDOINICIAL, -w(2), k1, KH2O);
    y(2) = holder(1);
    x(2) = holder(2);
    
    
    ERR = (w).*(w + x + y + H3OInicial)./(k2*(x-w)) - 1;
    d = 0;
    
    %   Ciclo para encontrar uma solução
    while ((abs(y(1)-y(2)) > 10^(-14))&&(abs((y(1)-y(2))/(y(1)+y(2))*2) > 10^(-5)))||...
          ((abs(w(1)-w(2)) > 10^(-14))&&(abs((w(1)-w(2))/(w(1)+w(2))*2) > 10^(-5)))||...
          ((abs(x(1)-x(2)) > 10^( -8))&&(abs((x(1)-x(2))/(x(1)+x(2))*2) > 10^(-6)))
        
        wMeio = (w(1)+w(2))/2;
        holder = AcidoFracoEX(OHInicial+z(c), H3OInicial+wMeio, ACIDOINICIAL, -wMeio, k1, KH2O);
        yMeio = holder(1);
        xMeio = holder(2);
        ERRMeio = (wMeio)*(wMeio + xMeio + yMeio + H3OInicial)./(k2*(xMeio-wMeio)) - 1;
        
        if ERRMeio < 0
            if ERR(1) < 0
                w(1)   = wMeio;
                x(1)   = xMeio;
                y(1)   = yMeio;
                ERR(1) = ERRMeio;
                FLAG   = 1;
            else
                w(2)   = wMeio;
                x(2)   = xMeio;
                y(2)   = yMeio;
                ERR(2) = ERRMeio;
                FLAG   = 2;
            end
        else
            if ERR(1) > 0
                w(1)   = wMeio;
                x(1)   = xMeio;
                y(1)   = yMeio;
                ERR(1) = ERRMeio;
                FLAG   = 1;
            else
                w(2)   = wMeio;
                x(2)   = xMeio;
                y(2)   = yMeio;
                ERR(2) = ERRMeio;
                FLAG   = 2;
            end
        end
        d = d+1;
        
    end
    
    H3OEquilibrio(:,c)    = y;
    
    RH2Equilibrio(:,c)    = ACIDOINICIAL - x;
    RHEquilibrio(:,c)     = x-w;
    REquilibrio(:,c)      = w;
    
    ph2(:,c) = 14+log(OHInicial+z(c)+y)/log(10);
    TEXT = [ num2str(c) ' de 201 feitos'];
    disp(TEXT)
    
end
hold on;
plot(z,ph2(1,:),'r')
plot(z,ph2(2,:),'c')
hold off;
disp("Pressione uma tecla para continuar Next")
pause;
hold on;
plot(z,abs(H3OEquilibrio(2,:)-H3OEquilibrio(1,:)),'black')
plot(z,abs(RH2Equilibrio(2,:)-RH2Equilibrio(1,:)),'r')
%plot(z,abs(RHEquilibrio(2,:)-RHEquilibrio(1,:)),'c')
plot(z,abs(REquilibrio(2,:)-REquilibrio(1,:)),'g')
hold off;
disp("Pressione uma tecla para continuar Next")
pause;
hold on;
plot(z,RH2Equilibrio(1,:),'r')
plot(z,RHEquilibrio(1,:),'c')
plot(z,REquilibrio(1,:),'g')
hold off;