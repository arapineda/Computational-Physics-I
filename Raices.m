% Raices: se encuentra la x donde f(x)=0.
% Se despliega una matriz X donde el primer rengl�n es el valor de x y el
% segundo rengl�n es la funci�n evaluada donde 

x1=linspace(0,2*pi,1000); % el vector x donde se define la funci�n a evaluar
y1=cos(x1); % La funci�n donde se busca el cero. Si se busca un minimo o un m�ximo, esta funci�n debe de ser la derivada, en este caso, se pueden quitar los comentarios de F para evaluar dicha funci�n donde se encuentra un minimo y un m�ximo.
%F=F % Funci�n a maximizar o minimizar.
[h,l]=size(y1);
i=1;
x=1;
ix=1;
while i<l;
    if y1(i)>0;
        while y1(i)>0 & i<l;
            i=i+1;
        end
        if i~=l
            x(1,ix)=x1(i);
            %x(2,ix)=F(i);
            ix=ix+1;
        end
    else y1(i)<0;
        while y1(i)<0 & i<l;
            i=i+1;
        end
        if i~=l
            x(1,ix)=x1(i);
            %x(2,ix)=F(i);
            ix=ix+1;
        end
    end
end
x % Se dexpliega x