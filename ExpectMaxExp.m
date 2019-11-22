function m = ExpectMaxExp(mu)
A = [-mu(1)];
alpha = [1];
temp  = size(mu);
n = temp(2); 

for i = 2:n
    B = [-mu(i)]; 
    betta = [1]; 
    [C, gamma] = MaxPH(A, B, alpha, betta);
    %Копирование параметров
    A = C; 
    alpha = gamma; 
end 
m = -alpha*inv(A)*ones(length(A), 1);
    
    
    