% a - ����������� �������� ����������, b - ������������ �������� ����������, 
%M - ����������� ������� 
%��������� ���� ��������� 

function space = states(a, b, M)
x = [a:1:b]; 
space= x; 
for i = [1:M-1]
space = cart(space, x);  
end
end
