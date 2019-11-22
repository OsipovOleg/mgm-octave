function P = MatrixServer(M, mu)
P = zeros(M^M, M^M); 
%������������� ��� ���������
S = states(1, M, M);
%�������� �� ���� ���������� 
for i = 1:M^M
  %�������� �� ����������� ��������� ��������
  for ni = 1:M
  %��� ������ ���������� ���� ����������� �������� � ������ ��������� �
  %��������� �����������
  s1 = S(:, i);
  s2=[];
  s2=s1; 
  s2(ni) = 1;
  
  %j = vectorfind(S,s2, 'c');
  j = find(ismember(S', s2', 'rows'));
  
  
  P(i,j) = P(i,j) +  mu(ni)*s1(ni)/(mu*s1); 

  end
 end
 
end
