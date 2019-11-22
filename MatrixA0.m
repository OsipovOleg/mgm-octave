function A0 = MatrixA0(M, mu)
A0 = zeros(M^M, M^M); 
%—генерировали все состо€ни€
S = states(1, M, M);
%ѕроходим по всем состо€ни€м 
for i = 1:M^M
  %ѕроходим по координатам состо€ни€ текущего
  for ni = 1:M
  %ƒл€ каждой координаты есть возможность перехода в другое состо€ние с
  %единичной координатой
  s1 = S(:, i);
  s2=[];
  s2=s1; 
  s2(ni) = 1;
  
  %j = vectorfind(S,s2, 'c');
  j = find(ismember(S', s2', 'rows'));
  
  A0(i,j) = A0(i,j) +  mu(ni)*s1(ni);

  end
 end
 
end
