function WT = ForkJoinMGM(Lambda, mu, M)
debug = 0; 
total_debug = 0;
%��������� ���������
S = states(1,M,M);
if debug
    disp('S');
    disp(S); 
end

A0 = MatrixA0(M, mu);
if(debug)
    disp('������� A0 = ');
    disp(size(A0));
    disp(A0);
    xlswrite('E:\Mess\A0.xls',A0);
    pause
end

A2 = MatrixA2(M, Lambda);
if(debug)
    disp('������� A2 = ');
    disp(size(A2));
    disp(A2);
    xlswrite('E:\Mess\A2.xls',A2);

    pause; 
end;


a1 = (A0 + A2)*ones(M^M, 1);
A1 = -diag(a1);
if(debug)
    disp('������� A1 = ');
    disp(size(A1));
    disp(A1);
    xlswrite('E:\Mess\A1.xls',A1);

    pause;
end



%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%������������ ������ Bxy 
S0 = states(0, M, M);
if(debug)
    disp('��������� ��������� S0 (n0=0) ������� ������������');
    for i = 1:(M+1)^M
        disp(i)
        disp(S0(:, i));
        pause
    end

end


temp = size(S0);
N = temp(2);
if(debug)
    disp('����� ���������');
    disp(N);
end 


%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%���������� ������� B00 (������ �������� ������)
B00 = zeros(N,N);

%������� ��� ����������� 
for i = 1:N
    %������ �� ������� ���������
    s1 = S0(:,i);
    %������� ��� ����������� ���������� 
    %(���� ���� ������� ����������, �� �������� �� ��� �� ������)
    zero_coords = find(s1==0); 
    if ~isempty(zero_coords)
        temp = size(zero_coords);
        
        d = temp(1); 
        s2 = s1; 
        s2(zero_coords) = d;
        
        
        %����� ������ ������ ���������
        j = find(ismember(S0', s2', 'rows'));
        B00(i,j) = B00(i,j) +  Lambda; 
        
        if(total_debug)
        s1
        i
        s2
        j
        disp('arrive')
        pause
        end
    end
end

%������� ��� ������������
for i = 1:N
    %������ �� ������� ���������
    s1 = S0(:,i);
    %������ �� ������ ����������
    for ni = 1:M
        if(s1(ni)>0)
            s2 = s1; 
            s2(ni) = 0; 
            %����� ������ ������ ���������
            j = find(ismember(S0', s2', 'rows'));
            B00(i,j) = B00(i,j) +  mu(ni)*s1(ni); 
            if(total_debug)
                    s1
                    i
                    s2
                    j
                disp(mu(ni)*s1(ni))
                disp('leave')
                pause
            end
        end
    end
end
%������� B00 ������������ ��� �� ���������, 
%���������� ��������� B01

%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

%���������� ������� B01 (������� � �������� ������ �� ������)
B01 = zeros(N,M^M);

%������� �������� ������ ��� ����������� ���
%������ ��������� �������
for i = 1:N
    %������ �� ������� ���������
    s1 = S0(:,i);
    zero_coords = find(s1==0); 
    if isempty(zero_coords)
        %������� � �� �� ��������� (�� ����� ��� ����� ������)
        %��������� ��������� �� S(n0>0)
        j = find(ismember(S', s1', 'rows'));
        B01(i,j) = B01(i,j) +  Lambda; 
    end
end

%������� B01 ������������ ��������� 


%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%����� ���� ��� ������������ ������� B01 
%����� ������������� ������� B00 (������������ ��������)
b00 = B01*ones(M^M,1) + B00*ones(N,1);
B00 = B00 - diag(b00); 

if(debug)
    disp('������� B00 = ')
    disp(B00);
    xlswrite('E:\Mess\B00.xls',B00);

    disp('������� B01 = ')
    disp(B01);
    xlswrite('E:\Mess\B01.xls',B01);
end


%S0
%pause

%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%������������ ������� B10 (������� � ������� ������ � �������)
B10 = zeros(M^M, N); 
%�������� �� ���� ���������� 
for i = 1:M^M
  %�������� �� ����������� ��������� ��������
  for ni = 1:M
  %��� ������ ���������� ���� ����������� �������� � ������ ��������� 
  s1 = S(:, i);
  s2=[];
  s2=s1; 
  s2(ni) = 1;
  j = find(ismember(S0', s2', 'rows'));

  B10(i,j) = B10(i,j)+ mu(ni)*s1(ni);
  end
end

 if(debug)
     disp('������� B10 = ');
     disp(B10);
     xlswrite('E:\Mess\B10.xls',B10);

 end 


%���������� ��������-��������������� ������ 
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%MGM Begin
eps = 10^(-20); 
R = A0*0; 
V = A2/A1;
W = A0/A1;
while(1)
    Rnext = -V - R*R*W; 
    if(norm(Rnext-R)<eps) 
        break; 
    end
    R = Rnext;
end

if(debug)
    disp('R = ');
    disp(R);
end

if(debug)
    disp('R^-1 = ');
    disp(inv(R));
end

%MGM End

Q = [
        B00, B01;
        B10, A1+ R*A0
    ];
 
if(debug)
    disp('Q = ');
    disp(Q);
end;


%Replacing the last equation with pi0_1 = 1
temp = size(Q);
Q(:, temp(1)) = [1; zeros(temp(1)-1, 1)];
%disp(Q);

b = [zeros(1, temp(1)-1) 1];

pi01 = b/Q;


temp = size(B00);
n = temp(1);
pi0 = pi01(1:n);
temp = size(pi01);
pi1 = pi01(n+1:temp(2));

if(debug)
    disp('pi0');
    disp(pi0);

    disp('pi1');
    disp(pi1);
end

%The normalization constant 
alpha = sum(pi0) + sum(pi1/(eye(size(R)) - R));

if(debug)
    disp('alpha');
    disp(alpha);
end

pi0 = pi0/alpha;
pi1 = pi1/alpha;

if(debug)
%After mormalization 
disp('pi0');
disp(pi0);

disp('pi1');
disp(pi1);



pi2 = pi1*R;
disp('pi2')
disp(pi2)

pi3 = pi2*R;
disp('pi3')
disp(pi3)

pi4 = pi3*R;
disp('pi4')
disp(pi4)


pi5 = pi4*R;
disp('pi5')
disp(pi5)


pi6 = pi5*R
disp('pi4')
disp(pi6)


disp('PI0')
disp(sum(pi0))


disp('PI1')
disp(sum(pi1))


disp('PI2')
disp(sum(pi2))

disp('PI3')
disp(sum(pi3))

disp('PI4')
disp(sum(pi4))

disp('PI5')
disp(sum(pi5))

disp('PI6')
disp(sum(pi6))

end


temp = size(R); 
disp('Queue length in the network '); 
queue_length = sum((pi1 - pi1*(R - 2*eye(size(R)))*R/((R-eye(size(R)))^2)));

%queue_length = sum(pi1*(eye(size(R))-R)^(-2));


W = queue_length/Lambda

if(debug)
    MMkappa(Lambda, mu, M)
end









%**************** ���������� �.�. ������������ ���������� 
%�������� �� ��������������� �������������� ��������� ������� 
V2=0;
for i = 1:N
    %������ �� ������� ���������
    s0 = S0(:,i);
    d = sum(s0==0); 
    if(d>0)
        vars = mu(s0==0)*d; 
        m = ExpectMaxExp(vars); 
        V2 = V2 + pi0(i)*m; 
    end
end

A = MatrixServer(M,mu); 

%�.�. ������������ ������������ ��� ������� ���������� � ������� ���
%������������� ��������� 
v = []; 
for i = 1:M^M
    s1 = S(:, i);
    vi = 0;
    for ni = 1:M
        alphai = mu(ni)*s1(ni)/(mu*s1); 
        vi = vi + alphai/mu(ni); 
    end
    v = [v; vi]; 
end 

 
V1 = pi1*sylvester(inv(R), -A, inv(R))*A*v
eig(inv(R))
eig(-A)
%pause;


%��� ������ ������� � ����������� �������� 
for i = 1:N
    %������ �� ������� ���������
    s0 = S0(:,i);
    d = sum(s0==0); 
    if(d==0)
        vi = 0;
        for ni = 1:M
            alphai = mu(ni)*s0(ni)/(mu*s0); 
            vi = vi + alphai/mu(ni); 
        end
        j = find(ismember(S0', s0', 'rows'));
        prob = pi0(j); 
        V1 = V1 + vi*prob; 

    end
end 
    
disp('V1 V2');     
V1
V2
V = V1 + V2;

disp('T = W + V'); 
T = W+V;
disp(T);

WT = [W,T]
disp('heterogeneous case')


end

