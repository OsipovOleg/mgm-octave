delete all
clc


TableW = []; 
TableT = [];
TablePSI = [];
TableLAMBDA = []; 


PSI0 = 0.1; 
PSI1 = 0.9;


Mmin =2; 
Mmax=4; 

MU = [1, 1, 1, 1 ]
for M = Mmin:Mmax
    mu = MU(1:M)
    %����� ����� PSI(i)��� ���������� ����� �����
    sample = 100; 
    %������������� ��������� ������
    LAMBDA= []; 
    %�-�� ������������� �������
    PSI = []
    %������������ ���������� ���������� � �������
    W = [];
    %������������ ���������� ���������� � ������� ������������
    T =[]; 
    
    %��� ��������� ������������� ��������� ������
    lambda = PSI0*sum(mu); 
    h = (sum(mu)*PSI1-lambda)/(sample-1);

    %���������� ����� ����� 
    for i = [1:sample]
        temp = lambda+h*(i-1);
        LAMBDA = [LAMBDA, temp]; 
        PSI = [PSI, temp/(sum(mu))]; 
        
        %���������� ������������� ���������� � ������� � � �������
        wt = ForkJoinMGM(temp, mu, M);
        W=[W, wt(1)]; 
        T=[T, wt(2)]; 
    end

    
    grid on;
    hold on;
    plot (PSI, T);

   
        
    TableW = [TableW; [M, W]];
    TableT = [TableT; [M, T]];
    TablePSI = [TablePSI; [M, PSI]];
    TableLAMBDA = [TableLAMBDA; [M, LAMBDA]];
end;

%���������� � ����
    xlswrite('E:\Mess\TableW.xls',TableW);
    xlswrite('E:\Mess\TableT.xls',TableT);
    xlswrite('E:\Mess\TablePSI.xls',TablePSI);
    xlswrite('E:\Mess\TableLAMBDA.xls',TableLAMBDA);
    
    
    
    for i=Mmin:Mmax
        j = i - Mmin + 1; 
        x = (TablePSI(j, 2:sample+1))';
        y = (TableT(j, 2:sample+1))'; 
        dlmwrite(strcat(strcat('E:\Mess\Example', num2str(i)), '.csv'),[x,y],' ');
    end

    
    
    

  


