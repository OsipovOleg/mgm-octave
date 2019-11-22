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
    %Число точек PSI(i)для построения одной линии
    sample = 100; 
    %Интенсивности входящего потока
    LAMBDA= []; 
    %К-ты использования системы
    PSI = []
    %Длительности пребывания требований в очереди
    W = [];
    %Длительности пребывания требований в системе обслуживания
    T =[]; 
    
    %Шаг изменения интенсивности входящего потока
    lambda = PSI0*sum(mu); 
    h = (sum(mu)*PSI1-lambda)/(sample-1);

    %Построение одной линии 
    for i = [1:sample]
        temp = lambda+h*(i-1);
        LAMBDA = [LAMBDA, temp]; 
        PSI = [PSI, temp/(sum(mu))]; 
        
        %Вычисление длиетльностей пребывания в очереди и в системе
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

%Сохранение в файл
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

    
    
    

  


