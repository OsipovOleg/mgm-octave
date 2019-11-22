function MMkappa(lambda, mu, kappa)

    disp('M/M/kappa (start)')
    psi = lambda/(kappa*mu); 
    temp = 0; 
    for n = 0:kappa-1
        temp= temp + ((kappa*psi)^n)/factorial(n); 
    end
    p0= 1/(((kappa*psi)^kappa)/(factorial(kappa)*(1-psi)) + temp); 
    disp('p0 = ');
    disp (p0);
    
    b = p0*(kappa^kappa)*(psi^(kappa+1))/...
        (factorial(kappa)*(1-psi)*(1-psi));
    
    w = b/lambda
    disp('M/M/kappa (end)')
   
end
