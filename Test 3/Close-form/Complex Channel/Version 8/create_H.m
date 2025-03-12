function H = create_H(K,N,var,scattering_num1,scattering_num2) 
    H = zeros(K,N);

    for k = 1:K 
        L_m = unifrnd(scattering_num1,scattering_num2);
        sum = 0;
            
        for i = 1:L_m 
            alpha_m = sqrt(var/2)*(randn(1,1)+1i*randn(1,1));
            teta_tm = unifrnd(-pi/2,pi/2);
            sum = sum + alpha_m*cal_a_t(teta_tm,N);
        end
      
        H(k,:) = sqrt(N/L_m)*sum;
    end
end