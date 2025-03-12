function a_t = cal_a_t(teta_tm,N) 
    a_t = [];

    for j = 1:N 
        a_t = [a_t exp(1i*2*pi*(j-1)*1/2*sin(teta_tm))];
    end
end