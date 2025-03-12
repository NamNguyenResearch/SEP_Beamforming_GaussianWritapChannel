function s_hat=ML_Detection_QPSK(r,H,W)

s_candidate_set=[]; % Define the candidate set
dist_set=[]; % Define the distance set

for b1_candidate=[-1-j, -1+j, 1-j, 1+j] % Candidate of antenna 1 includes {-1-j, -1+j, 1-j, 1+j} for QPSK modulation
    for b2_candidate=[-1-j, -1+j, 1-j, 1+j] % Candidate of antenna 2 includes {-1-j, -1+j, 1-j, 1+j} for QPSK modulation
        s_candidate=[b1_candidate b2_candidate]'; % Transimitted signal candidate
        r_candidate=H*W*s_candidate; % Received signal candidate
        dist= sum(abs(r-r_candidate).^2); % Distance between recieved signal and received signal candidate
                
        s_candidate_set=[s_candidate_set s_candidate]; %Add s_candidate as a new column of s_candidate_set
        dist_set=[dist_set dist]; %Add dist as an new element of dist_se   
    end
end

[A B]=min(dist_set); % Find the received signal candidate has the minimum distance with received signal
s_hat=s_candidate_set(:,B); % Detect that corresponding transmitted signal candidate as the transmitted signal 

