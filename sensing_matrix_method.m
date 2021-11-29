% Function to recover the output matrix from the input sigal and the
% sensing matrix. 
function recovered_percentage = sensing_matrix_method(U, signal, noisy_signal, s, n)
    % Output signal
    y = U*noisy_signal;
    % Orthogonal matching pursuit algorithm
    recovered_signal = algo_omp(s, U, y);
    recovered_count = 0;
    total_signal_count = 0;
    for i = 1:n
        if signal(i,1)~=0 
            total_signal_count = total_signal_count + 1;
            if 30*abs(snr(signal(i,1), recovered_signal(i,1)))+20<100
                recovered_count = recovered_count+1;
            end
        end
    end
    recovered_percentage = (recovered_count/total_signal_count);
end