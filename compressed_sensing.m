function cs = compressed_sensing(n,m, type, initial_sparsity, final_sparsity, step_size, figure_num)
    disp("Figure number "+string(figure_num) +": " +type+ " matrix");
    gauss_recovery = [];
    rdft_recovery = [];
    bch_recovery = [];

    if type=="Singer"   
        bch_row = 840;
        other_row = 820;
        singer_recovery = [];   
        [singer, unused] = generate_singer(n, 11, 1.5, 0.5);
        singer = singer(1:other_row,:);
        plot_x = 1:1:final_sparsity;
        p_ary = 29;
    else 
        bch_row = 124;
        other_row = 132;
        singer_recovery = [];
        [singer, unused] = generate_singer(n, 11, 1.5, 0.5);
        singer = singer(1:other_row,:);
        plot_x = 1:1:final_sparsity;
        p_ary = 5;
    end
    
    % Constructing singer and BCH matrix without varying sparsity as it
    % takes ages to compute if iterated
%     bch = gen_bch_matrix(n, initial_sparsity, p_ary);
%     bch = bch(1:bch_row,:);
    % Test for different sparsity levels
    for s = initial_sparsity:step_size:final_sparsity
        disp("Running tests for sparsity of " + s);
        % Setting the number of experiments for each sparsity order
        % This is 1000 in the paper, but that would take a lot of computational
        % power and lots of time. Reducing it to 10 for this project
        iterations = 5;
    
        % Creating arrays to store percentage recovery for different sensing
        % matrices
        gaussian_avg_recovery = [];
        rdft_avg_recovery = [];
        bch_avg_recovery = [];
        singer_avg_recovery = [];
    
        for tests = 1:iterations
    
            % ******************Input signal******************
            % Generating a sparse signal with s randomized values.
            sel = randperm(n); sel = sel(1:s);
            signal = zeros(n,m); signal(sel)=1;
            % Randomization of the signs and values
            signal = signal.*sign(randn(n,1)).*(1-.5*rand(n,1));
            p = other_row;
            P0 = mean(signal.^2);
            SNR = 10^(30/10);
            noisy_signal = signal + sqrt(P0/SNR)*randn(size(signal));
    
    
            % ****************Complex Gaussian****************
            % Complex Gaussian
            gaussian = 1/sqrt(2)*(rand(p, n) +1i*rand(p,n));
            % normalization
            gaussian = gaussian ./ repmat( sqrt(sum(gaussian.^2)), [p 1] );
            % Recovery
            gaussian_percentage = sensing_matrix_method(gaussian, signal, noisy_signal, s, n);
            gaussian_avg_recovery(tests) = gaussian_percentage;
    
    
            % **********************RDFT**********************
            % Creating a nxn DFT matrix
            rdft = dftmtx(n);
            % Randomly choosing m rows from the DFT matrix to form the RDFT
            % matrix
            k = randperm(n);
            rdft = rdft(k(1:p),:);
            % Recovery
            rdft_percentage = sensing_matrix_method(rdft, signal, noisy_signal, s, n);
            rdft_avg_recovery(tests) = rdft_percentage;

            % **********************BCH***********************
            bch = gen_bch_matrix(n, s, p_ary);
            bch = bch(1:bch_row,:);
            bch_percentage = sensing_matrix_method(bch, signal, noisy_signal, s, n);
            bch_avg_recovery(tests) = bch_percentage;

            % *****************Singer matrix******************
            singer_percentage = sensing_matrix_method(singer, signal, noisy_signal, s, n);
            singer_avg_recovery(tests) = singer_percentage;
        end
        gauss_recovery(s) = mean(gaussian_avg_recovery);
        rdft_recovery(s) = mean(rdft_avg_recovery);
        bch_recovery(s) = mean(bch_avg_recovery);
        singer_recovery(s) = mean(singer_avg_recovery);
        disp("The recovery percentage for gaussian is : " + gauss_recovery(s));
        disp("The recovery percentage for rdft is : " + rdft_recovery(s));
        disp("The recovery percentage for bch is : " + bch_recovery(s));
        disp("The recovery percentage for singer is : " + singer_recovery(s));
    end
    % addressing logical array of nonzero elements
    isNZ=(~gauss_recovery==0);
    isNZ2=(~rdft_recovery==0);
    isNZ3=(~bch_recovery==0);
    isNZ4=(~singer_recovery==0);
    figure(figure_num);
    % Plotting gauss recovery
    plot(plot_x(isNZ), gauss_recovery(isNZ), '-*', 'Color', [1, 0.4, 0.7],'LineWidth',1.5);
    % Plotting rdft recovery
    hold on;
    plot(plot_x(isNZ2), rdft_recovery(isNZ2),  '-x', 'Color', [0, 1, 0],'LineWidth',1.5);
    % Plotting bch recovery
    hold on;
    plot(plot_x(isNZ3), bch_recovery(isNZ3), '-s', 'Color', [1, 0, 0],'LineWidth', 1.5);
    % Plotting singer recovery
    hold on;
    plot(plot_x(isNZ4), singer_recovery(isNZ4), '-^', 'Color', [0, 0, 1],'LineWidth', 1.5);
    % Formatting the plot
    xlim([initial_sparsity final_sparsity])
    ylim([0 1])
    xlabel('Sparsity k')
    ylabel('Recovery percentage')
    legend('Complex Gaussian', 'RDFT', string(p_ary) + '-ary BCH', type)
    grid on
    title('Perfect recovery percentage of noiseless' + string(n) + '× 1 signals')
end