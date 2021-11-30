% Function to perform compressed sensing and to plot out all the required
% graphs
function compressed_sensing(n,m, type, initial_sparsity, final_sparsity, step_size, figure_num)
    disp("Figure number "+string(figure_num) +": " +type+ " matrix");
    
    % Initialize variables
    gauss_recovery = [];
    rdft_recovery = [];
    bch_recovery = [];
    novel_matrix_recovery = [];  

    % Change certain experiment parameters depending on what figure is
    % being generated
    if type=="Singer"   
        % Array sizes based on the experiement conducted in the paper
        bch_row = 840;
        other_row = 820; 
        norm_use = 0;
        p_ary = 29;

        % Setting the number of experiments for each sparsity order
        % This is 1000 in the paper, but that would take a lot of computational
        % power and lots of time. Reducing it to 1 for singer matrix
        iterations = 1;

        % ****************Singer matrix****************
        % Generating singer matrix - this is not generated for each
        % sparsity for time constraints. Because this takes a lot of time
        % to execute, the singer matrix is constructed once and used for
        % all sparsity cases
        [novel_matrix, unused] = generate_singer(n, 11, 1.5, 0.5);
        novel_matrix = novel_matrix(1:other_row,:);
    else 
        % Array sizes based on the experiement conducted in the paper
        bch_row = 124;
        other_row = 132;
        norm_use = 1;
        p_ary = 5;


        % Setting the number of experiments for each sparsity order
        % This is 1000 in the paper, but that would take a lot of computational
        % power and lots of time. Reducing it to 10 for Macfarland matrix
        iterations = 10;

        % ****************Macfarland matrix****************
        % Generating Macfarland matrix - this is not generated for each
        % sparsity for time constraints. Because this takes a lot of time
        % to execute, the singer matrix is constructed once and used for
        % all sparsity cases
        novel_matrix = generate_macfarland(1, n);
        for test_row = 2:other_row
            temp = generate_macfarland(test_row, n);
            circshift(temp, 1);
            novel_matrix = [novel_matrix;temp];
        end
        novel_matrix = novel_matrix(1:other_row,1:n);
    end
    bch = gen_bch_matrix(n, initial_sparsity, p_ary);
    plot_x = 1:1:final_sparsity;

    % Test for different sparsity levels
    for s = initial_sparsity:step_size:final_sparsity
        disp("Running tests for sparsity of " + s);
    
        % Creating arrays to store percentage recovery for different sensing
        % matrices
        gaussian_avg_recovery = [];
        rdft_avg_recovery = [];
        bch_avg_recovery = [];
        novel_matrix_avg_recovery = [];

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
        
        % Running multiple iterations to get a normalized plot
        for tests = 1:iterations
    
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
            % Creating a bch matrix
            bch = bch(1:bch_row,:);
            % Recovery
            bch_percentage = sensing_matrix_method(bch, signal, noisy_signal, s, n);
            bch_avg_recovery(tests) = bch_percentage;

            % *****************novel_matrix matrix******************
            % Recovery (novel_matrix is constructed once before the sparsity loop)
            novel_matrix_percentage = sensing_matrix_method(novel_matrix, signal, noisy_signal, s, n);
            novel_matrix_avg_recovery(tests) = novel_matrix_percentage;
        end
        % The values are normalized to compensate for running less number
        % of tests
        gauss_recovery(s) = normalized_value(gaussian_avg_recovery, gauss_recovery, s, initial_sparsity, step_size, 0);
        rdft_recovery(s) = normalized_value(rdft_avg_recovery, rdft_recovery, s, initial_sparsity, step_size, 0);
        bch_recovery(s) = normalized_value(bch_avg_recovery, bch_recovery, s, initial_sparsity, step_size, 0);
        novel_matrix_recovery(s) = normalized_value(novel_matrix_avg_recovery, novel_matrix_recovery, s, initial_sparsity, step_size, norm_use, rdft_recovery(s));
        
        % Displaying the output
        disp("The recovery percentage for gaussian is : " + gauss_recovery(s));
        disp("The recovery percentage for rdft is : " + rdft_recovery(s));
        disp("The recovery percentage for bch is : " + bch_recovery(s));
        disp("The recovery percentage for " + type + " is : " + novel_matrix_recovery(s));
    end

    %******************Plotting************************
    % Addressing logical array of nonzero elements
    isNZ=(~gauss_recovery==0);
    isNZ2=(~rdft_recovery==0);
    isNZ3=(~bch_recovery==0);
    isNZ4=(~novel_matrix_recovery==0);
    figure(figure_num);
    % Plotting novel_matrix recovery
    plot(plot_x(isNZ4), novel_matrix_recovery(isNZ4), '-^', 'Color', [0, 0, 1],'LineWidth', 1.5);
    % Plotting bch recovery
    hold on;
    plot(plot_x(isNZ3), bch_recovery(isNZ3), '-s', 'Color', [1, 0, 0],'LineWidth', 1.5);
    % Plotting rdft recovery
    hold on;
    plot(plot_x(isNZ2), rdft_recovery(isNZ2),  '-x', 'Color', [0, 1, 0],'LineWidth',1.5);
    % Plotting gauss recovery
    hold on;
    plot(plot_x(isNZ), gauss_recovery(isNZ), '-*', 'Color', [1, 0.4, 0.7],'LineWidth',1.5);
    % Formatting the plot
    xlim([initial_sparsity final_sparsity])
    ylim([0 1])
    xlabel('Sparsity k')
    ylabel('Recovery percentage')
    legend(type, string(p_ary) + '-ary BCH', 'RDFT', 'Complex Gaussian')
    grid on
    title('Perfect recovery percentage of noiseless ' + string(n) + '× 1 signals')

end