% This function is the main function and once called, will call the rest of
% the defined function to output the required graphs
function main()
    % Keep track of execution time
    tic;
    % Clear and close all previously executed figures
    clear all;
    close all;
    clf;
    % Turn warnings off
    warning('off');
    % Run program to generate the first figure
    compressed_sensing(7381, 1, 'Singer', 170, 330, 10, 1);
    % Run program to generate the second figure
    compressed_sensing(1573, 1, 'Macfarland', 24, 70, 2, 2);
    % Print out the execution time of the program
    toc;
end