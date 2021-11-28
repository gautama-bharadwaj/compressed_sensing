function main = main()
    tic;
    clf;
    warning('off');
%     compressed_sensing(7381, 1, 'Singer', 170, 330, 10, 1);
    compressed_sensing(1573, 1, 'Macfarland', 24, 70, 2, 1);
    toc;
end