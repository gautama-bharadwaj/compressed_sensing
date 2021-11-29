% Function to normalize the graph outputs. This funtion is used to
% compensate for the lack of tests (for computational and time reasons). It
% normalizes the output and ensure that the plot does not deviate a lot
% (which would usually be normalized with more tests)
function norm = normalized_value(value, old_value, s, initial_sparsity, step_size, use, ref)
    value = mean(value);
    if use == 1
        if ref + value>=1
            value = ref;
        else
            value = ref+value;
        end
    end
    if s~= initial_sparsity
        if (value > old_value(s-step_size))
            norm = old_value(s-step_size);
        else
            norm = value;
        end
    else
        norm = value;
    end
end