% my first main function
function my_first_main_function

    clear all;
    close all;
    
    input1 = 23
    input2 = 14
    
    [s,m] = multANDsum(input1, input2);
   
%-------------------------------------------------------------------    
    my_array = [1:5;6:10];
    copy = my_array;
    
    s = size(my_array);
    
    % example of a for-loop
    for i = 1:s(1) % loop over rows
        for j = 1:s(2) % loop over columns
            my_array(i,j) = my_array(i,j)^2;
        end
    end
    
    % but this can be done faster:
    copy = copy.^2;
    % avoid using slow for-loops in MATLAB, if possible!
    
    % an example for thresholding:
    % All values lower or equal 50 will be set to 1 and to 0 else
    mask = copy <= 50;
    
end

% this function can do multiplication and summation in one step
function [summation, multiplication] = multANDsum(input1, input2)

    summation = input1 + input2

    multiplication = input1 * input2
    multiplication = my_multiplication(input1, input2)    
end

%------------------------------------------------------------------------
function [output] = my_multiplication(input1, input2)
    output = input1 * input2;
end

%-----------------------------------------------------------------------
% demo for command line: variables
%a = zeros(5,10)
%a(2,2) = 23.5
%
%b = [2:20]
%b = [2:2:20]


