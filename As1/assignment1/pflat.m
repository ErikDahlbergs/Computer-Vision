function [output] = pflat(input)
%   Computer Exercise 1:
%   Divides the homogeneous coordinates with their last entry for points of any dimensionality
output = input ./ repmat(input(end,:), size(input, 1), 1);
end