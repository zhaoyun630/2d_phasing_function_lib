% this piece of code is used to create an obj and support

% written by Yun Zhao, July 23, 2014.
% Good good study, Day day up ^_^


disp('please give the image file: \n');
imagename = input(prompt,'s');

A = imread('test.jpg');
B = rgb2gray(A);
