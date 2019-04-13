clc
clear all

%
% Only needed because the task explicitly asks for a comparison with the
% matlab fft() function
%

n = 2;
times = [];
sizes = [];

%warm matlab up...
inputArr = rand(1, n);
start = tic();
fft(inputArr);
duration = toc(start);

f = fopen('matlab.csv', 'w');
while duration <= 1.0
    n = n * 2;
    inputArr = rand(1, n);
    
    start = tic();
    for i = 1:5
        fft(inputArr);
    end
    duration = toc(start) / 5;
    fprintf(f, '%i,%f\n', n, duration);
end
fclose(f);
