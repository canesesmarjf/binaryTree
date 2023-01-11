% This script read the output produced during Ex_5 and shows the result of producing a normally distributed set

clear all
close all
clc

% Read data produced and saved by C++ code:
home_dir = cd;
cd ..
cd output_files

fileName = 'data.txt';
set = load(fileName);

% Create analytic PDF:
x = linspace(-4,4,1e3);
y = normalDistPDF(x,0,8/10);

% Create histogram:
figure('color','w')
hold on
box on
h(1) = histogram(set,100,'Normalization','pdf');
h(2) = plot(x,y,'r','LineWidth',2);

legendText{1} = ['C++ armadillo data'];
legendText{2} = ['Theoretical PDF'];

title('Normally distributed data produced from C++ code')
hL = legend(h,legendText);
hL.Interpreter = 'Latex';
hL.FontSize = 13;

cd(home_dir);

