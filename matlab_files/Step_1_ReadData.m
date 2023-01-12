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

for qq = 1:3
    try
        fileName = ['result_',num2str(qq-1),'.txt'];
        xq{qq} = load(fileName);
    catch
        disp(['File note found: ',fileName]);
        continue;
    end
end

cd(home_dir);

% Create analytic PDF:
x = linspace(-4,4,1e3);
y = normalDistPDF(x,0,8/10);

% Create histogram:
if 0
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
end

% Plot raw data with queries from C++ binary search tree:
figure('color','w')
hold on
box on
hq(1) = plot(set,'k.');
grid on;

for qq = 1:3
    try
        ix     = xq{qq}+1;
        set_ix = set(ix);
    
        hq(qq+1) = plot(ix,set_ix,'ro');
        plot(ix,set_ix,'ro');
    catch
        disp(['Data not found']);
        continue;
    end
end

