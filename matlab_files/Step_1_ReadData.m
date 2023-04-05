% This script read the output produced during Ex_5 and shows the result of producing a normally distributed dataset

% Computational procedure:
% - Binary tree with bounds from -4 to +4 is created
% - Dataset is composed of 1E4 points distributed according to a gaussian
% centered at zero and standard deviation of 1/10th the size of domain.
% - Maximum depth of tree is k = 7 which leads to 128 cells.
% - To produce the data output, three queries are performed.
% These consist of asking the user for 3 axial position and the result
% return the node containing that point.
% - From that node, we can extract all the indices of the points inside that
% node and save to file

clear all
close all
clc

saveFig = 1;

% Get data used as input to binary tree:
% =========================================================================
fileName = '../output_files/data.txt';
dataset = load(fileName);

% Get output data produced by binary tree during queries:
% =========================================================================
fileName = dir('../output_files/result_*');
for qq = 1:numel(fileName)
    try
        xq{qq} = load(['../output_files/',fileName(qq).name]);
    catch
        disp(['File note found: ',fileName(qq).name]);
        continue;
    end
end

% Create histogram:
% =========================================================================
if 0
    % Create analytic PDF:
    x = linspace(-4,4,1e3);
    y = normalDistPDF(x,0,8/10);

    figure('color','w')
    hold on
    box on
    h(1) = histogram(dataset,100,'Normalization','pdf');
    h(2) = plot(x,y,'r','LineWidth',2);
    
    legendText{1} = ['C++ armadillo data'];
    legendText{2} = ['Theoretical PDF'];
    
    title('Normally distributed data produced from C++ code')
    hL = legend(h,legendText);
    hL.Interpreter = 'Latex';
    hL.FontSize = 13;
end

%% Plot results from C++ binary tree search:
% =========================================================================
figure('color','w')
hold on
box on
grid on;

markerColor = {'k.','r.','bl.','g.','m.','c.'};

% Plot the input data:
hin(1) = plot(dataset,'k.');

% Plot the data found using the binary tree:
for qq = 1:numel(// ================================================================================================================
)
    try
        ix     = xq{qq}+1;
        dataset_ix = dataset(ix);
        
        val = circshift(markerColor,qq-1);
        hq(qq) = plot(ix,dataset_ix,val{1});
    catch
        disp(['Data not found']);
        continue;
    end
end

% hL = legend([hin(1),hq(1)],'Input data','from binary search');
% set(hL,'interpreter','latex')
set(gca,'fontSize',14)
ylim([-1,+1]*2)

% Save figure:
% =========================================================================
if saveFig
    folderName = '../output_files/';
    figureName = 'Step_1_find_results';
    
    % PDF figure:
    exportgraphics(gcf,[folderName,figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[folderName,figureName,'.tiff'],'Resolution',600) 
end
