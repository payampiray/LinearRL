function out = def(action)

st = dbstack('-completenames');
sciencedir = fileparts(fileparts(fileparts(fileparts(st(1).file))));
% [~,rootname] = fileparts(fileparts(fileparts(st(1).file)));
rootname = 'linear_rl';
pipedir = fullfile(sciencedir,'xgit',rootname,'shared');
bigpipedir = fullfile(sciencedir,'xgit',rootname,'shared');
codedir = fullfile(sciencedir,'xgit',rootname,'sgared');
tempdir = fullfile(sciencedir,'xloc',rootname);
figdir  = fullfile(sciencedir,'xgit',rootname,'figs');

switch action
    case 'figdir'
        out = figdir;
        
    case 'pipedir'
    out = pipedir;
    case 'bigpipedir'
    out = bigpipedir;    
    case 'tempdir'
    out = tempdir;
%     case 'thisdir'
%         fname = st(2).name;
%         fdir = fullfile(tempdir,fname);
%         if ~exist(fdir,'dir')
%             mkdir(fdir);
%         end
%         out = fdir;
        
    case 'addpath'
        addpath(fullfile(codedir,'tools'));
        addpath(fullfile(codedir,'altmany-export_fig'));
        out = fullfile(codedir,'tools');
        
        
    case 'unp'
        out = 'Unpredictability';
    case 'vol'
        out = 'Volatility';
    case 'lr'
        out = 'Learning rate';
    case 'col1'
        out = [1 .2 .2];
    case 'col'
        out = [0    0.4470    0.7410; 0.8500    0.3250    0.0980];
        out = [.2 0 0; .8 .4 .1; 1 .2 .2] ;
%         out = [.8 .4 .1; 1 .2 .2] ;
    case 'alf'
        out = .6;
    case 'fs'
        out = 12;
    case 'fsy'
        out = 14;
    case 'fsA'
        out = 18;
    case 'fn'
        out = 'Calibri';
        out = 'Helvetica';
    case 'xsA'        
        out = -.15;
    case 'ysA'
        out = 1.1;
    case 'abc'
        out = 'abcdefghijkl';        
        
    otherwise
        error('!');
end
end