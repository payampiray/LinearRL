function out = def(action)

codedir = pwd;

switch action
    case 'pipedir'
    out = codedir;
    case 'bigpipedir'
    out = codedir;    
    case 'tempdir'
    out = codedir;
        
    case 'addpath'
        addpath(fullfile(codedir,'tools'));        
        
    case 'col1'
        out = [1 .2 .2];
    case 'col'
        out = [.2 0 0; .8 .4 .1; 1 .2 .2] ;
    case 'alf'
        out = .6;
    case 'fs'
        out = 12;
    case 'fsy'
        out = 14;
    case 'fsA'
        out = 18;
    case 'fn'
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