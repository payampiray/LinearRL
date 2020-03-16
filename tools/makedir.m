function [success]=makedir(fullpath)
if(~exist(fullpath,'dir')), mkdir(fullpath); end;    
success = true;
end