function write_txt_file(formatstr,info_str, formatdata,info_data,logpath,file_name)
% write_txt_file - creates a .txt file with desired information as
% specified by "info_str)  info_str,
% function calls :
%       WRITE_TXT_FILE(LOGPATH,FILENAME,INFO_STR)
%           LOGPATH = pathname (as a string) to the folder where the file is to be saved, 
%           ex: 'C:\Desktop'
%           FILENAME = name of the file (as a string), ex: 'test.txt'
%           INFO_STR = information that you want to be written in the file
%           ex: a variable or a string such as 'testing file'

    % adding '\' between pathname and filename
    if(logpath(end)=='\')
          file = strcat(logpath,file_name); 
    else
    file = strcat(logpath,'\',file_name);  
    end
    % open file
    if exist(file)
        fid = fopen(file,'wt');
    else
        fid = fopen(file,'wt');
    end
    if fid==-1
        errdlg('File cannot be created!');
    end
    % print info in the file
%     fprintf(fid,'**************************************************** \n');
%     fprintf(fid,'\n');
fprintf(fid, formatstr,info_str);
 fprintf(fid,'\n');


[m,n]=size(info_data);


for i=1:m
    fprintf(fid, formatdata,info_data(i,:)); 
    fprintf(fid,'\n');
end
%     fprintf(fid,'\n');  info_str,
%     fprintf(fid,'\n');
%     fprintf(fid,'**************************************************** \n');
     fclose(fid);
end  % func end