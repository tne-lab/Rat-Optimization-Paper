function full_tbl = load_opt_data(file_dir, func)
    addpath(genpath('C:\Users\Evan\Documents\GitHub\eelib'))
    if nargin < 2
        func = @read_ss_pyb_pr;
    end

    data_files = dir(file_dir);
    data_files = data_files(3:end);
    ss_files = cell(length(data_files),1);
    nrow = 0;
    
    for i=1:length(ss_files)
        ss_files{i} = func(fullfile(data_files(i).folder, data_files(i).name),false);
        nrow = nrow + length(ss_files{i}.data);
    end
    
    T = struct2table(ss_files{1}.data);
    types = varfun(@class, T, 'OutputFormat', 'cell');
    full_tbl = table('Size',[nrow,width(T)+4],'VariableTypes',horzcat({'string','string','string','double'},types),'VariableNames',["Rat","Session","Protocol","SessionNum",T.Properties.VariableNames]);
    
    %%
    M = containers.Map('KeyType','char','ValueType','double');
    row = 1;
    for i=1:length(ss_files)
        T = struct2table(ss_files{i}.data);
        if isfield(ss_files{i}, 'subject_config')
            full_tbl.Opt(row:row+height(T)-1) = ss_files{i}.subject_config.amps;
        else
            full_tbl.Opt(row:row+height(T)-1) = "";
        end
        full_tbl(row:row+height(T)-1,5:end-1)=T;
        full_tbl.Rat(row:row+height(T)-1) = ss_files{i}.subject;
        full_tbl.Session(row:row+height(T)-1) =num2str(i);
        full_tbl.Protocol(row:row+height(T)-1) = ss_files{i}.protocol;
        if isKey(M, ss_files{i}.subject)
            M(ss_files{i}.subject) = M(ss_files{i}.subject) + 1;
        else
            M(ss_files{i}.subject) = 1;
        end
        full_tbl.SessionNum(row:row+height(T)-1) = M(ss_files{i}.subject);
        row = row + height(T);
    end
end