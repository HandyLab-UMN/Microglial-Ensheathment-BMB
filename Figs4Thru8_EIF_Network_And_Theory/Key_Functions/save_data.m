function [] = save_data(data_dir,random_seed, rates_trial, ...
    spkMatrix,params,stim_num)


file_name = sprintf('EIF_sen_%.2f_%d.mat',params.s_en_vec(stim_num),random_seed);
name_full = strcat(data_dir,file_name);
save(name_full,'spkMatrix','rates_trial','params','-v7.3')

end

