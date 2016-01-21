
n_seq=5;
n_data_per_seq = 11;
dream_data_new = [];

ts1 = cell(n_seq, 1);
ts2 = cell(n_seq, 1);
for i=1:n_seq
    dream_start = (i-1)*n_data_per_seq+1;
    dream_end = i*n_data_per_seq;
    ts1{i}=timeseries(dream_data(dream_start:dream_end,:), 1:n_data_per_seq);
    ts2{i}=resample(ts1{i}, 1:0.5:n_data_per_seq);
    dream_data_new = [dream_data_new; [ [1:0.5:n_data_per_seq]' ts2{i}.data] ];
end

