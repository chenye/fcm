n_seq = 1;

alg_name = {'ACOR', 'DE', 'PSO', 'RCGA'};
test_name = {'n10-1', 'n10-2', 'n10-3', 'n10-4', 'n10-5', 'n50-1', 'n50-2', 'n50-3', 'n50-4', 'n50-5'};
result = cell(4, 10);
avg_data_error = zeros(4,10);
avg_model_error = zeros(4,10);
avg_auc = zeros(4,10);
avg_aupr = zeros(4,10);
for i=1:4
    for j=1:10
        filename = sprintf('result-s-%s-para1-np100-sp0.00-test-%s-seq-%d-10.txt', alg_name{i}, test_name{j}, n_seq);
        result{i,j}=importdata(filename);
        avg_data_error(i,j) = mean(result{i,j}(:,1));
        avg_model_error(i,j) = mean(result{i,j}(:,2));
        avg_auc(i,j) = mean(result{i,j}(:,6));
        avg_aupr(i,j) = mean(result{i,j}(:,7));
    end
end

figure;
subplot(411);
bar(avg_data_error');
subplot(412);
bar(avg_model_error');
subplot(413);
bar(avg_auc');
subplot(414);
bar(avg_aupr');